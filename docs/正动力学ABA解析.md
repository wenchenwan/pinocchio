# 正动力学 ABA 解析：`algorithm/aba.hxx`

> **ABA**（Articulated Body Algorithm，铰接体算法，Featherstone）解决**正动力学**问题：
> 给定力矩 $(q,\dot q,\tau)$，求加速度 $\ddot q$。$O(n)$ 三趟递推，是 $\ddot q=M^{-1}(\tau-b)$ 的高效等价实现。
> 本文讲透 [src/algorithm/aba.hxx](../include/pinocchio/src/algorithm/aba.hxx)：铰接体惯量的物理意义、三趟递推、
> $U/D/UD^{-1}$ 三件套怎么把 $M^{-1}$ 隐式解出来。
> 配套：[逆动力学RNEA解析.md](逆动力学RNEA解析.md)（对偶力传播、重力技巧）、[质量矩阵CRBA解析.md](质量矩阵CRBA解析.md)、[空间代数运算解析.md](空间代数运算解析.md)。

---

## 目录
- [1. 正动力学问题与两条路线](#1-正动力学问题与两条路线)
- [2. 核心概念：铰接体惯量（Articulated-Body Inertia）](#2-核心概念铰接体惯量articulated-body-inertia)
- [3. 三趟递推骨架](#3-三趟递推骨架)
- [4. Pass 1 前向：运动学 + 惯量/偏置力初始化](#4-pass-1-前向运动学--惯量偏置力初始化)
- [5. Pass 2 后向：铰接体惯量装配 + $U/D/UD^{-1}$](#5-pass-2-后向铰接体惯量装配--udud-1)
- [6. Pass 3 前向：解出加速度](#6-pass-3-前向解出加速度)
- [7. 为什么这样就等价于 $\ddot q=M^{-1}(\tau-b)$](#7-为什么这样就等价于-ddot-qm-1tau-b)
- [8. `computeMinverse`：直接求 $M^{-1}$](#8-computeminverse直接求-m-1)
- [9. 复杂度与约定（WORLD/LOCAL）](#9-复杂度与约定worldlocal)
- [10. 一句话总结](#10-一句话总结)

---

## 1. 正动力学问题与两条路线

**已知** $q,\dot q$ 和关节力矩 $\tau$，**求**加速度：

$$
\boxed{\,\ddot q=M(q)^{-1}\big(\tau-\underbrace{C\dot q-g}_{-b(q,\dot q)}\big)=M^{-1}(\tau-b)\,}
$$

两条实现路线：

| 路线 | 步骤 | 复杂度 |
|---|---|---|
| **朴素** | CRBA 求 $M$ → Cholesky 分解 → RNEA 求 $b$ → 解 $M\ddot q=\tau-b$ | $O(n^2\sim n^3)$ |
| **ABA** | 铰接体三趟递推，隐式解出 $\ddot q$，**从不显式组装 $M$** | $O(n)$ |

ABA 的精髓：用**铰接体惯量**沿树传播，把"求逆"分解成每个关节局部的 $1\times1$/小块求逆，避免了整机 $n\times n$ 求逆。对高自由度机器人（人形 30+ DoF）优势巨大。

---

## 2. 核心概念：铰接体惯量（Articulated-Body Inertia）

普通刚体惯量 $I_i$ 描述"单个连杆"。而**铰接体惯量** $I_i^A$ 描述的是：**连杆 $i$ 连同它下游整个子树（通过关节铰接）作为一个整体，对施加在 $i$ 上的力表现出的等效惯量**。

直觉：你推一节机械臂的某个连杆，它"感觉有多重"不只取决于该连杆本身，还取决于挂在它下面的所有连杆——但下游关节是**可自由转动的**（受 $\tau$ 驱动），所以子树不是刚性附着，而是"铰接"地跟随。$I_i^A$ 精确刻画了这种"带可动关节的等效惯量"。

铰接体的力-加速度关系（含偏置力 $p_i^A$）：

$$
f_i = I_i^A\,a_i + p_i^A
$$

- $I_i^A$（代码 `data.oYaba[i]`，6×6）：铰接体惯量。
- $p_i^A$（偏置力 bias force，藏在 `data.of[i]` 的传播里）：即使 $a_i=0$，因子树内部速度耦合和已知力矩产生的力。

Pass 2 的核心就是**递推地把子树惯量"塌缩"进父连杆**，得到每个 $I_i^A$。

---

## 3. 三趟递推骨架

```cpp
abaWorldConvention(model, data, q, v, tau) {
  data.oa_gf[0] = -model.gravity;          // 重力技巧（同 RNEA §4）

  for (i=1; i<njoints; ++i) Pass1::run(...);   // 前向：运动学 + I^A/f 初值
  for (i=njoints-1; i>0; --i) Pass2::run(...);  // 后向：装配 I^A，算 U,D⁻¹,UDinv
  for (i=1; i<njoints; ++i) Pass3::run(...);   // 前向：解 ddq
  return data.ddq;
}
```

注意与 RNEA 的**两趟**不同，ABA 是**前-后-前三趟**：多出的第三趟是因为加速度要"先知道父的加速度才能解自己"，而父的加速度又依赖已塌缩的子树惯量——必须先后向装配惯量，再前向解加速度。

> 这里展示的是 **WORLD 约定**版本（`AbaWorldConvention*`）：所有量都通过 `oMi[i].act(...)` 变换到世界系再算，数值更稳、便于复用 `data.J`。另有 LOCAL 约定版本，数学等价。

---

## 4. Pass 1 前向：运动学 + 惯量/偏置力初始化

```cpp
jmodel.calc(jdata, q, v);
data.liMi[i] = model.jointPlacements[i] * jdata.M();
data.oMi[i]  = (parent>0) ? data.oMi[parent]*data.liMi[i] : data.liMi[i];

jmodel.jointCols(data.J) = data.oMi[i].act(jdata.S());   // 世界系运动子空间列

ov = data.oMi[i].act(jdata.v());                          // 世界系速度
if (parent>0) ov += data.ov[parent];

data.oa_gf[i] = data.oMi[i].act(jdata.c());               // 偏置加速度
if (parent>0) data.oa_gf[i] += (data.ov[parent] ^ ov);    // 速度耦合项

data.oinertias[i] = data.oYcrb[i] = data.oMi[i].act(model.inertias[i]);
data.oYaba[i] = data.oYcrb[i].matrix();                   // ★ I^A 初值 = 本体惯量
data.oh[i] = data.oYcrb[i] * ov;                          // 世界系动量
data.of[i] = ov.cross(data.oh[i]);                        // ★ 偏置力初值 = v ×* (I v)
```

- 前半段是 [正向运动学](正向运动学解析.md) 的速度/加速度递推，只是全部搬到**世界系**（`ov, oa_gf`）。
- **两个关键初值**：
  - `oYaba[i] = 本体惯量 I_i`：铰接体惯量的**种子**，Pass 2 会把子树惯量累加进来。
  - `of[i] = ov ×* (I ov)`：偏置力的种子，就是 RNEA 里的陀螺力项（此刻还没有 $\tau$ 和子树贡献）。
- `oa_gf` 含 $-g$（根初值 $-g$ 传播下来），同 RNEA 重力技巧。

---

## 5. Pass 2 后向：铰接体惯量装配 + $U/D/UD^{-1}$

这是 ABA 的心脏。逆序遍历，对每个关节做"局部消元 + 塌缩给父"：

```cpp
Ia = data.oYaba[i];                                      // 当前铰接体惯量（含已并入的子树）
Jcols = jmodel.jointCols(data.J);                        // 世界系 S_i

// ① 关节力残差 u = τ_i − Sᵀ f_i
u_i = jmodel.jointVelocitySelector(data.u);
u_i -= Jcols.transpose() * data.of[i].toVector();

// ② 三件套
jdata.U()    = Ia * Jcols;                               // U = I^A S       (6×nv)
jdata.StU()  = Jcols.transpose() * jdata.U();            // D = Sᵀ I^A S     (nv×nv)
jdata.StU().diagonal() += armature;                      // + 转子惯量
matrix_inversion(jdata.StU(), jdata.Dinv());             // D⁻¹
jdata.UDinv() = jdata.U() * jdata.Dinv();                // U D⁻¹           (6×nv)

// ③ 把"消掉本关节自由度后"的铰接体惯量与偏置力塌缩给父
if (parent > 0) {
  Ia -= jdata.UDinv() * jdata.U().transpose();           // I^A_new = I^A − U D⁻¹ Uᵀ
  fi += Ia * data.oa_gf[i] + jdata.UDinv() * u_i;        // 偏置力更新
  data.oYaba[parent] += Ia;                              // 累加惯量给父
  data.of[parent]    += fi;                              // 累加偏置力给父
}
```

### 5.1 $U$、$D$、$D^{-1}$ 是什么

- **$D_i=S_i^\top I_i^A S_i$**（代码 `StU`，$n_{v_i}\times n_{v_i}$）：铰接体惯量在**本关节运动子空间**上的投影，即"沿关节 $i$ 自由度方向的等效惯量"。对 1-DoF 关节就是个标量。加上 `armature`（转子反射惯量）改善条件数。
- **$U_i=I_i^A S_i$**（6×nv）：惯量作用在子空间上的耦合矩阵，桥接"关节加速度"与"6D 空间力"。
- **$D_i^{-1}$**：正是这一步把"求逆"**局部化**到关节维度——整机 $M^{-1}$ 被分解成一串小 $D_i^{-1}$。

### 5.2 惯量塌缩公式 $I^A\!-\!UD^{-1}U^\top$（Schur 补！）

$$
\boxed{\,I_i^{A,\text{proj}} = I_i^A - U_i D_i^{-1} U_i^\top = I_i^A - I_i^A S_i(S_i^\top I_i^A S_i)^{-1}S_i^\top I_i^A\,}
$$

这是一个 **Schur 补 / 投影**：从铰接体惯量里"减掉关节 $i$ 能自由让开的那部分"。物理含义：因为关节 $i$ 可自由转动（在 $\tau_i$ 驱动下），子树对父连杆的等效惯量要扣除沿 $S_i$ 方向"被关节吸收"的分量。剩下的 $I_i^{A,\text{proj}}$ 才是**通过关节 $i$ 传给父连杆**的有效惯量，累加进 `oYaba[parent]`。

偏置力同理更新后累加给父。逆序保证处理 $i$ 时其所有子树已塌缩完毕，所以 `oYaba[i]` 进来时已是"$i$ + 全下游"的铰接体惯量。

---

## 6. Pass 3 前向：解出加速度

```cpp
data.oa_gf[i] += data.oa_gf[parent];                     // 父加速度已知，加进来
ddq_i = jdata.Dinv() * u_i - jdata.UDinv().transpose() * data.oa_gf[i].toVector();
data.oa_gf[i] += J_cols * ddq_i;                         // 用解出的 ddq 更新本体加速度
data.oa[i] = data.oa_gf[i] + model.gravity;              // 还原真实加速度（去掉 −g）
```

**关节加速度求解公式**：

$$
\boxed{\,\ddot q_i = D_i^{-1}u_i - (U_iD_i^{-1})^\top a_{\lambda(i)} = D_i^{-1}\big(u_i - U_i^\top a_{\lambda}\big)\,}
$$

- 正序遍历：处理 $i$ 时父加速度 $a_\lambda$ 已解出。
- 直觉：本关节加速度 = "自己的力残差 $u_i$ 除以等效惯量 $D_i$" 减去 "父连杆加速度通过耦合 $U_i$ 拽着走的部分"。
- 解出 $\ddot q_i$ 后立刻更新本体空间加速度 $a_i=a_\lambda+S_i\ddot q_i$（`oa_gf[i] += J_cols*ddq`），供子节点使用。
- 最后 `oa[i] = oa_gf[i] + gravity` 把重力技巧引入的 $-g$ 还原，得到**真实**空间加速度。

---

## 7. 为什么这样就等价于 $\ddot q=M^{-1}(\tau-b)$

ABA 本质是对 KKT/牛顿-欧拉方程组做**稀疏 $LDL^\top$ 消元**，只是沿运动树的拓扑顺序进行：

- Pass 2 的 $I^A\!-\!UD^{-1}U^\top$ 是**块 Gauss 消元**——把子节点自由度从系统里消去，等价于对 $M$ 做 Cholesky 时的 Schur 补。$D_i$ 恰好对应 [Cholesky 分解](multibody子系统解析.md) 里 $M=UDU^\top$ 的对角块 $D$。
- Pass 3 是**回代**（back-substitution），解出各 $\ddot q_i$。
- 因运动树的稀疏性（`parents[i]<i`，每个关节只耦合支撑链），消元只在父子间进行，故 $O(n)$ 而非稠密的 $O(n^3)$。

所以 ABA = "针对刚体树稀疏结构特化的 $M^{-1}(\tau-b)$ 求解器"，与朴素路线数学恒等，但快一个数量级。

---

## 8. `computeMinverse`：直接求 $M^{-1}$

`aba.hxx` 还提供 `computeMinverse(q)`：**不组装 $M$ 就直接得到 $M^{-1}$**。复用 ABA 的 $U/D^{-1}/UD^{-1}$ 三件套：把 Pass 3 的回代对**每个自由度**跑一遍（相当于解 $M^{-1}$ 的每一列），$O(n^2)$。比"CRBA 求 $M$ + 求逆 $O(n^3)$"更快，用于需要 $M^{-1}$ 的场合（如某些约束动力学、ABA 导数 $\partial\ddot q/\partial\tau=M^{-1}$）。

---

## 9. 复杂度与约定（WORLD/LOCAL）

- **复杂度** $O(n)$：三趟线性遍历，每关节一次小块求逆 $D_i^{-1}$（1-DoF 关节是标量除法）。
- **WORLD vs LOCAL 约定**：本文的 `AbaWorldConvention*` 把所有量变换到世界系（数值稳定、复用 `data.J`、便于导数）；LOCAL 版在各关节局部系递推，少几次 `oMi.act`。二者结果相同。默认接口 `aba(...)` 会选其一。
- **不支持 mimic 关节**（`assert(model.check(MimicChecker()))`）——mimic 破坏了每关节独立求逆的结构。
- 与 RNEA 共享重力技巧、共享 `oYcrb`/`J` 等中间量（故 `computeABADerivatives` 能顺带复用）。

---

## 10. 一句话总结

ABA = **铰接体惯量三趟 $O(n)$ 递推**：

$$
\underbrace{\text{根→叶}}_{\text{FK}+I^A,p^A\ \text{初值}}\ \Rightarrow\
\underbrace{\text{叶→根}}_{D=S^\top I^A S,\ I^A\mathrel{-}=UD^{-1}U^\top\ \text{塌缩给父}}\ \Rightarrow\
\underbrace{\text{根→叶}}_{\ddot q_i=D_i^{-1}(u_i-U_i^\top a_\lambda)}
$$

- **铰接体惯量** $I^A$ = "连杆 + 可动子树"的等效惯量，靠 Schur 补 $I^A-UD^{-1}U^\top$ 沿树塌缩；
- $D_i^{-1}$ 把整机求逆**局部化**到关节维度，这就是 ABA 达到 $O(n)$ 的关键；
- 数学上等价于对 $M\ddot q=\tau-b$ 做**沿树稀疏 $LDL^\top$ 消元**，比朴素 $M^{-1}$ 快一个量级。
