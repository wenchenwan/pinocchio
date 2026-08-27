# 质量矩阵 CRBA 解析：`algorithm/crba.hxx`

> **CRBA**（Composite Rigid Body Algorithm，组合刚体算法，Featherstone 第 6 章）计算**关节空间惯性矩阵** $M(q)$（质量矩阵）。
> 本文讲透 [src/algorithm/crba.hxx](../include/pinocchio/src/algorithm/crba.hxx)：组合刚体惯量的意义、
> $M_{ij}=S_i^\top(\sum I)S_j$ 的推导、为什么只算上三角、如何借世界系与力集实现。
> 配套：[逆动力学RNEA解析.md](逆动力学RNEA解析.md)、[正动力学ABA解析.md](正动力学ABA解析.md)、
> [multibody 子系统解析 §3.5（M 的稀疏 Cholesky）](multibody子系统解析.md)、[空间代数运算解析.md](空间代数运算解析.md)。

---

## 目录
- [1. 质量矩阵是什么](#1-质量矩阵是什么)
- [2. 核心概念：组合刚体惯量（Composite Rigid Body Inertia）](#2-核心概念组合刚体惯量composite-rigid-body-inertia)
- [3. 关键公式：$M_{ij}=S_i^\top Y_i^c S_j$](#3-关键公式m_ijs_itop-y_ic-s_j)
- [4. 两趟递推骨架](#4-两趟递推骨架)
- [5. Pass 1 前向：运动学 + 世界系惯量](#5-pass-1-前向运动学--世界系惯量)
- [6. Pass 2 后向：组合惯量累加 + 填 $M$](#6-pass-2-后向组合惯量累加--填-m)
- [7. 为什么只填上三角、`nvSubtree` 的作用](#7-为什么只填上三角nvsubtree-的作用)
- [8. 与质心动量矩阵 $A_g$、ABA、Cholesky 的关系](#8-与质心动量矩阵-a_gabacholesky-的关系)
- [9. 复杂度](#9-复杂度)
- [10. 一句话总结](#10-一句话总结)

---

## 1. 质量矩阵是什么

$M(q)$ 是拉格朗日动力学 $M\ddot q+b=\tau$ 里的**关节空间惯性矩阵**（$n_v\times n_v$，对称正定）：

- 动能 $T=\tfrac12\dot q^\top M(q)\dot q$；
- $M_{ij}$ = "关节 $j$ 加速时，关节 $i$ 感受到的惯性耦合力矩"；
- 正动力学要 $M^{-1}$，最优控制/碰撞/操作空间控制处处要它。

CRBA 用 $O(nd)$ 直接组装出 $M$ 的上三角（$d$=树深度），是求 $M$ 的标准方法。

---

## 2. 核心概念：组合刚体惯量（Composite Rigid Body Inertia）

**组合刚体惯量** $Y_i^c$（代码 `data.oYcrb[i]`）= 以连杆 $i$ 为根的**整个子树，当作一个刚体焊死**后的合惯量：

$$
\boxed{\,Y_i^c=\sum_{j\in\text{subtree}(i)} I_j\,}
$$

（各 $I_j$ 需变换到同一坐标系再相加。）

与 ABA 的**铰接体惯量** $I^A$ 对照：

| | 组合刚体惯量 $Y^c$（CRBA） | 铰接体惯量 $I^A$（ABA） |
|---|---|---|
| 子树关节视为 | **焊死（刚性）** | **可自由转动（铰接）** |
| 累加方式 | 简单相加 $\sum I_j$ | Schur 补塌缩 $I^A-UD^{-1}U^\top$ |
| 用途 | 组装 $M$ | 隐式求 $M^{-1}$ |

CRBA 的"焊死"假设正是它简单的原因：$M_{ij}$ 描述的是**惯性耦合**，等价于把两关节间的连杆看成刚性时的惯量投影。

---

## 3. 关键公式：$M_{ij}=S_i^\top Y_i^c S_j$

设关节 $i$ 是关节 $j$ 的祖先（$i\preceq j$，即 $i$ 在 $j$ 的支撑链上）。质量矩阵元：

$$
\boxed{\,M_{ij}=S_i^\top\,Y_{\max(i,j)}^c\,S_j\,}
$$

推导直觉：$M_{ij}$ = 关节 $j$ 单位加速度在关节 $i$ 轴上引起的力矩。关节 $j$ 加速 → 带动**子树 $j$**（组合惯量 $Y_j^c$）产生空间力 $Y_j^c S_j$ → 这个力沿支撑链传到祖先关节 $i$，在 $S_i$ 上的投影就是 $S_i^\top(Y_j^c S_j)$。因为 $j$ 是后代，$\max(i,j)=j$，用 $Y_j^c$。

代码里对应两步（Pass 2）：

```cpp
motionSet::inertiaAction(data.oYcrb[i], J_cols, Ag_cols);   // F = Y_i^c · S_i  （力集）
data.M.block(...) = J_cols.transpose() * data.Ag.middleCols(...);  // M块 = S^T · F
```

- **第一步** $F_i=Y_i^c S_i$：组合惯量作用在运动子空间上，得到一批**空间力**（力集 force set，见 [multibody §12.3](multibody子系统解析.md)）。存进 `data.Ag`（此列块也正是质心动量矩阵的列，见 §8）。
- **第二步** $M[\ldots]=S^\top F$：力集对支撑链上各关节子空间的投影，一次性填出 $M$ 的一整条（本关节行 × 子树列）。

---

## 4. 两趟递推骨架

```cpp
crbaWorldConvention(model, data, q) {
  for (i=1; i<njoints; ++i)  Pass1::run(...);   // 前向：FK + 世界系惯量
  for (i=njoints-1; i>0; --i) Pass2::run(...);   // 后向：组合惯量累加 + 填 M
  // （含 mimic 关节时有额外的 MimicStep 补齐）
}
```

- **Pass 1 正序**：算 `oMi`、把 $S_i$ 和 $I_i$ 搬到世界系。
- **Pass 2 逆序**：子树组合惯量自叶向根累加，同时用 $M_{ij}=S^\top Y^c S$ 填矩阵。
- 只需 $q$（$M$ 只依赖位形，不依赖速度）。

---

## 5. Pass 1 前向：运动学 + 世界系惯量

```cpp
jmodel.calc(jdata, q);
data.liMi[i] = model.jointPlacements[i] * jdata.M();
data.oMi[i]  = (parent>0) ? data.oMi[parent]*data.liMi[i] : data.liMi[i];

jmodel.jointExtendedModelCols(data.J) = data.oMi[i].act(jdata.S());  // 世界系 S_i → data.J
data.oYcrb[i] = data.oMi[i].act(model.inertias[i]);                  // 世界系本体惯量（组合惯量种子）
```

- `data.J` 的列存**世界系运动子空间** $\!^oS_i$（这也是整机雅可比！CRBA 顺带把它算好）。
- `oYcrb[i]` 初值 = 本体惯量搬到世界系，作为组合惯量 $Y_i^c$ 的**种子**（Pass 2 累加子树）。
- 全程 WORLD 约定：所有惯量、子空间在世界系表达，累加 $\sum I_j$ 时无需再做坐标变换（同系直接相加）。

---

## 6. Pass 2 后向：组合惯量累加 + 填 $M$

```cpp
Ag_cols = jmodel.jointCols(data.Ag);
J_cols  = jmodel.jointExtendedModelCols(data.J);           // 世界系 S_i

// ① 力集 F_i = Y_i^c · S_i
motionSet::inertiaAction(data.oYcrb[i], J_cols, Ag_cols);

// ② 填 M 的一条：本关节行 × 整个子树列
data.M.block(idx_v, idx_v, nv_i, nvSubtree[i]).noalias()
    = J_cols.transpose() * data.Ag.middleCols(idx_v, nvSubtree[i]);

// ③ 组合惯量塌缩给父（简单相加！）
data.oYcrb[parent] += data.oYcrb[i];
```

### 6.1 ③ 组合惯量累加 —— 与 ABA 的关键区别

$$
Y_\lambda^c \mathrel{+}= Y_i^c\qquad(\texttt{oYcrb[parent] += oYcrb[i]})
$$

**直接相加**，没有 Schur 补！因为组合刚体把子树当**焊死的刚体**，惯量就是简单叠加。逆序保证处理 $i$ 时子树已全部并入 `oYcrb[i]`，故它此刻就是完整的 $Y_i^c=\sum_{j\in\text{subtree}(i)}I_j$。对比 ABA 的 `oYaba[parent] += Ia`（Ia 是塌缩后的），这里是纯累加——这正是 CRBA 比 ABA 简单的地方。

### 6.2 ② 一次填一条（行×子树列）

`data.M.block(idx_v, idx_v, nv_i, nvSubtree[i])`：从对角块开始，向右填 `nvSubtree[i]` 列——即关节 $i$ 这一（组）行里，**列属于 $i$ 的子树**的那些元素 $M_{i,j}$（$j\in\text{subtree}(i)$）。这些正是"$i$ 与它所有后代"的惯性耦合，用力集 $S_i^\top F_{\text{子树}}$ 一次算完。

---

## 7. 为什么只填上三角、`nvSubtree` 的作用

- **对称性**：$M=M^\top$，只需算上三角，另一半复制即可（省一半计算）。
- **`data.nvSubtree[i]`**：以 $i$ 为根的子树的总速度维数。因 DFS 编号，子树的关节索引是**连续区间** `[idx_v(i), idx_v(i)+nvSubtree[i])`，所以"$i$ 行里属于子树的列"恰好是一段连续块 → 用 `middleCols(idx_v, nvSubtree[i])` 一次切出，无需逐元素跳。
- 于是 CRBA 天然只碰**祖先-后代对**（$M_{ij}\ne0$ 仅当 $i,j$ 有支撑关系）。**非支撑关系的关节对 $M_{ij}=0$** ——这就是质量矩阵的**分支诱导稀疏**，也是后续 [Cholesky 分解](multibody子系统解析.md) 沿 `_fromRow` 消元树高效的根源。见 [雅可比稀疏模式解析.md](雅可比稀疏模式解析.md)。

> **mimic 关节**：破坏"子树连续区间"假设，故有专门的 `CrbaWorldConventionMimicStep` 在主循环后补齐被驱动关节对矩阵的贡献（代码 §110 起的三步注释）。普通模型不触发。

---

## 8. 与质心动量矩阵 $A_g$、ABA、Cholesky 的关系

- **$A_g$（质心动量矩阵）**：注意 Pass 2 的力集存进了 **`data.Ag`**——这不是巧合。$Y_i^c S_i$ 既是填 $M$ 的中间力集，其在质心系的表达正是**质心动量矩阵 CMM** 的列（$h_g=A_g\dot q$）。所以 CRBA 和 `ccrba` 共享同一套组合惯量与力集，见 [算法总览 §E](算法总览.md)。
- **ABA**：都用世界系组合惯量作种子（`oYcrb`），但 ABA 用 Schur 补塌缩（铰接），CRBA 用直接相加（焊死）。
- **Cholesky**：CRBA 产出的 $M$（上三角 + 分支稀疏）正好喂给 `cholesky::decompose` 做 $M=UDU^\top$，两者的稀疏结构一致（都源于支撑链），配合起来 $O(nd)$ 解 $M x=y$。
- **正动力学**：$\ddot q=M^{-1}(\tau-b)$ 里，CRBA 求 $M$、RNEA 求 $b$、Cholesky 求解——这是 ABA 之外的另一条路。

---

## 9. 复杂度

- **CRBA** $O(nd)$：$d$ 为树平均深度。链式结构 $d\approx n$ → $O(n^2)$；分支多则 $d$ 小。填 $M$ 的每条只碰子树列，总代价正比于 $\sum_i n_{v_i}\cdot\text{nvSubtree}[i]$。
- 组装 $M$ 后若要 $M^{-1}$：Cholesky $O(nd)$ 解方程，或 `computeMinverse`（ABA 路线）$O(n^2)$ 直接求逆。
- 只依赖 $q$，可在位形固定时缓存复用。

---

## 10. 一句话总结

CRBA = **组合刚体惯量两趟 $O(nd)$ 递推**：

$$
\underbrace{\text{根→叶}}_{\text{FK}+\ ^oS_i,\ ^oI_i}\quad\Rightarrow\quad
\underbrace{\text{叶→根}}_{Y_\lambda^c\mathrel{+}=Y_i^c\ (\text{直接相加}),\quad M_{ij}=S_i^\top(Y_i^cS_j)}
$$

- **组合刚体惯量** $Y^c=\sum I_j$ = 子树"焊死"成一个刚体的合惯量，**直接累加**（对比 ABA 的 Schur 补塌缩）；
- **$M_{ij}=S_i^\top Y^c S_j$**：子树惯量作用在子空间得力集 $F=Y^cS$，再投影 $S^\top F$；
- 只填**上三角**、只碰**祖先-后代对**，天然给出质量矩阵的分支诱导稀疏，直通 Cholesky 与质心动量矩阵。
