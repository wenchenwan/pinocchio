# 逆运动学（IK）解析：`examples/inverse-kinematics.cpp` 逐行数学

> 本文把 [examples/inverse-kinematics.cpp](../examples/inverse-kinematics.cpp) 的**阻尼最小二乘 IK**从数学上讲透：
> 误差怎么定义、雅可比怎么修正、为什么解出来是"速度"、$n_q\neq n_v$ 怎么办、q/v/twist 到底谁是谁。
> 配套：[空间代数运算解析.md §7.7](空间代数运算解析.md)（log6/Jlog6/Jlog3 的推导）、[multibody子系统解析.md](multibody子系统解析.md)（关节/雅可比/李群）。

---

## 目录
- [1. 问题定义](#1-问题定义)
- [2. 算法总骨架（高斯-牛顿 on SE(3)）](#2-算法总骨架高斯-牛顿-on-se3)
- [3. 逐步解析](#3-逐步解析)
  - [3.1 `iMd = oMi.actInv(oMdes)`：局部系相对位姿](#31-imd--omiactinvomdes局部系相对位姿)
  - [3.2 `err = log6(iMd)`：流形上的误差](#32-err--log6imd流形上的误差)
  - [3.3 `computeJointJacobian`：LOCAL 右扰动雅可比](#33-computejointjacobianlocal-右扰动雅可比)
  - [3.4 `Jlog6(iMd.inverse())`：为什么修正、为什么传逆](#34-jlog6imdinverse为什么修正为什么传逆)
  - [3.5 阻尼最小二乘：为什么解出来是"速度"](#35-阻尼最小二乘为什么解出来是速度)
  - [3.6 `integrate`：在流形上更新](#36-integrate在流形上更新)
- [4. 当 $n_q\neq n_v$：雅可比如何计算](#4-当-n_qneq-n_v雅可比如何计算)
- [5. q / v / 速度旋量：三个概念务必分清](#5-q--v--速度旋量三个概念务必分清)
- [6. 一句话总结](#6-一句话总结)

---

## 1. 问题定义

求关节配置 $q$，使第 6 个关节的位姿等于目标：

$$M_6(q)={}^0M_6(q)=M_{\text{des}}\quad(\texttt{oMi[6]}=\texttt{oMdes})$$

$M_6(q)$ 是 $q$ 的高度非线性函数（一串 SE(3) 复合，见 [multibody 解析 §13](multibody子系统解析.md)），无闭式解，故迭代。

---

## 2. 算法总骨架（高斯-牛顿 on SE(3)）

$$
\begin{aligned}
&\textbf{repeat:}\\
&\quad M_6(q) \leftarrow \text{forwardKinematics}(q) &&\text{(SE3 递推)}\\
&\quad e \leftarrow \log_6\!\big(M_6^{-1}M_{\text{des}}\big) &&\text{(流形误差, }\in\mathbb{R}^6)\\
&\quad \textbf{if } \|e\|<\varepsilon:\ \textbf{break}\\
&\quad \tilde J \leftarrow -\,\mathrm{Jlog6}(\text{iMd}^{-1})\,J &&\text{(误差对 }q\text{ 的雅可比)}\\
&\quad v \leftarrow -\tilde J^\top(\tilde J\tilde J^\top+\lambda^2 I)^{-1}e &&\text{(阻尼最小二乘步)}\\
&\quad q \leftarrow q \oplus v\,\Delta t &&\text{(流形积分)}
\end{aligned}
$$

本质：把非线性方程 $e(q)=0$ 反复**线性化 → 解线性最小二乘 → 在流形上走一步**。$\lambda$ 让它成为 Levenberg–Marquardt（稳健化）。三处都"活在正确的几何空间"：误差用 $\log_6$、雅可比乘 $\mathrm{Jlog6}$、更新用 `integrate`。

---

## 3. 逐步解析

### 3.1 `iMd = oMi.actInv(oMdes)`：局部系相对位姿

```cpp
const pinocchio::SE3 iMd = data.oMi[JOINT_ID].actInv(oMdes);
```

`A.actInv(B) = A⁻¹·B`（不显式求逆，用 $R^\top,-R^\top p$ 拼，见 [空间代数解析 §3.2/§3.3](空间代数运算解析.md)）。于是

$$\texttt{iMd}={}^0M_6^{-1}\cdot{}^0M_{\text{des}}={}^6M_{\text{des}}$$

名字 `iMd`=${}^iM_d$ = **"站在关节 6 自身坐标系看，目标还差多少"**。到位时 $\texttt{iMd}=I$。选局部系是为了和后面的 `log6`、LOCAL 雅可比**坐标系对齐**。

### 3.2 `err = log6(iMd)`：流形上的误差

```cpp
err = pinocchio::log6(iMd).toVector();   // ∈ R⁶
```

位姿住在 SE(3) 流形，两个位姿不能相减。正确的"差"是相对变换取对数：$e=\log_6(\texttt{iMd})\in\mathfrak{se}(3)\cong\mathbb{R}^6$（前 3 线误差、后 3 角误差）。**收敛判据** $\|e\|<\varepsilon \iff \texttt{iMd}\approx I$，因为 $\log_6(I)=0$。数学见 [空间代数解析 §7.3/§7.7](空间代数运算解析.md)。

### 3.3 `computeJointJacobian`：LOCAL 右扰动雅可比

```cpp
pinocchio::computeJointJacobian(model, data, q, JOINT_ID, J);  // J ∈ R^{6×nv}
```

这个 5 参重载返回 **LOCAL（关节自身系）雅可比**，对应**右扰动**：

$$M_6(q+\delta q)\approx M_6\exp(J\delta q)\ \Longleftrightarrow\ {}^6v_6=J\dot q$$

即 $J\dot q$ 是 $M_6$ 的 body twist（右平移化速度）。**列数是 $n_v$（切空间维度），不是 $n_q$**（见 §4）。要世界系版本用 `getJointJacobian(..., WORLD/LOCAL_WORLD_ALIGNED, ...)`。

### 3.4 `Jlog6(iMd.inverse())`：为什么修正、为什么传逆

```cpp
pinocchio::Jlog6(iMd.inverse(), Jlog);
J = -Jlog * J;                          // J̃ = ∂e/∂q
```

**为什么要 Jlog6**：误差 $e=\log_6(\texttt{iMd})$ 是 iMd 的**非线性**函数。关节动 → （$J$）iMd 动 → （$\mathrm{Jlog6}$）log 读数动。链式法则：

$$\frac{\partial e}{\partial q}=\frac{\partial\log_6}{\partial\,\texttt{iMd}}\cdot\frac{\partial\,\texttt{iMd}}{\partial q}=\mathrm{Jlog6}\cdot J$$

不乘 Jlog6，靠近目标处偏差不大（$\mathrm{Jlog6}\to I_6$），但**远离目标时非线性显著、步子会偏**。直觉：类似 1D 牛顿法 $x\leftarrow x-f(x)/f'(x)$ 里那个 $f'$——Jlog6 是"误差读数对物理运动的敏感度/换算系数"。

**为什么传逆 + 负号**：`Jlog6` 是**右雅可比**（右扰动约定），代码可验证——其旋转块 $\mathrm{Jlog3}=J_r^{-1}=I+\tfrac12[\omega]_\times+\dots$，`addSkew(0.5*w,…)` 的 $+\tfrac12$ 即右。但关节右扰动经 $M_6^{-1}$ 反转后落到 iMd **左侧**（$\exp(-J\delta q)\cdot\texttt{iMd}$），对不上。改传 $\texttt{iMd}^{-1}=M_{\text{des}}^{-1}M_6$ 让扰动回到右边，再用 $\log(M^{-1})=-\log(M)$ 补负号：

$$\frac{\partial e}{\partial q}=-\,\mathrm{Jlog6}(\texttt{iMd}^{-1})\,J$$

冲突**不在**"左右雅可比"本身（关节 Jacobian 与 Jlog6 都是"右"），而在 $M_6^{-1}$ 把右扰动搬到了 iMd 左侧。**完整推导见 [空间代数解析 §7.7.4/§7.7.5](空间代数运算解析.md)**。记 $\tilde J=-\mathrm{Jlog6}(\texttt{iMd}^{-1})J=\partial e/\partial q$。

### 3.5 阻尼最小二乘：为什么解出来是"速度"

```cpp
JJt = J*Jᵀ;  JJt.diagonal() += damp;                  // (J̃J̃ᵀ + λ²I)，6×6
v = -Jᵀ · JJt.ldlt().solve(err);                       // v ∈ R^{nv}
```

**这是什么**：线性化 $e+\tilde J\,\delta q\approx0$ 后求 $\tilde J\,\delta q=-e$ 的**阻尼最小二乘/最小范数解**：

$$v=-\tilde J^\top(\tilde J\tilde J^\top+\lambda^2 I)^{-1}e\ =\ \arg\min_{v}\ \|\tilde J v+e\|^2+\lambda^2\|v\|^2$$

**为什么落在关节速度空间**——跟着维度走：

| 部件 | 维度 | 空间搬运 |
|------|------|----------|
| `err` $e$ | $\mathbb{R}^6$ | 任务空间（6D 位姿误差） |
| $(\tilde J\tilde J^\top+\lambda^2 I)^{-1}$ | $6\times6$ | 任务 → 任务 |
| $\tilde J^\top$ | $n_v\times6$ | **任务 → 关节速度空间** |
| **结果 `v`** | $\mathbb{R}^{n_v}$ | **关节速度空间** |

最后乘的 $\tilde J^\top$（$n_v\times6$，雅可比转置=对偶反传）把 6D 任务误差**回拉**到 $n_v$ 维关节空间，所以结果天生是关节速度向量。

**为什么叫"速度"而非"位移"**：下一行 `q = integrate(model, q, v*DT)` 把它当"速度 × 时间"用。整个 IK 等价于积分一条**连续牛顿流 / 梯度流**：

$$\dot q=-\tilde J^\top(\tilde J\tilde J^\top+\lambda^2 I)^{-1}e(q)$$

沿这个速度场走会让 $e(q)$ 单调降到 0；每次循环 = 对它做**一步欧拉积分**，步长 `DT=0.1`。在这个人造动力系统里 `v` 就是字面的速度。`DT=1` 且无阻尼时退化成整步高斯-牛顿。

**阻尼 $\lambda^2 I$ 的作用**：奇异构型下 $\tilde J\tilde J^\top$ 病态，纯伪逆会爆；加 $\lambda^2 I$（`damp=1e-6`）保证正定可逆、步长有界；`ldlt().solve` 用 $LDL^\top$ 分解稳健求解。用 $\tilde J\tilde J^\top$（$6\times6$）而非 $\tilde J^\top\tilde J$（$n_v\times n_v$）更省，且给"最小范数解"。

### 3.6 `integrate`：在流形上更新

```cpp
q = pinocchio::integrate(model, q, v * DT);   // q ← q ⊕ v·DT
```

$v$ 是切空间量，更新**不能** $q+v\Delta t$（会把四元数/浮动基踢出流形）。`integrate` 沿每个关节各自流形走测地线：欧氏关节就是 $q+v\Delta t$，球副/浮动基走 $q\cdot\exp(v\Delta t)$（四元数版 exp3/exp6），自动维持归一化。机制见 [multibody 解析 §10](multibody子系统解析.md)。

---

## 4. 当 $n_q\neq n_v$：雅可比如何计算

**结论：关节雅可比永远是 $6\times n_v$，列数是 $n_v$（切空间维度），从不是 $n_q$。**

雅可比定义 $v_{\text{spatial}}=J(q)\,v$ 里的 $v$ 是切空间速度（$n_v$ 维），$n_q$ 维的配置向量（如四元数 4 个数）**不作为雅可比的列**。$n_q>n_v$ 只是**存储表示的冗余**（4 个数存 3 自由度姿态），真实自由度/雅可比/动力学全是 $n_v$ 维。

| 关节 | $n_q$ | $n_v$ | 雅可比列块 | 运动子空间 S |
|------|------|------|-----------|-------------|
| 球副 Spherical | 4（四元数） | 3（$\omega$） | $6\times3$ | $[0_3;I_3]$ |
| 浮动基 FreeFlyer | 7（平移+四元数） | 6 | $6\times6$ | $I_6$ |
| 无界转动 RUBX | 2（$\cos,\sin$） | 1 | $6\times1$ | 单位角速度列 |

**代码怎么处理不等**（靠 [multibody 解析 §2.2](multibody子系统解析.md) 的三套索引各司其职）：

```cpp
q_joint = jmodel.jointConfigSelector(q);   // 用 (idx_q, nq) 读配置：球副取 4 个四元数
calc(jdata, q_joint);                       // 用 q 算关节变换 M 和 S(6×nv)
jmodel.jointCols(J) = oMi.act(jdata.S);     // 用 (idx_v, nv) 摆雅可比列：球副填 3 列
```

`jointConfigSelector` 用 `idx_q/nq` 读配置；`jointCols/jointVelocitySelector` 用 `idx_v/nv` 摆雅可比/速度。**天然分开。**

**$\dot q$（配置向量真实导数，$n_q$ 维）去哪了**：它 $\neq v$、不进雅可比，二者由配置相关的**切映射（tangent map）** $T(q)\in\mathbb{R}^{n_q\times n_v}$ 联系：$\dot q_{\text{config}}=T(q)\,v$。四元数版：

$$\dot{\mathbf q}=\tfrac12\,\mathbf q\otimes\begin{bmatrix}\omega\\0\end{bmatrix}=\tfrac12\,\Omega(\mathbf q)\,\omega,\quad \Omega(\mathbf q)\in\mathbb{R}^{4\times3}$$

（Pinocchio 里 `integrateCoeffWiseJacobian`/`tangentMap`）。它保证四元数速度切于单位球面。正因 $\dot q\neq v$，才必须 `integrate` 而非相加。**所以含球副/浮动基的机器人，这段 IK 代码一字不用改**：$J$ 天然 $6\times n_v$、$v$ 天然 $n_v$ 维、`integrate` 天然处理 $n_q\neq n_v$。

---

## 5. q / v / 速度旋量：三个概念务必分清

以典型 6 轴纯转动工业臂（如 UR5，6 个 revolute，$n_q=n_v=6$）为例，三样东西都"6 维"却完全不同：

$$
\underbrace{q}_{\text{关节角}[\text{rad}]}\ \xrightarrow{\ \frac{d}{dt}\ }\ \underbrace{v=\dot q}_{\text{关节角速度}[\text{rad/s}]}\ \xrightarrow{\ \times J(q)\ }\ \underbrace{v_{\text{spatial}}=Jv}_{\textbf{速度旋量 twist}\in\mathfrak{se}(3)}
$$

| 符号 | 是什么 | 维度 | 空间 |
|------|--------|------|------|
| `q` | 关节角度 $[\theta_1,\dots,\theta_6]$ | $\mathbb{R}^6$ | 配置空间 |
| `v` | 关节角速度 $\dot q$ | $\mathbb{R}^6$ | **关节速度空间（切空间）** |
| 速度旋量 twist | 末端 6D 空间速度 $(v_{\text{lin}};\omega)$ | $\mathbb{R}^6$ | $\mathfrak{se}(3)$，**= $J\cdot v$ 的输出** |

- ✅ **`q` 是关节配置**（6 个角度）。
- ❌ **`v` 不是速度旋量**——它是**关节速度** $\dot q$（6 个关节角速度），住在关节速度空间。速度旋量（twist）是 $J\cdot v$ 的**输出**，不是 `v` 本身。

**易混点：输入 `v` vs `data.v[i]`**。`forwardKinematics(model,data,q,v)` 里输入的 `v` 是关节速度 $\dot q$；算出的 `data.v[i]` 才是 twist（第 $i$ 关节的空间速度旋量，$\in\mathfrak{se}(3)$），由递推 `data.v[i]=liMi.actInv(data.v[parent]) + S_i·v_i` 累积（见 [multibody 解析 §13](multibody子系统解析.md)）。

纯转动臂因所有关节欧氏，$n_q=n_v=6$、`integrate` 退化成加法；但概念上 `v` 仍是"关节速度"而非"旋量"。

---

## 6. 一句话总结

> IK = 在 SE(3) 上跑**高斯-牛顿**：误差用 $\log_6$ 量在流形上（$e=\log_6(M_6^{-1}M_{\text{des}})$）、雅可比用 $\mathrm{Jlog6}$ 修正到 log 空间（传 `iMd.inverse()`+负号把右扰动掰正）、解**阻尼最小二乘**得到一个住在 $n_v$ 关节速度空间的步长 $v$（因末乘 $\tilde J^\top$ 回拉，又被 `integrate(q,v·DT)` 当速度×时间用，本质是连续牛顿流的一步欧拉）、最后在**流形上积分**更新 $q$。雅可比、速度、力矩永远是 $n_v$ 维，$n_q$ 只用于存配置；二者由切映射 $T(q)$ 相连，体现在 `integrate` 而非雅可比里。
