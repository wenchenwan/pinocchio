# 质量矩阵 CRBA 解析：`algorithm/crba.hxx`

> **CRBA**（Composite Rigid Body Algorithm，组合刚体算法，Featherstone 第 6 章）计算**关节空间惯性矩阵** $M(q)$（质量矩阵，Featherstone 记作 $H$）。
> 本文**逐行贴 [src/algorithm/crba.hxx](../include/pinocchio/src/algorithm/crba.hxx)** 讲透其两种实现约定：
> 组合刚体惯量的意义、$M_{ij}=S_i^\top Y^c S_j$ 的推导、**LOCAL（经典 Featherstone）vs WORLD（Pinocchio 优化）** 两条实现路线、
> 上三角与 `nvSubtree`、紧凑树要求、armature、mimic 补齐。
> 配套：[逆动力学RNEA解析.md](逆动力学RNEA解析.md)、[正动力学ABA解析.md](正动力学ABA解析.md)、
> [质心与质心动量解析.md](质心与质心动量解析.md)、[multibody 子系统解析 §3.5（M 的稀疏 Cholesky）](multibody子系统解析.md)、[空间代数运算解析.md](空间代数运算解析.md)。

---

## 目录

- [1. 质量矩阵是什么](#1-质量矩阵是什么)
- [2. 核心概念：组合刚体惯量](#2-核心概念组合刚体惯量composite-rigid-body-inertia)
- [3. 关键公式：$M_{ij}=S_i^\top Y^c S_j$](#3-关键公式m_ijs_itop-yc-s_j)
- [4. 顶层入口与两种约定](#4-顶层入口与两种约定)
- [5. 涉及的 `Data` 数据结构](#5-涉及的-data-数据结构)
- [6. LOCAL 约定实现（经典 Featherstone 体坐标版）](#6-local-约定实现经典-featherstone-体坐标版)
- [7. WORLD 约定实现（Pinocchio 优化版）](#7-world-约定实现pinocchio-优化版)
- [8. 两种约定对比：差异到底在哪](#8-两种约定对比差异到底在哪)
- [9. 上三角、`nvSubtree` 与紧凑树要求（`CRBAChecker`）](#9-上三角nvsubtree-与紧凑树要求crbachecker)
- [10. `armature`：转子惯量](#10-armature转子惯量)
- [11. mimic 关节：`MimicStep` 补齐](#11-mimic-关节mimicstep-补齐)
- [12. 与 $A_g$、ABA、Cholesky 的关系](#12-与-a_gabacholesky-的关系)
- [13. 复杂度](#13-复杂度)
- [14. 一句话总结](#14-一句话总结)

---

## 1. 质量矩阵是什么

$M(q)$ 是拉格朗日动力学 $M\ddot q+b=\tau$ 里的**关节空间惯性矩阵**（$n_v\times n_v$，对称正定）：

- 动能 $T=\tfrac12\dot q^\top M(q)\dot q$；
- $M_{ij}$ = "关节 $j$ 加速时，关节 $i$ 感受到的惯性耦合力矩"；
- 正动力学要 $M^{-1}$，最优控制/碰撞/操作空间控制处处要它。

CRBA 用 $O(nd)$ 直接组装出 $M$ 的上三角（$d$=树深度），是求 $M$ 的标准方法。

---

## 2. 核心概念：组合刚体惯量（Composite Rigid Body Inertia）

**组合刚体惯量** $Y_i^c$ = 以连杆 $i$ 为根的**整个子树，当作一个刚体焊死**后的合惯量：

$$
\boxed{\,Y_i^c=\sum_{j\in\text{subtree}(i)} I_j\,}\qquad\text{（各 } I_j\text{ 需变换到同一坐标系再相加）}
$$

递推形式（实现用这个）：$Y_i^c = I_i + \sum_{c\in\text{children}(i)} Y_c^c$。依据是空间惯量的**可加性**：同一坐标系下刚性连接的刚体，惯量直接相加。

与 ABA 的**铰接体惯量** $I^A$ 对照：

|              | 组合刚体惯量 $Y^c$（CRBA） | 铰接体惯量 $I^A$（ABA） |
| ------------ | ------------------------- | ------------------------------- |
| 子树关节视为 | **焊死（刚性）**          | **可自由转动（铰接）**          |
| 累加方式     | 简单相加 $\sum I_j$        | Schur 补塌缩 $I^A-UD^{-1}U^\top$ |
| 用途         | 组装 $M$                   | 隐式求 $M^{-1}$                  |

CRBA 的"焊死"假设正是它简单的原因：$M_{ij}$ 描述的是**惯性耦合**，等价于把两关节间的连杆看成刚性时的惯量投影。

---

## 3. 关键公式：$M_{ij}=S_i^\top Y^c S_j$

设关节 $i$ 是关节 $j$ 的祖先（$i$ 在 $j$ 的支撑链上）。质量矩阵元（Featherstone Eq. 6.16）：

$$
\boxed{\,M_{ij}=S_i^\top\,Y_{\max(i,j)}^c\,S_j\,}
$$

**推导直觉**：$M_{ij}$ = 关节 $j$ 单位加速度在关节 $i$ 轴上引起的力矩。关节 $j$ 加速 → 带动**子树 $j$**（组合惯量 $Y_j^c$）产生空间力 $Y_j^c S_j$ → 这个力沿支撑链传到祖先关节 $i$，在 $S_i$ 上投影即 $S_i^\top(Y_j^c S_j)$。因为 $j$ 是后代、子树更小，交集取小者，故用 $Y_{\max(i,j)}^c=Y_j^c$。

> 完整的动能推导（从 $T=\tfrac12\sum_k v_k^\top I_k v_k$ 逼出此式）见 Featherstone §6.2，或本仓库配套的 CRBA 章节笔记。这里聚焦**代码实现**。

实现拆成两步，两种约定都遵循这个骨架：

$$
\underbrace{F=Y^c S}_{\text{①惯量作用→力集}}\qquad
\underbrace{M[\ldots]=S^\top F}_{\text{②力集投影→矩阵块}}
$$

---

## 4. 顶层入口与两种约定

`crba(...)` 是一个**分派器**，按 `Convention` 选两种数学等价、但数据布局不同的实现（[crba.hxx:608](../include/pinocchio/src/algorithm/crba.hxx#L608)）：

```cpp
const MatrixXs & crba(model, data, q, const Convention convention) {
  switch (convention) {
  case Convention::LOCAL: return impl::crbaLocalConvention(model, data, q);  // 经典 Featherstone
  case Convention::WORLD: return impl::crbaWorldConvention(model, data, q);  // Pinocchio 优化
  default: throw std::invalid_argument("Bad convention.");
  }
}
```

| 约定 | 惯量/力集表达在 | 数据 | 特点 |
|---|---|---|---|
| **LOCAL** | 各连杆**自身体坐标系** | `Ycrb[i]`, `Fcrb[i]` | 教科书 Featherstone；累加时每步做 `liMi.act` 变换 |
| **WORLD** | **世界系** | `oYcrb[i]`, `Ag` | 前向先全变到世界系；后向累加**免变换**；顺带产出质心动量矩阵 |

两者算出**完全相同**的 $M$。下面分别精讲。

---

## 5. 涉及的 `Data` 数据结构

| 字段 | 约定 | 类型 | 含义 |
|---|---|---|---|
| `data.M` | 共用 | $n_v\times n_v$ | 输出：质量矩阵（只填上三角 + 对角） |
| `data.liMi[i]` | 共用 | SE3 | 父→子相对位姿（Pass1 填） |
| `data.Ycrb[i]` | LOCAL | Inertia | 组合刚体惯量，**在 body $i$ 系** |
| `data.Fcrb[i]` | LOCAL | $6\times n_v$ | 力集矩阵：子树各列 = $Y^cS$ 变换到 body $i$ 系 |
| `data.oYcrb[i]` | WORLD | Inertia | 组合刚体惯量，**在世界系** |
| `data.J` | WORLD | $6\times n_v$ | 世界系运动子空间 ${}^oX_iS_i$（也是整机雅可比）|
| `data.Ag` | WORLD | $6\times n_v$ | 力集 $Y^cS$，即**质心动量矩阵**列（见 §12）|

---

## 6. LOCAL 约定实现（经典 Featherstone 体坐标版）

这是与 Featherstone 教材 Table 6.2、式 6.17/6.19/6.20 **逐行对应**的实现，最能看清算法本质。

### 6.0 驱动函数 `crbaLocalConvention` 骨架

顶层函数（[crba.hxx:460](../include/pinocchio/src/algorithm/crba.hxx#L460)）本身很短——三趟循环 + 一行 armature，每趟背后是下面 §6.1–6.4 的数学：

```cpp
const MatrixXs & crbaLocalConvention(model, data, q) {
  assert(model.check(data));  PINOCCHIO_CHECK_ARGUMENT_SIZE(q.size(), model.nq);

  for (i=1; i<njoints; ++i)          Pass1::run(...);   // §6.1 前向：FK + 播种 Ycrb
  for (i=njoints-1; i>0; --i)        Pass2::run(...);   // §6.2 后向：Featherstone 主算法
  for (i=0; i<mimicking_joints; ++i) Pass3::run(...);   // §11 mimic 补齐

  data.M.diagonal() += model.armature;                  // §10 转子惯量
  return data.M;                                        // 只填上三角+对角
}
```

- 返回的 `data.M` **只填上三角 + 对角**，下三角留空（用时镜像）。
- 三趟都用编译期访问者按关节类型零开销分派（见 [正向运动学 §4](正向运动学解析.md)）。

### 6.1 Pass 1 前向：只做 FK + 播种本体惯量

```cpp
// CrbaLocalConventionForwardStep
jmodel.calc(jdata, q);
data.liMi[i] = model.jointPlacements[i] * jdata.M();   // 父→子相对位姿
data.Ycrb[i] = model.inertias[i];                      // ★ 种子 = 本体惯量（本系，不变换）
```

注意 **LOCAL 约定的 `Ycrb[i]` 种子无需任何坐标变换**——因为惯量本就存在 body $i$ 自己的坐标系里。这是与 WORLD 约定（§7.1 要 `oMi.act`）的第一个差别。

### 6.2 Pass 2 后向：Featherstone 四行

源码顶部的注释块（[crba.hxx:282](../include/pinocchio/src/algorithm/crba.hxx#L282)）就是 Featherstone 算法本身：

```
F[1:6,i]      = Y*S                    ← ①本关节的力集列
M[i,SUBTREE]  = S'*F[1:6,SUBTREE]      ← ②填 M 的这一行块
if li>0:
  Yli += liXi Yi                       ← ③组合惯量塌缩给父
  F[1:6,SUBTREE] = liXi F[1:6,SUBTREE] ← ④力集变换到父系，供父使用
```

对应代码及**每步的数学公式**（右侧标注 Featherstone 式号）：

```cpp
// ① F[:,i] = Y^c_i · S_i  —— 组合惯量作用于运动子空间，得空间力
jmodel.jointCols(data.Fcrb[i]) = data.Ycrb[i] * jdata.S();

// ② M[i, SUBTREE] = S_i^T · F[:,SUBTREE]  —— 一次填这一行块
data.M.block(idx_v, idx_v, nv, nvSubtree[i]).noalias()
    = jdata.S().transpose() * data.Fcrb[i].middleCols(idx_v, nvSubtree[i]);

if (parent > 0) {
  // ③ Y^c_parent += liXi · Y^c_i  —— 组合惯量变换到父系再相加
  data.Ycrb[parent] += data.liMi[i].act(data.Ycrb[i]);

  // ④ F_parent[SUBTREE] = liXi · F_i[SUBTREE]  —— 力集整体搬到父系
  Block jF = data.Fcrb[parent].middleCols(idx_v, nvSubtree[i]);
  Block iF = data.Fcrb[i].middleCols(idx_v, nvSubtree[i]);
  forceSet::se3Action(data.liMi[i], iF, jF);
}
```

- **①** ${}^iF_i=Y_i^c\,S_i$（Featherstone 式 6.19 种子）——"让整棵子树 $\nu(i)$ 沿 $S_i$ 方向获单位加速度所需的空间力"。$Y_i^c$ 此刻已完整（子树在更高编号，已由 ③ 累加进来）。
- **②** $M_{ij}=S_i^\top\,{}^iF_j=S_i^\top\,{}^iX_j^{*}\,Y_j^c\,S_j\ (j\in\text{subtree}(i))$（式 6.18/6.20 的 $j\in\nu(i)$ 情形，$i$ 祖先、$j$ 后代）。因 DFS 编号 $\text{idx\_v}(i)\le\text{idx\_v}(j)$，落在**上三角**；列 $=i$ 即对角 $M_{ii}=S_i^\top Y_i^cS_i$。
- **③** $Y_\lambda^c\mathrel{+}={}^{\lambda}X_i^{*}\,Y_i^c\,{}^iX_\lambda$（式 6.17 一项）。`SE3::act` 作用在 Inertia 上做**相合变换** $I\mapsto{}^\lambda X_i^{*}I\,{}^iX_\lambda$（子系惯量→父系）。**纯累加带一次变换，无 Schur 补**。
- **④** ${}^{\lambda}F_k={}^{\lambda}X_i^{*}\,{}^iF_k\ (k\in\text{subtree}(i))$（式 6.19 递推）。把子树力集批量搬到父系存入 `Fcrb[parent]`，供父的 ② 使用。

### 6.3 关键：`Fcrb[i]` 累积的是什么

处理关节 $i$ 时，`data.Fcrb[i]` 的**子树各列**已经被 $i$ 的后代在之前的迭代里（通过第 ④ 步 `forceSet::se3Action`）填好——每个后代 $k$ 的力集 $Y_k^cS_k$ 已被逐级变换到 **body $i$ 的坐标系**，即 Featherstone 的 ${}^iF_k={}^iX_k\,Y_k^cS_k$。于是第 ② 步：

$$
M_{i,k}=S_i^\top\,{}^iF_k=S_i^\top\,{}^iX_k\,Y_k^c\,S_k\qquad(k\in\text{subtree}(i))
$$

一次矩阵乘 $S_i^\top\cdot\text{Fcrb}[i]_{[\text{subtree}]}$ 就把关节 $i$ 这一（组）行、所有子树列一次算完。

### 6.4 力沿树上传（第 ④ 步的意义）

`forceSet::se3Action(liMi[i], iF, jF)` 把 $i$ 的子树力集从 body $i$ 系变换到父系，写入 `Fcrb[parent]`。这样父被处理时，`Fcrb[parent]` 已含 $i$ 整棵子树的力（在父系表达）。**这正是 Featherstone "力沿链向根传递、每级做一次 $X^*$"** 的实现——组合惯量 `Ycrb`（③）与力集 `Fcrb`（④）**双双向上传播**。

> `forceSet::se3Action` 是[空间力集的批量对偶变换](multibody子系统解析.md)（wrench 用 `act`，见 [§12.3 ForceSet](multibody子系统解析.md) 与 [RNEA §5.2](逆动力学RNEA解析.md) 力/运动对偶）。

### 6.5 代码 ↔ Featherstone 公式对照

`crbaLocalConvention` 就是 Featherstone 体坐标 CRBA（Table 6.2）的直译，逐行对上号：

| 代码 | Featherstone 式 | 数学 |
|---|---|---|
| `Ycrb[i] = inertias[i]` | 6.12/6.13 初值 | $Y_i^c\leftarrow I_i$ |
| `Fcrb[i][:,i] = Ycrb[i]*S` | 6.19 种子 | ${}^iF_i=Y_i^cS_i$ |
| `M[i,SUB] = Sᵀ·Fcrb[i][SUB]` | 6.18 / 6.20 | $M_{ij}=S_i^\top\,{}^iX_j^{*}Y_j^cS_j$ |
| `Ycrb[λ] += liMi.act(Ycrb[i])` | 6.17 | $Y_\lambda^c\mathrel{+}={}^\lambda X_i^{*}Y_i^c\,{}^iX_\lambda$ |
| `forceSet::se3Action(liMi,iF,jF)` | 6.19 递推 | ${}^\lambda F_k={}^\lambda X_i^{*}\,{}^iF_k$ |

要点串起来：前向播种本体惯量，后向沿树把**组合惯量 $Y^c$ 与力集 $F$ 同时向根传播**（③④），每到一个关节用 $M_{ij}=S_i^\top\,{}^iF_j$ 一次填出它整行的子树块（②），只碰祖先-后代对（上三角）；末尾补 armature、mimic。$O(1)$ 摊到每条树边，总 $O(nd)$。

> **对照 WORLD 约定**（§7）：LOCAL 的③④每步都要 `liMi.act`/`se3Action` 变换（因量在各自体系）；WORLD 把变换提前到 Pass 1，后向累加 `oYcrb[父]+=oYcrb[i]` 免变换。两者数学等价，见 §8。

---

## 7. WORLD 约定实现（Pinocchio 优化版）

WORLD 约定把所有量在 Pass 1 一次性搬到世界系，换来 Pass 2 累加**免坐标变换**，并顺带产出质心动量矩阵。

### 7.1 Pass 1 前向：FK + 全变到世界系

```cpp
// CrbaWorldConventionForwardStep
data.liMi[i] = model.jointPlacements[i] * jdata.M();
data.oMi[i]  = (parent>0) ? data.oMi[parent]*data.liMi[i] : data.liMi[i];

jmodel.jointExtendedModelCols(data.J) = data.oMi[i].act(jdata.S());  // ★ 世界系 S_i → data.J
data.oYcrb[i] = data.oMi[i].act(model.inertias[i]);                  // ★ 世界系本体惯量（种子）
```

- `data.J` 存**世界系运动子空间** ${}^oX_iS_i$（这也是整机雅可比）。
- `oYcrb[i]` = 本体惯量 `oMi.act` 到世界系，作组合惯量种子。
- 代价：比 LOCAL 多了每关节一次 `oMi.act`（惯量 + 子空间）。

### 7.2 Pass 2 后向：累加免变换

```cpp
// CrbaWorldConventionBackwardStep
Ag_cols = jmodel.jointCols(data.Ag);
J_cols  = jmodel.jointExtendedModelCols(data.J);          // 世界系 S_i

// ① 力集 F_i = Y^c_i · S_i（世界系）
motionSet::inertiaAction(data.oYcrb[i], J_cols, Ag_cols);  // Ag 列 = oYcrb[i]·J_cols

// ② 填 M 的一条：本关节行 × 子树列
data.M.block(idx_v, idx_v, nv, nvSubtree[i]).noalias()
    = J_cols.transpose() * data.Ag.middleCols(idx_v, nvSubtree[i]);

// ③ 组合惯量塌缩给父 —— 直接相加，无 liMi.act！
data.oYcrb[parent] += data.oYcrb[i];
```

**核心差别在第 ③ 步**：

$$
Y_\lambda^c \mathrel{+}= Y_i^c\qquad(\texttt{oYcrb[parent] += oYcrb[i]})
$$

**没有 `liMi.act` 变换！** 因为 `oYcrb[父]` 和 `oYcrb[子]` 都已在**同一个世界系**里表达，直接相加即可。对比 LOCAL 的第 ③ 步 `Ycrb[parent] += liMi[i].act(Ycrb[i])` 需要一次惯量变换——这就是 WORLD 约定"前向多花、后向省回"的本质。力集也不必像 LOCAL 那样逐级上传（`data.Ag` 的列直接就是世界系的，第 ② 步用 `data.Ag.middleCols` 取子树列即可）。

### 7.3 收尾：产出质心动量矩阵

WORLD 约定跑完后，`crbaWorldConvention` **顺手**从 `oYcrb[0]`（整机组合惯量）抽出质心，并把 `data.Ag` 平移到质心系（[crba.hxx:545](../include/pinocchio/src/algorithm/crba.hxx#L545)）：

```cpp
data.mass[0] = data.oYcrb[0].mass();
data.com[0]  = data.oYcrb[0].lever();                     // 整机质心
for (i=0; i<nv; ++i)
  Ag_ang.col(i) += Ag_lin.col(i).cross(data.com[0]);      // Ag 平移到质心 → CMM
```

于是 CRBA（WORLD）**一次遍历同时得到 $M$ 和质心动量矩阵 $A_g$**——LOCAL 约定不产出 $A_g$。详见 §12 与 [质心与质心动量解析.md](质心与质心动量解析.md)。

---

## 8. 两种约定对比：差异到底在哪

| 步骤 | LOCAL（Featherstone） | WORLD（Pinocchio） |
|---|---|---|
| Pass 1 惯量种子 | `Ycrb[i]=inertias[i]`（本系，免变换）| `oYcrb[i]=oMi.act(inertias[i])`（多一次变换）|
| 子空间 | `jdata.S()`（本系） | `oMi.act(jdata.S())`→`data.J`（世界系）|
| Pass 2 力集 | `Fcrb[i]=Ycrb[i]*S`，需逐级上传（④）| `Ag=oYcrb[i]*J`，列即世界系，不上传 |
| Pass 2 惯量累加 | `Ycrb[父]+=liMi.act(Ycrb[i])`（**有变换**）| `oYcrb[父]+=oYcrb[i]`（**无变换**）|
| 副产物 | 无 | **质心动量矩阵 $A_g$**、质心 `com[0]` |
| 出处 | 教科书标准 | 与雅可比/CoM/ccrba 共享世界系中间量 |

**取舍**：WORLD 把"每次累加的坐标变换"提前摊到 Pass 1 的一次性变换里，且复用了世界系 `data.J`（雅可比也要）、附赠 $A_g$，在 Pinocchio 的整体算法生态里更划算——所以很多上层默认走 WORLD。LOCAL 更贴教材、在只要 $M$ 且不需要世界系量时略省。

---

## 9. 上三角、`nvSubtree` 与紧凑树要求（`CRBAChecker`）

### 9.1 只填上三角

$M=M^\top$，两种约定都只算**上三角 + 对角**（`data.M.block(idx_v, idx_v, nv, nvSubtree[i])` 从对角块向右填），另一半用时再镜像。省一半计算。

### 9.2 `nvSubtree[i]`：子树是连续列区间

`data.M.block(idx_v, idx_v, nv, nvSubtree[i])`：从对角块开始，向右填 `nvSubtree[i]` 列——即关节 $i$ 行里"列属于 $i$ 子树"的元素 $M_{i,j}$（$j\in\text{subtree}(i)$）。**这一步能用连续 `middleCols` 切出，全靠子树在速度索引里是一段连续区间** `[idx_v(i), idx_v(i)+nvSubtree[i])`。

### 9.3 紧凑树要求：`CRBAChecker`

上面的"子树 = 连续区间"不是白来的，它要求树是**紧凑编号**的。`CRBAChecker::checkModel_impl`（[crba.hxx:579](../include/pinocchio/src/algorithm/crba.hxx#L579)）就检查这个不变量：

```cpp
// 对每个 i，区间 i+1..n-1 能切成 [i+1..k]（全是 i 的后代）和 [k+1..n-1]（全不是）
for (i=1; i<njoints-1; ++i) {
  k = i+1;
  while (isDescendant(model, k, i)) ++k;      // 后代必须紧邻 i 之后连续排布
  for (; k<njoints; ++k)
    if (isDescendant(model, k, i)) return false;  // 后代之后不能再冒出后代
}
```

即"**一个节点的所有后代必须在 `parents` 数组里紧跟其后连续存放**"。DFS 建树天然满足（见 [multibody §2.4](multibody子系统解析.md)）。满足它，`nvSubtree[i]` 才能当连续块用；否则 CRBA 的块操作会取错列。

### 9.4 分支诱导稀疏

因为只碰**祖先-后代对**，**非支撑关系的关节对 $M_{ij}=0$**（结构性精确为零，如人形左右臂之间）。这就是质量矩阵的**分支诱导稀疏**，也是后续 [Cholesky 分解](multibody子系统解析.md) 沿 `_fromRow` 消元树高效的根源。见 [雅可比稀疏模式解析.md](雅可比稀疏模式解析.md)。

---

## 10. `armature`：转子惯量

两种约定收尾都加一行（[crba.hxx:494](../include/pinocchio/src/algorithm/crba.hxx#L494)）：

```cpp
data.M.diagonal() += model.armature;
```

`armature` 是电机转子经减速比折算的**反射惯量**（reflected inertia），在关节空间是纯对角项，直接加到 $M$ 的对角线。它让有效质量矩阵变成 $M+\text{diag}(\text{armature})$，改善数值条件、贴合带减速器的真实关节（同 [RNEA §8.2](逆动力学RNEA解析.md)）。

---

## 11. mimic 关节：`MimicStep` 补齐

### 11.1 为什么需要单独一趟

mimic（仿从）关节的 `idx_v` **指向被仿（主动）关节的槽位**（见 [RNEA §7.1](逆动力学RNEA解析.md)）。主循环的后向趟（Pass 2）对 mimic 关节做了**截断**——只做组合惯量/力集的向上传播，**不填 M 的行块**。否则它会用 mimic 的 `idx_v`（= 主动关节槽）去写 M，把主动关节自己的数据**覆盖掉**。所以 mimic 对 $M$ 的贡献**推迟到 Pass 3**（`Crba*ConventionMimicStep`）逐个补回（[crba.hxx:344](../include/pinocchio/src/algorithm/crba.hxx#L344) 三步注释）：

```cpp
for (i=0; i<model.mimicking_joints.size(); ++i)
  Pass3::run(joints[mimicking_joints[i]], ..., i);
```

普通模型 `mimicking_joints` 为空，Pass 3 空转。下面以 LOCAL 约定（[crba.hxx:387](../include/pinocchio/src/algorithm/crba.hxx#L387)）为例，拆成四部分。

### 11.2 ① 算 mimic 的力集 `F = Y^c_mimic · S`

```cpp
auto & F = data.Fcrb[mimicking_id];
jmodel.jointCols(F) = data.Ycrb[mimicking_id] * jdata.S();
```

- `Ycrb[mimicking_id]` = 以 mimic 为根的**子树组合惯量**（Pass 2 逆序已把后代累加进来）；
- `jdata.S()` 是**带耦合倍率 $s$** 的运动子空间 $s\cdot S$（mimic 的 S 是 `ScaledJointMotionSubspace`，见 [§6](#6-local-约定实现经典-featherstone-体坐标版) 与 [RNEA §7.2](逆动力学RNEA解析.md)）；
- 得 ${}^{\text{mim}}F=Y^c_{\text{mim}}\cdot(sS)$。

### 11.3 ② 自身贡献 → 加到主动关节对角块

```cpp
jmodel.jointBlock(data.M).noalias() += jdata.S().transpose() * jmodel.jointCols(F);
```

`jointBlock(M)` 取 mimic 的 `(idx_v, idx_v)` 块，而 mimic 的 `idx_v` = 主动关节槽，故往 $M[\text{prim},\text{prim}]$ 累加：

$$
M_{\text{prim},\text{prim}}\mathrel{+}=(sS)^\top Y^c_{\text{mim}}(sS)=s^2\,S^\top Y^c_{\text{mim}}S
$$

正是 [RNEA §7.4](逆动力学RNEA解析.md) "mimic 惯性按 $s^2$ 反射到主动关节" 的质量矩阵版。

### 11.4 ③ 与下游真实 DoF 子树的耦合

```cpp
JointIndex sub_mimic_id = data.mimic_subtree_joint[mims_id];
if (sub_mimic_id != 0)
  jmodel.jointRows(data.M).middleCols(idx_vs[sub_mimic_id], nvSubtree[sub_mimic_id])
      += jdata.S().transpose() * F.middleCols(idx_vs[sub_mimic_id], nvSubtree[sub_mimic_id]);
```

`mimic_subtree_joint[mims_id]` = mimic 子树里**第一个有独立 DoF 的关节**（[data.hxx:861](../include/pinocchio/src/multibody/data.hxx#L861)），没有则为 0。这些后代**有自己的 M 列**（不共享主动槽），故它们与 mimic 的耦合要单独补 $M_{\text{prim},\,\text{sub}}\mathrel{+}=S^\top F_{[\text{sub 列}]}$。

### 11.5 ④ 沿支撑链往根：补 mimic 与各祖先的耦合

```cpp
const auto & supports = model.supports[mimicking_id];   // 根 → mimic 路径
size_t j = supports.size() - 2;                         // 从 mimic 的父开始（size-1 是自己）
for (; j > 0; j--) {
  const JointIndex & li = supports[j];        // 当前祖先
  const JointIndex & i  = supports[j + 1];    // 路径上的子（提供 liMi[i]）

  forceSet::se3Action(data.liMi[i], jmodel.jointCols(F), jmodel.jointCols(F));  // F 原地上移一级
  int parent_idx = model.idx_vs[li];
  applyConstraintOnForceVisitor<SETTO>(data.joints[li], jmodel.jointCols(F), temp_JAG.noalias());
  // temp_JAG = S_li^T · F  ——  祖先 li 子空间对 mimic 力集的投影 = 耦合块

  if (parent_idx >= mims_idx_v) {                          // 祖先槽 ≥ mimic(主动)槽
    data.M.block(mims_idx_v, parent_idx, ...) += temp_JAG; //   填 (mim行, 祖先列)
    if (mims_idx_v == parent_idx)                          //   祖先与 mimic 同主动槽 → 计两次
      data.M.block(mims_idx_v, parent_idx, ...) += temp_JAG;
  } else {                                                 // 祖先槽 < mimic 槽
    data.M.block(parent_idx, mims_idx_v, ...) += temp_JAG; //   填 (祖先行, mim列)
  }
}
```

- **`forceSet::se3Action(liMi[i], F, F)`**：`src==dst`，把 mimic 力集**原地**从 body $i$ 系变换到父系 $li$（逐列余伴随，视图就地写入见 [C++技巧 §18](C++语法技巧.md)）。每转一级得 ${}^{li}F_{\text{mim}}={}^{li}X_{\text{mim}}^{*}Y^c_{\text{mim}}(sS)$。
- **`applyConstraintOnForceVisitor<SETTO>`**（[joint-basic-visitors.hxx:1617](../include/pinocchio/src/multibody/joint/joint-basic-visitors.hxx#L1617)）对祖先 `li` 算 `temp_JAG = S_li^T · F`，即耦合 $M_{li,\text{mim}}$。

**上三角 index 定位（`min/max` 逻辑）**——普通关节里祖先 `idx_v` **总比后代小**，$M[\text{祖先},\text{后代}]$ 天然上三角；但 mimic 的 `idx_v` = 主动关节槽，主动关节在树里位置与 mimic 不同，祖先槽可能**大也可能小**。故比较 `parent_idx` 与 `mims_idx_v`，把耦合块放到"**较小者做行**"的位置，**始终落上三角**（下三角靠对称补）。`mims_idx_v == parent_idx` 时祖先**共享同一主动槽**（祖先就是被仿关节、或另一同主动的 mimic），贡献**计两次**。

### 11.6 小结与 WORLD 差异

本质：**mimic 没有自己的行列，它的所有惯性耦合都要手动映射/反射到主动关节的槽位再累加**（`idx_v` 映射 + 倍率 $s$/$s^2$），与 RNEA 里 mimic 处理同源。WORLD 约定的 MimicStep 结构相同，还额外把 mimic 对**质心动量矩阵**的贡献加回 `Ag`（[crba.hxx:217](../include/pinocchio/src/algorithm/crba.hxx#L217)）。

---

## 12. 与 $A_g$、ABA、Cholesky 的关系

- **$A_g$（质心动量矩阵）**：WORLD 约定 Pass 2 的力集存进 **`data.Ag`** 并非巧合——$Y_i^cS_i$ 既是填 $M$ 的中间力集，其在质心系的表达正是 CMM 的列（$h_g=A_g\dot q$）。所以 CRBA(WORLD) 与 [`ccrba`](质心与质心动量解析.md) 共享同一套组合惯量与力集，$A_g$ 是 CRBA 的免费副产物。
- **ABA**：都用组合惯量作种子，但 ABA 用 Schur 补塌缩（铰接体），CRBA 用直接相加（焊死）。见 [正动力学ABA解析.md](正动力学ABA解析.md)。
- **Cholesky**：CRBA 产出的 $M$（上三角 + 分支稀疏）正好喂给 `cholesky::decompose` 做 $M=UDU^\top$，稀疏结构一致（都源于支撑链），配合 $O(nd)$ 解 $Mx=y$。见 [multibody §3.5](multibody子系统解析.md)。
- **正动力学**：$\ddot q=M^{-1}(\tau-b)$ 里，CRBA 求 $M$、RNEA 求 $b$、Cholesky 求解——这是 ABA 之外的另一条路。

---

## 13. 复杂度

- **CRBA** $O(nd)$：$d$ 为树平均深度。链式结构 $d\approx n$ → $O(n^2)$；分支多则 $d$ 小。填 $M$ 每条只碰子树列，总代价正比于 $\sum_i n_{v_i}\cdot\text{nvSubtree}[i]$。
- 组装 $M$ 后若要 $M^{-1}$：Cholesky $O(nd)$ 解方程，或 `computeMinverse`（ABA 路线）$O(n^2)$ 直接求逆。
- 只依赖 $q$，可在位形固定时缓存复用。

---

## 14. 一句话总结

CRBA = **组合刚体惯量两趟 $O(nd)$ 递推**，Pinocchio 提供 LOCAL / WORLD 两种等价实现：

$$
\underbrace{\text{根→叶}}_{\text{FK}+\text{播种}Y^c}\ \Rightarrow\
\underbrace{\text{叶→根}}_{F=Y^cS,\ M_{ij}=S_i^\top F,\ Y_\lambda^c\mathrel{+}=Y_i^c}
$$

- **组合刚体惯量** $Y^c=\sum I_j$ = 子树"焊死"成一个刚体的合惯量，**直接累加**（对比 ABA 的 Schur 补塌缩）；
- **$M_{ij}=S_i^\top Y^c S_j$**：子树惯量作用于子空间得力集 $F=Y^cS$，再投影 $S^\top F$；
- **LOCAL** 在各自体系递推、累加需 `liMi.act`（经典 Featherstone）；**WORLD** 先全变世界系、累加免变换、附赠 $A_g$（Pinocchio 优化）；
- 只填**上三角**、只碰**祖先-后代对**（要求紧凑树 `CRBAChecker`），天然给出分支诱导稀疏，直通 Cholesky 与质心动量矩阵；末尾补 `armature`，mimic 由 `MimicStep` 补齐。
