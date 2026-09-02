# 质量矩阵 CRBA 解析：`algorithm/crba.hxx`

> **CRBA**（Composite Rigid Body Algorithm，组合刚体算法，Featherstone 第 6 章）计算**关节空间惯性矩阵** $M(q)$（质量矩阵，Featherstone 记作 $H$）。
> 本文**逐行贴 [src/algorithm/crba.hxx](../include/pinocchio/src/algorithm/crba.hxx)** 讲透其两种实现约定：
> 组合刚体惯量的意义、$M_{ij}=S_i^\top Y^c S_j$ 的推导、**LOCAL（经典 Featherstone）vs WORLD（Pinocchio 优化）** 两条实现路线、
> 上三角与 `nvSubtree`、紧凑树要求、armature、mimic 补齐。
> **§6 / §7 按代码执行顺序逐行读完两个驱动函数**（驱动循环 → Pass 1 → Pass 2 普通/mimic 重载 → Pass 3 → 收尾），
> 并拆到底层原语：访问者分派、索引访问器族、`nvSubtree` 的写入范围、`Inertia::__mult__`、`SE3::act(Inertia)`、
> `forceSet::se3Action`、`motionSet::inertiaAction` 的真实算式，另附五关节分叉树的完整走查与耗时/内存实测。
> 配套：[逆动力学RNEA解析.md](逆动力学RNEA解析.md)、[正动力学ABA解析.md](正动力学ABA解析.md)、
> [质心与质心动量解析.md](质心与质心动量解析.md)、[multibody 子系统解析 §3.5（M 的稀疏 Cholesky）](multibody子系统解析.md)、[空间代数运算解析.md](空间代数运算解析.md)。

---

## 目录

- [1. 质量矩阵是什么](#1-质量矩阵是什么)
- [2. 核心概念：组合刚体惯量](#2-核心概念组合刚体惯量composite-rigid-body-inertia)
- [3. 关键公式：$M_{ij}=S_i^\top Y^c S_j$](#3-关键公式m_ijs_itop-yc-s_j)
- [4. 顶层入口与两种约定](#4-顶层入口与两种约定)
- [5. 涉及的 `Data` 数据结构](#5-涉及的-data-数据结构)
  - [5.1 共享机制①：访问者如何拿到具体关节类型](#51-三趟-pass-共享的机制之一访问者如何拿到具体关节类型)
  - [5.2 共享机制②：索引访问器族（含 mimic 重写）](#52-三趟-pass-共享的机制之二索引访问器族)
  - [5.3 共享机制③：`nvSubtree` 决定往 M 的哪块写](#53-三趟-pass-共享的机制之三nvsubtree-决定往-m-的哪块写)
- [6. LOCAL 约定：逐行读完 `crbaLocalConvention`](#6-local-约定实现逐行读完-crbalocalconvention)
  - [6.1 驱动函数：三趟循环](#61-驱动函数三趟循环)
  - [6.2 Pass 1](#62-pass-1crbalocalconventionforwardstep)
  - [6.3 Pass 2（普通关节）](#63-pass-2普通关节crbalocalconventionbackwardstepalgo_impl)
  - [6.4 在五关节分叉树上完整走一遍](#64-在五关节分叉树上完整走一遍)
  - [6.5 Pass 2（mimic 关节）](#65-pass-2mimic-关节同名函数的另一个重载)
  - [6.6 Pass 3](#66-pass-3crbalocalconventionmimicstep)
- [7. WORLD 约定：逐行读完 `crbaWorldConvention`](#7-world-约定实现逐行读完-crbaworldconvention)
  - [7.1 驱动函数：多了首尾两段](#71-驱动函数多了首尾两段)
  - [7.2 Pass 1](#72-pass-1crbaworldconventionforwardstep)
  - [7.3 Pass 2（普通关节）：从四步减到三步](#73-pass-2普通关节从四步减到三步)
  - [7.4 Pass 2（mimic 关节）：只剩一行](#74-pass-2mimic-关节只剩一行)
  - [7.5 Pass 3](#75-pass-3crbaworldconventionmimicstep)
  - [7.6 收尾那 6 行：质心动量矩阵](#76-收尾那-6-行质心动量矩阵是怎么顺带算出来的)
- [8. 两种约定对比：差异到底在哪](#8-两种约定对比差异到底在哪)
  - [8.1 实测数据：耗时、内存与一致性](#81-实测数据耗时内存与一致性)
- [9. 上三角、`nvSubtree` 与紧凑树要求（`CRBAChecker`）](#9-上三角nvsubtree-与紧凑树要求crbachecker)
- [10. `armature`：转子惯量](#10-armature转子惯量)
- [11. mimic 关节：为什么要拆成"截断 + 补账"](#11-mimic-关节为什么要拆成截断--补账)
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

### 5.1 三趟 Pass 共享的机制之一：访问者如何拿到具体关节类型

三趟 Pass 都写成 `fusion::JointUnaryVisitorBase` 的派生类，调用形式统一是：

```cpp
Pass2::run(model.joints[i], data.joints[i], typename Pass2::ArgsType(model, data));
```

`model.joints[i]` 的静态类型是 `JointModel`——一个 **`boost::variant`**，运行期才知道里面装的是 `JointModelRX` 还是 `JointModelFreeFlyer`。`run` 内部做的事是（[joint-unary-visitor.hxx:221](../include/pinocchio/src/multibody/visitor/joint-unary-visitor.hxx#L221)）：

```cpp
return boost::apply_visitor(visitor, jmodel);       // 按 variant 的 which() 分派

template<typename JointModelDerived>
ReturnType operator()(const JointModelBase<JointModelDerived> & jmodel) const
{
  return bf::invoke(
    &JointVisitorDerived::template algo<JointModelDerived>,   // ← 取具体类型的 algo 实例
    bf::append(boost::ref(jmodel.derived()),
               boost::ref(boost::get<JointDataDerived>(jdata)), args));
}
```

**关键后果**：`algo<JointModel>` 会为**每一种关节类型各编译一份**。于是在 `algo` 体内：

- `JointModel::NV` 是**编译期常量**（旋转关节 = 1，free-flyer = 6，`JointModelComposite`/`JointModelMimic` = `Eigen::Dynamic`）；
- `jdata.S()` 的类型是该关节特有的运动子空间。**很多关节的子空间是不含任何数据成员的空类型**——
  实测 `sizeof(JointDataRX::Constraint_t) == 1`、`sizeof(JointDataFreeFlyer::Constraint_t) == 1`
  （C++ 空类占 1 字节）。旋转关节的 `S` 只是"第 4/5/6 行取 1"的编译期标记，`S^T * F` 会被优化成
  **一次取行**，不产生任何乘法；free-flyer 的 `S` 是 $6\times6$ 单位阵，同样无需存储。
  （反例：`JointModelRevoluteUnaligned` 的轴向是运行期数据，其 `S` 才真的占内存。）
- 所有 `middleCols` / `block` 都走 `SizeDepType<NV>`，定长时返回**固定尺寸的 Eigen 块表达式**，编译器可完全展开循环。

这就是 Pinocchio "运行期决定自由度、却仍享受 Eigen 定长优化"的落地点——**一次 `apply_visitor` 的间接跳转，换来函数体内部全部定长**。

> 注意两趟的 `run` 签名不同：LOCAL 的 Pass 2 传 `(model.joints[i], data.joints[i], args)`（需要 `jdata.S()`），
> **WORLD 的 Pass 2 只传 `(model.joints[i], args)`**——因为它的子空间早在 Pass 1 就写进 `data.J` 了，
> 后向趟不再需要 `jdata`。这是 WORLD 约定省下的一笔隐性开销。

### 5.2 三趟 Pass 共享的机制之二：索引访问器族

CRBA 全程不手写 `idx_v`/`nv` 的加减，而是用 `JointModelBase` 的一族访问器。它们的**实际展开**（[joint-model-base.hxx:412](../include/pinocchio/src/multibody/joint/joint-model-base.hxx#L412)）：

| 访问器 | 展开成 | 起点 | 宽度 |
|---|---|---|---|
| `jmodel.jointCols(A)` | `SizeDepType<NV>::middleCols(A, idx_v(), nv())` | `idx_v` | `nv` |
| `jmodel.jointExtendedModelCols(A)` | `SizeDepType<NVExtended>::middleCols(A, idx_vExtended(), nvExtended())` | `idx_vExtended` | `nvExtended` |
| `jmodel.jointRows(A)` | `middleRows(A, idx_v(), nv())` | `idx_v` | `nv` |
| `jmodel.jointBlock(M)` | `block(M, idx_v(), idx_v(), nv(), nv())` | 对角 | `nv × nv` |

**两套索引的分工**（普通模型下两者完全相同，只有 mimic 模型才分叉）：

| | 索引 | 用在 |
|---|---|---|
| **扩展索引** `idx_vExtended` / `nvExtended` | 每个关节（含 mimic）**各占一列** | `data.J`（世界系运动子空间） |
| **真实索引** `idx_v` / `nv` | mimic **不占列**，与主动关节**共用** | `data.M`、`data.Ag`、`data.Fcrb` |

所以 WORLD 的 Pass 1 写 `jointExtendedModelCols(data.J)`、Pass 2 写 `jointCols(data.Ag)`，**不是笔误**——前者按扩展列布局，后者按真实列布局。

**`JointModelMimic` 重写了这三个访问器**（[joint-mimic.hxx:815](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L815)）：

```cpp
jointCols_impl(A)  { return SizeDepType<NV>::middleCols(A, Base::i_v, m_nvExtended); }
jointRows_impl(A)  { return SizeDepType<NV>::middleRows(A, Base::i_v, m_nvExtended); }
jointBlock_impl(M) { return SizeDepType<NV,NV>::block(M, Base::i_v, Base::i_v,
                                                      m_nvExtended, m_nvExtended); }
```

即 **起点取 `i_v`（= 主动关节的列号），宽度取 `m_nvExtended`（= 自己的维度）**。必须重写，因为 mimic 的 `nv()` 恒为 0（[joint-mimic.hxx:618](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L618)），沿用基类版本会得到**宽度为 0 的空视图**，Pass 3 什么也写不进去。

> 实测（12 个 mimic 测试用例、共 20 个 mimic 关节）：
> `nv()==0` 命中 20/20；`idx_v` 落在某个 `nv>0` 关节的列块内 20/20——**列别名 100% 成立**。

---

### 5.3 三趟 Pass 共享的机制之三：`nvSubtree` 决定"往 M 的哪块写"

这是读懂 Pass 2 唯一那行 `data.M.block(...)` 的前提，单独讲清楚。

#### 它是怎么算出来的

`Data` 构造时一次性算好（[data.hxx:851](../include/pinocchio/src/multibody/data.hxx#L851)）：

```cpp
for (int i = model.njoints - 1; i >= 0; --i)   // 逆序 = 从叶到根
{
  const Index parent = model.parents[i];
  nvSubtree[i]      += model.joints[i].nv();   // 先把自己的自由度算上
  nvSubtree[parent] += nvSubtree[i];           // 再整个甩给父节点
}
```

逆序保证处理到 `i` 时它所有孩子都已把自己的总数汇给它。结果：

$$ \texttt{nvSubtree}[i] \;=\; \sum_{j\,\in\,\text{subtree}(i)} n_v(j) \qquad(\text{含 } i \text{ 自己}) $$

#### 它为什么能当"连续列区间"用

`nvSubtree[i]` 只是一个**数量**。Pass 2 却直接拿它当**列宽**：

```cpp
data.M.block(jmodel.idx_v(), jmodel.idx_v(), jmodel.nv(), data.nvSubtree[i])
//            ↑行起点        ↑列起点         ↑行高        ↑列宽
```

这隐含了一个强假设：**关节 $i$ 子树里所有关节的速度索引，恰好构成从 `idx_v(i)` 开始的一段连续区间**

$$ \bigl\{\,\texttt{idx\_v}(j)\ :\ j\in\text{subtree}(i)\,\bigr\} \;=\; \bigl[\ \texttt{idx\_v}(i),\ \ \texttt{idx\_v}(i)+\texttt{nvSubtree}[i]\ \bigr) $$

这个假设由**紧凑树编号**保证（DFS 建树天然满足），并由 `CRBAChecker` 在 `model.check(data)` 里验证——见 §9.3。不满足时 `block(...)` 会静默取错列。

#### 实测：五关节分叉树

```
     0 ──j1── j2 ── j3
            └─ j4 ── j5
```

| id | 关节 | 父 | nv | idx_v | subtree | nvSubtree | Pass 2 写入的 M 块 |
|---|---|---|---|---|---|---|---|
| 1 | j1 | 0 | 1 | 0 | {1,2,3,4,5} | **5** | `M[0, 0..4]` |
| 2 | j2 | 1 | 1 | 1 | {2,3} | **2** | `M[1, 1..2]` |
| 3 | j3 | 2 | 1 | 2 | {3} | **1** | `M[2, 2..2]` |
| 4 | j4 | 1 | 1 | 3 | {4,5} | **2** | `M[3, 3..4]` |
| 5 | j5 | 4 | 1 | 4 | {5} | **1** | `M[4, 4..4]` |

把 `data.M` 预置成 `NaN` 再跑一次 `crba`，看哪些格子真的被写过（`W` = 被写，`.` = 仍是 NaN）：

```
      0 1 2 3 4
  0   W W W W W      ← j1 是所有人的祖先，整行写满
  1   . W W . .      ← j2 只碰 {j2,j3}
  2   . . W . .
  3   . . . W W      ← j4 只碰 {j4,j5}
  4   . . . . W
```

**每个关节写且只写"自己这行 × 自己子树那几列"**，五块拼起来正好覆盖上三角的非零部分。

`M(2,4)`（j3 与 j5）**始终是 NaN**——它们分属左右两支，互不为祖先后代，CRBA **根本没访问过这个格子**。这就是分支诱导稀疏（§9.4）在代码层面的样子。

> ⚠️ **实用提醒**：上面这张图说明 `data.M` 里"没被写"的格子保留的是**调用前的旧值**。
> `Data` 构造时 `M` 初始化为零，所以通常读出来是 0；但如果你自己往 `data.M` 写过东西，
> CRBA **不会替你清零**那些结构零位置。

#### 它顺带说明了复杂度

Pass 2 写入的 M 元素总数就是 $\sum_i \texttt{nvSubtree}[i]$。实测 30 关节浮动基人形（$n_v=34$）：

| 量 | 值 |
|---|---|
| $\sum_i \texttt{nvSubtree}[i]$（实际写入的元素数） | **152** |
| $n_v(n_v{+}1)/2$（稠密上三角元素数） | 595 |

只写了 **26%**——剩下 74% 是结构零，一次都没碰。这就是 CRBA 比"$\sum_i J_i^\top I_i J_i$ 硬算"快的根本原因。

---

## 6. LOCAL 约定实现：逐行读完 `crbaLocalConvention`

本节按**代码的实际执行顺序**走一遍：先看驱动函数，再依次拆 Pass 1 / Pass 2 / Pass 3，最后在一棵具体的树上完整走一遍数。

### 6.1 驱动函数：三趟循环

[crba.hxx:464](../include/pinocchio/src/algorithm/crba.hxx#L464)，去掉模板参数后的全貌：

```cpp
const MatrixXs & crbaLocalConvention(const Model & model, Data & data, const VectorXs & q)
{
  assert(model.check(data) && "data is not consistent with model.");        // ← 校验紧凑树，见 §9.3
  PINOCCHIO_CHECK_ARGUMENT_SIZE(q.size(), model.nq, "...");

  // Pass 1：根 → 叶
  typedef CrbaLocalConventionForwardStep<...> Pass1;
  for (JointIndex i = 1; i < (JointIndex)model.njoints; ++i)
    Pass1::run(model.joints[i], data.joints[i], typename Pass1::ArgsType(model, data, q));

  // Pass 2：叶 → 根
  typedef CrbaLocalConventionBackwardStep<...> Pass2;
  for (JointIndex i = (JointIndex)(model.njoints - 1); i > 0; --i)
    Pass2::run(model.joints[i], data.joints[i], typename Pass2::ArgsType(model, data));

  // Pass 3：只遍历 mimic 关节（普通模型下这个循环一次都不进）
  typedef CrbaLocalConventionMimicStep<...> Pass3;
  for (size_t i = 0; i < model.mimicking_joints.size(); i++)
    Pass3::run(model.joints[model.mimicking_joints[i]],
               data.joints[model.mimicking_joints[i]],
               typename Pass3::ArgsType(model, data, i));

  data.M.diagonal() += model.armature;    // ← §10
  return data.M;
}
```

三个循环的**方向**是算法的骨架：

| 趟 | 方向 | 循环变量 | 为什么必须是这个方向 |
|---|---|---|---|
| Pass 1 | `i = 1 → njoints-1` | 关节 id 升序 | 算 `oMi[i]` 要先有 `oMi[parent]`，父 id 恒小于子 id |
| Pass 2 | `i = njoints-1 → 1` | 关节 id 降序 | 用 `Ycrb[i]` 时它必须已经吸收了全部后代 |
| Pass 3 | 只跑 `mimicking_joints` | — | 补 Pass 2 对 mimic 故意跳过的部分 |

> **"父 id 恒小于子 id"** 由 `addJoint` 的构建顺序保证，`CRBAChecker` 会验证（§9.3）。
> 正因为如此，"降序遍历 id" 就等价于 "从叶到根"，不需要真的做后序 DFS。

### 6.2 Pass 1：`CrbaLocalConventionForwardStep`

[crba.hxx:237](../include/pinocchio/src/algorithm/crba.hxx#L237)，函数体只有三行：

```cpp
const JointIndex & i = jmodel.id();

jmodel.calc(jdata.derived(), q.derived());               // ①
data.liMi[i] = model.jointPlacements[i] * jdata.M();     // ②
data.Ycrb[i] = model.inertias[i];                        // ③
```

**① `jmodel.calc(jdata, q)`**
先用 `jointConfigSelector(q)` 从整段 $q$ 里切出本关节那几个分量（`SizeDepType<NQ>::segment(q, idx_q(), nq())`，零拷贝视图），算出两样东西塞进 `jdata`：

- `jdata.M()`——关节自身的位形变换 ${}^{J_i}M_i(q_i)$。绕 X 轴旋转关节就是 $(\mathrm{Rot}_x(\theta),\,0)$。
- `jdata.S()`——运动子空间 $S_i$（$6\times n_v(i)$）。

> `S` 常常**不占内存**：实测 `sizeof(JointDataRX::Constraint_t) == 1`、`sizeof(JointDataFreeFlyer::Constraint_t) == 1`（C++ 空类占 1 字节）。
> 旋转关节的 `S` 只是"第 4 行取 1"的编译期标记，后面 `S^T * F` 会被优化成**一次取行**，不产生乘法。
> （反例：`JointModelRevoluteUnaligned` 的轴向是运行期数据，它的 `S` 才真占内存。）

**② `liMi[i] = jointPlacements[i] * jdata.M()`**
把 URDF 里写死的安装位姿和关节转动串起来，得到父连杆到子连杆的相对位姿：

$$ {}^{\lambda(i)}M_i \;=\; \underbrace{{}^{\lambda(i)}M_{J_i}}_{\texttt{jointPlacements[i]}}\;\cdot\;\underbrace{{}^{J_i}M_i(q_i)}_{\texttt{jdata.M()}} $$

**③ `Ycrb[i] = model.inertias[i]`**
组合惯量的**种子**——纯拷贝 10 个数（$m$、质心 $c$、绕质心转动惯量 $\bar I_C$ 的 6 个分量），**不做任何坐标变换**。因为 `model.inertias[i]` 本来就存在连杆 $i$ 自己的坐标系里，而 LOCAL 约定的 `Ycrb[i]` 要的正是这个系。

> **注意 Pass 1 完全没碰 `data.M` 和 `data.Fcrb`。** 它只是把后面要用的 `liMi`、`Ycrb` 铺好。

### 6.3 Pass 2（普通关节）：`CrbaLocalConventionBackwardStep::algo_impl`

[crba.hxx:276](../include/pinocchio/src/algorithm/crba.hxx#L276)。这是整个 CRBA 的心脏，四步：

```cpp
const JointIndex & i = jmodel.id();

// ① 本关节的力集列：F[:,i] = Y^c_i · S_i
jmodel.jointCols(data.Fcrb[i]) = data.Ycrb[i] * jdata.S();

// ② 填 M：本关节这(几)行 × 自己子树那几列
data.M.block(jmodel.idx_v(), jmodel.idx_v(), jmodel.nv(), data.nvSubtree[i]).noalias() =
    jdata.S().transpose() * data.Fcrb[i].middleCols(jmodel.idx_v(), data.nvSubtree[i]);

const JointIndex & parent = model.parents[i];
if (parent > 0)
{
  // ③ 组合惯量塌缩给父
  data.Ycrb[parent] += data.liMi[i].act(data.Ycrb[i]);

  // ④ 把自己子树的力集整体搬到父的坐标系
  Block jF = data.Fcrb[parent].middleCols(jmodel.idx_v(), data.nvSubtree[i]);
  Block iF = data.Fcrb[i].middleCols(jmodel.idx_v(), data.nvSubtree[i]);
  forceSet::se3Action(data.liMi[i], iF, jF);
}
```

#### 先搞清楚 `Fcrb[i]` 里到底存着什么

这是最容易卡住的地方。`data.Fcrb` 是**每个关节一个 $6\times n_v$ 矩阵**（[data.hxx:713](../include/pinocchio/src/multibody/data.hxx#L713)），但任一时刻只有一小段列是有意义的。

**当 Pass 2 处理到关节 $i$ 时**（此时它所有后代都已处理完），`Fcrb[i]` 的第 `idx_v(i) .. idx_v(i)+nvSubtree[i]-1` 列存的是：

$$ \texttt{Fcrb}[i]\bigl[:,\ k\bigr] \;=\; {}^{i}X_{j(k)}^{*}\;Y^c_{j(k)}\,S_{j(k)} \qquad k\in\text{子树列区间} $$

读作：**"让子树里第 $k$ 个自由度获得单位加速度，需要施加的空间力——写在关节 $i$ 自己的坐标系里"**。

这些列是怎么进来的？

- **第 `idx_v(i)` 列**（自己那列）由**本次迭代的第 ① 步**刚写入；
- **其余列**由各个后代在**它们自己那次迭代的第 ④ 步**逐级搬上来——每上一级做一次 `forceSet::se3Action`，所以搬到 $i$ 时已经换算成 $i$ 的坐标系了。

#### 于是第 ② 步就是一次矩阵乘搞定一整行

$$ M_{i,k} \;=\; S_i^\top\,\texttt{Fcrb}[i][:,k] \;=\; S_i^\top\,{}^{i}X_{j(k)}^{*}\,Y^c_{j(k)}\,S_{j(k)} $$

这正是 §3 的 $M_{ij}=S_i^\top Y^c S_j$，只是把"变换到公共坐标系"这件事分摊到了第 ④ 步。因为子树列是**连续区间**（§5.3），`middleCols` 一刀切出来，一次 GEMM 填完关节 $i$ 的整个行块。

#### 第 ③ 步：`Ycrb[parent] += liMi[i].act(Ycrb[i])`

$$ Y^c_{\lambda} \mathrel{+}= {}^{\lambda}X_i^{*}\,Y^c_i\,{}^{i}X_{\lambda} $$

`InertiaTpl::se3Action_impl`（[inertia.hxx:1034](../include/pinocchio/src/spatial/inertia.hxx#L1034)）**从不组装 $6\times6$ 矩阵**，只有一行：

```cpp
return InertiaTpl(mass(), M.translation() + M.rotation() * lever(), inertia().rotate(M.rotation()));
```

$$ (m,\;c,\;\bar I_C)\;\longmapsto\;(m,\;t + Rc,\;R\,\bar I_C\,R^\top) $$

质量不变、质心按 SE3 搬、转动惯量只做旋转相合。平移不影响 $\bar I_C$ ——因为它是**关于质心**的惯量，平行轴定理那个 $-m\hat c\hat c$ 项被"显式存质心 $c$"这个表示吸收掉了。本该是两次 $6\times6$ 乘法（约 400 次浮点乘），压成一次 $3\times3$ 旋转相合。

`operator+=` 则是"质量相加 + 质心按质量加权平均 + 转动惯量平移到新质心后相加"。

**这一步是 CRBA 与 ABA 的分水岭**：CRBA 直接相加（子关节视为焊死），ABA 要做 Schur 补（子关节可自由转动）。

#### 第 ④ 步：`forceSet::se3Action(liMi[i], iF, jF)`

把子树力集从 $i$ 系批量搬到父系，写进 `Fcrb[parent]` 的**同一段列**：

$$ {}^{\lambda}F_k = {}^{\lambda}X_i^{*}\;{}^{i}F_k,\qquad
\begin{bmatrix} f'\\ \tau'\end{bmatrix}=\begin{bmatrix} R & 0\\ \hat t R & R\end{bmatrix}\begin{bmatrix} f\\ \tau\end{bmatrix} $$

实现是**逐列**调 `SE3::act(Force)`（[act-on-set.hxx:186](../include/pinocchio/src/spatial/act-on-set.hxx#L186)）：

```cpp
f.linear()  = m.rotation() * linear();
f.angular() = m.rotation() * angular();
f.angular() += m.translation().cross(f.linear());
```

源码在这里留了句实话：

> *Specialized implementation of block action, using colwise operation. It is empirically **much faster** than the true block operation, although I do not understand why.*

即**逐列变换比"先组装 $6\times6$ 的 $X^*$ 再一次矩阵乘"更快**——组装矩阵的开销加上缓存不友好，盖过了 GEMM 的向量化收益。

`se3Action` 的模板参数 `Op ∈ {SETTO, ADDTO, RMTO}` 决定写入方式是 `=` / `+=` / `-=`。这里用默认的 `SETTO`（**覆盖**）——待会儿 mimic 那版就是把它换成 `ADDTO`，其余不变。

### 6.4 在五关节分叉树上完整走一遍

```
     0 ──j1── j2 ── j3            idx_v:  j1=0  j2=1  j3=2  j4=3  j5=4
            └─ j4 ── j5           nvSubtree: 5, 2, 1, 2, 1
```

Pass 2 按 `i = 5,4,3,2,1` 执行：

| 迭代 | ① 写入 `Fcrb[i]` 的列 | ② 写入 `M` | ③ 惯量塌缩 | ④ 力集上搬 |
|---|---|---|---|---|
| `i=5` | `Fcrb[5][:,4] = Ycrb[5]·S₅` | `M[4, 4]` | `Ycrb[4] += ᵃX·Ycrb[5]` | `Fcrb[4][:,4] ← Fcrb[5][:,4]` |
| `i=4` | `Fcrb[4][:,3] = Ycrb[4]·S₄` | `M[3, 3..4]` | `Ycrb[1] += ᵃX·Ycrb[4]` | `Fcrb[1][:,3..4] ← Fcrb[4][:,3..4]` |
| `i=3` | `Fcrb[3][:,2] = Ycrb[3]·S₃` | `M[2, 2]` | `Ycrb[2] += ᵃX·Ycrb[3]` | `Fcrb[2][:,2] ← Fcrb[3][:,2]` |
| `i=2` | `Fcrb[2][:,1] = Ycrb[2]·S₂` | `M[1, 1..2]` | `Ycrb[1] += ᵃX·Ycrb[2]` | `Fcrb[1][:,1..2] ← Fcrb[2][:,1..2]` |
| `i=1` | `Fcrb[1][:,0] = Ycrb[1]·S₁` | `M[0, 0..4]` | `parent==0`，跳过 | 跳过 |

几个值得盯住的点：

- **`i=4` 的第 ② 步**读的是 `Fcrb[4]` 的第 3、4 两列：第 3 列是本次 ① 刚写的，第 4 列是 `i=5` 迭代的 ④ 搬上来的。
- **`i=1` 的第 ② 步**一次填满 `M[0, 0..4]`：第 0 列本次写，第 1–2 列来自 j2 分支（`i=2` 的 ④），第 3–4 列来自 j4 分支（`i=4` 的 ④）。此时 `Ycrb[1]` 已经吸收了整棵树。
- **两条支链从不相遇**：j2 分支只往 `Fcrb[·]` 的第 1–2 列写，j4 分支只往第 3–4 列写，两段列区间不相交。所以 `M(2,4)`（j3 × j5）**永远不会被写**——结构零（§5.3 的实测图里它一直是 `NaN`）。

### 6.5 Pass 2（mimic 关节）：同名函数的另一个重载

[crba.hxx:318](../include/pinocchio/src/algorithm/crba.hxx#L318)。注意这**不是 if 分支**，而是 C++ 的**函数重载**——`algo` 里那句 `algo_impl(jmodel, jdata, model, data)` 在编译期就按 `jmodel` 的具体类型选中了不同的函数体，运行期零开销：

```cpp
static void algo_impl(const JointModelBase<JointModelMimic> & jmodel,
                      JointDataBase<...> & /* jdata 根本没用 */,
                      const Model & model, Data & data)
{
  const JointIndex & i = jmodel.id();
  const JointIndex & parent = model.parents[i];
  if (parent > 0)
  {
    data.Ycrb[parent] += data.liMi[i].act(data.Ycrb[i]);        // ③ 照做

    typedef typename Data::Matrix6x & Matrix;                   // ← 注意是【引用】
    Matrix jF = data.Fcrb[parent];
    Matrix iF = data.Fcrb[i];
    forceSet::se3Action<ADDTO>(data.liMi[i], iF, jF);           // ④ 改成 ADDTO，且整块搬
  }
}
```

和普通版逐条对比：

| 步骤 | 普通关节 | mimic 关节 |
|---|---|---|
| ① `Fcrb[i] = Ycrb[i]·S` | 做 | **不做** |
| ② 写 `data.M` | 做 | **不做** |
| ③ `Ycrb[parent] +=` | 做 | 做 |
| ④ 力集上搬 | `SETTO`，只搬 `middleCols(idx_v, nvSubtree[i])` | **`ADDTO`**，搬**整个** $6\times n_v$ |

**为什么 ① ② 被砍掉？** 因为 mimic 关节的 `idx_v` **指向被模仿（主动）关节的列**，而 `nv()` 恒为 0。若照常执行第 ②，那句 `=` 会**覆盖掉主动关节刚算好的 M 条目**。所以只能先跳过，欠的账推迟到 Pass 3 用 `+=` 补（§11）。

**为什么第 ④ 从 `SETTO` 变 `ADDTO`、还搬整块？** 源码注释讲得很清楚：

> *Since we don't just copy columns, we need to use ADDTO method, to avoid ecrasing already filled columns, in case of parallel arms, it is important to not ecrase already filled data from the other part of the tree.*

mimic 的列不再是一段干净的连续子树区间（它和主动关节共用列），没法像普通关节那样切一刀；用 `ADDTO` 整块搬才不会抹掉并联支链已经传上来的力。

### 6.6 Pass 3：`CrbaLocalConventionMimicStep`

[crba.hxx:391](../include/pinocchio/src/algorithm/crba.hxx#L391)。同样用重载：普通关节命中的是**空函数体**，只有 mimic 走真正的实现。

```cpp
const JointIndex mimicking_id = jmodel.id();
auto & F = data.Fcrb[mimicking_id];

// ① 补上 Pass 2 没做的第 ①：算 mimic 的力集
jmodel.jointCols(F) = data.Ycrb[mimicking_id] * jdata.S();

// ② 自身贡献 → 累加到主动关节的对角块
jmodel.jointBlock(data.M).noalias() += jdata.S().transpose() * jmodel.jointCols(F);

// ③ 与下游"有独立 DoF 的子树"的耦合
JointIndex sub_mimic_id = data.mimic_subtree_joint[mims_id];
if (sub_mimic_id != 0)
  jmodel.jointRows(data.M)
        .middleCols(model.idx_vs[sub_mimic_id], data.nvSubtree[sub_mimic_id]).noalias()
      += jdata.S().transpose()
       * F.middleCols(model.idx_vs[sub_mimic_id], data.nvSubtree[sub_mimic_id]);

// ④ 沿支撑链往根走，补 mimic 与每个祖先的耦合
const auto & supports = model.supports[mimicking_id];    // = {0, ..., λ(mim), mim}
size_t j = supports.size() - 2;                          // 从 mimic 的父开始
for (; j > 0; j--)                                       // j=0 是 universe，不处理
{
  const JointIndex & li = supports[j];        // 当前祖先
  const JointIndex & i  = supports[j + 1];    // 路径上更深一级（提供 liMi[i]）

  forceSet::se3Action(data.liMi[i], jmodel.jointCols(F), jmodel.jointCols(F));  // 原地上移一级
  applyConstraintOnForceVisitor<SETTO>(data.joints[li], jmodel.jointCols(F), temp_JAG.noalias());
  //  temp_JAG = S_li^T · F   ← 祖先子空间对 mimic 力集的投影

  int parent_idx = model.idx_vs[li];
  if (parent_idx >= mims_idx_v) {                                       // 祖先列号更大 → 填行
    data.M.block(mims_idx_v, parent_idx, mims_nvExtended, model.nvExtendeds[li]) += temp_JAG;
    if (mims_idx_v == parent_idx)                                       // 同一列 → 再加一次
      data.M.block(mims_idx_v, parent_idx, mims_nvExtended, model.nvExtendeds[li]) += temp_JAG;
  } else {                                                              // 祖先列号更小 → 填列
    data.M.block(parent_idx, mims_idx_v, model.nvExtendeds[li], mims_nvExtended) += temp_JAG;
  }
}
```

四处需要解释的实现细节：

**(a) `jmodel.jointCols` / `jointBlock` / `jointRows` 在 mimic 上是被重写过的。**
基类版本用 `(idx_v(), nv())`，而 mimic 的 `nv()==0`，会得到**宽度为 0 的空视图**，什么也写不进去。`JointModelMimic` 重写成 `(i_v, m_nvExtended)`——**起点取主动关节的列号，宽度取自己的维度**（[joint-mimic.hxx:815](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L815)）。所以第 ② 步那句 `jointBlock(data.M) += ...` 实际落在主动关节的对角块上。

**(b) `forceSet::se3Action(liMi[i], F_cols, F_cols)` 是原地写。**
源和目标是同一个视图——因为每列的变换只依赖该列自己，逐列覆盖是安全的。每转一级，$F$ 就从 $i$ 系换算到 $li$ 系。

**(c) 为什么要用 `applyConstraintOnForceVisitor` 而不是直接写 `S^T·F`？**
因为循环里访问的是**别的关节** `data.joints[li]`，它的具体类型在这个函数里是未知的（`data.joints[li]` 是 variant）。所以要再套一层访问者做运行期分派，它做的事就一句（[joint-basic-visitors.hxx:1617](../include/pinocchio/src/multibody/joint/joint-basic-visitors.hxx#L1617)）：

```cpp
const_cast<ExpressionType &>(R) = jdata.S().transpose() * F;   // Op == SETTO
```

**(d) 那个 `if (parent_idx >= mims_idx_v)` 的行列翻转。**
普通关节里祖先的 `idx_v` **总是**比后代小，$M[\text{祖先},\text{后代}]$ 天然落在上三角。但 mimic 的 `idx_v` 是主动关节的列号，主动关节在树里的位置和 mimic 不同，**祖先的列号可能比它大也可能比它小**。所以要比大小，把耦合块放到"列号小的做行"的位置，保证始终落在上三角。

`mims_idx_v == parent_idx` 时说明这个祖先**和 mimic 共用同一列**（祖先就是被模仿的那个关节，或另一个同主动的 mimic），此时 $M_{ps}$ 与 $M_{sp}$ 折叠到同一格，所以**加两次**——数学上的来历见 §11.2。

**(e) `mims_id` 是序号，不是关节 id。** 驱动里 `Pass3::run(..., i)` 传的是**循环下标** `i`（[crba.hxx:485](../include/pinocchio/src/algorithm/crba.hxx#L485)），即"这是第几个 mimic"。它专门用来索引**和 `model.mimicking_joints` 平行的数组** `data.mimic_subtree_joint[mims_id]`（长度 = mimic 个数，非关节数）。真正的关节 id 是 `mimicking_id = jmodel.id() = model.mimicking_joints[mims_id]`，用于 `subtrees`/`supports`/`idx_vs`/`Fcrb` 等按关节 id 排的量。二者别混。

**(f) `sub_mimic_id` 是 mimic 的下游后代，不是被模仿关节。** `data.mimic_subtree_joint[mims_id]` = 该 mimic **下游第一个"真实自由度"后代的关节 id**（官方原话 "first non mimic child"，[data.hxx:860](../include/pinocchio/src/multibody/data.hxx#L860)），没有则存哨兵 **0**（关节 0 是 universe，永不为后代，零歧义）。**它在下游（id 更大），与上游的被模仿主动关节方向相反**。之所以要单存这个句柄：mimic 自己的 `idx_v` 已被劫持指向主动关节，无法用它定位"mimic 自己下游子树"的列，只能靠 `sub_mimic_id` 拿到下游列的起点 `idx_vs[sub_mimic_id]` + 宽度 `nvSubtree[sub_mimic_id]`。

> **第 ③ 步的物理意义**：填的是 $M[\text{主动关节},\ \text{mimic 下游后代}]$——**mimic 与它下游那截胳膊本是刚性串联、有惯性耦合的一对；但 mimic 没有独立账户，这份"mimic ↔ 下游"耦合按齿轮比 $r$ 折算、记到主动关节行上**。即"推下游关节时主动电机会感到的那份力矩"（反之亦然）。这里只带**一个** $r$（一端是 mimic、一端是真关节），对比第 ② 步两端都是 mimic 得 $r^2$。这份账 Pass 2 漏算了（mimic 那行被截断），故 Pass 3 补。

> ⚠️ **已知假设 / 局限：mimic 下游须为单链。** 这步用的宽度是 `nvSubtree[sub_mimic_id]`（**第一棵**下游真子树），不是 `nvSubtree[mimicking_id]`（全部下游）。当 mimic 下游是**单棵子树**（单链，所有 `MimicTestCases` 与典型用法都如此）时二者相等、覆盖完整；但若 mimic **直接分叉成多棵下游真子树**，只有第一棵被填，其余分支的耦合会被漏掉——属**未测试/未支持的边界**。验证法：建"mimic 直接接两条链分叉"的模型，比对 `crba` 与 $G^\top M_{\text{full}}G$。

---

## 7. WORLD 约定实现：逐行读完 `crbaWorldConvention`

同样按代码顺序走。WORLD 和 LOCAL 的**骨架完全一样**（三趟循环、同样的方向、同样的重载分派），差别只在"量表达在哪个坐标系"，以及由此带来的**一步消失、一步新增**。

### 7.1 驱动函数：多了首尾两段

[crba.hxx:508](../include/pinocchio/src/algorithm/crba.hxx#L508)：

```cpp
const MatrixXs & crbaWorldConvention(const Model & model, Data & data, const VectorXs & q)
{
  assert(model.check(data) && "data is not consistent with model.");
  PINOCCHIO_CHECK_ARGUMENT_SIZE(q.size(), model.nq, "...");

  data.oYcrb[0].setZero();                                    // ★ LOCAL 没有这句

  typedef CrbaWorldConventionForwardStep<...> Pass1;
  for (JointIndex i = 1; i < (JointIndex)(model.njoints); ++i)
    Pass1::run(model.joints[i], data.joints[i], typename Pass1::ArgsType(model, data, q));

  typedef CrbaWorldConventionBackwardStep<...> Pass2;
  for (JointIndex i = (JointIndex)(model.njoints - 1); i > 0; --i)
    Pass2::run(model.joints[i], typename Pass2::ArgsType(model, data));   // ★ 不传 jdata

  typedef CrbaWorldConventionMimicStep<...> Pass3;
  for (size_t i = 0; i < model.mimicking_joints.size(); i++)
    Pass3::run(model.joints[model.mimicking_joints[i]], typename Pass3::ArgsType(model, data, i));

  data.M.diagonal() += model.armature;

  // ★ 收尾：把 data.Ag 变成质心动量矩阵（LOCAL 没有这一段）
  data.mass[0] = data.oYcrb[0].mass();
  data.com[0]  = data.oYcrb[0].lever();
  const Block3x Ag_lin = data.Ag.template middleRows<3>(Force::LINEAR);
  Block3x       Ag_ang = data.Ag.template middleRows<3>(Force::ANGULAR);
  for (long i = 0; i < model.nv; ++i)
    Ag_ang.col(i) += Ag_lin.col(i).cross(data.com[0]);

  return data.M;
}
```

三处 ★ 标记的差异：

1. **开头 `data.oYcrb[0].setZero()`**——槽 0（universe）要当累加器用，见 §7.6。
2. **Pass 2 的 `run` 不传 `data.joints[i]`**——因为它不需要 `jdata.S()`：运动子空间在 Pass 1 就已经算好并存进 `data.J` 了。少传一个 variant，少一次 `boost::get`。
3. **结尾多了 6 行**——顺带产出质心动量矩阵，见 §7.6。

### 7.2 Pass 1：`CrbaWorldConventionForwardStep`

[crba.hxx:33](../include/pinocchio/src/algorithm/crba.hxx#L33)，比 LOCAL 多两行：

```cpp
const JointIndex & i = jmodel.id();
jmodel.calc(jdata.derived(), q.derived());                          // ① 同 LOCAL

data.liMi[i] = model.jointPlacements[i] * jdata.M();                // ② 同 LOCAL

const JointIndex & parent = model.parents[i];                       // ③ 新增：累积到世界系
if (parent > 0) data.oMi[i] = data.oMi[parent] * data.liMi[i];
else            data.oMi[i] = data.liMi[i];

jmodel.jointExtendedModelCols(data.J) = data.oMi[i].act(jdata.S()); // ④ 新增：S 搬到世界系
data.oYcrb[i] = data.oMi[i].act(model.inertias[i]);                 // ⑤ 惯量也搬到世界系
```

**③ `oMi[i] = oMi[parent] * liMi[i]`**——就是标准的前向运动学。这一步本来就有（`forwardKinematics` 也这么写），只是 LOCAL 版 CRBA 用不到所以省了。

**④ `jointExtendedModelCols(data.J) = oMi[i].act(jdata.S())`**——把运动子空间搬到世界系：${}^oX_i S_i$。

> 这里用的是 **`jointExtendedModelCols`**（`idx_vExtended` / `nvExtended`），不是 `jointCols`。
> 因为 `data.J` 按**扩展列布局**——每个关节（包括 mimic）各占一列；而 `data.M`、`data.Ag` 按**真实列布局**——mimic 与主动关节共用列。普通模型下两者完全一致，只有 mimic 模型才分叉（§5.2）。
>
> 顺带一提：写完之后 `data.J` 就是**整机雅可比**，`getJointJacobian` 可以直接取用——WORLD 约定的又一个副产物。

**⑤ `oYcrb[i] = oMi[i].act(model.inertias[i])`**——组合惯量种子，但这次要做一次 §6.3 讲的 10 参变换，搬到世界系。这是 WORLD 比 LOCAL 多花的两次 `act` 之一。

### 7.3 Pass 2（普通关节）：从四步减到三步

[crba.hxx:78](../include/pinocchio/src/algorithm/crba.hxx#L78)。和 LOCAL 的四步（§6.3）对着看：

```cpp
const JointIndex & i = jmodel.id();

// ① 力集：Ag 的本关节列 = oYcrb[i] · J 的本关节列
ColsBlock Ag_cols = jmodel.jointCols(data.Ag);                    // 真实列布局
ColsBlock J_cols  = jmodel.jointExtendedModelCols(data.J);        // 扩展列布局
motionSet::inertiaAction(data.oYcrb[i], J_cols, Ag_cols);

// ② 填 M：与 LOCAL 结构完全相同，只是数据源换成 data.Ag
data.M.block(jmodel.idx_v(), jmodel.idx_v(), jmodel.nv(), data.nvSubtree[i]).noalias() =
    J_cols.transpose() * data.Ag.middleCols(jmodel.idx_v(), data.nvSubtree[i]);

// ③ 组合惯量塌缩给父 —— 【没有任何坐标变换】
const JointIndex & parent = model.parents[i];
data.oYcrb[parent] += data.oYcrb[i];

// ④ 【整步消失】
```

#### 消失的第 ④ 步是 WORLD 的核心收益

LOCAL 需要第 ④ 步，是因为 `Fcrb[i]` 的各列写在**关节 $i$ 自己的坐标系**里，父关节要用就必须先搬过去，每层搬一次。

WORLD 里，`data.Ag` 是**一个全局共享的 $6\times n_v$ 矩阵**（不是每关节一个），每一列从写下去那一刻起就在**世界系**——所有关节看到的是同一个系，**不需要搬**。所以第 ② 步直接 `data.Ag.middleCols(idx_v, nvSubtree[i])` 取子树列就行。

同理第 ③ 步的 `oYcrb[parent] += oYcrb[i]` **连 `.act()` 都不用**，纯 10 参加法。

> **代价量化**（30 关节浮动基人形，$n_v=34$）：LOCAL 的第 ④ 步累计要搬 $\sum_{i:\lambda(i)>0}\texttt{nvSubtree}[i] = \mathbf{118}$ 列；
> WORLD 全部省掉，换成 Pass 1 里 $29$ 次额外的 `act`。实测净结果 WORLD 反而略慢 5%（§8.1）——
> **省下的没有多花的多**。WORLD 的真正价值在副产物（`data.J`、`data.Ag`、`com[0]`），不在速度。

#### `motionSet::inertiaAction` 的底层

[act-on-set.hxx:646](../include/pinocchio/src/spatial/act-on-set.hxx#L646)，同样是**逐列**：

```cpp
for (int col = 0; col < jV.cols(); ++col)
  motionSet::inertiaAction<Op>(I, iV.col(col), jV.col(col));
```

单列落到 `InertiaTpl::__mult__`（[inertia.hxx:890](../include/pinocchio/src/spatial/inertia.hxx#L890)）——同样不组装 $6\times6$：

```cpp
f.linear().noalias() = mass() * (v.linear() - lever().cross(v.angular()));
Symmetric3::rhsMult(inertia(), v.angular(), f.angular());     // τ = Ī_C · ω
f.angular() += lever().cross(f.linear());                     // τ += c × f
```

即用 $f = m(v - c\times\omega)$、$\tau = \bar I_C\omega + c\times f$ 代替

$$
\begin{bmatrix} f\\ \tau\end{bmatrix}
=\begin{bmatrix} m\mathbb 1 & -m\hat c\\[2pt] m\hat c & \bar I_C - m\hat c\hat c\end{bmatrix}
\begin{bmatrix} v\\ \omega\end{bmatrix}
$$

展开 $c\times f = m\hat c v - m\hat c\hat c\,\omega$ 可见两者逐项相等，但乘法从 36 次降到约 12 次 + 2 次叉积，且不需要 $6\times6$ 临时内存。

`Op` 同样是编译期参数：这里默认 `SETTO`；Pass 3 给 mimic 往 `Ag` 补贡献时用 `ADDTO`。

### 7.4 Pass 2（mimic 关节）：只剩一行

[crba.hxx:99](../include/pinocchio/src/algorithm/crba.hxx#L99)。比 LOCAL 的 mimic 重载还要简单——整个函数体就一句：

```cpp
static void algo_impl(const JointModelBase<JointModelMimic> & jmodel,
                      const Model & model, Data & data)
{
  const JointIndex & i = jmodel.id();
  const JointIndex & parent = model.parents[i];
  data.oYcrb[parent] += data.oYcrb[i];      // 就这一句
}
```

对照普通版：第 ① 步（写 `Ag` 列）、第 ② 步（写 `M`）**全部跳过**，只保留惯量累加。

LOCAL 的 mimic 重载还得处理力集上搬（`ADDTO` 整块搬），WORLD 这里连这个都不需要——因为根本没有"力集上搬"这一步。**这是 WORLD 约定在 mimic 处理上的额外简化。**

### 7.5 Pass 3：`CrbaWorldConventionMimicStep`

[crba.hxx:148](../include/pinocchio/src/algorithm/crba.hxx#L148)。结构与 LOCAL 版（§6.6）一一对应，但每一步都因为"全在世界系"而变简单：

```cpp
JointIndex mimicking_id = jmodel.id();

// ① 算 mimic 的力集（临时变量，不占用 data.Ag）
auto J_cols = jmodel.jointExtendedModelCols(data.J);          // 世界系 S（带耦合倍率）
Matrix<Scalar,6,Dynamic,...> Ag_sec(6, jmodel.nvExtended());
motionSet::inertiaAction(data.oYcrb[mimicking_id], J_cols, Ag_sec);

// ② 自身贡献 → 主动关节对角块
data.M.block(jmodel.idx_v(), jmodel.idx_v(), jmodel.nvExtended(), jmodel.nvExtended()).noalias()
    += J_cols.transpose() * Ag_sec;

// ③ 与下游有独立 DoF 的子树的耦合
JointIndex sub_mimic_id = data.mimic_subtree_joint[mims_id];
if (sub_mimic_id != 0)
  jmodel.jointRows(data.M).middleCols(model.idx_vs[sub_mimic_id], data.nvSubtree[sub_mimic_id])
      .noalias() += J_cols.transpose()
                  * data.Ag.middleCols(model.idx_vs[sub_mimic_id], data.nvSubtree[sub_mimic_id]);

// ④ 沿支撑链往根
for (size_t j = supports.size() - 2; j > 0; j--)
{
  int sup_idx_v         = model.idx_vs[supports[j]];
  int sup_idx_vExtended = model.idx_vExtendeds[supports[j]];
  int sup_nvExtended    = model.nvExtendeds[supports[j]];

  temp_JAG.noalias() = data.J.middleCols(sup_idx_vExtended, sup_nvExtended).transpose() * Ag_sec;

  if (sup_idx_v >= mims_idx_v) { ... += temp_JAG;  if (相等) 再加一次; }
  else                         { ... += temp_JAG; }     // 行列翻转，同 §6.6(d)
}

// ⑤ mimic 对质心动量矩阵的贡献，补回主动关节那列
motionSet::inertiaAction<ADDTO>(data.oYcrb[mimicking_id], J_cols, jmodel.jointCols(data.Ag));
```

**和 LOCAL 版 Pass 3 的三处关键区别：**

| | LOCAL（§6.6） | WORLD |
|---|---|---|
| 力集存哪 | 借用 `data.Fcrb[mimicking_id]` | **临时变量 `Ag_sec`**，不污染 `data.Ag` |
| 沿支撑链上溯 | 每级做一次 `forceSet::se3Action`**原地变换** `F` | **不变换**，直接读 `data.J` 对应祖先的列 |
| 取祖先子空间 | 需要 `applyConstraintOnForceVisitor` 运行期分派 | 直接 `data.J.middleCols(...)`，**无需访问者** |

第二、三行就是 WORLD 的威力：**每个关节的世界系子空间在 Pass 1 就都算好躺在 `data.J` 里了**，所以循环里既不用逐级变换力集，也不用为了拿 `S_li` 再套一层 variant 分派——一次 `middleCols` 取出来转置相乘即可。

**第 ⑤ 步是 LOCAL 完全没有的**：mimic 也会影响质心动量，得把 $Y^c_{\text{mim}}\cdot S_{\text{mim}}$ 用 `ADDTO` 累加回主动关节在 `data.Ag` 里的那一列，否则 §7.6 算出的 $A_g$ 会漏掉 mimic 的贡献。

### 7.6 收尾那 6 行：质心动量矩阵是怎么"顺带"算出来的

```cpp
data.mass[0] = data.oYcrb[0].mass();
data.com[0]  = data.oYcrb[0].lever();

const Block3x Ag_lin = data.Ag.template middleRows<3>(Force::LINEAR);
Block3x       Ag_ang = data.Ag.template middleRows<3>(Force::ANGULAR);
for (long i = 0; i < model.nv; ++i)
  Ag_ang.col(i) += Ag_lin.col(i).cross(data.com[0]);
```

**第一问：`oYcrb[0]` 为什么是整机惯量？**

回看 Pass 2 的第 ③ 步 `data.oYcrb[parent] += data.oYcrb[i]`——它对**每个** $i\ge1$ 都执行，包括 `parent == 0` 的情形（普通关节版本里没有 `if (parent > 0)` 的保护，这点和 LOCAL 不同）。所以循环跑到 $i=1$ 结束时，整棵树都累加进了槽 0：

$$ Y^c_0=\sum_{j=1}^{n}{}^oX_j^{*}\,I_j\,{}^jX_o $$

这也解释了驱动函数开头那句 `data.oYcrb[0].setZero()`——**槽 0 是累加器，每次调用前必须清零**，否则会把上次的结果叠上去。

于是 `.mass()` 直接给出总质量，`.lever()` 直接给出**世界系下的整机质心** $c$。

> LOCAL 约定没有这一步，因为各 `Ycrb[i]` 表达在各自体坐标系，**不能直接相加**。

**第二问：为什么只改 angular 行？**

Pass 2 写进 `data.Ag` 的列是 $Y^c_i\,{}^oX_iS_i$——一个**以世界原点 $O$ 为参考点**的空间动量。质心动量矩阵要求参考点在**质心 $G$**。两者差一个**纯平移的余伴随**：

$$
{}^GX_O^{*}=\begin{bmatrix}\mathbb 1 & 0\\ -\hat c & \mathbb 1\end{bmatrix}
\qquad\Longrightarrow\qquad
\begin{cases} f_G = f_O & \text{线性分量不变}\\[2pt] \tau_G = \tau_O - c\times f_O\end{cases}
$$

线性分量不变 → 代码只碰 `Ag_ang`。而代码写的是 `+= f × c`，注意

$$ f\times c=-\,c\times f $$

与公式一致。**这就是那一行 `cross` 的全部来历**：把动量的参考点从世界原点搬到质心。

**第三问：产物是什么？**

`data.Ag` 至此就是**质心动量矩阵** $A_g$，满足 $h_G=A_g\dot q$，且前三行满足 $A_g^{\text{lin}}\dot q = m\,v_{\text{com}}$。

> **实测**（30 关节浮动基人形）：
> `||Ag(crba WORLD) − Ag(ccrba)|| = 0.000e+00`、`||com − com(ccrba)|| = 0.000e+00`——**逐位相同**。
> 既要 $M$ 又要 $A_g$ 时，一次 `crba(..., Convention::WORLD)` 就够，不必再调 `ccrba`。

---

## 8. 两种约定对比：差异到底在哪

把 §6、§7 逐行读到的差异汇总成一张表：

| | LOCAL（§6） | WORLD（§7） |
|---|---|---|
| 驱动函数开头 | 无 | `data.oYcrb[0].setZero()`（槽 0 当累加器） |
| Pass 2 的 `run` | 传 `(jmodel, jdata, args)` | **不传 `jdata`**（子空间已在 `data.J`） |
| Pass 1 惯量种子 | `Ycrb[i] = inertias[i]`（纯拷贝） | `oYcrb[i] = oMi.act(inertias[i])`（多一次变换） |
| Pass 1 子空间 | 只放在 `jdata.S()`（本系） | 额外写进 `data.J`（世界系，多一次变换） |
| Pass 2 力集容器 | `Fcrb[i]`——**每关节一个** $6\times n_v$ | `data.Ag`——**全局共享一个** $6\times n_v$ |
| Pass 2 步数 | **四步**（含力集逐级上搬） | **三步**（第 ④ 步整个消失） |
| Pass 2 惯量累加 | `Ycrb[父] += liMi.act(Ycrb[i])`，且有 `if (parent>0)` | `oYcrb[父] += oYcrb[i]`，**无变换、无保护** |
| Pass 2 的 mimic 重载 | 保留惯量累加 + `ADDTO` 整块搬力集 | **只保留惯量累加一行** |
| Pass 3 沿支撑链 | 逐级 `forceSet::se3Action` 原地变换 + 访问者取 `S_li` | 直接读 `data.J` 的祖先列，**两者都不需要** |
| Pass 3 步数 | 四步 | **五步**（多一步把 mimic 贡献补回 `data.Ag`） |
| 收尾 | 只加 armature | 加 armature + **产出 $A_g$、`com[0]`、`mass[0]`** |

**该怎么选？**

- **只要 $M$**：选 LOCAL。实测它反而**略快 5%**（§8.1）——WORLD 在 Pass 1 多花的两次 `act`，没有被 Pass 2 省下的 118 次列变换赚回来。而且 LOCAL 更贴 Featherstone 教材，读代码更直观。
- **还要 $A_g$ / 质心 / 雅可比**：选 WORLD。它把 `data.J`（整机雅可比）、`data.Ag`（质心动量矩阵）、`com[0]`、`mass[0]` 一并算出来，且与 `ccrba` 的结果**逐位相同**。Pinocchio 上层很多算法默认走 WORLD，正是为了复用这些中间量。
- **内存敏感**：选 WORLD。`Fcrb` 是 $O(n\cdot n_v)$ 的（每关节一个 $6\times n_v$），实测 47.8 KB vs WORLD 的 3.2 KB，差一个数量级。

---

### 8.1 实测数据：耗时、内存与一致性

> 测量条件：30 关节浮动基人形（$n_q=35$、$n_v=34$），`-O3 -DNDEBUG`，N=2000 次随机位形取平均。

**耗时**

| 调用 | 单次耗时 |
|---|---|
| `crba(..., Convention::LOCAL)` | **2.09 µs** |
| `crba(..., Convention::WORLD)` | **2.20 µs** |

WORLD **略慢约 5%**——Pass 1 多的两次 `act`（§7.4）没有被 Pass 2 省下的列变换完全抵消。
所以"WORLD 更快"是个误解：**如果只要 $M$，LOCAL 反而略优**；WORLD 的价值在于**顺带产出 $A_g$ 与 `com[0]`**。

**中间量内存**

| 约定 | 中间量 | 大小 |
|---|---|---|
| LOCAL | `Fcrb`：`njoints × 6 × nv × 8B` | **47.8 KB** |
| WORLD | `J` + `Ag`：`2 × 6 × nv × 8B` | **3.2 KB** |
| 共用 | `Ycrb`/`oYcrb`：`njoints × 10 × 8B` | 2.3 KB |

`Fcrb` 是 $O(n\cdot n_v)=O(n^2)$ 的，[data.hxx:713](../include/pinocchio/src/multibody/data.hxx#L713) 里为**每个关节**都分配了一整个 $6\times n_v$ 矩阵。WORLD 约定只需两个 $6\times n_v$，**内存少一个数量级**。

**一致性**

$$ \lVert M_{\text{LOCAL}} - M_{\text{WORLD}}\rVert = 6.27\times10^{-14} $$

> 另有一个未深究的观察：同一环境下单独调 `ccrba` 实测约 **610 µs**（可稳定复现，比 `crba(WORLD)` 慢约 280×）。
> 因此在本构建里，需要 $A_g$ 时用 `crba(..., WORLD)` 明显更划算。这个差距的成因我没有追查，暂记于此。

---

## 9. 上三角、`nvSubtree` 与紧凑树要求（`CRBAChecker`）

### 9.1 只填上三角——精确说法

$M=M^\top$，两种约定都只沿"对角块向右"填：

```cpp
data.M.block(jmodel.idx_v(), jmodel.idx_v(), jmodel.nv(), data.nvSubtree[i])
```

**但"只填上三角"这句话需要限定。** 实测（浮动基人形）：

```
||M 的严格下三角||               = 1.442e+01     ← 不是 0！
  其中 落在关节自身 nv×nv 对角块内 = 1.442e+01
  其中 跨关节的严格下三角          = 0.000e+00    ← 精确为 0
```

原因很直白：上面那个 `block` 的**起点在对角线上、高度是 `nv`**，所以对 $n_v>1$ 的关节（如 free-flyer 的 $6\times6$ 块），它会把**自己那个对角方块整块写满**，自然包含该块的下三角。

于是准确的表述是：

> **关节与关节之间只填上三角；每个关节自身的 $n_v\times n_v$ 对角块是完整（对称）填好的。**

对使用者的影响：`M.selfadjointView<Eigen::Upper>()` 或转置镜像**依然是正确且必需的**做法——对角块被重复写一遍不影响结果（它本来就对称）。

### 9.2 `nvSubtree[i]`：子树是连续列区间

已在 [§5.3](#53-三趟-pass-共享的机制之三nvsubtree-决定往-m-的哪块写) 完整讲过（含构造代码、连续区间假设、实测写入图与工作量统计），此处不再重复。
一句话：`data.M.block(idx_v, idx_v, nv, nvSubtree[i])` 能成立，全靠"子树的速度索引是一段连续区间"，而这由下面 §9.3 的紧凑树要求保证。

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

## 11. mimic 关节：为什么要拆成"截断 + 补账"

§6.5/§6.6（LOCAL）与 §7.4/§7.5（WORLD）已经把代码逐行读完了。本节补上**数学上的理由**——为什么非得这么拆。

### 11.1 mimic 关节在索引上是什么

`JointModelMimic` 的三个关键性质（[joint-mimic.hxx:615](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L615)）：

```cpp
inline int nq_impl()        const { return 0; }              // 不贡献位置维度
inline int nv_impl()        const { return 0; }              // 不贡献自由度
inline int nvExtended_impl()const { return m_nvExtended; }   // 但在扩展空间仍占一列
```

而 `setIndexes_impl` **故意把传进来的 q、v 参数注释掉不用**，另设 `setMimicIndexes` 让 `idx_q` / `idx_v` 指向**被模仿的主动关节**。源码注释写得很直白：

> *idx_q and idx_v should remain pointing to the mimicked joint.*

所以一句话：**mimic 关节的 `idx_v` 是主动关节的列号——它和主动关节共用 `data.M` 的同一行同一列。**

> **实测**（12 个 mimic 测试用例、共 20 个 mimic 关节）：
> `nv()==0` 命中 **20/20**；`idx_v` 落在某个 `nv>0` 关节的列块内 **20/20**——列别名 100% 成立。

### 11.2 数学上：mimic 把 M 变成一个合同变换

设扩展速度与真实速度的映射为 $\dot q_{\text{ext}} = G\,\dot q$（$G\in\mathbb R^{n_{v\text{Ext}}\times n_v}$，mimic 关节按传动比 $r$ 跟随主动关节）。动能不变给出

$$\boxed{\;M \;=\; G^\top M_{\text{ext}}\, G\;}$$

对某个被共享的自由度（主动 $p$、从动 $s$、传动比 $r$），$G$ 的这一列有两个非零元（$1$ 和 $r$），展开得

$$
M_{aa}
\;=\;\underbrace{M_{pp}}_{\text{Pass 2 已算}}
\;+\;\underbrace{2r\,M_{ps}}_{\text{Pass 3 的 ④}}
\;+\;\underbrace{r^2 M_{ss}}_{\text{Pass 3 的 ②}}
$$

**这三项正好对应代码里的三处：**

- $M_{pp}$——Pass 2 处理主动关节 $p$ 时正常写入。
- $r^2M_{ss}$——Pass 3 第 ② 步 `jointBlock(data.M) += S^T·F`。传动比 $r$ 已经被吸收进 `ScaledJointMotionSubspace`（[joint-mimic.hxx:183](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L183) 的 `S = m_scaling_factor * m_constraint.matrix_impl()`），左右各一个 $S$ 于是自动得到 $r^2$。
- $2r\,M_{ps}$——Pass 3 第 ④ 步支撑链循环里那段"**看起来很怪**"的代码：

  ```cpp
  if (mims_idx_v == sup_idx_v)
    data.M.block(...) += temp_JAG;      // 同一个 temp_JAG 再加一遍
  ```

  交叉项 $M_{ps}$ 与 $M_{sp}$ 对称，折叠到同一个上三角格子里就是**加两次**，因子 2 由此而来。

#### 多个 mimic 跟随同一主动关节

若主动关节 $p$ 被 $m$ 个 mimic 同时跟随（$s_1,\dots,s_m$，传动比 $r_1,\dots,r_m$），$G$ 主动那一列变成**所有 mimic 的和**：

$$g_p = e_p + \sum_{k=1}^{m} r_k\,e_{s_k}$$

代进二次型 $M_{pp}^{\text{red}}=g_p^\top M_{\text{ext}}\,g_p$，用 $e_a^\top M e_b=M_{ab}$ 与对称性展开，比单 mimic **多出一类项**：

$$
\boxed{\,M_{pp}^{\text{red}}
=\underbrace{M_{pp}}_{\text{主动自己}}
+\underbrace{2\sum_k r_k M_{p,s_k}}_{\text{主动–mimic 交叉}}
+\underbrace{\sum_k r_k^2 M_{s_k s_k}}_{\text{各 mimic 反射}}
+\underbrace{2\!\sum_{k<l} r_k r_l\,M_{s_k s_l}}_{\textbf{mimic–mimic 交叉（新增）}}\,}
$$

**新增的 $2\sum_{k<l}r_kr_l M_{s_ks_l}$** 是两个 mimic 通过共享主动 DoF 产生的惯性交叉耦合。由分支诱导稀疏（§9.4），$M_{s_ks_l}\ne0$ **仅当一个 mimic 落在另一个的子树里（嵌套）**；两个 mimic 分处不同分支时该项为 0，退化为"各自独立反射"。

**代码如何无特判地实现**：Pass 3 遍历 `model.mimicking_joints`（**所有** mimic），逐个处理、全程 `+=`；因每个 mimic 的 `idx_v` 都指向同一主动槽，贡献自然汇总。所有交叉项的**因子 2** 都由那句双加产生——注释说得很准：

```cpp
// check if support is either a mimic with the same primary or jmodel primary
if (mims_idx_v == sup_idx_v)     // 支撑链祖先的列号 == 主动槽
  data.M.block(...) += temp_JAG; // 再加一遍
```

处理 $s_k$ 沿支撑链上行时，祖先 `sup_idx_v == 主动槽` 恰有两种情形，正好对上两类因子 2：

| 支撑链祖先是 | 贡献项 | 因子 2 来源 |
|---|---|---|
| **主动关节 $p$ 本身** | $2r_k M_{p,s_k}$ | $p$ 的 `idx_v` = 主动槽 |
| **另一个同主动的 mimic $s_l$** | $2r_k r_l M_{s_k s_l}$ | $s_l$ 的 `idx_v` 也 = 主动槽 |

即注释里"support 是**同主动的 mimic** 或**主动关节本身**就加两次"正是为多 mimic 而写：它同时补上了两类交叉项的对称第二份。§11.4 的 case 3/case 7（双 mimic）正是对此的数值验证。

### 11.3 于是"截断"是必然的

如果 Pass 2 对 mimic 照常执行

```cpp
data.M.block(jmodel.idx_v(), jmodel.idx_v(), jmodel.nv(), data.nvSubtree[i]) = ...;
//                ↑ 其实是主动关节的列号          ↑ 而且是 "=" 覆盖赋值
```

它会拿 mimic 的数据**覆盖掉主动关节刚算好的 $M_{pp}$**——而正确答案是三项**相加**。所以只能：

1. **Pass 2 跳过写 M**（截断），保留主动关节自己算的 $M_{pp}$；
2. **Pass 3 用 `+=` 把 $r^2M_{ss}$ 和 $2rM_{ps}$ 补上去**。

同样的道理，LOCAL 的 Pass 2 里力集上搬要从 `SETTO` 改成 `ADDTO`（§6.5）——`Fcrb` 的列也是别名的，覆盖会抹掉并联支链已经传上来的力。

### 11.4 数值验证

拿 `MimicTestCases` 的三个用例，对比 `crba(model_mimic)` 与 $G^\top M_{\text{full}}G$：

| 用例 | 说明 | $\lVert M_{\text{LOCAL}}-G^\top M_{\text{full}}G\rVert$ | $\lVert M_{\text{WORLD}}-G^\top M_{\text{full}}G\rVert$ |
|---|---|---|---|
| case 0 | 单 mimic，父子相邻 | $1.64\times10^{-14}$ | $1.82\times10^{-13}$ |
| case 3 | 双 mimic，跨左右支 | $7.47\times10^{-14}$ | $1.50\times10^{-13}$ |
| case 7 | 串联双 mimic，末端 | $4.50\times10^{-14}$ | $1.88\times10^{-13}$ |

两种约定也互相吻合（$\sim10^{-13}$）。**合同变换关系在数值上成立**，说明"截断 + 补账"这套拆法是对的。

### 11.5 一张对照表

| | 普通关节 | mimic 关节 |
|---|---|---|
| `nv()` | $\ge 1$ | **0** |
| `nvExtended()` | $=n_v$ | $\ge1$（`data.J` 里独占一列） |
| `idx_v` | 自己的列 | **主动关节的列（别名）** |
| `jointCols` / `jointBlock` / `jointRows` | 基类版 `(idx_v, nv)` | **重写为** `(i_v, nvExtended)` |
| Pass 2 写 `data.M` | ✅ `=` 覆盖 | ❌ **完全跳过** |
| Pass 2 力集上搬（仅 LOCAL） | `SETTO`，切子树列 | **`ADDTO`**，整块搬 |
| Pass 2 累加 `Ycrb`/`oYcrb` | ✅ | ✅ |
| Pass 3 | 空函数体 | 四（LOCAL）/ 五（WORLD）步补账，全用 `+=` |

**一句话记忆**：mimic 关节在 Pass 2 里"**只传惯量、不写矩阵**"，因为它的 `idx_v` 会踩到主动关节的格子；账留到 Pass 3 用累加方式补上。

---

## 12. 与 $A_g$、ABA、Cholesky 的关系

- **$A_g$（质心动量矩阵）**：WORLD 约定 Pass 2 的力集存进 **`data.Ag`** 并非巧合——$Y_i^cS_i$ 既是填 $M$ 的中间力集，其在质心系的表达正是 CMM 的列（$h_g=A_g\dot q$）。所以 CRBA(WORLD) 与 [`ccrba`](质心与质心动量解析.md) 共享同一套组合惯量与力集，$A_g$ 是 CRBA 的免费副产物。
- **ABA**：都用组合惯量作种子，但 ABA 用 Schur 补塌缩（铰接体），CRBA 用直接相加（焊死）。见 [正动力学ABA解析.md](正动力学ABA解析.md)。
- **Cholesky**：CRBA 产出的 $M$（上三角 + 分支稀疏）正好喂给 `cholesky::decompose` 做 $M=UDU^\top$，稀疏结构一致（都源于支撑链），配合 $O(nd)$ 解 $Mx=y$。见 [multibody §3.5](multibody子系统解析.md)。
- **正动力学**：$\ddot q=M^{-1}(\tau-b)$ 里，CRBA 求 $M$、RNEA 求 $b$、Cholesky 求解——这是 ABA 之外的另一条路。

---

## 13. 复杂度

- **CRBA** $O(nd)$：$d$ 为树平均深度。链式结构 $d\approx n$ → $O(n^2)$；分支多则 $d$ 小。
  Pass 2 写入的 $M$ 元素总数恰为 $\sum_i \texttt{nvSubtree}[i]$——实测 30 关节浮动基人形（$n_v=34$）为 **152**，
  而稠密上三角有 $n_v(n_v{+}1)/2=595$ 个元素，**只碰了 26%**，其余是结构零（§5.3）。
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
- **关节之间只填上三角**（关节自身的 $n_v\times n_v$ 对角块是完整的，§9.1）、只碰**祖先-后代对**
  （要求紧凑树 `CRBAChecker`），天然给出分支诱导稀疏，直通 Cholesky 与质心动量矩阵；
- 末尾补 `armature`；mimic 关节在 Pass 2 被**截断**（只传惯量、不写矩阵），由 Pass 3 用 `+=` 补账（§11）。
