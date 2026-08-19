# Pinocchio `multibody/` 子系统解析

> 本文件系统解析 `include/pinocchio/multibody/`（声明）与 `include/pinocchio/src/multibody/`（模板实现）下**全部文件**，
> 讲清每个文件的作用、内部实现与涉及的数学。配套阅读：
> [源码解析.md](源码解析.md)（模板机制/访问者三件套）、
> [空间代数运算解析.md](空间代数运算解析.md)（SE3/Motion/Force/Inertia/exp-log 的数学）、
> [类结构图.md](类结构图.md)（UML）、[C++语法技巧.md](C++语法技巧.md)、
> [逆运动学IK解析.md](逆运动学IK解析.md)（IK 例子逐行数学：误差/Jlog6/阻尼最小二乘/$n_q\neq n_v$）。

---

## 目录

- [0. 这个子系统在 Pinocchio 里的位置](#0-这个子系统在-pinocchio-里的位置)
- [1. 文件地图（53 个实现文件）](#1-文件地图53-个实现文件)
- [2. `ModelTpl`：运动学树的静态定义](#2-modeltpl运动学树的静态定义)
- [3. `DataTpl`：算法工作区](#3-datatpl算法工作区)
- [4. `FrameTpl` / `ModelItem`：坐标系与树节点基类](#4-frametpl--modelitem坐标系与树节点基类)
- [5. 关节系统总览：Model/Data 分离 + CRTP + variant 双层](#5-关节系统总览modeldata-分离--crtp--variant-双层)
- [6. `JointMotionSubspace`（运动子空间 S）](#6-jointmotionsubspace运动子空间-s)
- [7. `calc` / `calc_aba`：关节的核心计算语义](#7-calc--calc_aba关节的核心计算语义)
- [8. 全关节清单（NQ/NV/流形/S）](#8-全关节清单nqnv流形s)
- [9. 复合关节：Composite / Mimic / Unaligned / Unbounded](#9-复合关节composite--mimic--unaligned--unbounded)
- [10. 李群系统：位形流形的加减法](#10-李群系统位形流形的加减法)
- [11. 访问者系统：把 variant 变成可调用的算法](#11-访问者系统把-variant-变成可调用的算法)
- [12. 其余文件：pool / sample-models / force-set](#12-其余文件pool--sample-models--force-set)
- [13. 一次 `forwardKinematics` 如何串起整个子系统](#13-一次-forwardkinematics-如何串起整个子系统)
- [14. 速查表](#14-速查表)

---

## 0. 这个子系统在 Pinocchio 里的位置

Pinocchio 分三大层：

1. **`spatial/`**：单个刚体的空间代数（SE3、Motion、Force、Inertia、exp/log）——见 [空间代数运算解析.md](空间代数运算解析.md)。
2. **`multibody/`（本文）**：把很多刚体用关节连成一棵**运动学树**，并定义树的**静态结构**（`Model`）、**计算缓存**（`Data`）、**关节**（`Joint*`）、**位形流形**（`LieGroup`）、**类型擦除分派**（`visitor`）。
3. **`algorithm/`**：在 `Model`+`Data` 上跑 RNEA/CRBA/ABA/雅可比/碰撞等具体算法。

一句话：**`multibody/` 是"机器人是什么"的数据结构层，`algorithm/` 是"机器人怎么算"的算法层**。所有算法都以 `const Model &` + `Data &` 为输入输出。

---

## 1. 文件地图（53 个实现文件）

> 三文件结构：`multibody/xxx.hpp`（对外声明/聚合）→ `src/multibody/xxx.hxx`（模板实现，`#pragma once` + IWYU private）→ `src/*.cpp`（显式实例化）。下面按子系统列出 `.hxx`。

| 子系统 | 文件 | 作用 |
|--------|------|------|
| **前置声明** | `fwd.hxx` | `ModelTpl`/`DataTpl`/`FrameTpl` 前置声明、`Index` 系列 typedef、`ReferenceFrame`(WORLD/LOCAL/LWA)、`KinematicLevel`、`Convention` 枚举 |
| **核心结构** | `model.hxx` | `ModelTpl`：运动学树的静态定义 + `addJoint`/`addFrame`/`cast`/`createData` |
| | `data.hxx` | `DataTpl`：所有算法的中间量缓存（~150 个字段） |
| | `frame.hxx` | `FrameTpl` + `FrameType` 枚举（OP_FRAME/JOINT/FIXED_JOINT/BODY/SENSOR） |
| | `model-item.hxx` | `ModelItem`：`Frame` 等树节点的公共基类（name/parentJoint/parentFrame/placement） |
| | `force-set.hxx` | 力集合/运动集合的批量列操作辅助 |
| **关节基类** | `joint/joint-model-base.hxx` | `JointModelBase<Derived>`（CRTP）+ 一大套 `PINOCCHIO_JOINT_*` 宏、索引访问器、段/列/块选择器 |
| | `joint/joint-data-base.hxx` | `JointDataBase<Derived>`（CRTP）+ 访问器宏 |
| | `joint/joint-collection.hxx` | `JointCollectionDefaultTpl`：关节"菜单" → `JointModelVariant`/`JointDataVariant` |
| | `joint/joint-generic.hxx` | `JointModelTpl`/`JointDataTpl`：变体包装器（同时继承 CRTP 基类和 variant） |
| | `joint/joint-basic-visitors.hxx` | `calc`/`calc_aba`/`nv`/`nq`/`cast` 等对 variant 的自由函数入口 |
| | `joint/joint-common-operations.hxx` | 关节共用的小工具（如默认 `calc_aba` 内积） |
| | `joint/fwd.hxx` | 关节前置声明、`MAX_JOINT_NV=6` |
| **运动子空间** | `joint-motion-subspace-base.hxx` | `JointMotionSubspaceBase<Derived>`（CRTP，约束矩阵 S 的接口） |
| | `joint-motion-subspace-generic.hxx` | `JointMotionSubspaceTpl`：稠密 6×Dim 的通用 S 实现 |
| **具体关节** | `joint/joint-revolute.hxx` | 定轴转动 RX/RY/RZ（原型） |
| | `joint/joint-revolute-unaligned.hxx` | 任意轴转动 |
| | `joint/joint-revolute-unbounded{,-unaligned}.hxx` | 无界转动（用 cos/sin 存 q，SO(2) 流形） |
| | `joint/joint-prismatic{,-unaligned}.hxx` | 移动副 PX/PY/PZ / 任意轴 |
| | `joint/joint-helical{,-unaligned}.hxx` | 螺旋副（转+移耦合） |
| | `joint/joint-spherical.hxx` | 球副（四元数，SO(3)） |
| | `joint/joint-spherical-ZYX.hxx` | 球副（ZYX 欧拉角参数化） |
| | `joint/joint-translation.hxx` | 3D 平移 |
| | `joint/joint-planar.hxx` | 平面副（x,y,θ，SE(2)） |
| | `joint/joint-free-flyer.hxx` | 自由浮动基（SE(3)，浮动基机器人根关节） |
| | `joint/joint-ellipsoid.hxx` | 椭球约束关节 |
| | `joint/joint-universal.hxx` | 万向节（2 DoF） |
| | `joint/joint-composite.hxx` | 复合关节：把多个关节串成一个 |
| | `joint/joint-mimic.hxx` | 镜像关节：q 线性绑定到另一个关节 |
| **李群** | `liegroup/liegroup-base.hxx` | `LieGroupBase<Derived>`（CRTP）：integrate/difference/雅可比等 API |
| | `liegroup/vector-space.hxx` | ℝⁿ（欧氏，默认流形） |
| | `liegroup/special-orthogonal.hxx` | SO(2)/SO(3) |
| | `liegroup/special-euclidean.hxx` | SE(2)/SE(3) |
| | `liegroup/cartesian-product{,-variant,-variant-fwd}.hxx` | 群的笛卡尔积（整机位形空间由此拼成） |
| | `liegroup/liegroup-map.hxx` | **关节类型 → 流形**的映射表 |
| | `liegroup/liegroup-generic.hxx` | 变体化的 `LieGroupGenericTpl` |
| | `liegroup/liegroup-collection.hxx` | 李群"菜单" → variant |
| | `liegroup/liegroup-variant-visitors.hxx` | 对李群 variant 的访问者入口 |
| | `liegroup/liegroup-algo.hxx` | 沿整棵树逐关节应用李群操作（integrate/difference 整机版） |
| | `liegroup/liegroup-joint.hxx` | 关节与其李群的桥接 |
| | `liegroup/fwd.hxx` | 前置声明 |
| **访问者** | `visitor/fusion.hxx` | `bf::append`：把参数包压进 Boost.Fusion 序列 |
| | `visitor/joint-unary-visitor.hxx` | 一元访问者（对单个关节 variant 分派） |
| | `visitor/joint-binary-visitor.hxx` | 二元访问者（对两个关节 variant 分派，如惯量复合） |
| **并行/采样** | `pool/model.hxx`、`pool/geometry.hxx`、`pool/fwd.hxx` | 多线程并行用的 (Model,Data) 池 |
| | `sample-models.hxx` | 生成测试模型（humanoid/manipulator） |

---

## 2. `ModelTpl`：运动学树的静态定义

`src/multibody/model.hxx`。`Model = ModelTpl<double,0,JointCollectionDefaultTpl>`。它是**只读的机器人定义**——一旦 `addJoint`/`addFrame` 完成，算法阶段不再改动它。

### 2.1 类型与继承

```cpp
struct ModelTpl
: serialization::Serializable<ModelTpl>   // 序列化
, NumericalBase<ModelTpl>                 // 提供 Scalar typedef
, ModelEntity<ModelTpl>                   // CRTP 实体标记
```

`traits<ModelTpl>` 暴露 `Scalar`/`Options`/`Data`/`JointCollection`。注意 **`Data` 类型是从 `Model` 反查出来的**（`traits<ModelTpl>::Data`），所以 `Model` 与 `Data` 严格配对。

### 2.2 三套维度与索引体系（关键）

这是理解 Pinocchio 数据布局的钥匙。每个关节在几个不同维度的**全局向量**里各占一段：

| 空间 | 总维度 | 每关节段长 / 起点 | 含义 |
|------|--------|-------------------|------|
| **位形空间** $q$ | `nq` | `nqs[i]` / `idx_qs[i]` | 位形向量维度。可 > nv（球副 nq=4、nv=3） |
| **速度/切空间** $v$ | `nv` | `nvs[i]` / `idx_vs[i]` | 广义速度、力矩、加速度维度 |
| **扩展速度空间** | `nvExtended` | `nvExtendeds[i]` / `idx_vExtendeds[i]` | 雅可比列空间，**展开 mimic 后**的速度维度 |
| 树节点计数 | `njoints` / `nbodies` / `nframes` | — | 关节数（含 universe）/ 刚体数 / 坐标系数 |

声明处：`nq/nv/nvExtended` 在 `model.hxx:93-100`，三对索引表在 `model.hxx:120-136`。

#### 2.2.1 `nqs[i]` / `idx_qs[i]` 的读法

这两张表就是**关节 id → 位形向量 $q$ 中的切片**的查找表：

- `nqs[i]`：第 $i$ 个关节在位形空间里占几维，即 $\dim\mathcal{Q}_i$
- `idx_qs[i]`：这一段在全局 $q$ 里的起始下标

$$q_i = q\big[\;\texttt{idx\_qs}[i]\;:\;\texttt{idx\_qs}[i]+\texttt{nqs}[i]\;\big]\in\mathbb{R}^{\texttt{nqs}[i]}$$

**源码里从不手写这个 segment**，而是用 `joint-model-base.hxx:330` 的 `jointConfigSelector`，
它内部就是 `segment(idx_q(), nq())`；速度侧对应 `jointVelocitySelector`（`segment(idx_v(), nv())`）。
矩阵版还有 `jointCols` / `jointRows` / `jointBlock` 三组，语义同理。

#### 2.2.2 为什么必须存这两张表（而不是用 `i` 直接索引 `q`）

因为 **$\dim q \neq \dim v$**。位形空间是流形 $\mathcal{Q}=\prod_i\mathcal{Q}_i$，
速度活在它的切空间 $T_q\mathcal{Q}$（李代数），两者维数不同：

| 关节类型 | `nqs[i]` | `nvs[i]` | 位形的参数化 |
|---|---|---|---|
| Revolute / Prismatic | 1 | 1 | $\theta$ |
| RevoluteUnbounded（连续转） | 2 | 1 | $(\cos\theta,\sin\theta)$ |
| Spherical | 4 | 3 | 单位四元数 $\mathbf{q}\in S^3$ |
| FreeFlyer（浮动基） | 7 | 6 | $(p,\mathbf{q})\in SE(3)$ |
| Universe（`i=0`） | 0 | 0 | — |

所以关节 $i$ 在 $q$ 里的偏移和在 $v$ 里的偏移是**两套完全不同的累加**，必须分别记：
`idx_qs/nqs` 走 $q$，`idx_vs/nvs` 走 $v$。

这也是 `integrate` 存在的原因——不能写 $q \leftarrow q + v\,\delta t$，只能

$$q \leftarrow q \oplus (v\,\delta t),\qquad \oplus:\mathcal{Q}\times\mathbb{R}^{n_v}\to\mathcal{Q}$$

而 $\oplus$ 的实现（`liegroup-algo.hxx`）正是按 `nqs[i]/nvs[i]` 逐关节分派到各自李群的 `exp` 上（见 §10）。

#### 2.2.3 `njoints` 与 `nq` 的关系

**两者之间没有固定关系式**——`njoints` 数的是"树上有几个铰"，`nq` 数的是"位形向量有多长"。
唯一严格成立的是求和式：

$$n_q=\sum_{i=1}^{\text{njoints}-1}\texttt{nqs}[i]$$

求和从 **1** 开始，因为下标 0 是 universe（`nqs[0] = 0`）。即「真实关节数 = `njoints - 1`」。
换句话说：**`njoints` 是计数，`nq` 是加权和**，权重就是每个关节的 `nqs[i]`。

本仓库自带样例模型实测：

```
manipulator : njoints=7   nq=6   nv=6
humanoid    : njoints=30  nq=35  nv=34
```

humanoid 的构成（按关节类型拆开）：

| 关节类型 | 个数 | 每个 `nq` | 对 `nq` 的贡献 |
|---|---|---|---|
| universe（id=0） | 1 | 0 | 0 |
| `JointModelFreeFlyer` | 1 | 7 | 7 |
| 各类 revolute | 28 | 1 | 28 |
| **合计** | **30** | | **35** |

三种典型情形：

| 情形 | 关系 | 例子 |
|---|---|---|
| 纯单自由度关节链（工业臂） | $n_q = \text{njoints}-1$ | manipulator：$7-1=6$ |
| 带浮动基（人形/四足） | $n_q = \text{njoints}+5$ | humanoid：$30+5=35$ |
| 含 mimic 关节 | $n_q < \text{njoints}-1$ | mimic 的 `nqs[i] = 0`，白占一个 id |

浮动基那条为什么是 $+5$：FreeFlyer 一个关节贡献 7 维而不是 1 维，
比"每关节 1 维"的基线多出 6，于是 $(\text{njoints}-1)+6$。

上下界：$0 \le n_q \le 7\,(\text{njoints}-1)$，上界取在全是 FreeFlyer 时，下界取在全是 mimic 时。

> **三个容易混的点**
> 1. `njoints` **不是自由度数**。要自由度用 `nv`，要状态向量长度用 `nq`，两个都不等于 `njoints`。
> 2. URDF 里的 `<joint>` 标签数 **≠** `njoints`：`fixed` 类型关节在解析时被折叠，
>    不进 `joints` 而是变成一个 `Frame`。40 个 `<joint>` 的 URDF 可能只建出 20 个关节。
> 3. `nbodies` 与 `njoints` 通常相等但语义不同：前者数刚体，后者数铰。
>    用 `appendBodyToJoint` 往同一关节挂多个刚体时两者分离。

#### 2.2.4 索引表如何建立：前缀和

`addJoint`（`model.hxx:866-898`）里同步维护三套表：

```cpp
jmodel.setIndexes(joint_id, nq, nv, nvExtended);  // 把当前累加值当作本关节的起始偏移
...
nq += joint_nq;              // 先给自己发号，再累加
nqs.push_back(joint_nq);
idx_qs.push_back(joint_idx_q);
```

关键点：`setIndexes` 传进去的 `nq` 是**添加本关节之前**的总维数，所以它天然就是本关节的 `idx_q`。
也就是前缀和

$$\texttt{idx\_qs}[i]=\sum_{k<i}\texttt{nqs}[k]$$

因此**布局顺序 = 关节添加顺序**；而 Pinocchio 保证 `parents[i] < i`（父节点先添加），
于是 $q$ 里的分段天然按拓扑序排列，`idx_qs` 严格递增。
`setIndexes` 同时把全局起点**写回关节对象自身**——这就是 `jmodel.idx_q()` 能从全局 `q` 里
取出自己那几维的机制。

#### 2.2.5 三个容易踩的陷阱

1. **下标 0 永远是 universe**。构造函数里 `idx_qs(1,0), nqs(1,0), idx_vs(1,0), nvs(1,0)`
   （`model.hxx:264-269`）是给 `id=0` 的 universe 关节占位的，它不消耗任何自由度，全为 0。
   真实关节从 1 开始——这是读 Pinocchio 时最常见的偏移错误。
2. **`nqs[i]` 不是子树的维度**，只是这一个关节自己的。子树用 `subtrees[i]` / `data.nvSubtree[i]`。
3. **`Data` 里的量按 `idx_v` 排布，不是 `idx_q`**：`J`、`M`、`tau`、`nle` 全是 $n_v$ 维/列的
   （`model.hxx:959-966` 建稀疏模式时用的就是 `idx_vs/nvs`）。只有 $q$ 本身，以及
   `lowerPositionLimit/upperPositionLimit/positionLimitMargin` 用 `idx_q`。

#### 2.2.6 `nvExtended` 概览

`nvExtended`（`model.hxx:99-100`，文档里叫 *jacobian space*）是第三套索引空间：

- `nv` = **独立**自由度个数（$\dot q$、$M$、$\tau$ 的维数）
- `nvExtended` = **树里所有关节各自贡献的运动子空间列数之和**，不管它们是否独立

每个具体关节的 `traits` 里 `NVExtended` 与 `NV` 取**同一个编译期常量**
（如 `joint-revolute.hxx:641` 的 `NVExtended = 1`、`joint-free-flyer.hxx:158` 的 `= 6`），
默认的 `nvExtended_impl()` 直接返回它（`joint-model-base.hxx:167-170`）。
所以**普通模型里 `nvExtended == nv`，三套索引完全重合，这套机制零开销地消失**。
唯一重写它的是 **`JointModelMimic`**。实测（Baxter，含两对 mimic 夹爪）：

| 模型 | `nq` | `nv` | `nvExtended` |
|---|---|---|---|
| `buildModelFromUrdf(path)`（忽略 mimic） | 19 | 19 | 19 |
| `buildModelFromUrdf(path, mimic=True)` | 17 | **17** | **19** |

即 mimic 模型少了 2 个独立自由度，但内部仍需按 19 列展开来做递推，
再用传动矩阵 $G$ 投影回 17 维。**完整机制见 §9.2**。

### 2.3 全部数据字段（按用途分组）

`ModelTpl` 有 40+ 个公开字段。按用途归类如下（可当速查表用）：

**① 维度与计数**

| 字段 | 类型 | 含义 |
|---|---|---|
| `nq` / `nv` / `nvExtended` | `int` | 三套总维度（见上） |
| `njoints` | `int` | 关节数，**含 universe(0)** |
| `nbodies` | `int` | 刚体数 |
| `nframes` | `int` | 坐标系数 |
| `name` | `std::string` | 模型名 |

**② 逐关节的静态属性**（长度均为 `njoints`）

| 字段 | 类型 | 含义 |
|---|---|---|
| `joints` | `JointModelVector` | 关节模型对象（variant，见 §5） |
| `jointPlacements[i]` | `SE3` | 关节 $i$ 相对**父关节系**的固定位姿 ${}^{\lambda(i)}M_i$ |
| `inertias[i]` | `Inertia` | 关节 $i$ 所驮刚体的空间惯量（`appendBodyToJoint` 累加而来） |
| `names[i]` | `std::string` | 关节名（`getJointId` 的反查表） |
| `idx_qs` / `nqs` | `vector<int>` | 关节 $i$ 在 $q$ 中的起点与长度 |
| `idx_vs` / `nvs` | `vector<int>` | 关节 $i$ 在 $v$ 中的起点与长度 |
| `idx_vExtendeds` / `nvExtendeds` | `vector<int>` | 扩展速度空间中的起点与长度 |

**③ 树拓扑**（这组是所有 $O(n)$ 算法的稀疏性来源）

| 字段 | 含义 | 实测示例（样例人形，关节 5 = `rleg_elbow_joint`） |
|---|---|---|
| `parents[i]` | 父关节，`universe`=0 为根 | `parents[5] = 4` |
| `children[i]` | 子关节列表 | `children[5] = [6]` |
| `supports[j]` | 从根到 $j$ 的**路径**（含两端） | `supports[5] = [0,1,2,3,4,5]` |
| `subtrees[j]` | 以 $j$ 为根的**整棵子树**（含自身） | `subtrees[5] = [5,6,7]` |
| `sparsity_pattern_vector[i]` | 关节 $i$ 的雅可比**哪些列非零**（布尔） | — |
| `span_indexes_vector[i]` | 同上，索引列表形式 | — |

> **怎么记**：`supports` 往**上**看（我依赖谁）→ 雅可比列稀疏性；
> `subtrees` 往**下**看（谁依赖我）→ CRBA 复合惯量的累加范围（见 [空间代数 §6.3](空间代数运算解析.md)）。

**④ 关节限位与驱动参数**（长度 `nv` 或 `nq`）

| 字段 | 维度 | 含义 |
|---|---|---|
| `lowerPositionLimit` / `upperPositionLimit` | `nq` | 位置限位。`randomConfiguration` 在此范围采样 |
| `positionLimitMargin` | `nq` | 位置限位的**提前激活缓冲带**（详见下方专栏，⚠️ 仅正向动力学用） |
| `lowerVelocityLimit` / `upperVelocityLimit` | `nv` | 速度限位 |
| `lowerEffortLimit` / `upperEffortLimit` | `nv` | 力矩限位 |
| `lowerDryFrictionLimit` / `upperDryFrictionLimit` | `nv` | 库仑（干）摩擦上下限 |
| `damping` | `nv` | 粘滞阻尼系数 |
| `armature` | `nv` | **电枢惯量**：加到 $M(q)$ 对角上的等效转子惯量 |
| `rotorInertia` / `rotorGearRatio` | `nv` | 转子惯量与减速比（`armature` 的来源） |

> ⚠️ **浮动基的位置限位默认是 $\pm\infty$**，所以 `randomConfiguration` 在人形上会失败——
> 必须先手工设 `model.lowerPositionLimit[:7] = -1` 之类（见 `examples/run-algo-in-parallel.py`）。

> **`armature` 的实际意义**：谐波减速器的转子折算到关节侧的惯量为 $n^2 I_{\text{rotor}}$（$n$ 为减速比），
> 对高减速比关节这一项可能与连杆惯量同量级。忽略它会让模型显著偏软。

#### 2.3.1 专栏：`positionLimitMargin` 是什么、只在正向动力学用

> ⚠️ **备注：此字段只在正向动力学（约束/接触前向动力学，即仿真）中使用；逆动力学（RNEA）用不到。后续读到"约束前向动力学"算法时再回看本节。**

**它是什么**：`nq` 维向量，每个配置坐标一个值，是关节限位约束的**提前激活缓冲带**——决定关节离限位还有多远时，就把限位约束**提前打开**。由 `addJoint` 写入，**默认 0**。消费者是 [`constraints/joint-limit-constraint.hxx`](../include/pinocchio/src/constraints/joint-limit-constraint.hxx) 的 `JointLimitConstraintModel`，激活判据：

$$\text{下限激活：}\ q_i-\text{lower}_i\le \text{margin}_i,\qquad \text{上限激活：}\ \text{upper}_i-q_i\le \text{margin}_i$$

**为什么只属于正向动力学**：

| | 输入→输出 | 限位有没有用 |
|---|---|---|
| 逆动力学 RNEA | $(q,v,a)\to\tau$ | ❌ 运动已给定，只反算力矩，不需要"阻止越界" |
| 正向动力学（约束） | $(q,v,\tau)\to a$ | ✅ 运动未知；若 $\tau$ 会把关节推过限位，需**约束力**顶住 |

关节限位是**单边（不等式）约束**（像一堵"墙"：没到墙边约束不存在，接近/触墙才产生约束力）。"要不要产生约束力"只有在**解算运动**（正向动力学）时才有意义，故 margin 只出现在约束前向动力学 / 仿真里。

**典型用法（仿真主循环）**：

```cpp
model.positionLimitMargin = VectorXd::Constant(model.nq, 0.05);   // 设 0.05 缓冲带
JointLimitConstraintModel limit_cm(model, activable_joints);      // 读入 lower/upper/margin
auto limit_cd = limit_cm.createData();

for (每个仿真步) {
  tau = controller(q, v);
  limit_cm.makeSelectionFilteredByLimitProximity(q);  // ★ 按 margin 选出接近限位的关节进激活集
  a   = 约束前向动力学(model, data, q, v, tau, limit_cm, limit_cd); // 激活约束产生约束力顶住边界
  v  += a*dt;  q = integrate(model, q, v*dt);
}
```

**margin 的直观效果**：`margin=0` → 到限位才激活，一步可能冲过再弹回（抖动/穿透）；`margin>0` → 到墙之前提前介入、平滑减速停在边界内。它是可调的**安全提前量**（每坐标独立）。约束动力学的 `Data` 字段见 [§3.8](#38-辨识回归量与约束动力学)。

**⑤ mimic（耦合关节）**

| 字段 | 含义 |
|---|---|
| `mimicking_joints` | 跟随者关节 ID 列表（实测 Baxter：`[10, 19]`） |
| `mimicked_joints` | 被跟随者 ID 列表（实测：`[9, 18]`，与上一一对应） |
| `mimic_joint_supports` | mimic 关节的支撑路径 |

**⑥ 其余**

| 字段 | 含义 |
|---|---|
| `frames` | 所有 Frame（见 §4） |
| `referenceConfigurations` | 命名配置字典（SRDF 的 `half_sitting` 等） |
| `gravity` | 重力，`Motion` 类型，默认线性部分 $(0,0,-9.81)$ |
| `gravity981` | 静态常量 $(0,0,-9.81)$ |

> **重力为什么是 `Motion` 而不是 `Vector3`**：RNEA 把根节点加速度初始化为 $a_0 = -g$，
> 让重力**沿递推自动传遍全身**（见 [GUIDE §4.3.2](../PINOCCHIO_GUIDE.md#432-第一趟正向传播运动学根--叶)）。
> 存成 `Motion` 才能直接参与空间加速度递推。把 `gravity` 置零即可做"零重力"仿真。

### 2.4 全部方法

**① 构造与赋值**

| 方法 | 说明 |
|---|---|
| `ModelTpl()` | 默认构造：建立 `universe`（关节 0），`nq=nv=0`, `njoints=1` |
| `ModelTpl(const ModelTpl&)` | 拷贝 |
| `ModelTpl(const ModelTpl<S,O,OtherCollection>&)` | 跨**关节集合**转换 |
| `operator=` ×2 | 同上两种赋值 |
| `operator==` / `!=` | 逐字段比较（`mimic_dynamics.py` 用它验证手工 mimic 与 URDF 解析结果一致） |
| `cast<NewScalar>()` | 换标量类型，autodiff/多精度的入口（见 [GUIDE §11.2](../PINOCCHIO_GUIDE.md#112-tpl-模板--context-默认标量)） |

**② 建树**

```cpp
JointIndex addJoint(JointIndex parent, const JointModelBase<D>& jmodel,
                    const SE3& joint_placement, const std::string& name, ...);
```

**7 个重载**，参数逐级增加（是否给力矩/速度/位置限位、摩擦、阻尼），最终都**汇聚到同一个 13 参总实现**（见 §2.5 逐步解析）。

> ⚠️ **必须按深度优先顺序添加**：`parents[i] < i` 是所有递推算法的前提
> （正向遍历 `for i=1..njoints` 时父节点必已算完）。传入未注册的 `parent` 会破坏这个不变量。

### 2.5 `addJoint` 主实现逐步解析（建树的核心）

所有 `addJoint` 重载最终调用同一个 **13 参总实现**（[model.hxx:802](../include/pinocchio/src/multibody/model.hxx#L802)）。它做的事：**把一个关节追加进树，并同步维护 Model 里几十个并行数组与拓扑结构**。理解它 = 理解 §2.3 那些字段是怎么长出来的。

**13 个参数的含义**（用户常问的那一长串）：

| 参数 | 维度 | 作用 |
|------|------|------|
| `parent` | — | 父关节 id（必须已注册，保证 `parent < 新 id`） |
| `joint_model` | — | 关节模型（决定 nq/nv、运动子空间 S） |
| `joint_placement` | SE3 | 关节相对父关节系的固定位姿 ${}^{\lambda(i)}M_i$ |
| `joint_name` | — | 关节名（`getJointId` 反查用） |
| `min_effort`/`max_effort` | nv | 力矩（广义力）下/上限 |
| `min_velocity`/`max_velocity` | nv | 速度下/上限 |
| `min_config`/`max_config` | nq | 位置下/上限 |
| `config_limit_margin` | nq | 限位提前激活缓冲带（见 §2.3.1，仅正向动力学用） |
| `min_joint_friction`/`max_joint_friction` | nv | 干摩擦下/上限 |
| `joint_damping` | nv | 粘滞阻尼 |

**十个执行阶段**（源码已加对应中文注释）：

| 阶段 | 做什么 | 关键点 |
|------|--------|--------|
| **① 校验** | 断言 4 个并行数组长度 == `njoints`；`PINOCCHIO_CHECK` 所有限位向量尺寸、`min≤max`、`parent` 合法 | 不变量：`joints/inertias/parents/jointPlacements` 始终同步等长 |
| **② 分配 id + 拷贝 + 写回索引** | `joint_id = njoints++`；`joints.push_back(拷贝)`；`jmodel.setIndexes(id, nq, nv, nvExtended)` | **把"本关节在全局 q/v 里的起点 = 当前累计 nq/nv"写回关节自身**——这就是 `jmodel.idx_q()` 能定位的机制 |
| **③ 读回尺寸** | 从 `jmodel` 取 `joint_nq/idx_q/joint_nv/idx_v/...` | 供后续累加与稀疏模式用 |
| **④ push 树属性** | `inertias`(先 Zero)、`parents`、`children`、`children[parent]`(反向登记)、`jointPlacements`、`names` 各追加一格 | 惯量置零，等 `appendBodyToJoint` 再填 |
| **⑤ 累加全局维度** | `nq += joint_nq` 等；push `nqs/idx_qs/nvs/idx_vs/nvExtendeds/idx_vExtendeds` | 全局 nq/nv/nvExtended 在此增长 |
| **⑥ resize + 写限位/驱动** | 每个限位/驱动向量 `conservativeResize(nv 或 nq)` 保留旧值扩容，再用 `jointVelocitySelector`/`jointConfigSelector` 精确写入本关节段 | fixed 关节（nq==0）跳过。armature/rotorInertia 置零、rotorGearRatio 置一 |
| **⑦ subtrees** | 本关节子树先只含自己；`addJointIndexToParentSubtrees` 沿 parents 上溯，把 id 加进**每个祖先**的子树 | CRBA 累加范围的来源 |
| **⑧ supports** | 继承父的"根→父"路径，末尾接自己 = **根→本关节** 的完整路径 | 雅可比列稀疏性的来源 |
| **⑨ 雅可比稀疏模式** | 先把已有布尔向量扩容到新 nv 并清零尾部；再由 supports 路径构建本关节的**非零列集合** `extended_support`（= 各祖先的 v 段 + 自己的 v 段），写进 `span_indexes_vector`（索引）与 `sparsity_pattern_vector`（布尔掩码） | 含义：末端关节速度由"根到它路径上所有关节的速度"决定，故这些列非零 |
| **⑩ mimic 记账** | `mimic_joint_supports` 继承父路径；若本关节是 mimic 类型，追加自己并登记 `mimicking_joints`(跟随者)→`mimicked_joints`(被跟随者) | 见 [§9.2](#92-mimicjoint-mimichxx) 与 `nvExtended` |

**一句话**：`addJoint` = 校验 → 把关节拷进树并把全局起点写回关节自身 → 同步累加维度与所有并行数组 → 用 supports/subtrees/稀疏模式把树拓扑固化下来（供 $O(n)$ 算法用）→ 处理 mimic 耦合。**它是 §2.3 所有字段的唯一构造者。**

> ⚠️ **源码里发现一处疑似上游 BUG**（[model.hxx](../include/pinocchio/src/multibody/model.hxx) 阶段⑥）：`lowerVelocityLimit` 被赋 `max_velocity`，但按同类项模式（`lowerEffortLimit=min_effort`、`lowerPositionLimit=min_config`、`lowerDryFrictionLimit=min_joint_friction`）应为 `min_velocity`。已在源码加注释标记，**未改动行为**（改动会影响下游，需上游确认）。

**③ 其它建树辅助方法**

| 方法 | 作用 |
|------|------|
| `appendBodyToJoint(joint_id, Y, placement)` | 把刚体空间惯量 `Y`（换算到关节系后）累加进 `inertias[joint_id]` |
| `addJointFrame(joint_id, prev_frame)` | 给关节挂一个同名 `JOINT` 型 Frame（冗余但便于查询） |
| `addBodyFrame(name, parentJoint, placement, parentFrame)` | 加一个 `BODY` 型 Frame |
| `addFrame(frame, append_inertia)` | 注册任意 Frame；若带惯量且 `append_inertia`，惯量并入父关节 |
| `createData()` | 按当前 Model 分配配对的 `Data` |
| `getJointId/getBodyId/getFrameId` `existJointName/...` | 按名字反查索引（线性扫 `names`/`frames`） |

**④ `cast<NewScalar>()` / `operator=` / `operator==`**

- `cast`：逐字段 `.cast<>()`——标量数组、`joints[k].cast()`、`inertias[k].cast()`、`frames[k].cast()`、`referenceConfigurations` 逐条转换，产出 `ModelTpl<NewScalar,…>`（autodiff/多精度入口，见 [源码解析.md 条目 2](源码解析.md)）。
- `operator=`（跨关节集合模板版）：逐字段拷贝，`joints` 用 `push_back` 逐个转换。
- `operator==`：逐字段比较，含 `referenceConfigurations` 逐键逐向量、各限位向量尺寸+数值。

| 其余建树方法 | 说明 |
|---|---|
| `appendBodyToJoint(joint_id, Y, placement)` | 把刚体惯量 $Y$ 变换到关节系后**累加**到 `inertias[joint_id]`（多个几何体可叠加到同一关节） |
| `addJointFrame(joint_index, previous_frame_index)` | 为关节自动注册一个 `JOINT` 类型的 Frame |
| `addBodyFrame(...)` | 注册 `BODY` 类型 Frame |
| `addFrame(frame, append_inertia=true)` | 注册任意 Frame；若带惯量且 `append_inertia=true`，惯量并入父关节 |
| `addJointIndexToParentSubtrees(joint_id)` | 内部维护 `subtrees` 用 |

**③ 名称查询**（URDF 名 ↔ 内部索引）

| 方法 | 说明 |
|---|---|
| `getJointId(name)` / `existJointName(name)` | 关节名 → ID。**查不到返回 `njoints`（不抛异常）**，故批量转换时应先 `exist*` 判断（`build-reduced-model.py` 就是这么写的） |
| `getFrameId(name, type)` / `existFrame(name, type)` | Frame 名 → ID，可按类型过滤 |
| `getBodyId(name)` / `existBodyName(name)` | 刚体名 → ID |

> ⚠️ URDF 里的**固定关节**在 Pinocchio 中不是关节而是 Frame（见 §4），
> 所以很多"关节名"必须用 `getFrameId` 而非 `getJointId` 才查得到。

**④ 校验与工具**

| 方法 | 说明 |
|---|---|
| `createData()` | 按当前 Model 分配配套 `Data`。**Model 改了必须重建 Data** |
| `check()` / `check(data)` / `check(checker)` | 模型自洽性校验；带 `Data` 的版本检查二者是否匹配 |
| `hasConfigurationLimit()` | 返回长度 `nq` 的 `vector<bool>`：每个位形分量是否**有**限位（浮动基的四元数分量为 false） |
| `hasConfigurationLimitInTangent()` | 同上但按 `nv` 维 |
| `getChildJoints()` | 取叶子关节列表 |

---

## 3. `DataTpl`：算法工作区

`src/multibody/data.hxx`。如果说 `Model` 是"机器人的图纸"，`Data` 就是**"算草稿纸"**——~150 个预分配的字段，让 RNEA/CRBA/ABA 等算法**零动态分配**地反复运行。字段按算法族分组：

#### 命名约定（读懂 Data 的第一把钥匙）

`Data` 字段名遵循几条固定前后缀规则，认得它们就能猜出大半字段的含义：

| 记号 | 含义 | 例 |
|------|------|-----|
| 前缀 `o` | **世界系（origin）表达** | `v[i]` 是 LOCAL 速度，`ov[i]` 是同一速度在世界系的表达 |
| 前缀 `li` | local relative to parent | `liMi[i]` = ${}^{\lambda(i)}M_i$ |
| 后缀 `_fromRow` | 按**自由度行**（而非关节）重排的索引 | `parents_fromRow` |
| 后缀 `_augmented` | mimic 展开后的增广版本 | `joints_augmented` |
| `d` 前缀 / `_d*` | 时间导数或偏导 | `dJ`、`dtau_dq` |
| `crb` | Composite Rigid Body | `Ycrb`、`oYcrb` |
| `_in` | 算法入参的缓存副本 | `q_in`、`v_in` |

> **为什么同一物理量要存 LOCAL 和 WORLD 两份**：RNEA 的递推在 LOCAL 系最省
> （关节子空间 $S$ 是常向量），但导数与质心量在世界系表达更方便。二者用伴随变换互转
> （见 [空间代数 §4](空间代数运算解析.md)），预存两份是**空间换时间**。

### 3.1 运动学量（每关节一个，`std::vector` 长度 = njoints）

| 字段 | 含义 | 由谁写入 |
|---|---|---|
| `oMi[i]` | 关节 $i$ 在世界系的位姿 ${}^0M_i$ | `forwardKinematics` 主输出 |
| `liMi[i]` | 关节 $i$ 相对父关节的位姿 ${}^{\lambda(i)}M_i$（**含**关节自由度产生的变换） | `forwardKinematics` |
| `v[i]` / `ov[i]` | 空间速度（LOCAL / WORLD） | `forwardKinematics(q,v)` |
| `a[i]` / `oa[i]` | 空间加速度（LOCAL / WORLD） | `forwardKinematics(q,v,a)` |
| `a_gf[i]` / `oa_gf[i]` | **含重力**的加速度（gf = gravity field） | RNEA。构造时 `a_gf[0] = -model.gravity` |
| `oa_drift` | 漂移加速度 $\dot J v$（零加速度下的 $a$） | 约束动力学 |
| `oMf[i]` | Frame 在世界系的位姿 | `updateFramePlacements` |
| `joints` / `joints_augmented` | 各关节的 `JointData`（缓存 $S$、$M_J$ 等） | `jmodel.calc()` |
| `q_in` / `v_in` / `a_in` / `tau_in` | 入参缓存 | 各算法入口 |

> `a_gf[0] = -model.gravity` 这一行是**重力注入的全部机密**：把根节点加速度设为 $-g$，
> 递推自然把重力效应传遍全身，无需在每个连杆单独加重力项。

### 3.2 力与动量
- `f[i]`/`of[i]`：作用在关节 i 上的空间力（RNEA 反向递推）。
- `h[i]`/`oh[i]`：关节 i 的空间动量 $\mathbf{h}=\mathbf{I}v$。
- `hg`/`dhg`/`Ig`：质心处总动量/其导数/合惯量（质心动量矩阵用）。

### 3.3 质量矩阵与逆动力学
- `M`：关节空间惯量矩阵（CRBA 输出，`nv×nv`，对称）。
- `Minv`：其逆（`RowMatrixXs`）。
- `C`：科氏/离心矩阵。
- `nle`/`g`/`tau`：非线性效应项 / 重力项 / 关节力矩。
- `Ycrb[i]`/`oYcrb[i]`：复合刚体惯量（子树累加，CRBA 核心，见 [空间代数解析 §6.3](空间代数运算解析.md)）；`dYcrb`/`doYcrb` 其导数。

### 3.4 ABA（前向动力学）

| 字段 | 含义 |
|---|---|
| `Yaba[i]` / `oYaba[i]` | **铰接体惯量** $I^A_i$（子关节自由响应后的等效惯量，见 [GUIDE §4.4.2](../PINOCCHIO_GUIDE.md#442-核心概念关节化体惯量)） |
| `u[i]` | 关节 $i$ 的偏置力项 |
| `U` / `D` / `Dinv` | Schur 补分解量：$U=I^A S$、$D=S^\top I^A S$、`Dinv`$=D^{-1}$ |
| `SDinv` / `UDinv` / `IS` | 递推中间矩阵（避免重复乘法） |
| `ddq` | **关节加速度输出** |
| `B` / `vxI` / `Ivx` | $v\times^*I$、$I\times v$ 等惯量时间导数（解析梯度用） |

> `Yaba` 与 `Ycrb` 的区别正是 ABA 与 CRBA 的分野：**`crb` 假设子关节锁死**（对应 $M$），
> **`aba` 假设子关节自由响应**（对应 $M^{-1}$）。

### 3.5 稀疏 Cholesky（$M = U D U^\top$）

| 字段 | 含义 |
|---|---|
| `U` | 单位上三角因子 |
| `D` / `Dinv` | 对角块及其逆 |
| `tmp` | 求解时的工作向量 |
| `parents_fromRow[k]` | 自由度 $k$ 的"父自由度" |
| `nvSubtree_fromRow[k]` | 从行 $k$ 起子树占的自由度数 |
| `supports_fromRow[k]` | 行 $k$ 的支撑集 |
| `start_idx_v_fromRow` / `end_idx_v_fromRow` | 行区间边界 |
| `mimic_parents_fromRow` / `non_mimic_parents_fromRow` | mimic 场景下的分支版本 |

**为什么要 `_fromRow` 这一套**：Model 里的 `parents`/`subtrees` 是**按关节**组织的，
而稀疏 Cholesky 要**按自由度行**遍历（一个球关节占 3 行）。构造 `Data` 时由
`computeParents_fromRow` / `computeSupports_fromRow` / `computeNvSubtree` 一次性把
关节级拓扑"展开"到自由度级，之后分解直接按行跳转，无需反复换算。

### 3.6 解析导数（RNEA/ABA 对 q,v,τ 的偏导）
- `dtau_dq`/`dtau_dv`：逆动力学导数。
- `ddq_dq`/`ddq_dv`/`ddq_dtau`：前向动力学导数。
- `dVdq`/`dAdq`/`dAdv`/`dFdq`... ：中间量。
- `J`/`dJ`/`ddJ`：关节雅可比及其时间导数；`kinematic_hessians`：二阶（`Tensor3x`）。

### 3.7 质心与能量
- `com[i]`/`vcom[i]`/`acom[i]`/`mass[i]`：各子树质心位置/速度/加速度/质量。
- `Jcom`：质心雅可比。
- `kinetic_energy`/`potential_energy`/`mechanical_energy`。

### 3.8 辨识回归量与约束动力学
- `staticRegressor`/`bodyRegressor`/`jointTorqueRegressor`：惯量辨识回归矩阵（见 [空间代数解析 §6.10](空间代数运算解析.md)）。
- `JMinvJt`/`lambda_c`/`impulse_c`/`osim`（操作空间惯量逆）/`KA`/`LA`/`lA`：约束/接触动力学（见 [floating-base-and-contact-dynamics.md](floating-base-and-contact-dynamics.md)）。
- **关节限位约束**也在此类算法里生效：`Model.positionLimitMargin` 决定限位约束的提前激活（**仅正向动力学用**，见 [§2.3.1 专栏](#231-专栏positionlimitmargin-是什么只在正向动力学用)）。后续读约束前向动力学源码时对照该专栏。

> 记忆法：**`Data` 的每个字段几乎都精确对应某个算法的某个中间量**。看到陌生字段，去 `algorithm/` 里 grep 它的名字即可定位用途。

---

## 4. `FrameTpl` / `ModelItem`：坐标系与树节点基类

### 4.1 `ModelItem`（`model-item.hxx`，79 行）

`Frame`（及未来其他"挂在树上的东西"）的公共基类，只有四个字段：

| 字段 | 类型 | 含义 |
|---|---|---|
| `name` | `std::string` | 名称 |
| `parentJoint` | `JointIndex` | **挂在哪个关节上**（决定它随谁运动） |
| `parentFrame` | `FrameIndex` | 父坐标系（主要供 URDF 层级记录/第三方使用，不参与运算） |
| `placement` | `SE3` | 相对**父关节系**的固定位姿 ${}^iM_f$（常量，不随 $q$ 变） |

抽出这个基类是为了让未来的传感器、附着体等复用同一套"挂载"语义。

### 4.2 `FrameTpl`（`frame.hxx`，223 行）

在 `ModelItem` 基础上增加两个字段：

| 字段 | 含义 |
|---|---|
| `type` | `FrameType` 枚举，见下 |
| `inertia` | 该 Frame 携带的惯量（`addFrame(frame, append_inertia=true)` 时会并入父关节） |

**方法**（都很轻）：

| 方法 | 说明 |
|---|---|
| `FrameTpl()` | 默认构造 |
| `FrameTpl(name, parentJoint, placement, type, inertia)` 等 3 个重载 | 常规构造 |
| `operator==` / `!=` | 支持**跨标量类型**比较（模板参数 `S2,O2`） |
| `cast<NewScalar>()` | 换标量类型 |

### 4.3 `FrameType`：位标志枚举

```cpp
enum FrameType {
  OP_FRAME    = 0x1 << 0,   // 1  用户运行时自定义的操作坐标系
  JOINT       = 0x1 << 1,   // 2  关节坐标系（关节的镜像）
  FIXED_JOINT = 0x1 << 2,   // 4  URDF 中被折叠掉的固定关节
  BODY        = 0x1 << 3,   // 8  连杆的惯量/视觉/碰撞坐标系
  SENSOR      = 0x1 << 4    // 16 传感器坐标系
};
```

**用位标志而非顺序枚举，是为了支持按位或的组合查询**：

```cpp
model.getFrameId("tool0", FrameType(JOINT | BODY));   // 在两类里找
```

> 实测（样例人形）：70 个 Frame 的类型分布为 `BODY:40, JOINT:29, FIXED_JOINT:1`；
> 按 `JOINT|BODY` 组合查询能正确命中。`getFrameId` 的 `type` 参数默认是**全部类型的并集**。

### 4.4 Frame vs Joint：最常见的困惑

| | **Joint** | **Frame** |
|---|---|---|
| 是不是树节点 | ✅ 运动学树的**真实节点** | ❌ 只是**挂件** |
| 有无自由度 | 有（贡献 nq/nv） | **无** |
| 参与递推吗 | 参与全部动力学递推 | **不参与** |
| 位姿从哪来 | `data.oMi[i]`，递推算出 | `data.oMf[f] = oMi[parentJoint] * placement`，**事后**乘出来 |
| 数量级 | 少（实测人形 30） | 多（实测人形 70） |

**关键机制**：URDF 里的**固定关节不会进入运动学树**，而是被降级成 `FIXED_JOINT` 类型的 Frame。
这就是"URDF 里几十个 link + 固定 joint"被压缩成"少数活动关节 + 一堆 Frame"的原因，
也是很多人 `getJointId("某个URDF关节名")` 查不到的根源——**得用 `getFrameId`**。

> Frame 层的位姿更新 API（`framesForwardKinematics` / `updateFramePlacement(s)`）
> 及其两个易错陷阱，见 [GUIDE §4.1.2–4.1.3](../PINOCCHIO_GUIDE.md#412-frame-位姿更新的三个-api)。

---

## 5. 关节系统总览：Model/Data 分离 + CRTP + variant 双层

关节是本子系统最精巧的部分。它要同时满足三个矛盾需求：
**(a)** 每种关节数学不同（转/移/球/浮动…）；
**(b)** 一棵树里混装不同关节，要能存进一个 `std::vector`；
**(c)** 递推热循环里不能有虚函数开销。
解法是**四个正交设计**叠加。

### 5.1 Model 对象 vs Data 对象（职责分离）

每种关节都成对出现：

| | `JointModelXxx` | `JointDataXxx` |
|---|---|---|
| 存什么 | **静态参数**：轴向、索引 `idx_q/idx_v`、id | **每次 calc 的结果**：$S$、$M_J$、$v_J$、$c_J$、$U$、$D^{-1}$ |
| 存在哪 | `model.joints[i]` | `data.joints[i]` |
| 可变性 | 只读（算法阶段） | 每步被 `calc` 覆写 |
| 线程 | 多线程共享 | **每线程一份** |

这与 `Model`/`Data` 的整体分离是同一思想的下沉（见 [GUIDE §3](../PINOCCHIO_GUIDE.md)）。

### 5.2 CRTP 基类（静态多态，零开销）

```
JointModelBase<JointModelRX>  ←  JointModelRX
JointDataBase <JointDataRX>   ←  JointDataRX
```

基类提供**统一 API**（`calc`、`nq()`、`idx_q()`、段选择器…），派生类提供 `*_impl`。
全程 `static_cast` 静态分派，**可内联、无虚表**（原理见 [GUIDE §11.3](../PINOCCHIO_GUIDE.md#113-crtp--boostfusion-访问者)）。

### 5.3 `JointCollectionDefaultTpl` → `boost::variant`（类型擦除）

`joint/joint-collection.hxx` 列出库支持的**全部**关节类型，打包成：

```cpp
typedef boost::variant<JointModelRX, JointModelRY, ..., JointModelFreeFlyer, ...> JointModelVariant;
```

于是 `model.joints` 可以是 `std::vector<JointModelVariant>`——**异构关节存进同一个容器**。
代价是访问时需要 `boost::apply_visitor` 分派（见 §11）。

### 5.4 `JointModelTpl`（`joint-generic.hxx`）：变体包装器

它**同时**继承 `JointModelBase<JointModelTpl>` 和持有 `JointModelVariant`：

```cpp
struct JointModelTpl : JointModelBase<JointModelTpl>, JointCollection::JointModelVariant
```

作用是把"变体"重新包装成一个**满足 CRTP 接口**的类型，让上层代码
（`model.joints[i].calc(...)`）写起来和具体关节一模一样，内部自动做 variant 分派。

> 这四层叠起来的效果：**写代码像多态、跑起来像手写 switch**。

### 5.5 `JointModelBase` 完整 API

**① 核心计算**（转发到派生类的 `calc_impl` / 自由函数）

| 方法 | 说明 |
|---|---|
| `calc(data, q)` | 零阶：只算位形相关量（$M_J$、$S$） |
| `calc(data, Blank(), v)` | 一阶但**跳过位形**：`Blank` 是占位标签，表示"只更新速度" |
| `calc(data, q, v)` | 一阶：位形 + 速度（算出 $v_J$、$c_J$） |
| `calc_aba(data, armature, I, update_I)` | ABA 专用，见 §7 |

**② 维度与索引**

| 方法 | 说明 |
|---|---|
| `nq()` / `nv()` / `nvExtended()` | 本关节占的维度 |
| `idx_q()` / `idx_v()` / `idx_vExtended()` | 在**全局**向量中的起始下标 |
| `id()` | 关节在树中的 ID |
| `setIndexes(id, q, v[, vExtended])` | 由 `Model::addJoint` 调用，把全局起点写进关节自身 |

**③ 段/块选择器**（本类最实用的一组）

这些函数把"从全局向量里取出属于本关节的那几维"封装成一次调用，**返回视图、零拷贝**
（视图机制见 [GUIDE §11.5](../PINOCCHIO_GUIDE.md#115-eigen-视图类型零拷贝的实现基础)）：

| 选择器 | 从什么里取 | 取出什么 |
|---|---|---|
| `jointConfigSelector(q)` | 全局 $q$（nq 维） | 本关节的 $q_i$ 段 |
| `jointVelocitySelector(v)` | 全局 $v$（nv 维） | 本关节的 $v_i$ 段 |
| `JointMappedConfigSelector(q)` | 同上 | mimic 映射后的版本 |
| `JointMappedVelocitySelector(v)` | 同上 | 同上 |
| `jointCols(J)` | $6\times n_v$ 矩阵 | 本关节对应的**列块**（雅可比装配用） |
| `jointExtendedModelCols(J)` | 同上 | 扩展速度空间的列块 |
| `jointRows(A)` | $n_v\times k$ 矩阵 | 本关节对应的**行块** |
| `jointBlock(M)` | $n_v\times n_v$ 矩阵 | 本关节的**对角方块**（质量矩阵装配用） |

每个都有 const / 非 const 两个重载。CRBA 里那句
`data.M.block(idx_v(), idx_v(), nv(), nvSubtree[i])` 就是这类操作的手写展开。

**④ 其余**

| 方法 | 说明 |
|---|---|
| `lieGroup()` | 返回本关节对应的**李群对象**（见 §10.4 `LieGroupMap`） |
| `hasConfigurationLimit()` / `...InTangent()` | 各分量是否有限位（浮动基四元数分量为 false） |
| `shortname()` / `classname()` | 类型名字符串，如 `"JointModelRX"`（`disp` 与调试用） |
| `cast<NewScalar>()` | 换标量 |
| `operator==` / `isEqual` / `hasSameIndexes` | 比较；`hasSameIndexes` 只比索引不比参数 |
| `disp(os)` | 打印 |

### 5.6 `JointDataBase` 完整 API

`JointData` 是 `calc` 的**输出容器**，接口很短但每一项都对应递推公式里的一个符号：

| 访问器 | 数学符号 | 含义 |
|---|---|---|
| `S()` | $S_i$ | **运动子空间**（$6\times n_v$），见 §6 |
| `M()` | $M_{J_i}(q_i)$ | 关节自由度产生的 SE(3) 变换 |
| `v()` | $v_{J_i} = S_i\dot q_i$ | 关节自身贡献的空间速度 |
| `c()` | $c_{J_i}$ | **速度积偏置项**（bias），$\dot S_i\dot q_i$ 类项 |
| `U()` / `Dinv()` / `UDinv()` | $I^AS$、$(S^\top I^AS)^{-1}$ | ABA 的 Schur 补分解量 |
| `StU()` | $S^\top U$ | 同上中间量 |

其余：`shortname()`、`disp()`、`operator==`、`isEqual()`。

> **对照递推公式记忆**（见 [GUIDE §4.3.2](../PINOCCHIO_GUIDE.md#432-第一趟正向传播运动学根--叶)）：
> $\nu_i = {}^iX_{\lambda(i)}\nu_{\lambda(i)} + \underbrace{S_i\dot q_i}_{\texttt{v()}}$，
> $a_i = {}^iX_{\lambda(i)}a_{\lambda(i)} + S_i\ddot q_i + \underbrace{\nu_i\times S_i\dot q_i}_{\text{与 }\texttt{c()}\text{ 相关}}$。
> **`JointData` 的每个字段都在公式里有对应项**。

---

## 6. `JointMotionSubspace`（运动子空间 S）

### 6.1 定义

$S_i \in \mathbb{R}^{6\times n_{v_i}}$ 把关节速度映射为空间速度：

$$\nu_{J_i} = S_i\,\dot q_i$$

对偶地，力矩由虚功原理投影得到：

$$\tau_i = S_i^\top f_i$$

**它是"这个关节允许什么方向的运动"的矩阵表述**，也是关节类型之间唯一的本质差异——
上层 FK/RNEA/ABA/CRBA 的递推公式对所有关节**长得完全一样**，差异全被 $S$ 吸收。

### 6.2 各关节的 $S$

| 关节 | $n_v$ | $S$（按 Pinocchio 的 $[v;\omega]$ 顺序） |
|---|---|---|
| `JointModelRZ` | 1 | $[0,0,0,\ 0,0,1]^\top$ |
| `JointModelPX` | 1 | $[1,0,0,\ 0,0,0]^\top$ |
| `JointModelSpherical` | 3 | $\begin{bmatrix}0_{3\times3}\\ \mathbb{1}_3\end{bmatrix}$ |
| `JointModelTranslation` | 3 | $\begin{bmatrix}\mathbb{1}_3\\ 0_{3\times3}\end{bmatrix}$ |
| `JointModelFreeFlyer` | 6 | $\mathbb{1}_6$ |
| `JointModelRevoluteUnaligned` | 1 | $[0_3;\ \hat n]$，$n$ 为任意单位轴 |

### 6.3 为什么要为每种关节写专门的 $S$ 类型

因为**大多数 $S$ 都极度稀疏**。若统一用稠密 $6\times n_v$ 矩阵，
$I S$、$S^\top I S$ 这些运算会白白做大量乘 0。

Pinocchio 的做法是给每类关节一个**专用 $S$ 类型**（如 `JointMotionSubspaceRevoluteTpl`），
并重载其上的运算，使乘法在编译期退化成"取一列 / 取一个元素"：

| 表达式 | 稠密写法 | 定轴转动关节的实际开销 |
|---|---|---|
| $S^\top f$ | $6n_v$ 次乘加 | **取 $f$ 的一个分量** |
| $I S$ | $36 n_v$ 次乘加 | **取 $I$ 的一列** |
| $S^\top I S$ | 更多 | **取 $I$ 的一个对角元** |

这套"把轴编码进类型"的手法与 [空间代数 §2](空间代数运算解析.md) 的 `CartesianAxis`/`SpatialAxis` 一脉相承，
是 Pinocchio 单关节开销极低的微观原因。

### 6.4 `JointMotionSubspaceBase` 与通用实现

- `joint-motion-subspace-base.hxx`（191 行）：CRTP 基类，规定 $S$ 必须提供
  `matrix()`、`motionAction()`、`se3Action()`、以及与 Force/Inertia 相乘的接口。
- `joint-motion-subspace-generic.hxx`：`JointMotionSubspaceTpl` —— **稠密**的通用 $S$，
  用于复合关节等无法静态确定稀疏模式的场合，或运行时才定维度的情形。

---

## 7. `calc` / `calc_aba`：关节的核心计算语义

### 7.1 `calc`：填充 `JointData`

每种关节必须实现三个重载（签名见 §5.5）：

```cpp
void calc(JointDataDerived& data, const ConfigVector& qs) const;              // 零阶
void calc(JointDataDerived& data, const Blank, const TangentVector& vs) const;// 只更新速度
void calc(JointDataDerived& data, const ConfigVector& qs,
                                   const TangentVector& vs) const;            // 一阶
```

**做的事**：从**全局** `q`/`v` 里取出本关节那几维（用段选择器），
算出 $M_J(q_i)$、$S_i$、$v_{J_i}$、$c_{J_i}$，写进 `data`。

以定轴转动关节为例（`joint-revolute.hxx`）：

```cpp
data.joint_q[0] = qs[idx_q()];       // 从全局 q 取出本关节的角度
SINCOS(data.joint_q[0], &sa, &ca);   // 一次算出 sin/cos
data.M.setValues(sa, ca);            // 填充绕轴转 q 的 SE3
```

`Blank` 重载的存在意义：让通用代码统一写成 `calc(jdata, q_或_Blank, v)`，
传 `Blank` 即表示"跳过位形、只更新速度"——它是个空结构体标签（见 `src/fwd.hxx`）。

### 7.2 `calc_aba`：铰接体惯量的关节内积

```cpp
void calc_aba(JointDataDerived& data, const VectorLike& armature,
              const Matrix6Like& I, const bool update_I) const;
```

对应 ABA 第二趟递推的 Schur 补（见 [GUIDE §4.4.3](../PINOCCHIO_GUIDE.md#443-递推公式的推导)）：

$$U = I^A S,\qquad D = S^\top I^A S + \text{armature},\qquad
I^A_{\text{父}} \mathrel{-}= U D^{-1} U^\top$$

**定轴转动关节的实现**（`joint-revolute.hxx`）把这套公式压到了极致：

```cpp
data.U = I.col(Inertia::ANGULAR + axis);                       // I·S → 取一列，0 次乘法
data.Dinv[0] = 1 / (I(ANGULAR+axis, ANGULAR+axis) + armature[0]);  // SᵀIS → 取一个对角元
data.UDinv.noalias() = data.U * data.Dinv[0];
if (update_I)
  I.const_cast_derived().noalias() -= data.UDinv * data.U.transpose();  // Schur 补
```

因为 $S = e_{\text{ANGULAR}+axis}$ 是**单位列向量**，$IS$ 就是取 $I$ 的第 $(\text{ANGULAR}+axis)$ 列，
$S^\top IS$ 就是取那个对角元 —— §6.3 所说的收益在此具体兑现。

**`armature` 参数**：加在 $D$ 的对角上，即电枢（转子）惯量。
高减速比关节的 $n^2I_{\text{rotor}}$ 可能与连杆惯量同量级，忽略会让模型明显偏软。

**`update_I` 参数**：是否就地把 Schur 补扣除写回 `I`。ABA 反向递推时为 `true`
（结果要传给父关节）；只想取 $U,D$ 时传 `false`。

---

## 8. 全关节清单（NQ/NV/流形/S）

下表的 `nq`/`nv` **均由本仓库实测得出**（构造各关节对象后读 `nq()`/`nv()`）：

| 关节类型 | `nq` | `nv` | 位形流形 | 说明 |
|---|---|---|---|---|
| `JointModelRX/RY/RZ` | 1 | 1 | $\mathbb{R}$ | 定轴转动。$S$ 为单位列，开销最低 |
| `JointModelRevoluteUnaligned` | 1 | 1 | $\mathbb{R}$ | **任意轴**转动，轴向运行时给定 |
| `JointModelRUBX/RUBY/RUBZ` | **2** | 1 | $SO(2)$ | **无界**转动，$q$ 存 $(\cos\theta,\sin\theta)$ |
| `JointModelRevoluteUnboundedUnaligned` | **2** | 1 | $SO(2)$ | 无界 + 任意轴 |
| `JointModelPX/PY/PZ` | 1 | 1 | $\mathbb{R}$ | 定轴移动 |
| `JointModelPrismaticUnaligned` | 1 | 1 | $\mathbb{R}$ | 任意轴移动 |
| `JointModelHX/HY/HZ` | 1 | 1 | $\mathbb{R}$ | **螺旋副**：转与移按导程耦合 |
| `JointModelHelicalUnaligned` | 1 | 1 | $\mathbb{R}$ | 任意轴螺旋 |
| `JointModelSpherical` | **4** | 3 | $SO(3)$ | 球副，$q$ 存**四元数** |
| `JointModelSphericalZYX` | **3** | 3 | $\mathbb{R}^3$ | 球副，$q$ 存 **Z-Y-X 欧拉角**（有万向锁！） |
| `JointModelTranslation` | 3 | 3 | $\mathbb{R}^3$ | 三维平移 |
| `JointModelPlanar` | **4** | 3 | $SE(2)$ | 平面副，$q$ 存 $(x,y,\cos\theta,\sin\theta)$ |
| `JointModelFreeFlyer` | **7** | 6 | $SE(3)$ | **浮动基**，$q$ 存 $(位置_3,四元数_4)$ |
| `JointModelUniversal` | 2 | 2 | $\mathbb{R}^2$ | 万向节（两正交轴） |
| `JointModelEllipsoid` | — | — | — | 椭球约束关节 |
| `JointModelComposite` | 累加 | 累加 | 笛卡尔积 | 见 §9.1 |
| `JointModelMimic` | 0 | 0 | — | 见 §9.2 |

### 8.1 三组 $n_q > n_v$ 的关节及其原因

这是 $n_q \neq n_v$ 的**全部来源**，值得逐一看清：

| 关节 | $n_q$ vs $n_v$ | 多出来的那一维用来干嘛 |
|---|---|---|
| `RUB*`（无界转动） | 2 vs 1 | 存 $(\cos\theta,\sin\theta)$ 而非 $\theta$ —— 这样转过 $\pm\pi$ 不会跳变，可无限旋转（轮子、传送带） |
| `Spherical` | 4 vs 3 | 四元数的**归一化约束**吃掉一维 |
| `Planar` | 4 vs 3 | 同上，姿态部分用 $(\cos,\sin)$ |
| `FreeFlyer` | 7 vs 6 | 同上，四元数 4 维对应 3 维角速度 |

> **反例**：`SphericalZYX` 是 $n_q = n_v = 3$，因为它用**欧拉角**参数化——代价是存在**万向锁**
> （$\text{pitch}=\pm90°$ 时丢一个自由度）。`Spherical` 用四元数则无此问题。
> 需要大范围姿态运动时优先选 `Spherical`。

> 这也解释了为何 `q += v*dt` 是错的、必须用 `integrate()`：对上述四类关节，
> $q$ 与 $v$ 维度都不同，逐元素相加根本无从谈起（见 §10）。

---

## 9. 复合关节：Composite / Mimic / Unaligned / Unbounded

### 9.1 Composite（`joint-composite.hxx`）

把**多个关节串联**成一个逻辑关节：$n_q$、$n_v$ 为各子关节之和，
$S$ 为各子关节 $S$ 经相对变换后拼成的 $6\times n_v$ 稠密矩阵，位形流形是各子流形的**笛卡尔积**。

**用途**：URDF 里用几个单自由度关节串出的复合机构（如 RPY 三连转），
合成一个关节后减少树节点数；也用于自定义关节原型。

**代价**：$S$ 变稠密，失去 §6.3 的稀疏红利，
故仅在确有需要时使用。其 `calc` 内部对子关节逐个调用并做变换复合。

### 9.2 Mimic（`joint-mimic.hxx`）与 `nvExtended` 的完整机制

把一个关节的位形**线性绑定**到另一个关节：

$$q_m = k\,q_p + b,\qquad \dot q_m = k\,\dot q_p$$

（$k=$ `scaling`，$b=$ `offset`，下标 $m$ = mimicking，$p$ = mimicked/primary。）

Model 里用 `mimicking_joints` / `mimicked_joints` 记录绑定关系
（实测 Baxter：`mimicking=[10,19]`、`mimicked=[9,18]`）。
**典型场景**：夹爪两指齿轮耦合、并联传动、差速驱动。

#### 9.2.1 三个维度的取值

`JointModelMimic` 的三个维度（`joint-mimic.hxx:614-623`）：

```cpp
inline int nq_impl()         const { return 0; }            // 不占 q
inline int nv_impl()         const { return 0; }            // 不占 v
inline int nvExtended_impl() const { return m_nvExtended; } // 但占雅可比一列
```

物理含义：mimic 关节**不是**一个新自由度（所以 `nq=nv=0`），
但它**确确实实是树上一个会转的铰**，有自己的 $S_m$、自己的 $^0X_m$，
前向递推必须给这一列留位置。这块位置就是 `idx_vExtended`。

配套地，`setMimicIndexes`（`joint-mimic.hxx:651-661`）把 `idx_q/idx_v` 直接指向**被模仿关节**
的那一段——注意 `setIndexes_impl` 刻意**不覆盖** `i_q/i_v`，只写 `i_vExtended`。
于是 `calc` 里 `qs.segment(Base::i_q, m_nqExtended)` 读到的就是主关节的位形。

#### 9.2.2 例子：耦合夹爪

| id | 关节 | `nq` | `nv` | `nvExtended` | `idx_q` | `idx_v` | `idx_vExtended` |
|---|---|---|---|---|---|---|---|
| 0 | universe | 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | finger_left (revolute) | 1 | 1 | 1 | 0 | 0 | 0 |
| 2 | finger_right (mimic of 1, $k=-1$) | 0 | 0 | 1 | **0** | **0** | **1** |

$\Rightarrow$ `model.nq = 1`，`model.nv = 1`，`model.nvExtended = 2`。

注意 id=2 那行：`idx_q/idx_v` 回指到 0（读同一个 $q$ 分量），`idx_vExtended` 是自己独有的 1。

#### 9.2.3 数学表述

设扩展雅可比 $J_{\text{ext}}\in\mathbb{R}^{6\times n_v^{\text{ext}}}$，每一列是一个关节的运动子空间
在世界系下的表达：

$$J_{\text{ext}}=\begin{bmatrix} {}^{0}X_{1}S_{1} & \cdots & {}^{0}X_{i}S_{i}\end{bmatrix}$$

这正是 `data.hxx:361-369` 注释里写的东西，也是那里特别强调
"This Jacobian has no special meaning" 的原因。

耦合关系是一个常值线性映射 $G\in\mathbb{R}^{n_v^{\text{ext}}\times n_v}$：

$$\dot q_{\text{ext}} = G\,\dot q,\qquad J = J_{\text{ext}}\,G$$

夹爪例子里 $G=\begin{bmatrix}1\\-1\end{bmatrix}$。动力学侧对应
$\tau_{\text{mimic}}=G^\top\tau_{\text{full}}$、$M_{\text{mimic}}=G^\top M_{\text{full}}G$
（见 `examples/mimic_dynamics.py`，该例还用 `model == model2` 验证手工
`transformJointIntoMimic` 与 URDF 解析结果一致）。

**关键实现细节**：Pinocchio 不显式存 $G$。缩放系数 $k$ 被提前折进了运动子空间本身——
`ScaledJointMotionSubspaceTpl::matrix_impl` 直接返回 $k\,S$（`joint-mimic.hxx:182`）。
所以 $J_{\text{ext}}$ 的第 2 列已经是 $-\,{}^{0}X_2S_2$，$G$ 退化成纯 0/1 的"求和模板"，
折叠时就是**列相加**。

#### 9.2.4 落地：`Data.J` 是 `nvExtended` 列

```cpp
, J(Matrix6x::Zero(6, model.nvExtended))     // data.hxx:728
, dJ(Matrix6x::Zero(6, model.nvExtended))
, ddJ(Matrix6x::Zero(6, model.nvExtended))
```

`getJointJacobian` 的输入输出维度差异写得非常明白（`jacobian.hxx:188-191`）：

```cpp
PINOCCHIO_CHECK_ARGUMENT_SIZE(Jin.cols(),  model.nvExtended);  // 内部缓冲
PINOCCHIO_CHECK_ARGUMENT_SIZE(Jout.cols(), model.nv);          // 交给用户的
```

折叠动作是 `jacobian.hxx:219-235` 的**两趟循环**：

```cpp
// 第一趟：非 mimic 的支撑链，赋值
for (jExtended = colRef; jExtended >= 0;
     jExtended = data.non_mimic_parents_fromRow[jExtended])
    v_out = v_in;      // v_out 绑定 Jout.col(idx_vExtended_to_idx_v_fromRow[jExtended])

// 第二趟：mimic 链，累加到被模仿列上
for (jExtended = colRefMimicPass; jExtended >= 0;
     jExtended = data.mimic_parents_fromRow[jExtended])
    v_out += v_in;     // 同一个输出列，+=
```

三张查找表各司其职（`data.hxx:322-351`，在 `DataTpl` 构造函数 `data.hxx:891-950` 里填好）：

| 表 | 含义 |
|---|---|
| `idx_vExtended_to_idx_v_fromRow[e] = v` | 扩展列 $e$ 该加到输出列 $v$ 上，即 $G$ 的稀疏表示 |
| `non_mimic_parents_fromRow[e]` | 支撑链上前一个**非 mimic** 列 |
| `mimic_parents_fromRow[e]` | 支撑链上前一个 **mimic** 列 |

后两张是把原来的 `parents_fromRow` 一分为二：一条支撑链拆成两条子链，
才能一趟 `=` 一趟 `+=`；否则先赋值的会被后来的覆盖掉。
`check-data.hxx:84-142` 校验了这三张表的长度都是 `nvExtended`，
且无 mimic 时 `parents_fromRow` 与两者一致。

CRBA 同理（`crba.hxx:186-211`）：用 `idx_vExtendeds` 取 $J_{\text{ext}}$ 的列算力，
但写回 `data.M` 时用 `idx_vs`——因为 $M\in\mathbb{R}^{n_v\times n_v}$ 始终是压缩后的。

#### 9.2.5 为什么要设计这第三套索引

因为**前向递推必须是局部的**。`forwardKinematics` / `computeJointJacobians` 沿树走一遍，
每个关节只知道自己的 $^0X_i$ 和 $S_i$，不知道自己有没有被别人 mimic、被谁 mimic。
如果直接往 $n_v$ 列的矩阵里写，两个耦合关节会抢同一列，就得在递推内部做条件判断和累加，
破坏 pass 的规整性与 SIMD 友好性。

`nvExtended` 的做法是：**递推阶段一律按"每个铰一列"无脑写，耦合的归约推迟到取用阶段**
（`getJointJacobian` 的那两趟循环）。代价是 $J$ 多占几列内存，
收益是核心 RNEA/CRBA/ABA 的 pass 结构对 mimic 完全无感知。

#### 9.2.6 三句话总结

1. `nv` 是**求解器看到的**维数（$M$、$\tau$、$\dot q$）；`nvExtended` 是**几何递推看到的**维数（$J$、$\dot J$、$\ddot J$ 的列数）。
2. 没有 mimic 关节时两者恒等，`idx_vs == idx_vExtendeds`，这套机制零开销地消失。
3. 遇到 `data.J` 时永远别假设它有 `nv` 列——要拿有物理意义的雅可比，
   必须走 `getJointJacobian` / `getFrameJacobian`，它们负责把 $J_{\text{ext}}$ 乘上 $G$。

### 9.3 Unaligned 系列

`RevoluteUnaligned` / `PrismaticUnaligned` / `HelicalUnaligned`：
轴向不是 X/Y/Z 而是**运行时给定的任意单位向量** $n$。

**代价**：定轴版本的 $S$ 是编译期已知的单位列（$IS$ = 取一列）；
Unaligned 版本的 $S$ 依赖运行时的 $n$，只能做真正的矩阵向量乘。
**所以 URDF 里能对齐坐标轴时，让轴与 X/Y/Z 对齐是有性能意义的**。

### 9.4 Unbounded 系列

`RUBX/RUBY/RUBZ` / `RevoluteUnboundedUnaligned`：位形存 $(\cos\theta,\sin\theta)$（$n_q=2$），
流形是 $SO(2)$ 而非 $\mathbb{R}$。

**解决什么问题**：普通转动关节的 $\theta$ 是实数，转到 $\pm\pi$ 会遇到限位或数值跳变；
无界版本在圆周上积分，可无限旋转且无跳变。**轮子、转台、连续旋转关节应当用它**。

### 9.5 Ellipsoid（`joint-ellipsoid.hxx`）

椭球约束关节，用于把运动限制在椭球面上的特殊机构建模（较少用）。

---

## 10. 李群系统：位形流形的加减法

### 10.1 为什么需要它

§8.1 说明了 $n_q\neq n_v$ 的四类来源。
一旦维度都对不上，`q + v*dt` 就无从谈起；即便维度相同（如 `SphericalZYX`），
逐元素相加也会离开流形（欧拉角相加不等于旋转复合）。

**李群系统的职责就是提供流形上正确的"加减法"**：

$$q \oplus v \;\equiv\; \texttt{integrate}(q, v), \qquad
q_1 \ominus q_0 \;\equiv\; \texttt{difference}(q_0, q_1)$$

> **实测**（样例人形，34 DoF）：
> `difference(q0, integrate(q0,v)) − v` 残差 $4.6\times10^{-16}$；
> `integrate(q0, difference(q0,q1))` 与 $q_1$ 相差 $6.3\times10^{-16}$ —— 二者严格互逆。

### 10.2 `LieGroupBase<Derived>` 完整 API（`liegroup-base.hxx`）

**① 基本运算**

| 方法 | 数学 | 说明 |
|---|---|---|
| `integrate(q, v, qout)` | $q\oplus v$ | 流形上"加法"，内部走各群的 exp |
| `difference(q0, q1, vout)` | $q_1\ominus q_0$ | 流形上"减法"，内部走 log |
| `interpolate(q0, q1, u, qout)` | 测地线插值 | $u{=}0\to q_0$、$u{=}1\to q_1$（实测端点残差 $0$ / $6.3\times10^{-16}$） |
| `distance(q0, q1)` / `squaredDistance` | $\lVert q_1\ominus q_0\rVert$ | 实测与 `‖difference‖` 一致（$1.8\times10^{-15}$） |
| `random(qout)` / `randomConfiguration(lb, ub, qout)` | — | 流形上均匀采样 / 限位内采样 |
| `normalize(qout)` | — | 投影回流形（四元数除以模长） |
| `isNormalized(q, prec)` | — | 是否在流形上 |
| `isSameConfiguration(q0, q1, prec)` | — | 两位形是否等价（注意四元数双覆盖） |

> ⚠️ **Python 绑定的 `normalize` 不就地修改**：C++ 的 `normalize(qout)` 写入参数，
> 而 `pin.normalize(model, q)` **返回**归一化后的新数组、原 `q` 不变。实测：
> 修改四元数使 $\lvert q\rvert=1.05$ 后调用，原数组仍为 1.05，返回值才是 1.0。
> 写成 `q = pin.normalize(model, q)` 才有效。

**② 一阶导数（轨迹优化必需）**

| 方法 | 含义 |
|---|---|
| `dIntegrate(q, v, J, arg)` | $\partial(q\oplus v)/\partial q$ 或 $/\partial v$，由 `arg` 选（`ARG0`/`ARG1`） |
| `dIntegrate_dq` / `dIntegrate_dv` | 上面两个的显式命名版本 |
| `dDifference(q0, q1, J, arg)` | $\partial(q_1\ominus q_0)/\partial q_0$ 或 $/\partial q_1$ |
| `integrateCoeffWiseJacobian(q, J)` | $\partial q/\partial(\text{参数})$ 的逐系数版（$n_q\times n_v$） |

> 实测 `dIntegrate` 返回 $34\times34$（即 $n_v\times n_v$）。
> 这些是 DDP/iLQR 在浮动基上做状态扰动时的必需件（见 [GUIDE §4.8.3](../PINOCCHIO_GUIDE.md#483-与-ddp-的状态空间矩阵的关系)）。

**③ 雅可比传输（transport）**

| 方法 | 用途 |
|---|---|
| `dIntegrateTransport(q, v, Jin, Jout, arg)` | 把一个定义在 $q$ 处切空间的雅可比，**搬运**到 $q\oplus v$ 处的切空间 |
| `dIntegrateTransport_dq` / `_dv` | 同上的具名版本 |

**为什么需要"传输"**：流形上不同点的切空间是**不同的向量空间**，
在 $q$ 处算出的梯度不能直接用在 $q\oplus v$ 处，必须先做平行移动。
这是李群优化与欧氏优化最容易被忽略的差别。

**④ 切空间映射**

| 方法 | 含义 |
|---|---|
| `tangentMap(q, TM)` | 切空间到位形空间增量的映射（$n_q\times n_v$） |
| `tangentMapProduct` / `tangentMapTransposeProduct` | 免构造矩阵的乘积形式 |

**⑤ 元信息**：`nq()`、`nv()`、`name()`、`neutral()`、`operator==`。

### 10.3 三种基础群 + 笛卡尔积

| 文件 | 群 | 用于哪些关节 |
|---|---|---|
| `vector-space.hxx` | $\mathbb{R}^n$ | 转动/移动/欧拉角球副等**普通**关节（默认流形） |
| `special-orthogonal.hxx` | $SO(2)$ / $SO(3)$ | 无界转动 / 球副 |
| `special-euclidean.hxx` | $SE(2)$ / $SE(3)$ | 平面副 / 浮动基 |
| `cartesian-product.hxx` | $G_1\times G_2\times\cdots$ | **整机位形空间** |

**整机流形是各关节流形的笛卡尔积**。例如人形：

$$\mathcal{Q} = \underbrace{SE(3)}_{\text{浮动基}} \times \underbrace{\mathbb{R}\times\cdots\times\mathbb{R}}_{\text{各转动关节}}$$

`cartesian-product-variant.hxx` 提供运行时可变长度的版本（关节数在编译期未知时用）。

### 10.4 `LieGroupMap`（`liegroup-map.hxx`）—— 关节→流形映射表（关键枢纽）

一个**编译期查表**结构，回答"某关节类型对应哪个流形"：

```cpp
template<typename JointModel> struct operation { typedef ... type; };
// JointModelRX          → VectorSpaceOperation<1>
// JointModelSpherical   → SpecialOrthogonalOperation<3>
// JointModelFreeFlyer   → SpecialEuclideanOperation<3>
// JointModelRUBX        → SpecialOrthogonalOperation<2>
```

`JointModelBase::lieGroup()` 就是查这张表返回对应的群对象。

> **它是整个设计的枢纽**：新增一种关节时，只要在这张表里登记它的流形，
> `integrate`/`difference`/`randomConfiguration` 等**整机级**操作就自动支持了新关节，
> 无需改动任何算法代码。

### 10.5 其余李群文件

| 文件 | 作用 |
|---|---|
| `liegroup-collection.hxx` | 李群"菜单" → `LieGroupGenericTpl` variant |
| `liegroup-generic.hxx` | 变体包装器（与 `JointModelTpl` 同构的设计） |
| `liegroup-variant-visitors.hxx` | 对李群 variant 的访问者入口 |
| `liegroup-algo.hxx` | **沿整棵树**逐关节应用李群操作 —— `pin.integrate(model,q,v)` 的实现落点 |
| `liegroup-joint.hxx` | 关节与其李群的桥接 |
| `fwd.hxx` | 前置声明 |

> 调用链：`pin::integrate(model,q,v)` → `liegroup-algo.hxx` 的树遍历
> → 每个关节 `jmodel.lieGroup()` 查 `LieGroupMap`
> → 调用该群的 `integrate` → 写回 `q` 的对应段。

---

## 11. 访问者系统：把 variant 变成可调用的算法

`model.joints[i]` 是 `boost::variant`，运行时可能是任意一种关节。
访问者负责**在不用虚函数的前提下**跳到正确的具体实现。

### 11.1 `JointUnaryVisitorBase`（`visitor/joint-unary-visitor.hxx`，328 行）

对**单个**关节 variant 分派。用法三要素（CRTP）：

```cpp
struct MyStep : public fusion::JointUnaryVisitorBase<MyStep>   // ① CRTP 自引用
{
  typedef boost::fusion::vector<const Model&, Data&, ...> ArgsType;  // ② 透传参数包

  template<typename JointModel>                                       // ③ 对每种关节各实例化一份
  static void algo(const JointModelBase<JointModel>& jmodel,
                   JointDataBase<typename JointModel::JointDataDerived>& jdata,
                   const Model& model, Data& data, ...)
  { /* 这里 jmodel 已是具体类型，静态分派、可内联 */ }
};
// 调用
MyStep::run(model.joints[i], data.joints[i], MyStep::ArgsType(model, data, ...));
```

**分派链路**：`run` → `boost::apply_visitor` → 编译期已生成的分支表 → `algo<JointModelRX>` →
`jmodel.calc(...)` → `joint-revolute.hxx` 的 `calc`。
运行时只是"选分支"，**没有虚表查找**。

> 有多个重载：带/不带 `JointData`、带/不带返回值。
> 手写一个访问者并与官方 `forwardKinematics` 比对的可运行例子见
> [`doc/architecture-tutorial/03-visitor/`](../doc/architecture-tutorial/03-visitor/main.cpp)。

### 11.2 `JointBinaryVisitorBase`（`joint-binary-visitor.hxx`，234 行）

对**两个**关节 variant 同时分派（需要 $N^2$ 个分支组合）。
用于涉及关节对的运算，如某些约束/耦合处理。

### 11.3 `visitor/fusion.hxx`

只做一件事：`bf::append` —— 把可变参数包压进 Boost.Fusion 序列，
使 `ArgsType` 能承载任意数量、任意类型的透传参数。

> **读算法的固定套路**：打开任意 `algorithm/*.hxx`，数一数里面有几个
> `struct ...Step : JointUnaryVisitorBase<...>`，就知道该算法有几趟递推 ——
> RNEA 两个（对应 [GUIDE §4.3](../PINOCCHIO_GUIDE.md#43-逆向动力学-rnea) 的正/反两趟）、
> ABA 三个、FK 一个。**公式与 Step 一一对应**。

---

## 12. 其余文件：pool / sample-models / force-set

### 12.1 `pool/model.hxx`、`pool/geometry.hxx`：并行计算的对象池

**要解决的问题**：`Model` 只读、可多线程共享，但 `Data` 是可写的，
每个线程必须有**自己的一份**。手工管理容易出错。

`ModelPoolTpl` 就是这个管理器：

| 方法 | 说明 |
|---|---|
| `ModelPoolTpl(model, pool_size)` | 构造：内部复制 `pool_size` 份 `Data` |
| `size()` / `resize(n)` | 池容量查询 / 调整（线程数变化时） |
| `getModel(i)` / `getModels()` | 取第 $i$ 个（或全部）`Model` |
| `getData(i)` / `getDatas()` | 取第 $i$ 个（或全部）`Data` —— **每线程用自己的索引** |
| `update(data)` | 用给定 `Data` 刷新池中所有副本 |

`GeometryPoolTpl` 是它对几何模型的对应物（并行碰撞检测用）。

**上层用法**见 [`examples/run-algo-in-parallel.py`](../examples/run-algo-in-parallel.py)：
`pin.ModelPool(model)` + `rneaInParallel(num_threads, pool, q, v, a, res)`，
对 $(n_q\times B)$ 的一批配置并行跑 RNEA。

### 12.2 `sample-models.hxx`：内置测试模型

| 函数 | 产出 |
|---|---|
| `buildModels::manipulator(model)` | 6-DoF 串联机械臂（固定基） |
| `buildModels::humanoid(model)` | 人形（浮动基） |
| `buildModels::humanoidRandom(model)` | 人形，惯量参数随机 |
| `buildModels::manipulatorGeometries(...)` 等 | 配套几何模型 |

**价值**：不依赖任何 URDF 文件即可跑通全部算法 —— 单元测试、基准测试、
以及本文档所有数值验证都用它（实测样例人形：`nq=35, nv=34, njoints=30, nframes=70`）。

### 12.3 `force-set.hxx`

对**一批** Force / Motion（存成 $6\times N$ 矩阵，每列一个）做批量变换的辅助。
与 [空间代数 §act-on-set](空间代数运算解析.md) 是同一套思路：
把公共部分（取 $R$、算 $\hat t$）提到循环外，让 Eigen 对整块矩阵向量化，
远快于逐列调用。雅可比的参考系转换即依赖它。

---

## 13. 一次 `forwardKinematics` 如何串起整个子系统

把上述所有部件连起来看一遍数据流（`algorithm/kinematics` 调 `multibody/` 的组件）：

```
输入：Model（静态树）+ q（位形向量）
for i = 1 .. njoints（深度优先）:
    jmodel = model.joints[i]           # JointModelVariant（类型擦除）
    jdata  = data.joints[i]            # JointDataVariant

    calc(jmodel, jdata, q)             # ── 访问者恢复真实关节类型（§11）
      └─ 具体 calc：取 q 的本关节段（jointConfigSelector，§5.2）
         填 jdata.M = 关节运动产生的 SE3（§7）
         填 jdata.S = 运动子空间（§6）

    data.liMi[i] = model.jointPlacements[i] * jdata.M    # 父→子固定位姿 · 关节运动
    data.oMi[i]  = data.oMi[parent] · data.liMi[i]       # 世界系位姿递推（SE3 复合，见空间代数解析 §3）

# 之后 updateFramePlacements：
for each frame f:
    data.oMf[f] = data.oMi[f.parentJoint] * f.placement  # §4
```

- **SE3 复合** `oMi[parent] * liMi` 的数学在 [空间代数解析 §3.1](空间代数运算解析.md)。
- 若要速度：`calc(jmodel,jdata,q,v)` 还会填 `jdata.v = S·v̇`，然后 `data.v[i] = liMi.actInv(data.v[parent]) + jdata.v`（速度在坐标系间用伴随变换传递，见 [空间代数解析 §4.1](空间代数运算解析.md)）。
- 若要 `pinocchio::integrate(model,q,v)`（在流形上走一步）：走 §10 的 `liegroup-algo`，逐关节查 `LieGroupMap` 调对应群 exp。

**这就是 `multibody/` 的全貌**：`Model` 提供静态树 → 访问者按真实类型调 `calc` → 关节用 `S`/`M` 把广义坐标变成空间运动 → SE3 递推填满 `Data` → 上层算法消费 `Data`；而 `LieGroup` 保证一切"位形加减"发生在正确的弯空间里。

---

## 14. 速查表

| 想找… | 去哪个文件 | 关键类型/函数 |
|-------|-----------|--------------|
| 机器人的树结构、维度 | `model.hxx` | `ModelTpl`，`nq/nv`，`parents/supports/subtrees` |
| `njoints` 和 `nq` 差在哪 | `model.hxx`（§2.2.3） | $n_q=\sum_{i\ge1}\texttt{nqs}[i]$；浮动基 $n_q=\text{njoints}+5$ |
| 某关节在 q / v 里占哪几维 | `model.hxx`（§2.2） | `nqs/idx_qs`、`nvs/idx_vs`，用 `jointConfigSelector`/`jointVelocitySelector` 取 |
| `data.J` 为什么列数不等于 `nv` | `model.hxx`+`jacobian.hxx`（§9.2） | `nvExtended`、`idx_vExtendeds`、`idx_vExtended_to_idx_v_fromRow` |
| 某算法的中间量存哪 | `data.hxx` | `DataTpl` 对应字段（grep 字段名到 algorithm/） |
| 坐标系 vs 关节 | `frame.hxx`,`model-item.hxx` | `FrameTpl`，`FrameType`，`oMf=oMi·placement` |
| 加新关节类型 | `joint/joint-xxx.hxx` + `joint-collection.hxx` | 仿 revolute 写 Model/Data/Constraint/Motion/Transform |
| 关节允许的运动方向 | `joint-motion-subspace-*.hxx` | $v=S\dot q$，$\tau=S^\top f$ |
| 关节怎么把 q 变成位姿 | 各关节 `calc` | `data.M`，`data.S`，`data.v` |
| ABA 消去关节自由度 | 各关节 `calc_aba` | $U=IS,\ D=S^\top IS,\ I\!\mathrel{-}=UD^{-1}U^\top$ |
| q 怎么积分/差分 | `liegroup/*.hxx` | `integrate=q⊕v`，`difference=q1⊖q0` |
| 某关节用哪个流形 | `liegroup-map.hxx` | `LieGroupMap::operation<Joint>` |
| SO(3)/SE(3) 的 exp/log/雅可比 | `special-{orthogonal,euclidean}.hxx` → `spatial/log.hxx` | `exp3/log3/Jlog3`，`exp6/log6/Jlog6`（[空间代数解析 §7](空间代数运算解析.md)） |
| variant 上怎么调具体关节 | `visitor/*.hxx` | `apply_visitor`，unary/binary visitor |
| 多线程批量算 | `pool/model.hxx` | `ModelPoolTpl` |

> 数学细节（SE3 复合、伴随/余伴随、exp/log、惯量复合与辨识）统一在 [空间代数运算解析.md](空间代数运算解析.md)；模板机制（CRTP、variant、fusion、cast、`::template`）统一在 [源码解析.md](源码解析.md) 与 [C++语法技巧.md](C++语法技巧.md)；继承关系图在 [类结构图.md](类结构图.md)。
