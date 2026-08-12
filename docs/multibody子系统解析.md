# Pinocchio `multibody/` 子系统解析

> 本文件系统解析 `include/pinocchio/multibody/`（声明）与 `include/pinocchio/src/multibody/`（模板实现）下**全部文件**，
> 讲清每个文件的作用、内部实现与涉及的数学。配套阅读：
> [源码解析.md](源码解析.md)（模板机制/访问者三件套）、
> [空间代数运算解析.md](空间代数运算解析.md)（SE3/Motion/Force/Inertia/exp-log 的数学）、
> [类结构图.md](类结构图.md)（UML）、[C++语法技巧.md](C++语法技巧.md)。

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

这是理解 Pinocchio 数据布局的钥匙。每个关节在三个不同维度的**全局向量**里各占一段：

| 空间 | 维度符号 | 每关节段长 / 起点 | 含义 |
|------|----------|-------------------|------|
| **位形空间** q | `nq` | `nqs[i]` / `idx_qs[i]` | 位形向量维度。可 > nv（如球副 nq=4 四元数、nv=3） |
| **速度/切空间** v | `nv` | `nvs[i]` / `idx_vs[i]` | 广义速度、力矩、加速度维度 |
| **雅可比扩展空间** | `nvExtended` | `nvExtendeds[i]` / `idx_vExtendeds[i]` | mimic 等关节雅可比列数 ≠ nv 时用 |

> **为什么 nq≠nv**：位形住在**弯曲流形**上（球副的姿态是 SO(3)，用 4 维四元数存但只有 3 个自由度），速度住在**切空间**（3 维）。二者的桥是李群 exp/log（见 §10）。

`addJoint` 里维护这些：`nq += joint_nq; idx_qs.push_back(...)` 等，并用 `jmodel.setIndexes(joint_id, nq, nv, nvExtended)` 把全局起点写回关节对象。

### 2.3 树结构字段

- `parents[i]`：关节 i 的父关节（`universe`=0 是根）。
- `children[i]`：子关节列表。
- `supports[j]`：从 `universe` 到 j 的路径（含两端）——雅可比稀疏性、CRBA 的基础。
- `subtrees[j]`：j 支撑的整棵子树（含 j 自己）——CRBA 复合惯量累加范围。
- `sparsity_pattern_vector[i]` / `span_indexes_vector[i]`：关节 i 对应雅可比**哪些列非零**（布尔向量 / 索引列表）。
- `jointPlacements[i]`：关节 i 相对父关节坐标系的固定位姿 `^{parent}M_i`（$X_T$，见 §13）。
- `inertias[i]`：关节 i 支撑的刚体空间惯量（`appendBodyToJoint` 累加进来）。
- `mimicking_joints`/`mimicked_joints`/`mimic_joint_supports`：镜像关节的绑定关系。

### 2.4 关键方法

- **`addJoint(parent, jmodel, placement, name, ...)`**：向树追加关节。**必须深度优先顺序**添加。内部：`njoints++` → `joints.push_back(JointModel(jmodel.derived()))` → `setIndexes` → 累加 nq/nv → `conservativeResize` 各种极限向量并写入 → 更新 `subtrees`/`supports`/稀疏模式 → 若是 mimic 记录绑定。有多个重载（是否给力/速/位极限、摩擦、阻尼），最终都汇聚到那个 13 参数的总实现。
- **`appendBodyToJoint(joint_id, Y, placement)`**：把刚体惯量 `Y`（换算到关节系后）加到 `inertias[joint_id]`。
- **`addFrame(frame, append_inertia)`**：注册操作坐标系；若带惯量且 `append_inertia`，惯量并入父关节。
- **`createData()`**：按当前 Model 分配一个匹配的 `Data`。
- **`cast<NewScalar>()`**：逐字段 `.cast<>()`，把整个模型换标量类型（float↔double↔ADScalar，autodiff 用；见 [源码解析.md](源码解析.md) 的 cast 条目）。
- **`check(checker)`**：模板化的模型合法性校验钩子。

`gravity`（`Motion`，默认 `(0,0,-9.81)` 线性部分）与静态成员 `gravity981` 也在此定义。

---

## 3. `DataTpl`：算法工作区

`src/multibody/data.hxx`。如果说 `Model` 是"机器人的图纸"，`Data` 就是**"算草稿纸"**——~150 个预分配的字段，让 RNEA/CRBA/ABA 等算法**零动态分配**地反复运行。字段按算法族分组：

### 3.1 运动学量（每关节一个，`std::vector` 长度 = njoints）
- `oMi[i]`：关节 i 在世界系的位姿 `^0M_i`（forwardKinematics 主输出）。
- `liMi[i]`：关节 i 相对父关节的位姿 `^{parent}M_i`（含关节运动）。
- `v[i]`/`a[i]`：关节 i 的空间速度/加速度（LOCAL）；`ov[i]`/`oa[i]`：世界系版本。
- `oMf[i]`：操作坐标系在世界系的位姿（updateFramePlacements 输出）。

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
- `Yaba[i]`/`oYaba[i]`：铰接体惯量（Articulated Body Inertia）。
- `u[i]`/`U`/`D`/`Dinv`：ABA 三趟递推的偏置力与 $U D^{-1} U^\top$ 分解量。
- `ddq`：关节加速度输出。

### 3.5 稀疏 Cholesky（`M = U D Uᵀ`）
- `U`（单位上三角）/`D`/`Dinv`/`tmp`；`parents_fromRow`/`nvSubtree_fromRow`/`supports_fromRow` 是**按行（自由度）而非按关节**重排的树结构，供稀疏分解按 DoF 遍历。

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

> 记忆法：**`Data` 的每个字段几乎都精确对应某个算法的某个中间量**。看到陌生字段，去 `algorithm/` 里 grep 它的名字即可定位用途。

---

## 4. `FrameTpl` / `ModelItem`：坐标系与树节点基类

### 4.1 `ModelItem`（`model-item.hxx`）
`Frame`（及未来的其他树挂件）的公共基类，四个字段：
`name`、`parentJoint`（挂在哪个关节）、`parentFrame`（父坐标系，多为文档/第三方用）、`placement`（相对父关节系的固定位姿 SE3）。

### 4.2 `FrameTpl`（`frame.hxx`）
在 `ModelItem` 基础上加 `FrameType type` 和 `Inertia inertia`。

**Frame vs Joint 的区别**（常见困惑，另见 [源码解析.md](源码解析.md) §3.1）：
- **Joint** 是运动学树的**真实节点**，有自由度、参与所有动力学递推；`oMi[i]` 是它的位姿。
- **Frame** 是**附着在某个关节上的固定标记**（`placement` 是常量），本身**无自由度**、不参与递推；它的位姿 `oMf = oMi[parentJoint] * placement` 是**事后**由关节位姿算出来的。
- `FrameType`：`OP_FRAME`（用户操作系）、`JOINT`（关节系的冗余镜像）、`FIXED_JOINT`（URDF 里被折叠掉的固定关节——**Pinocchio 不把固定关节放进运动学树，而是降级成 Frame**）、`BODY`（连杆的惯量/视觉/碰撞系）、`SENSOR`。

这正是 URDF 里"很多 link/固定 joint"被压缩成"少数活动关节 + 一堆 Frame"的机制。

---

## 5. 关节系统总览：Model/Data 分离 + CRTP + variant 双层

关节是本子系统最精巧的部分。它要同时满足三个矛盾需求：**(a)** 每种关节数学不同（转/移/球/浮动…）；**(b)** 一棵树里混装不同关节，要能存进一个 `std::vector`；**(c)** 递推热循环里不能有虚函数开销。解法是**四个正交设计**叠加：

### 5.1 Model 对象 vs Data 对象（职责分离）
每种关节都拆成两个类：
- **`JointModelXxx`**（住在 `Model.joints`）：**静态参数 + 计算逻辑**。存关节轴、索引 `idx_q/idx_v`，提供 `calc`/`calc_aba`。**无状态**（不随 q 变）。
- **`JointDataXxx`**（住在 `Data.joints`）：**计算缓存**。存当前 `q,v` 下算出的 `M`（关节变换）、`S`（运动子空间）、`v`（关节速度）、`c`（偏置）、ABA 量 `U/Dinv/UDinv`。

`calc(jdata, q)` 就是"给定 Model 的静态参数 + 输入 q，填充 Data 的缓存"。

### 5.2 CRTP 基类（静态多态，零开销）
`JointModelBase<Derived>`（`joint-model-base.hxx`）与 `JointDataBase<Derived>`。基类通过 `derived()` 静态转发到具体实现，没有虚函数。基类还提供一大套**段/列/块选择器**：
- `jointConfigSelector(q)` → q 里属于本关节的那 `nq` 段；
- `jointVelocitySelector(v)`、`jointCols(J)`、`jointBlock(M)` … 用 `SizeDepType<NV>` 在**编译期已知固定尺寸时**返回定长块（更快），动态时退化为运行时块。

`PINOCCHIO_JOINT_TYPEDEF_TEMPLATE` 宏把 `traits<Joint>` 里的 `Constraint_t/Transformation_t/Motion_t/Bias_t/U_t/D_t/...` 一次性拉进作用域，避免每个关节手写十几行 typedef。

### 5.3 `JointCollectionDefaultTpl` → `boost::variant`（类型擦除）
`joint-collection.hxx` 把所有具体关节类型列进一份"菜单"，末尾聚成：
```cpp
typedef boost::variant<JointModelRX, JointModelRY, ..., JointModelFreeFlyer, ...,
    boost::recursive_wrapper<JointModelComposite>,
    boost::recursive_wrapper<JointModelMimic>> JointModelVariant;
```
这样 `std::vector<JointModelVariant>` 就能混装任意关节。Composite/Mimic 用 `recursive_wrapper`（它们内部又含关节，类型递归）。

### 5.4 `JointModelTpl`（`joint-generic.hxx`）：变体包装器
`JointModelTpl` **同时**继承 `JointModelBase<JointModelTpl>`（CRTP，对外像个普通关节）和持有 `JointModelVariant`（对内是任意具体关节）。它把 `calc` 等调用通过 `apply_visitor` 转发到 variant 里真正的类型。`Model.joints` 的元素类型就是它。

> 完整机制（apply_visitor + fusion 打包 + `::template` 消歧）已在 [源码解析.md](源码解析.md) 详解，这里不重复。

---

## 6. `JointMotionSubspace`（运动子空间 S）

`joint-motion-subspace-{base,generic}.hxx`。这是关节**运动学的数学核心**。

### 6.1 定义
一个关节把 `nv` 维广义速度 $\dot q_j$ 映射成刚体的 6 维空间速度（twist）：

$$v_{\text{joint}} = S\,\dot q_j,\qquad S \in \mathbb{R}^{6\times n_v}$$

$S$ 就是**运动子空间矩阵**（constraint / motion subspace）。它的列张成"这个关节允许的瞬时运动方向"。对偶地，关节能传递的力/力矩是 $\tau_j = S^\top f$（`S.transpose() * force`）。

### 6.2 例子
- **绕 X 轴转动 RX**：$S = e_{\text{ANGULAR}+0}=[0,0,0,1,0,0]^\top$（只有绕 X 的角速度）。所以 `data.v.angularRate() = q̇`（见 §7 revolute 的 calc）。
- **沿 Z 移动 PZ**：$S=[0,0,1,0,0,0]^\top$。
- **球副**：$S=\begin{bmatrix}0_{3}\\ I_3\end{bmatrix}$（3 个角速度自由度）。
- **自由浮动基 FreeFlyer**：$S=I_6$（6 维全通）。

Pinocchio 为定轴关节写了**专用轻量 S 类型**（如 `JointMotionSubspaceRevoluteTpl`，根本不存那 6 个数，`S*v̇`/`Sᵀf` 直接取对应分量），只有通用/复合关节才用 `joint-motion-subspace-generic.hxx` 里那个**稠密 6×Dim** 的 `JointMotionSubspaceTpl`。

### 6.3 `JointMotionSubspaceTpl`（通用稠密版）实现要点
- 存 `DenseBase S`（`6×Dim`，`Dim` 可为 `Eigen::Dynamic`，带 `MaxDim` 上界避免堆分配）。
- `__mult__(v̇)` = `S*v̇` → `JointMotion`。
- 内嵌 `Transpose` 代理：`transpose()*f` = `Sᵀf`（力→关节力矩）、`transpose()*F`（对力集合批量）。
- `se3Action(M)` = $^{B}X_A$ 变换 S（把子空间搬到另一坐标系，`motionSet::se3Action`）。
- `motionAction(v)` = $v\times S$（用于偏置/科氏项）。
- `operator*(Inertia Y, S)` = $Y S$（`6×Dim`，CRBA/ABA 里 $I S$）。
- `StDiagonalMatrixSOperation` = $S^\top S$（约束正规化）。

这些操作全部走 `MotionAlgebraAction`/`SE3GroupAction`/`ConstraintForceOp` 等 traits 计算返回类型，保证表达式模板不产生临时。

---

## 7. `calc` / `calc_aba`：关节的核心计算语义

以 revolute（`joint-revolute.hxx`）为原型，`JointModelRevoluteTpl::calc` 有三个重载：

```cpp
// (1) 只更新位形相关：算关节变换 M(q)
void calc(JointData& data, q) {
  data.joint_q[0] = qs[idx_q()];
  SINCOS(data.joint_q[0], &sa, &ca);   // 一次同时算 sin/cos
  data.M.setValues(sa, ca);            // M = 绕轴转 q 的 SE3（只存 sin/cos，不建 3×3）
}
// (2) 只更新速度：v = S·q̇
void calc(JointData& data, Blank, v) {
  data.joint_v[0] = vs[idx_v()];
  data.v.angularRate() = data.joint_v[0]; // 因为 S=e_ANGULAR+axis
}
// (3) 位形+速度一起
void calc(JointData& data, q, v) { calc(data,q); ...v... }
```

`Blank`（`boost::blank`）作为**空标签**用于区分"只传 v"的重载（否则 `calc(data,v)` 和 `calc(data,q)` 签名撞车）——这是 tag dispatch，见 [C++语法技巧.md](C++语法技巧.md)。

`JointDataRevoluteTpl` 的字段正好对应递推所需：`joint_q,joint_v,S,M,v,c` + ABA 专用 `U,Dinv,UDinv,StU`。

### `calc_aba`：铰接体惯量的关节内积
```cpp
void calc_aba(data, armature, I /*6×6 铰接体惯量*/, update_I) {
  data.U = I.col(ANGULAR+axis);                        // U = I·S（S 是单位列，取一列即可）
  data.Dinv[0] = 1 / (I(ANGULAR+axis,ANGULAR+axis) + armature[0]); // D = SᵀIS + armature，求逆
  data.UDinv = data.U * data.Dinv[0];
  if (update_I) I -= data.UDinv * data.U.transpose();  // 把该关节自由度"投影消去"：I ← I − U D⁻¹ Uᵀ
}
```
这正是 ABA（Featherstone 前向动力学）里"从铰接体惯量中消去当前关节自由度"的一步。定轴关节因为 $S$ 是单位向量，$IS$ 退化成取矩阵一列、$S^\top I S$ 退化成取一个对角元——**这就是 Pinocchio 快的原因之一：为每种关节把 S 的稀疏结构写死**。`armature` 是电机转子等效惯量（加在对角上）。

---

## 8. 全关节清单（NQ/NV/流形/S）

| 关节 | 文件 | NQ | NV | 位形流形 | 运动子空间 S | 典型用途 |
|------|------|----|----|----------|-------------|----------|
| Revolute RX/RY/RZ | joint-revolute | 1 | 1 | ℝ(角度) | 单位角速度列 | 定轴铰链 |
| RevoluteUnaligned | joint-revolute-unaligned | 1 | 1 | ℝ | 任意单位轴的角速度 | 斜轴铰链 |
| RevoluteUnbounded RUBX/Y/Z | joint-revolute-unbounded | **2** | 1 | **SO(2)**（存 cosθ,sinθ） | 单位角速度 | 连续转动（无角度环绕问题） |
| RevoluteUnboundedUnaligned | joint-revolute-unbounded-unaligned | 2 | 1 | SO(2) | 任意轴角速度 | 同上，斜轴 |
| Prismatic PX/PY/PZ | joint-prismatic | 1 | 1 | ℝ | 单位线速度列 | 定轴滑动 |
| PrismaticUnaligned | joint-prismatic-unaligned | 1 | 1 | ℝ | 任意轴线速度 | 斜轴滑动 |
| Helical Hx/Hy/Hz | joint-helical | 1 | 1 | ℝ | 转+移耦合列（含螺距） | 螺旋副 |
| HelicalUnaligned | joint-helical-unaligned | 1 | 1 | ℝ | 任意轴螺旋 | 斜轴螺旋 |
| Spherical | joint-spherical | **4** | 3 | **SO(3)**（四元数） | $[0_3;I_3]$（3 角速度） | 球关节 |
| SphericalZYX | joint-spherical-ZYX | 3 | 3 | ℝ³(欧拉角) | 依赖欧拉角的 3×3 块 | 球关节（欧拉参数化） |
| Translation | joint-translation | 3 | 3 | ℝ³ | $[I_3;0_3]$（3 线速度） | 3D 平移台 |
| Planar | joint-planar | **4** | 3 | **SE(2)** | x,y 线速度 + z 角速度 | 平面移动底盘 |
| FreeFlyer | joint-free-flyer | **7** | 6 | **SE(3)**（平移+四元数） | $I_6$ | 浮动基根关节 |
| Ellipsoid | joint-ellipsoid | 3 | 3 | 椭球约束 | 依赖位形 | 椭球接触约束 |
| Universal | joint-universal | 2 | 2 | ℝ² | 两正交轴角速度 | 万向节 |
| Composite | joint-composite | 动态 | 动态 | 各子关节流形笛卡尔积 | 拼接 | 多自由度组合 |
| Mimic | joint-mimic | 动态 | 动态 | 继承被镜像关节 | 被镜像 S × 缩放 | 齿轮/耦合关节 |

> **NQ>NV 的关节**（RevoluteUnbounded、Spherical、Planar、FreeFlyer）就是位形住在弯流形、需要李群 exp/log 的那些——见 §10。

---

## 9. 复合关节：Composite / Mimic / Unaligned / Unbounded

### 9.1 Composite（`joint-composite.hxx`）
把若干关节**串联**成一个逻辑关节，NQ/NV 是各子关节之和（`Eigen::Dynamic`）。它的 `calc` 依次算各子关节的变换并复合 $M = M_1 M_2 \cdots$，S 由各子空间经 SE3 变换拼接。用于把"肩 3 转 + 肘…"打包，或造非标准关节。

### 9.2 Mimic（`joint-mimic.hxx`）
**镜像关节**：让 q 线性绑定到另一个（被镜像）关节：

$$q_{\text{mimic}} = s\cdot q_{\text{mimicked}} + o$$

（`s`=scaling、`o`=offset）。核心是 `ScaledJointMotionSubspaceTpl`：把被镜像关节的约束 $S$ 整体乘以缩放因子 $s$——`S_mimic = s·S_mimicked`，`Sᵀf = s·(S_mimickedᵀf)`。这样两个关节共享同一个自由度（齿轮组、平行四连杆）。`Model` 里用 `mimicking_joints`/`mimicked_joints`/`mimic_joint_supports` 记账，这也是 `nvExtended`（雅可比列数 ≠ nv）存在的原因。

### 9.3 Unaligned 系列
定轴版（RX/PZ…）把轴写死进模板参数 `axis`（0/1/2），S 的分量在编译期已知、最快。**Unaligned** 版把轴变成运行时 `Vector3 axis` 成员，`calc` 用 Rodrigues 现算旋转——灵活但略慢。

### 9.4 Unbounded 系列
普通 revolute 用**角度值**存 q，转过 ±π 会有环绕/插值问题。Unbounded 改用 **(cosθ, sinθ)** 两个数存 q（NQ=2），位形流形是 **SO(2)**，integrate/difference 走李群运算——适合需要连续多圈转动的轮子/连续关节。

---

## 10. 李群系统：位形流形的加减法

`liegroup/`。这是解决"**位形不能直接加减**"的子系统。

### 10.1 为什么需要它
广义速度积分 $q \leftarrow q + v\,\Delta t$ 只对**欧氏**自由度成立。对姿态（四元数）直接相加会**跳出流形**（不再是单位四元数）。正确做法是在**李群**上：

$$q \oplus v \;=\; q \cdot \exp(v)\quad(\text{integrate}),\qquad q_1 \ominus q_0 \;=\; \log(q_0^{-1} q_1)\quad(\text{difference})$$

于是 `integrate`/`difference` 成了所有涉及位形的算法（数值积分、IK、插值、随机采样、有限差分求导）的地基。

### 10.2 `LieGroupBase<Derived>`（`liegroup-base.hxx`）API
CRTP 基类，统一接口（节选）：

| 方法 | 数学 | 说明 |
|------|------|------|
| `integrate(q,v,qout)` | $q\oplus v$ | 沿切向量走一步 |
| `difference(q0,q1,d)` | $q_1\ominus q_0$ | 两位形之差（切向量） |
| `dIntegrate_dq/dv` | $\partial(q\oplus v)/\partial q,\partial v$ | 积分的雅可比（= Jexp 类） |
| `dDifference<ARG0/ARG1>` | $\partial(q_1\ominus q_0)/\partial q_{0/1}$ | 差分的雅可比（= Jlog 类） |
| `interpolate(q0,q1,u)` | $q_0\oplus(u\,(q_1\ominus q_0))$ | 测地线插值 |
| `dIntegrateTransport` | 平行移动 | 把切空间量沿积分搬运 |
| `randomConfiguration` / `normalize` / `isNormalized` | — | 流形上采样/投影/校验 |
| `squaredDistance` / `distance` | $\|q_1\ominus q_0\|^2$ | 测地距离 |
| `nq()` `nv()` `neutral()` `name()` | — | 维度、单位元、名字 |

### 10.3 三种基础群 + 笛卡尔积

**① `VectorSpaceOperationTpl<Dim>`（`vector-space.hxx`）—— ℝⁿ**
平凡群：`integrate` = $q+v$，`difference` = $q_1-q_0$，雅可比 = $I$。绝大多数关节（转/移/欧拉球…）默认用它。

**② `SpecialOrthogonalOperationTpl<N>`（`special-orthogonal.hxx`）—— SO(2)/SO(3)**
SO(3)：nq=4（四元数）、nv=3。核心直接调用 `spatial/` 的四元数 exp/log：
```cpp
integrate: quat_out = quat * exp3(v); firstOrderNormalize(quat_out);   // q·exp(v)
difference: d = log3(quat0.conjugate() * quat1);                       // log(q0⁻¹q1)
dDifference<ARG1>: Jlog3(R);   dDifference<ARG0>: -Jlog3(R)·Rᵀ
```
`neutral` = (0,0,0,1)。SO(2) 版用 (cosθ,sinθ) 存位形。这里的 `Jlog3`/`exp3` 正是 [空间代数解析 §7](空间代数运算解析.md) 里详解的那些。

**③ `SpecialEuclideanOperationTpl<N>`（`special-euclidean.hxx`）—— SE(2)/SE(3)**
SE(3)：nq=7（平移 3 + 四元数 4）、nv=6、`neutral`=(0,0,0,0,0,0,1)。
```cpp
difference: d = log6( SE3(quat0⁻¹·quat1, quat0⁻¹·(t1−t0)) ).toVector();  // 完整 SE3 log
integrate: 用 quaternion::exp6(v) 得到增量，再复合平移与姿态
dDifference<ARG1>: Jlog6(M);   dDifference<ARG0>: 复合 -R^T/skew 块后左乘 Jlog6
```
`R3crossSO3_t` 类型注释点明：SE(3) 位形空间可视为 ℝ³×SO(3)，但**integrate/difference 用的是真正的 SE(3) 测地线**（`exp6`/`log6`，平移与旋转耦合），不是分开的 ℝ³ 与 SO(3)。这与 [空间代数解析 §7.7](空间代数运算解析.md) 的 log6/Jlog6 完全对应。

**④ 笛卡尔积（`cartesian-product*.hxx`）**
整机位形空间 = 各关节流形的笛卡尔积 $\prod_i \mathcal{G}_i$。`integrate` 逐段调用各子群的 `integrate`。`liegroup-generic.hxx`/`liegroup-collection.hxx` 把这些群也做成 variant，`cartesian-product-variant` 支持运行时拼装。

### 10.4 `LieGroupMap`（`liegroup-map.hxx`）—— 关节→流形映射表（关键枢纽）
一张编译期查找表，决定"每种关节的 q 该用哪个群积分"：

| 关节 | 映射到的流形 |
|------|-------------|
| 默认（转/移/欧拉球/万向…） | `VectorSpaceOperationTpl<NQ>`（欧氏） |
| `JointModelSpherical` | `SpecialOrthogonalOperationTpl<3>`（SO(3)） |
| `JointModelFreeFlyer` | `SpecialEuclideanOperationTpl<3>`（SE(3)） |
| `JointModelPlanar` | `SpecialEuclideanOperationTpl<2>`（SE(2)） |
| `JointModelRevoluteUnbounded{,Unaligned}` | `SpecialOrthogonalOperationTpl<2>`（SO(2)） |
| Composite/Mimic/generic | 笛卡尔积 variant |

`liegroup-algo.hxx` 就是**沿整棵树对每个关节查这张表、调对应群的 integrate/difference**，从而实现整机版的 `pinocchio::integrate(model, q, v)`。

---

## 11. 访问者系统：把 variant 变成可调用的算法

`visitor/`。`Model.joints[i]` 是 `JointModelVariant`（类型被擦除），要在它上面调 `calc` 就得**恢复真实类型**——靠 Boost.Variant 的 `apply_visitor`。

- **`fusion.hxx`**：`bf::append(t, ts...)` 把一组参数 `push_front` 进 Boost.Fusion 序列。访问者需要把"除关节外的其它参数（data、q、v…）"打包成一个序列，随关节类型一起送进 visitor。
- **`joint-unary-visitor.hxx`**：一元访问者。对**单个**关节 variant 分派，如 `calc`、`nv`、`cast`。用户写一个带模板 `operator()(JointModelDerived&, args...)` 的 struct，框架用 `apply_visitor` 在运行时选中真实类型、编译期为每种类型实例化一份。
- **`joint-binary-visitor.hxx`**：二元访问者。对**两个**关节 variant 同时分派（笛卡尔积展开 N×N 种组合），用于需要两个关节类型的操作。
- **`joint-basic-visitors.hxx`**：把上面封装成好用的自由函数：`calc_zero_order(model,data,q)`、`calc_first_order(...)`、`calc_aba(...)`、`nv(jmodel)`、`cast(...)` 等——算法层直接调这些。

> 三件套（variant + apply_visitor + fusion 打包）的最小可运行 demo 与 `::template` 消歧细节见 [源码解析.md](源码解析.md)，此处不重复。

---

## 12. 其余文件：pool / sample-models / force-set

- **`pool/model.hxx`（`ModelPoolTpl`）**：并行计算池。持有 `vector<Model> m_models` 与 `vector<Data> m_datas`（每线程一份拷贝），配合 OpenMP（`omp_get_max_threads`）跑批量前向/碰撞。`update(data)` 广播一份 data 到全池。`pool/geometry.hxx` 是带碰撞几何的版本。
- **`sample-models.hxx`**：`buildModels::humanoid/manipulator(...)` 程序化生成测试机器人（单元测试、benchmark、教程用）。
- **`force-set.hxx`**：力集合/运动集合（`6×N` 矩阵）的批量列变换辅助，配合 `JointMotionSubspace` 的批量 `SᵀF`。

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
