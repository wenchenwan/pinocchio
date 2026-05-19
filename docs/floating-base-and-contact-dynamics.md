# Pinocchio 浮动基座动力学与接触力学文档

> 基于 Pinocchio 源码 (e:\03-自我学习\pinocchio) 整理，面向已有工业臂背景的读者，补齐浮动基座与接触力学知识。

---

## 目录

1. [核心概念对比：固定基座 vs 浮动基座](#1-核心概念对比)
2. [浮动基座动力学](#2-浮动基座动力学)
3. [接触力学基础](#3-接触力学基础)
4. [约束求解器](#4-约束求解器)
5. [冲量动力学（碰撞响应）](#5-冲量动力学)
6. [仿真环境搭建](#6-仿真环境搭建)
7. [完整仿真循环示例](#7-完整仿真循环示例)
8. [API 速查表](#8-api-速查表)

---

## 1. 核心概念对比

| 维度 | 工业臂（固定基座） | 移动机器人（浮动基座） |
|---|---|---|
| 基座 | 固定在地面/框架 | 自由漂浮，随机体运动 |
| 构型空间 | `R^n`，直接用关节角 | `SE(3) × R^n`，包含位姿 |
| `nq` | = 关节数 | = 关节数 + **7**（四元数位姿） |
| `nv` | = 关节数 | = 关节数 + **6**（速度旋量） |
| 欠驱动 | 否（每轴均有电机） | 是（基座 6 DOF 无驱动） |
| 主要挑战 | 逆运动学/动力学 | **接触力维持、全身运动规划** |

**关键不等式**：浮动基座时 `nq ≠ nv`，这是 Pinocchio 所有算法使用 `v`（速度）而非 `q`（构型）作为线性化基的根本原因。

---

## 2. 浮动基座动力学

### 2.1 Free-Flyer 关节

浮动基座在 Pinocchio 中用 **`JointModelFreeFlyer`** 表示，作为机器人运动链的根关节插入。

```python
import pinocchio as pin

# 加载 URDF 时指定浮动基座
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_path, mesh_dir,
    pin.JointModelFreeFlyer()   # <-- 插入 SE(3) 根关节
)
```

**参见源码**：[include/pinocchio/src/multibody/joint/joint-free-flyer.hxx](../include/pinocchio/src/multibody/joint/joint-free-flyer.hxx)

内部实现：`JointMotionSubspaceIdentityTpl`，运动子空间为 6×6 单位矩阵（6 自由度完全开放）。

### 2.2 构型空间与速度空间

```
构型向量 q  ∈ R^{nq}:  [x, y, z, qx, qy, qz, qw,  θ₁, θ₂, ..., θₙ]
                        └─── 基座位姿（SE3，7维）───┘  └── 关节角 ──┘

速度向量 v  ∈ R^{nv}:  [vx, vy, vz, ωx, ωy, ωz,   θ̇₁, θ̇₂, ..., θ̇ₙ]
                        └─── 基座速度旋量（6维）───┘  └── 关节速度 ──┘
```

**重要**：`nq = nv + 1`（四元数比角速度多1维），因此不能直接做 `q += v * dt`，必须用：

```python
q = pin.integrate(model, q, v * dt)   # 李群上的正确积分
```

### 2.3 选择矩阵 S（驱动分配）

浮动基座系统中，基座 6 个 DOF 无电机驱动，力矩向量结构：

```python
tau = np.zeros(model.nv)
# tau[:6]  = 0        # 浮动基座：无驱动（欠驱动部分）
# tau[6:]  = u        # 关节电机力矩

# 选择矩阵：将关节空间力矩映射到全广义坐标空间
S = np.zeros((model.nv - 6, model.nv))
S.T[6:, :] = np.eye(model.nv - 6)
```

**参见**：[examples/simulation-contact-dynamics.py:79-80](../examples/simulation-contact-dynamics.py#L79)

### 2.4 核心动力学方程

浮动基座多体系统的运动方程：

```
M(q) v̇ + C(q,v) v + g(q) = S^T τ + J^T(q) λ
J(q) v̇ + J̇(q,v) v = 0          （接触约束）
```

其中：
- `M(q)` — 广义质量矩阵（含惯量耦合项）
- `C(q,v)v` — 科里奥利/离心力项
- `g(q)` — 重力项
- `J(q)` — 接触约束雅可比矩阵
- `λ` — 接触力（拉格朗日乘子）

### 2.5 主要算法

| 算法 | 函数 | 复杂度 | 用途 |
|---|---|---|---|
| ABA（正向动力学） | `pin.aba(model, data, q, v, tau)` | O(n) | 无约束正向积分 |
| RNEA（逆向动力学） | `pin.rnea(model, data, q, v, a)` | O(n) | 计算所需力矩 |
| CRBA | `pin.crba(model, data, q)` | O(n²) | 计算质量矩阵 M |
| 约束动力学 | `pin.constraintDynamics(...)` | O(n + c²) | 有接触的正向动力学 |

---

## 3. 接触力学基础

### 3.1 接触类型

Pinocchio 定义了两种刚性接触类型（[contact-info.hxx:17-22](../include/pinocchio/src/constraints/contact-info.hxx#L17)）：

```python
pin.ContactType.CONTACT_3D   # 点接触：约束 3D 位置（法向 + 切向）
pin.ContactType.CONTACT_6D   # 面接触：约束完整 SE(3) 位姿（位置 + 姿态）
```

**选择准则**：
- 点接触（脚尖/轮子）→ `CONTACT_3D`（3个约束）
- 面接触（平足/法兰盘）→ `CONTACT_6D`（6个约束）

### 3.2 RigidConstraintModel（刚性约束模型）

这是接触力学的核心数据结构，描述**一个接触点的几何与约束信息**。

```python
# 源码位置：include/pinocchio/src/constraints/contact-info.hxx

# 构造函数签名（Python）：
contact_model = pin.RigidConstraintModel(
    contact_type,    # CONTACT_3D 或 CONTACT_6D
    model,           # 机器人模型
    joint1_id,       # 接触关节 ID
    joint1_placement,# 接触点在关节局部坐标系中的位姿（SE3）
    joint2_id=0,     # 参考关节 ID（0 = 世界坐标系）
    joint2_placement=SE3.Identity()  # 接触目标位姿（世界坐标系）
)
```

**关键属性**：

| 属性 | 含义 |
|---|---|
| `contact_model.joint1_id` | 接触侧关节索引 |
| `contact_model.joint1_placement` | 接触点在 joint1 局部系的位姿 |
| `contact_model.joint2_id` | 参考侧关节（通常为 0=世界系） |
| `contact_model.reference_frame` | 约束参考坐标系 |
| `contact_model.size()` | 约束维度（3 或 6） |

### 3.3 面接触示例（双足机器人）

```python
# 参见：examples/simulation-contact-dynamics.py:62-69

feet_name = ["left_sole_link", "right_sole_link"]
contact_models = []
contact_datas = []

for frame_name in feet_name:
    frame_id = model.getFrameId(frame_name)
    frame = model.frames[frame_id]

    contact_model = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_6D,     # 平足，约束全 6 DOF
        model,
        frame.parentJoint,              # 接触关节 ID
        frame.placement                 # 接触点在关节坐标系中的偏置
    )
    contact_models.append(contact_model)
    contact_datas.append(contact_model.createData())
```

### 3.4 点接触示例（四足机器人）

```python
# 参见：examples/anymal-simulation.py:31-38

for j, frame_id in enumerate(foot_frame_ids):
    contact_model = pinocchio.RigidConstraintModel(
        pinocchio.ContactType.CONTACT_3D,       # 点接触，约束 3D 位置
        robot.model,
        foot_joint_ids[j],
        robot.model.frames[frame_id].placement,
        0,                                       # 世界坐标系为参考
        data.oMf[frame_id],                     # 目标位姿（当前脚落地位置）
    )
    constraint_models.append(contact_model)
```

### 3.5 约束动力学的数学原理

#### 3.5.1 从拉格朗日力学到约束方程

对于无约束的浮动基座系统，运动方程由欧拉-拉格朗日方程给出：

```
M(q) q̈ + h(q, q̇) = τ

其中：
  M(q)      ∈ R^{nv×nv}  广义质量矩阵（对称正定）
  h(q, q̇)   ∈ R^{nv}     非线性项 = C(q,q̇)q̇ + g(q)
  τ         ∈ R^{nv}     广义力（关节力矩 + 外力映射）
  q̈         ∈ R^{nv}     广义加速度（求解目标）
```

**加入接触约束**后，每个接触点对系统施加一个**保持接触点不动**的几何约束。设接触点在笛卡尔空间的位置为 `φ(q)`，则：

```
φ(q) = 常数（接触点不动）

对时间求导：  J(q) q̇ = 0         （速度约束）
再次求导：    J(q) q̈ + J̇(q,q̇) q̇ = 0  （加速度约束）
```

定义 **约束漂移** `γ(q, q̇) = J̇(q,q̇) q̇`，加速度层约束写为：

```
J(q) q̈ = -γ(q, q̇)
```

接触力 `f ∈ R^{3 or 6}`（笛卡尔空间）通过雅可比矩阵映射到广义力空间：`τ_contact = J^T f`。

综合后，**含约束的运动方程**为：

```
M(q) q̈ + h(q, q̇) = τ + J^T(q) λ      …(1)
J(q) q̈ + γ(q, q̇) = 0                 …(2)
```

其中 `λ` 是拉格朗日乘子，物理含义是接触力（以广义坐标为基底的等效力）。

---

#### 3.5.2 最优化视角：M-范数最小化

Pinocchio 将约束动力学表述为一个**带等式约束的二次规划**（QP）问题（[constrained-dynamics.hpp:79-84](../include/pinocchio/algorithm/constrained-dynamics.hpp#L79)）：

```
min_{q̈}   ½ ‖q̈ - q̈_free‖²_{M(q)}

s.t.       J(q) q̈ + γ(q, q̇) = 0
```

**各项含义**：

- `‖·‖²_{M}` 表示以质量矩阵为度量的加权范数：`‖x‖²_M = x^T M x`
  - 这在物理上等价于**最小化附加动能**，即在满足接触约束的前提下，系统"尽量少改变"自由运动状态
- `q̈_free = M⁻¹(τ - h)` 是**无约束自由加速度**，即若地面突然消失时系统的加速度

这个目标函数的选择具有物理意义：接触力所做的功最小（无能量注入），系统只通过约束力被动改变运动方向。

---

#### 3.5.3 KKT 条件与线性系统求解

对上述 QP 建立拉格朗日函数并对 `q̈` 和 `λ` 求偏导，令其为零（KKT 一阶必要条件）：

```
∂L/∂q̈ = 0:    M(q̈ - q̈_free) + J^T λ = 0   →   M q̈ + J^T λ = M q̈_free = τ - h
∂L/∂λ = 0:    J q̈ + γ = 0
```

合并写成**鞍点线性系统（KKT 系统）**：

```
┌          ┐ ┌    ┐   ┌              ┐
│  M    J^T│ │ q̈  │   │ τ - h(q,q̇) │
│          │ │    │ = │              │
│  J    0  │ │ -λ │   │  -γ(q,q̇)   │
└          ┘ └    ┘   └              ┘

其中：
  左侧 (nv+nc) × (nv+nc) 块矩阵称为 KKT 矩阵
  nc = 总约束维度（所有接触点约束之和）
  h = C(q,q̇)q̇ + g(q)（由 RNEA 计算）
```

**解析求解过程**（Schur 补）：

```
第一步：从第一行解出  q̈ = M⁻¹(τ - h - J^T λ) = q̈_free - M⁻¹ J^T λ

第二步：代入第二行约束：
         J(q̈_free - M⁻¹ J^T λ) = -γ
         J M⁻¹ J^T λ = J q̈_free + γ

第三步：定义操作空间惯量矩阵（Operational Space Inertia Matrix）：
         Λ = (J M⁻¹ J^T)  ∈ R^{nc×nc}   （也称 "Delassus 矩阵"）
         求解：λ = Λ⁻¹ (J q̈_free + γ)

第四步：代回求加速度：
         q̈ = q̈_free - M⁻¹ J^T λ
```

`Λ = J M⁻¹ J^T` 是**接触空间的有效惯量**，它描述了在接触点施加单位冲量时，接触点产生的加速度响应——这与工业臂末端的操作空间惯量矩阵概念完全一致。

---

#### 3.5.4 约束漂移 γ 与 Baumgarte 稳定化

**γ 的物理来源**

理想情况下，从满足 `J q̇ = 0`（接触点速度为零）出发，对时间求导得：

```
J q̈ + J̇ q̇ = 0   →   γ = J̇ q̇
```

但实际仿真中由于数值误差，约束在位置和速度层都会产生漂移：

```
位置漂移：φ(q) ≠ 0    （接触点实际上已经"穿进"地面或离开地面）
速度漂移：J q̇ ≠ 0    （接触点速度不完全为零）
```

**Baumgarte 稳定化**：在 γ 中加入反馈修正项，让漂移以指数衰减方式收敛：

```
γ_baumgarte = J̇ q̇ + 2α (J q̇) + β² φ(q)

其中：
  α > 0  速度漂移增益（阻尼项）
  β > 0  位置漂移增益（刚度项）
  α = β  时对应临界阻尼（最快无超调收敛）
```

Pinocchio 中 `cm.calc(model, data, cd)` 会计算当前接触点的位置误差 `cd.c1Mc2`，这个误差就是用于 Baumgarte 修正的原始量。

---

#### 3.5.5 正则化与近端迭代（处理冗余约束）

当约束冗余时（如双足同时着地，两脚约束在垂直方向线性相关），`Λ = J M⁻¹ J^T` 会变得**奇异或病态**，直接求逆不稳定。

Pinocchio 使用**Tikhonov 正则化（近端算法）**：

```
正则化 KKT 系统：
┌            ┐ ┌    ┐   ┌              ┐
│  M    J^T  │ │ q̈  │   │ τ - h        │
│            │ │    │ = │              │
│  J   -μI  │ │ -λ │   │ -γ + μ λ_k  │
└            ┘ └    ┘   └              ┘

其中 μ > 0 是正则化参数（ProximalSettings.mu）
λ_k 是上一次迭代的接触力估计（近端点）
```

**迭代过程**（近端点算法）：

```
初始化：λ_0 = 0

第 k 次迭代：
  1. 求解正则化 KKT 系统得到 (q̈_k, λ_k+1)
  2. 检查残差：‖λ_k+1 - λ_k‖ < accuracy?
  3. 若收敛则停止，否则继续

μ 的作用：μ → 0 时趋近精确约束；μ 较大时约束允许轻微违反（软约束效果）
```

每次迭代的 Schur 补变为：

```
Λ_μ = J M⁻¹ J^T + μ I
```

加上 `μI` 后矩阵始终正定，Cholesky 分解不再奇异。这就是 `mu=1e-12` 能处理冗余约束的原因。

---

#### 3.5.6 Cholesky 分解求解 KKT 系统

Pinocchio 不直接对完整 KKT 矩阵做 Cholesky，而是利用多体系统的**稀疏树结构**，通过**Composite Rigid Body Algorithm (CRBA)** 高效分解 `M`，再基于 `M` 的因子化构建 `Λ`：

```
M 的 Cholesky 分解：M = L D L^T （利用运动链树结构，O(n) 复杂度）

Λ = J M⁻¹ J^T 的构建：
  对每列 j 的 J，解线性系统 M x = J_j^T，得到 x = M⁻¹ J_j^T
  然后 Λ_ij = J_i · x_j

Λ 的维度为 nc×nc（接触约束维度），通常远小于 nv×nv
→ 总体复杂度 O(n·nc + nc³)
```

这正是 `ConstraintCholeskyDecomposition` 所封装的逻辑（[constraint-cholesky.hpp](../include/pinocchio/algorithm/constraint-cholesky.hpp)）。

---

#### 3.5.7 约束力的物理解读

求解后，接触力 `λ`（`data.lambda_c`）的物理含义取决于接触类型：

```
CONTACT_3D（点接触）：
  λ = [fx, fy, fz]^T   在参考坐标系下的 3D 接触力
  单位：N（牛顿）

CONTACT_6D（面接触）：
  λ = [fx, fy, fz, mx, my, mz]^T   力旋量（力 + 力矩）
  单位：[N, N·m]
```

力学上的合法性检验（摩擦锥，需用户自行验证）：

```
法向力（防穿透）：fz ≥ 0
摩擦约束（Coulomb）：sqrt(fx² + fy²) ≤ μ_friction × fz
```

Pinocchio 的刚性接触假设不自动施加摩擦约束，如需摩擦锥约束需使用 `CoulombFrictionCone` 约束集（[constraints/sets/coulomb-friction-cone.hxx](../include/pinocchio/src/constraints/sets/coulomb-friction-cone.hxx)）。

---

### 3.6 constraintDynamics 完整调用流程

```
调用前（仅初始化一次，循环外）
  └── initConstraintDynamics(model, data, cms, cds)
        ├── 重建约束 Cholesky 结构（根据约束拓扑分配内存）
        ├── 分配 lambda_c, impulse_c, osim 等缓存矩阵
        └── 所有缓存清零

每步仿真循环内
  ├── Step 1: computeJointJacobians(model, data, q)
  │           计算所有关节的体雅可比矩阵 J（存入 data.J）
  │
  ├── Step 2: data.q_in = q（告知内部状态当前构型）
  │
  ├── Step 3: cm.calc(model, data, cd)   [对每个接触约束]
  │           ├── 计算接触点当前位姿 c1Mc2（相对误差 SE3）
  │           └── 用于构建约束漂移 γ（位置和速度误差）
  │
  └── Step 4: constraintDynamics(model, data, q, v, tau, cms, cds, prox)
              ├── 调用 RNEA 计算 h = C v + g（非线性项）
              ├── 调用 CRBA 计算 M 并做 Cholesky 分解
              ├── 构建 Λ = J M⁻¹ J^T（+ μI 正则化）
              ├── 近端迭代（max_iter 次）求解 KKT 系统
              ├── 输出 data.ddq = q̈（关节加速度）
              └── 输出 data.lambda_c = λ（接触力）
```

```python
# 完整调用示例（参见 examples/simulation-contact-dynamics.py:120-127）

prox_settings = pin.ProximalSettings(1e-12, 1e-12, 10)

# ── 初始化（仅调用一次）──────────────────────────
pin.initConstraintDynamics(model, data, contact_models, contact_datas)

# ── 仿真循环内 ────────────────────────────────────
pin.computeJointJacobians(model, data, q)   # Step 1
data.q_in = q                               # Step 2

for cm, cd in zip(contact_models, contact_datas):
    cm.calc(model, data, cd)                # Step 3：更新位置误差

a = pin.constraintDynamics(                 # Step 4：KKT 求解
    model, data, q, v, tau,
    contact_models, contact_datas,
    prox_settings
)

# 读取结果
q̈ = data.ddq         # 形状 (nv,)，关节加速度
λ = data.lambda_c    # 形状 (nc,)，接触力（拉格朗日乘子）
print(f"接触力范数：{np.linalg.norm(λ):.4f} N")
print(f"约束违反：{np.linalg.norm(J @ a + γ):.2e}")   # 应接近 0
```

---

## 4. 约束求解器

### 4.1 ProximalSettings 参数详解

近端求解器的参数控制**精度与速度的权衡**（[proximal.hxx:25-113](../include/pinocchio/src/algorithm/proximal.hxx#L25)）：

```python
prox_settings = pin.ProximalSettings(
    absolute_accuracy,  # 绝对收敛阈值：‖λ_{k+1} - λ_k‖ < 该值时停止
    relative_accuracy,  # 相对收敛阈值：‖λ_{k+1} - λ_k‖/‖λ_k‖ < 该值时停止
    mu,                 # Tikhonov 正则化参数：Λ_reg = Λ + μI
    max_iter            # 最大迭代次数上限
)
```

| 参数 | 典型值 | 物理含义 |
|---|---|---|
| `absolute_accuracy` | `1e-12` | 接触力估计的绝对精度（N） |
| `relative_accuracy` | `1e-12` | 接触力相对变化容忍度 |
| `mu` | `1e-12`（冗余约束） / `0`（满秩约束） | KKT 矩阵正则化强度 |
| `max_iter` | `1`（最快）/ `10~100`（精确） | 近端迭代上限 |

**读取求解状态**：

```python
a = pin.constraintDynamics(...)

print(prox_settings.iter)              # 实际迭代次数
print(prox_settings.absolute_residual) # 最终绝对残差 ‖λ_{k+1} - λ_k‖
print(prox_settings.relative_residual) # 最终相对残差
```

**调优策略**：

```
场景一：单一接触，约束满秩
  → mu=0, max_iter=1  （一次 Cholesky，最快）

场景二：双脚/四足同时接触，约束冗余
  → mu=1e-12, max_iter=1  （正则化后仍单次求解，通常够用）

场景三：精确力估计（如力控制、碰撞检测）
  → mu=1e-12, max_iter=50, accuracy=1e-10

场景四：柔顺控制（允许接触点轻微滑动）
  → mu=1e-6 ~ 1e-3, max_iter=1  （较大 μ 产生"软"约束效果）
```

### 4.2 ConstraintCholeskyDecomposition 内部结构

`ConstraintCholeskyDecomposition` 封装了对以下块结构矩阵的高效分解：

```
KKT 矩阵结构：
┌            ┐
│  M    J^T  │   ← nv 行
│            │
│  J   -μI  │   ← nc 行
└            ┘
  ↑ nv 列  ↑ nc 列

总规模：(nv + nc) × (nv + nc)
```

分解过程利用了 `M` 的**稀疏树状结构**（由运动链决定），避免了对完整稠密矩阵做 Cholesky：

```python
# 创建分解器（根据约束拓扑预分配内存）
kkt = pinocchio.ConstraintCholeskyDecomposition(model, constraint_models)

# 每步计算（更新当前 q 下的因子）
kkt.compute(model, data, constraint_models, constraint_datas, mu=0.0)

# 求解线性系统 KKT · z = rhs
#   rhs 的结构：[约束残差部分 (nc), 广义力部分 (nv)]
#   z   的结构：[对偶变量增量 (nc), 构型增量 (nv)]
rhs = np.concatenate([-constraint_value, generalized_force_error])
z = kkt.solve(rhs)

delta_lambda = z[:nc]   # 接触力修正量
delta_q = z[nc:]        # 构型修正量（用于逆运动学）
```

**两种用途对比**：

| 用途 | 求解变量 | 典型 rhs 结构 |
|---|---|---|
| 正向动力学（constraintDynamics） | `(q̈, λ)` | `[τ-h, -γ]` |
| 约束逆运动学（anymal-simulation） | `(dλ, dq)` | `[-约束误差, 任务误差]` |

---

## 5. 冲量动力学

处理**瞬时碰撞/冲击**事件，计算碰撞后速度。

### 5.1 物理背景：连续 vs 瞬时接触

| 特性 | 连续接触（constraintDynamics） | 瞬时碰撞（impulseDynamics） |
|---|---|---|
| 时间尺度 | 持续多步 | 瞬时（Δt → 0） |
| 输出量 | 加速度 `q̈` | 速度跳变 `Δq̇ = q̇⁺ - q̇⁻` |
| 接触力 | 有限大小（`data.lambda_c`，单位 N） | 冲量（`data.impulse_c`，单位 N·s） |
| 典型场景 | 站立、行走支撑相 | 脚着地瞬间、物体被抓取 |

### 5.2 数学推导

**从动量定理出发**（[impulse-dynamics.hpp:41-46](../include/pinocchio/algorithm/impulse-dynamics.hpp#L41)）：

碰撞过程中，接触力极大但时间极短，其积分（冲量）为有限值：

```
M(q) (q̇⁺ - q̇⁻) = J^T(q) p_c

其中 p_c ∈ R^{nc} 是碰撞冲量（impulse），即接触力对时间的积分
```

接触约束在速度层施加：

```
J(q) q̇⁺ = -ε J(q) q̇⁻

ε = 0：完全非弹性碰撞，碰后接触点速度为零（q̇⁺ 满足约束）
ε = 1：完全弹性碰撞，接触点速度完全反弹
ε ∈ (0,1)：部分弹性，Newton 恢复系数模型
```

**等价 QP 问题**：

```
min_{q̇⁺}   ½ ‖q̇⁺ - q̇⁻‖²_{M(q)}

s.t.        J(q) q̇⁺ = -ε J(q) q̇⁻
```

**解析解**（与 3.5.3 中 Schur 补推导完全对称）：

```
第一步：q̇⁺ = q̇⁻ + M⁻¹ J^T p_c

第二步：代入速度约束：
        J(q̇⁻ + M⁻¹ J^T p_c) = -ε J q̇⁻
        Λ p_c = -(1+ε) J q̇⁻
        p_c = -Λ⁻¹ (1+ε) J q̇⁻

第三步：代回：
        q̇⁺ = q̇⁻ - M⁻¹ J^T Λ⁻¹ (1+ε) J q̇⁻
```

注意：`ε = 0` 时，`q̇⁺` 满足 `J q̇⁺ = 0`（接触点速度为零，与连续接触初始条件一致）。

### 5.3 能量分析

碰撞过程的动能变化：

```
ΔKE = ½(q̇⁺)^T M q̇⁺ - ½(q̇⁻)^T M q̇⁻

完全非弹性（ε=0）：ΔKE ≤ 0，能量只能损失（满足热力学第二定律）
完全弹性（ε=1）：ΔKE = 0，动能守恒
```

### 5.4 调用方式

```python
# 完整碰撞处理流程

# 检测碰撞（此处假设已知碰撞发生）
r_coeff = 0.0                # 完全非弹性（最常用）
prox_settings = pin.ProximalSettings(1e-12, 1e-12, 10)

# 注意：contact_models 应仅包含正在碰撞的接触点
pin.impulseDynamics(
    model, data,
    q,               # 碰撞瞬间的构型
    v_before,        # 碰撞前速度（q̇⁻）
    contact_models,  # 新激活的接触约束
    contact_datas,
    r_coeff,         # 恢复系数 ε
    prox_settings
)

v_after = data.dq_after    # 碰撞后速度（q̇⁺），形状 (nv,)
impulse = data.impulse_c   # 碰撞冲量 p_c，形状 (nc,)，单位 N·s

# 典型仿真循环中的使用
v = data.dq_after.copy()   # 用碰撞后速度替换当前速度，继续积分
```

---

## 6. 仿真环境搭建

### 6.1 模型加载

```python
import pinocchio as pin
import numpy as np

# 方法 1：从 URDF 加载（固定基座）
model = pin.buildModelFromUrdf("robot.urdf")

# 方法 2：从 URDF 加载（浮动基座）
model, _, visual_model = pin.buildModelsFromUrdf(
    "robot.urdf", mesh_dir, pin.JointModelFreeFlyer()
)

# 方法 3：编程构建（简单模型）
model = pin.Model()
joint_id = model.addJoint(0, pin.JointModelRX(), pin.SE3.Identity(), "joint1")
model.appendBodyToJoint(joint_id, pin.Inertia.FromSphere(1.0, 0.1), pin.SE3.Identity())
```

### 6.2 数据结构

```python
data = model.createData()   # 算法计算结果的缓存结构

# 主要数据字段
data.ddq          # 关节加速度（constraintDynamics 的输出）
data.lambda_c     # 接触力（拉格朗日乘子）
data.dq_after     # 碰撞后速度（impulseDynamics 的输出）
data.M            # 质量矩阵（crba 后可用）
data.nle          # 非线性效应 C(q,v)v + g(q)（computeAllTerms 后）
data.com[0]       # 全局质心位置
data.oMi[i]       # 第 i 个关节的全局位姿
data.oMf[f]       # 第 f 个框架的全局位姿
```

### 6.3 可视化（MeshcatVisualizer）

```python
from pinocchio.visualize import MeshcatVisualizer

viz = MeshcatVisualizer(model, collision_model, visual_model)
viz.initViewer(open=True)    # open=True 自动打开浏览器
viz.loadViewerModel()
viz.display(q)               # 显示当前构型
```

MeshCat 在浏览器中通过 WebGL 渲染，访问 `http://127.0.0.1:7000/static/` 查看。

### 6.4 运动学计算

```python
# 正向运动学
pin.forwardKinematics(model, data, q)
pin.forwardKinematics(model, data, q, v)        # 同时计算速度
pin.forwardKinematics(model, data, q, v, a)     # 同时计算加速度

# 更新所有框架位姿
pin.framesForwardKinematics(model, data, q)

# 计算所有项（M, C, g, Jacobians 等）
pin.computeAllTerms(model, data, q, v)

# 获取框架雅可比矩阵（6 × nv）
J = pin.getFrameJacobian(
    model, data, frame_id,
    pin.ReferenceFrame.LOCAL_WORLD_ALIGNED  # 或 LOCAL / WORLD
)
```

---

## 7. 完整仿真循环示例

### 7.1 无接触仿真（单摆）

```python
# 参见：examples/simulation-pendulum.py

import pinocchio as pin
import numpy as np

# 构型、速度初始化
q = pin.neutral(model)
v = np.zeros(model.nv)
tau = np.zeros(model.nv)
dt = 0.01

for k in range(N):
    # 控制律（小阻尼）
    tau = -0.1 * v

    # ABA：正向动力学，O(n) 复杂度
    a = pin.aba(model, data, q, v, tau)

    # 半隐式欧拉积分
    v += a * dt
    q = pin.integrate(model, q, v * dt)   # 注意：必须用 integrate！

    viz.display(q)
```

### 7.2 有接触仿真（双足机器人）

```python
# 参见：examples/simulation-contact-dynamics.py

# ── 初始化（循环外）──────────────────────────────────────────
contact_models, contact_datas = [], []
for frame_name in feet_name:
    frame_id = model.getFrameId(frame_name)
    frame = model.frames[frame_id]
    cm = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_6D, model,
        frame.parentJoint, frame.placement
    )
    contact_models.append(cm)
    contact_datas.append(cm.createData())

pin.initConstraintDynamics(model, data_sim, contact_models, contact_datas)
prox_settings = pin.ProximalSettings(1e-12, 1e-12, 10)

q = q0.copy()
v = np.zeros(model.nv)

# ── 仿真循环 ─────────────────────────────────────────────────
while t <= T:
    # 1. 计算控制力矩（RNEA 静态补偿 + 姿态 PD）
    pin.computeJointJacobians(model, data_control, q)
    b = pin.rnea(model, data_control, q, v, np.zeros(model.nv))
    sol = np.linalg.lstsq(A.T, b, rcond=None)[0]
    tau = np.concatenate([np.zeros(6), sol[:model.nv - 6]])
    tau[6:] += -Kp * err_q[6:] - Kv * err_v[6:]

    # 2. 更新约束（计算约束残差）
    pin.computeJointJacobians(model, data_sim, q)
    data_sim.q_in = q
    for cm, cd in zip(contact_models, contact_datas):
        cm.calc(model, data_sim, cd)

    # 3. 约束正向动力学
    a = pin.constraintDynamics(
        model, data_sim, q, v, tau,
        contact_models, contact_datas, prox_settings
    )

    # 4. 半隐式欧拉积分
    v += a * dt
    q = pin.integrate(model, q, v * dt)

    viz.display(q)
    t += dt
```

### 7.3 约束下的逆运动学（重心控制）

```python
# 参见：examples/anymal-simulation.py:65-128

kkt = pinocchio.ConstraintCholeskyDecomposition(model, constraint_models)
mu = 0.0
y = np.ones(constraint_size)     # 对偶变量

for k in range(N):
    pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))
    pinocchio.computeJointJacobians(model, data, q)

    # 计算约束残差
    kkt.compute(model, data, constraint_models, constraint_datas, mu)
    constraint_value = np.concatenate([cd.c1Mc2.translation for cd in constraint_datas])

    # 重心误差
    com_err = data.com[0] - com_desired

    # 构建右端项并求解 KKT 系统
    rhs = np.concatenate([-constraint_value - y * mu,
                           kp * mass * com_err,
                           np.zeros(model.nv - 3)])
    dz = kkt.solve(rhs)

    dq = dz[constraint_size:]
    q = pinocchio.integrate(model, q, -alpha * dq)
```

---

## 8. API 速查表

### 8.1 模型构建

| 函数 | 说明 |
|---|---|
| `pin.buildModelFromUrdf(path)` | 加载固定基座 URDF |
| `pin.buildModelsFromUrdf(path, mesh, joint)` | 加载 URDF 并指定根关节 |
| `pin.JointModelFreeFlyer()` | 浮动基座关节（SE3，7/6 维） |
| `pin.JointModelRX/RY/RZ()` | 绕轴旋转关节 |
| `pin.JointModelPX/PY/PZ()` | 平移关节 |
| `model.createData()` | 创建算法数据缓存 |

### 8.2 运动学

| 函数 | 说明 |
|---|---|
| `pin.forwardKinematics(model, data, q[, v[, a]])` | 正向运动学 |
| `pin.framesForwardKinematics(model, data, q)` | 更新所有框架位姿 |
| `pin.computeJointJacobians(model, data, q)` | 计算所有关节雅可比 |
| `pin.getFrameJacobian(model, data, fid, ref)` | 获取指定框架雅可比 |
| `pin.integrate(model, q, dq)` | 李群上的构型积分 |
| `pin.difference(model, q1, q2)` | 构型差值（切空间） |

### 8.3 动力学

| 函数 | 返回 | 说明 |
|---|---|---|
| `pin.aba(model, data, q, v, tau)` | `data.ddq` | 关节加速度（无约束） |
| `pin.rnea(model, data, q, v, a)` | `data.tau` | 关节力矩（逆动力学） |
| `pin.crba(model, data, q)` | `data.M` | 质量矩阵 |
| `pin.computeAllTerms(model, data, q, v)` | 多项 | M, C, g, Jacobians 全计算 |

### 8.4 接触动力学

| 函数 | 说明 |
|---|---|
| `pin.RigidConstraintModel(type, model, j1, p1[, j2, p2])` | 创建接触约束模型 |
| `contact_model.createData()` | 创建约束数据 |
| `pin.initConstraintDynamics(model, data, cms, cds)` | 初始化约束动力学（仅一次） |
| `cm.calc(model, data, cd)` | 更新约束残差 |
| `pin.constraintDynamics(model, data, q, v, tau, cms, cds[, prox])` | 有接触正向动力学 |
| `pin.impulseDynamics(model, data, q, v⁻, cms, cds, ε, prox)` | 碰撞冲量动力学 |
| `data.lambda_c` | 接触力（向量，dim = 总约束维度） |
| `data.dq_after` | 碰撞后速度 |

### 8.5 ProximalSettings 参数

```python
# 典型值
prox = pin.ProximalSettings(
    accuracy = 1e-12,   # 收敛精度
    mu       = 1e-12,   # 正则化（冗余约束时用）
    max_iter = 10       # 迭代上限
)
# 单次求解（最快）
prox = pin.ProximalSettings(0, 0, 1)
```

---

## 附录：文件导航

| 主题 | 关键文件 |
|---|---|
| 浮动基座关节 | [include/pinocchio/src/multibody/joint/joint-free-flyer.hxx](../include/pinocchio/src/multibody/joint/joint-free-flyer.hxx) |
| 接触约束模型 | [include/pinocchio/src/constraints/contact-info.hxx](../include/pinocchio/src/constraints/contact-info.hxx) |
| 约束正向动力学 | [include/pinocchio/algorithm/constrained-dynamics.hpp](../include/pinocchio/algorithm/constrained-dynamics.hpp) |
| 冲量动力学 | [include/pinocchio/algorithm/impulse-dynamics.hpp](../include/pinocchio/algorithm/impulse-dynamics.hpp) |
| 近端求解器设置 | [include/pinocchio/src/algorithm/proximal.hxx](../include/pinocchio/src/algorithm/proximal.hxx) |
| Cholesky 约束求解 | [include/pinocchio/algorithm/constraint-cholesky.hpp](../include/pinocchio/algorithm/constraint-cholesky.hpp) |
| 双足仿真示例 | [examples/simulation-contact-dynamics.py](../examples/simulation-contact-dynamics.py) |
| 四足仿真示例 | [examples/anymal-simulation.py](../examples/anymal-simulation.py) |
| 单摆仿真示例 | [examples/simulation-pendulum.py](../examples/simulation-pendulum.py) |
