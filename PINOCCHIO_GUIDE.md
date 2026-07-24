# Pinocchio 技术指南

> 面向人形机器人方向的完整参考：数学原理 · 架构设计 · 全部 Demo 注释

---

## 目录

1. [项目定位](#1-项目定位)
2. [核心数学基础](#2-核心数学基础)
   - [2.1 李群 SE(3)](#21-李群-se3)
   - [2.2 空间速度与空间力](#22-空间速度与空间力)
   - [2.3 空间惯量](#23-空间惯量)
3. [Model / Data 架构](#3-model--data-架构)
4. [核心算法原理](#4-核心算法原理)
   - [4.1 正向运动学 FK](#41-正向运动学-fk)
   - [4.2 逆向运动学 IK](#42-逆向运动学-ik)
   - [4.3 逆向动力学 RNEA](#43-逆向动力学-rnea)
   - [4.4 正向动力学 ABA](#44-正向动力学-aba)
   - [4.5 质量矩阵 CRBA](#45-质量矩阵-crba)
   - [4.6 质心动量学](#46-质心动量学)
   - [4.7 接触约束动力学](#47-接触约束动力学)
   - [4.8 解析梯度](#48-解析梯度)
5. [Demo 注释：Python](#5-demo-注释python)
6. [Demo 注释：C++](#6-demo-注释c)
7. [人形机器人开发路径](#7-人形机器人开发路径)
8. [与工业机器人的对比](#8-与工业机器人的对比)
9. [关键头文件索引](#9-关键头文件索引)
10. [源码目录地图](#10-源码目录地图)
11. [源码架构三大机制](#11-源码架构三大机制)
    - [11.1 三层文件布局 .hpp / .hxx / .cpp](#111-三层文件布局-hpp--hxx--cpp)
    - [11.2 Tpl 模板 + context 默认标量](#112-tpl-模板--context-默认标量)
    - [11.3 CRTP + Boost.Fusion 访问者](#113-crtp--boostfusion-访问者)
12. [实战：追踪一个算法从 API 到实现](#12-实战追踪一个算法从-api-到实现)
13. [Python 绑定如何映射到 C++](#13-python-绑定如何映射到-c)
14. [源码阅读推荐顺序](#14-源码阅读推荐顺序)

---

## 1. 项目定位

**Pinocchio** 是 INRIA/LAAS-CNRS 开发的高性能刚体动力学 C++ 模板库（v3.6.0），实现了
Roy Featherstone 体系中最先进的多体动力学算法，并提供所有主要算法的**解析梯度**。

它是以下人形/腿足机器人框架的核心计算引擎：

| 框架 | 用途 |
|------|------|
| **Crocoddyl** | DDP/iLQR 轨迹优化 |
| **Humanoid Path Planner** | 人形运动规划 |
| **Stack-of-Tasks** | 层次全身控制 (WBC) |
| **MPC in Robotics** | 实时模型预测控制 |

**性能特点**：C++ 模板展开 + 缓存友好数据布局，单台笔记本（Intel i7 @ 2.4 GHz）可
在 1 ms 内完成百自由度机器人的动力学计算，支持 OpenMP 并行批量计算。

---

## 2. 核心数学基础

### 2.1 李群 SE(3)

刚体位姿不是普通向量，而属于**特殊欧氏群**（Special Euclidean Group）：

$$SE(3) = \left\{ T = \begin{bmatrix} R & p \\ 0 & 1 \end{bmatrix} \;\middle|\; R \in SO(3),\; p \in \mathbb{R}^3 \right\}$$

其中 $R$ 为旋转矩阵（满足 $R^\top R = I,\;\det R = 1$），$p$ 为平移向量。

**关键运算：**

| 运算 | 数学含义 | Pinocchio API |
|------|----------|---------------|
| $T_1 \cdot T_2$ | 变换复合 | `T1 * T2` |
| $T^{-1}$ | 逆变换 | `T.inverse()` |
| $T^{-1} \cdot T_{des}$ | 相对变换（误差） | `T.actInv(T_des)` |
| $\log(T) \in \mathfrak{se}(3)$ | 对数映射→切空间 | `pinocchio.log6(T)` |
| $\exp(\xi) \in SE(3)$ | 指数映射→群 | `pinocchio.exp6(xi)` |
| $q \oplus v\Delta t$ | 流形上积分 | `pinocchio.integrate(model, q, v*dt)` |
| $q_1 \ominus q_2$ | 流形上差分 | `pinocchio.difference(model, q1, q2)` |

**为什么需要流形积分？**  
人形机器人的浮动基旋转部分以四元数表示（单位球面 $S^3$），不能直接做 $q \mathrel{+}= \dot{q}\Delta t$。
必须在李群上做指数映射才能保持旋转的约束性。这是人形机器人与固定基工业机器人
最本质的数学差异。

**状态空间维度：**

$$n_q = n_{joints} + 7 \quad (\text{4 元数旋转 } + \text{ 3 位移})$$
$$n_v = n_{joints} + 6 \quad (\text{切空间维度，速度/加速度用此维度})$$

### 2.2 空间速度与空间力

Pinocchio 采用 Featherstone **空间代数（Spatial Algebra）** 统一描述运动与力：

**空间速度（Twist）**—— Motion 对象，$\mathbf{v} \in \mathbb{R}^6$：

$$\mathbf{v} = \begin{bmatrix} \omega \\ v \end{bmatrix}$$

$\omega \in \mathbb{R}^3$ 为角速度，$v \in \mathbb{R}^3$ 为线速度（在 LOCAL 坐标系下）。

**空间力（Wrench）**—— Force 对象，$\mathbf{f} \in \mathbb{R}^6$：

$$\mathbf{f} = \begin{bmatrix} \tau \\ f \end{bmatrix}$$

$\tau$ 为力矩，$f$ 为力。两者构成对偶关系，功率：

$$P = \mathbf{f}^\top \mathbf{v} = \tau^\top\omega + f^\top v$$

**坐标系约定（`ReferenceFrame`）**：

| 枚举值 | 含义 |
|--------|------|
| `LOCAL` | 量在关节局部系下表达 |
| `WORLD` | 量在世界系下表达 |
| `LOCAL_WORLD_ALIGNED` | 原点在关节处，轴与世界系对齐 |

### 2.3 空间惯量

每个刚体的惯量以 $6\times6$ 矩阵表示：

$$\mathbf{I} = \begin{bmatrix} I_c + m\,\hat{c}\hat{c}^\top & m\hat{c} \\ m\hat{c}^\top & mE \end{bmatrix}$$

其中 $I_c$ 为质心处惯量张量，$m$ 为质量，$c$ 为质心位置，$\hat{c}$ 为 $c$ 的反对称矩阵
（使得 $\hat{c}v = c\times v$）。

空间动量：$\mathbf{h} = \mathbf{I}\,\mathbf{v}$，空间力方程：$\mathbf{f} = \mathbf{I}\,\mathbf{a} + \mathbf{v}\times^*\mathbf{I}\,\mathbf{v}$。

---

## 3. Model / Data 架构

Pinocchio 强制分离不可变参数（`Model`）与可变中间结果（`Data`）：

```
Model（只读，多线程共享）
├── joints[]         关节类型列表（JointModelRX / FreeFlyer / ...）
├── inertias[]       各连杆空间惯量
├── jointPlacements[] 关节相对父连杆的固定变换
├── names[]          关节名称
├── nq / nv          配置/速度空间维度
└── referenceConfigurations  参考配置（half_sitting 等）

Data（可变，每线程独立一份）
├── oMi[]       世界系下各关节 SE(3) 变换
├── oMf[]       世界系下各 Frame 的 SE(3) 变换
├── v[]         各关节空间速度
├── a[]         各关节空间加速度
├── f[]         各关节空间力
├── tau         关节力矩（RNEA 输出）
├── ddq         关节加速度（ABA 输出）
├── M           质量矩阵（CRBA 输出）
├── Minv        质量矩阵逆（ABA 导数输出）
├── hg          质心动量
├── com[]       质心位置
├── lambda_c    接触力 Lagrange 乘子
└── ...
```

**标准用法模式：**

```python
# 固定基机器人
model = pin.buildModelFromUrdf("robot.urdf")

# 浮动基人形（必须加 FreeFlyer）
model = pin.buildModelFromUrdf("humanoid.urdf", pin.JointModelFreeFlyer())

data = model.createData()   # 每个线程各创建一份

# 算法调用：model 只读，data 存储中间结果
pin.forwardKinematics(model, data, q)
result = data.oMi[joint_id]
```

**多线程**：`ModelPool` 管理多个 `Data` 副本，支持 OpenMP 批量并行。

---

## 4. 核心算法原理

### 4.1 正向运动学 FK

**目标**：给定关节配置 $q$，计算所有连杆在世界系中的位姿 ${}^0T_i$。

**数学原理**：沿运动链递推变换：

$${}^0T_i = {}^0T_{p(i)} \cdot {}^{p(i)}T_i(q_i)$$

其中 $p(i)$ 为第 $i$ 关节的父关节，${}^{p(i)}T_i(q_i)$ 由 URDF 参数和关节类型（旋转/移动/球形）决定。

**API**：

```python
pin.forwardKinematics(model, data, q)          # 只更新 data.oMi[]
pin.forwardKinematics(model, data, q, v)       # 同时计算速度 data.v[]
pin.forwardKinematics(model, data, q, v, a)    # 同时计算加速度 data.a[]
pin.framesForwardKinematics(model, data, q)    # 更新 data.oMf[]（Frame 层）
pin.updateFramePlacement(model, data, frame_id)# 更新单个 Frame
```

**复杂度**：$O(n)$，$n$ 为关节数。

---

### 4.2 逆向运动学 IK

**目标**：给定末端目标位姿 $T_{des}$，求关节配置 $q$。

**方法**：带阻尼最小二乘的迭代 Jacobian 法（Levenberg-Marquardt 思路）。

**迭代步骤**：

**①** 计算 SE(3) 误差（在关节局部系下的对数映射）：

$$e_k = \log\!\left({}^iT_{des}^{-1}\right) \in \mathbb{R}^6$$

**②** 计算修正后的几何 Jacobian：

$$J_{eff} = -J_{\log_6}\!\left({}^iT_{des}^{-1}\right) \cdot {}^0J_i(q)$$

其中 $J_{\log_6}$ 是 $SE(3)$ 对数映射关于其参数的 Jacobian（通过 `Jlog6` 获得）。

**③** 阻尼最小二乘求速度：

$$\dot{q} = -J_{eff}^\top \left(J_{eff}J_{eff}^\top + \lambda^2 I_6\right)^{-1} e_k$$

阻尼系数 $\lambda$ 防止 Jacobian 奇异时解爆炸。

**④** 流形积分更新配置（非普通加法）：

$$q_{k+1} = q_k \oplus (\dot{q} \cdot \Delta t)$$

**收敛判据**：$\|e_k\|_2 < \varepsilon$（典型值 $10^{-4}$）。

---

### 4.3 逆向动力学 RNEA

**目标**：给定 $(q, \dot{q}, \ddot{q})$，求所需关节力矩 $\tau$。

**方程来源**：多体动力学方程（欧拉-拉格朗日方程）：

$$M(q)\ddot{q} + C(q,\dot{q})\dot{q} + g(q) = \tau + J_c^\top\lambda$$

无接触时 $\lambda=0$，RNEA 直接计算右侧 $\tau$。

**两趟递推，复杂度 $O(n)$**：

**第一趟（从根到叶，正向）** —— 计算各连杆速度与加速度：

$$\mathbf{v}_i = {}^iX_{p(i)}\,\mathbf{v}_{p(i)} + \mathbf{s}_i\dot{q}_i$$

$$\mathbf{a}_i = {}^iX_{p(i)}\,\mathbf{a}_{p(i)} + \mathbf{s}_i\ddot{q}_i + \mathbf{v}_i\!\times\!\mathbf{s}_i\dot{q}_i$$

**第二趟（从叶到根，反向）** —— 计算各连杆受力并向上传递：

$$\mathbf{f}_i = \mathbf{I}_i\mathbf{a}_i + \mathbf{v}_i\!\times^*\!\mathbf{I}_i\mathbf{v}_i - {}^i\mathbf{f}^{ext}_i + \sum_{j\in\text{child}(i)} {}^iX_j^*\,\mathbf{f}_j$$

**投影到关节轴**：

$$\tau_i = \mathbf{s}_i^\top\mathbf{f}_i$$

其中 $\mathbf{s}_i$ 为关节子空间（旋转关节为 $[0,0,0,0,0,1]^\top$ 等），${}^iX_j$ 为空间坐标变换算子，$\times^*$ 为空间力的伴随（反对称）运算。

---

### 4.4 正向动力学 ABA

**目标**：给定 $(q, \dot{q}, \tau)$，求关节加速度 $\ddot{q}$。

**原理**：直接求解 $\ddot{q} = M(q)^{-1}(\tau - h(q,\dot{q}))$ 但不显式构造 $M$ 也不做 $O(n^3)$ 矩阵求逆，
通过三趟 $O(n)$ 递推完成：

**关键量**：关节化体惯量（Articulated Body Inertia）$\mathbf{I}_i^A$，吸收子树贡献：

$$\mathbf{I}_i^A = \mathbf{I}_i + \sum_{j\in\text{child}(i)} {}^iX_j^*\!\left(\mathbf{I}_j^A - \frac{\mathbf{I}_j^A\mathbf{s}_j\mathbf{s}_j^\top\mathbf{I}_j^A}{\mathbf{s}_j^\top\mathbf{I}_j^A\mathbf{s}_j}\right){}^jX_i$$

括号内即 Schur 补，等效于"已知子关节会自由响应后"的有效惯量。

**三趟**：第一趟正向（速度/偏置力），第二趟反向（$\mathbf{I}^A$、偏置项），第三趟正向（加速度）。

---

### 4.5 质量矩阵 CRBA

**目标**：计算广义质量矩阵 $M(q)\in\mathbb{R}^{n_v\times n_v}$。

**复合刚体算法（Composite Rigid Body Algorithm）**：

$$M_{ij} = \mathbf{s}_i^\top\,\mathbf{I}_{\text{subtree}(j)}^c\,\mathbf{s}_j, \quad i \leq j$$

$\mathbf{I}_{\text{subtree}(j)}^c$ 为以 $j$ 为根的子树所有连杆的复合空间惯量，利用运动树的稀疏性只计算上三角。

$M$ 对称正定，可用 $LDL^\top$ 分解高效求逆，Pinocchio 提供 `cholesky` 模块专门处理。

---

### 4.6 质心动量学

这是**人形机器人控制的核心工具**，也是 Pinocchio 区别于其他库的特色之一。

**质心动量矩阵（CMM）** $A_g(q)\in\mathbb{R}^{6\times n_v}$：

$$\mathbf{h}_g = \begin{bmatrix} L_g \\ p \end{bmatrix} = A_g(q)\,\dot{q}$$

$L_g$：相对质心的系统角动量；$p = m\dot{c}$：系统线动量（$m$ 总质量，$c$ 质心）。

**整体牛顿-欧拉方程**（人形平衡控制的出发点）：

$$\dot{p} = \sum_k f_k^{contact} + mg$$

$$\dot{L}_g = \sum_k (p_k^{contact} - c) \times f_k^{contact} + \sum_k \tau_k^{contact}$$

这两个方程构成对接触力的线性约束，是 **ZMP 约束**、**接触力可行域（摩擦锥）** 的数学基础，
也是 MPC 质心轨迹优化的状态方程。

**API**：

```python
pin.computeCentroidalMomentum(model, data, q, v)
# → data.hg        : 质心动量 (6D)
# → data.dhg       : 质心动量时导数
# → data.com[0]    : 质心位置
# → data.vcom[0]   : 质心速度
# → data.Ag        : 质心动量矩阵 A_g (6 × nv)
```

---

### 4.7 接触约束动力学

**目标**：人形足底接触时，在约束下求 $\ddot{q}$ 和接触力 $\lambda$。

**最优化形式**：

$$\min_{\ddot{q}} \;\tfrac{1}{2}\|\ddot{q} - \ddot{q}_{free}\|_{M(q)}^2 \quad \text{s.t.}\quad J(q)\ddot{q} + \dot{J}(q,\dot{q})\dot{q} = 0$$

等价于求解 KKT（鞍点）系统：

$$\begin{bmatrix} M & J^\top \\ J & 0 \end{bmatrix} \begin{bmatrix} \ddot{q} \\ \lambda \end{bmatrix} = \begin{bmatrix} \tau - C\dot{q} - g \\ -\dot{J}\dot{q} \end{bmatrix}$$

$\lambda\in\mathbb{R}^{n_c}$ 即接触力（Lagrange 乘子），$J$ 为接触约束 Jacobian。

**高效求解**：Pinocchio 利用 **Contact Cholesky Decomposition** 充分利用运动树稀疏性，
比朴素高斯消元快数倍（具体实现见 `contact-cholesky.hpp`）。

**约束类型**：

| `ContactType` | 含义 | $n_c$ per contact |
|---------------|------|-------------------|
| `CONTACT_3D`  | 3D 点接触（位置约束）| 3 |
| `CONTACT_6D`  | 6D 面接触（SE(3) 约束）| 6 |

---

### 4.8 解析梯度

Pinocchio 提供所有主要算法的**解析偏导数**，这是 DDP/iLQR 等轨迹优化算法高效运行的基础。
相比数值差分，解析梯度更精确（不依赖步长选择）且通常快 10 倍以上。

**RNEA 梯度**（轨迹优化中 $\partial \tau / \partial (q,\dot{q},\ddot{q})$）：

$$\frac{\partial\tau}{\partial q},\quad \frac{\partial\tau}{\partial\dot{q}},\quad \frac{\partial\tau}{\partial\ddot{q}} = M(q)$$

**ABA 梯度**（动力学线性化，DDP 中的 $A_k, B_k$ 矩阵）：

$$\frac{\partial\ddot{q}}{\partial q} = A_k,\quad \frac{\partial\ddot{q}}{\partial\dot{q}},\quad \frac{\partial\ddot{q}}{\partial\tau} = M(q)^{-1} = B_k$$

**运动学梯度**（用于 Frame 速度/加速度对 $q, \dot{q}$ 的 Jacobian，WBC 必需）：

$$\frac{\partial v_i}{\partial q},\quad \frac{\partial a_i}{\partial q},\quad \frac{\partial a_i}{\partial\dot{q}},\quad \frac{\partial a_i}{\partial\ddot{q}}$$

---

## 5. Demo 注释：Python

### 5.1 `overview-simple.py` — 最小完整示例

```python
import pinocchio

# 创建内置 6-DOF 串联机械臂模型（无需 URDF）
model = pinocchio.buildSampleModelManipulator()
# 为该模型创建算法数据容器（存储中间计算结果）
data = model.createData()

# 生成零位配置 q=0（关节角全为0）
q = pinocchio.neutral(model)
# 速度与加速度置零
v = pinocchio.utils.zero(model.nv)
a = pinocchio.utils.zero(model.nv)

# 调用 RNEA：逆动力学，计算产生加速度 a 所需的关节力矩 τ
# 此处 a=0, v=0，故 τ 只反映重力补偿项 g(q)
tau = pinocchio.rnea(model, data, q, v, a)
print("tau = ", tau.T)
```

**涉及算法**：RNEA（逆动力学）。

---

### 5.2 `overview-urdf.py` — 从 URDF 加载并做正向运动学

```python
from pathlib import Path
from sys import argv
import pinocchio

pinocchio_model_dir = Path(__file__).parent.parent / "models"

# URDF 路径：默认使用 UR5，可命令行传入自定义路径
urdf_filename = (
    pinocchio_model_dir / "example-robot-data/robots/ur_description/urdf/ur5_robot.urdf"
    if len(argv) < 2
    else argv[1]
)

# 从 URDF 构建运动学/动力学模型
# 固定基机器人：不传第二个参数
# 浮动基（人形）：传入 pinocchio.JointModelFreeFlyer()
model = pinocchio.buildModelFromUrdf(urdf_filename)
print("model name: " + model.name)

# 创建算法数据容器，与 model 绑定
data = model.createData()

# 在配置空间中随机采样一个合法配置（自动满足关节限位）
q = pinocchio.randomConfiguration(model)
print(f"q: {q.T}")

# 正向运动学：递推计算所有关节在世界坐标系中的 SE(3) 变换
# 结果存入 data.oMi[i]
pinocchio.forwardKinematics(model, data, q)

# 打印每个关节的平移部分（世界系下的 3D 坐标）
for name, oMi in zip(model.names, data.oMi):
    print("{:<24} : {: .2f} {: .2f} {: .2f}".format(name, *oMi.translation.T.flat))
```

**涉及算法**：FK（正向运动学）。

---

### 5.3 `inverse-kinematics.py` — 6D 逆向运动学（SE3 误差）

```python
import numpy as np
import pinocchio
from numpy.linalg import norm, solve

model = pinocchio.buildSampleModelManipulator()
data = model.createData()

JOINT_ID = 6                                       # 末端关节 ID
# 目标位姿：旋转=单位阵（无旋转），平移=[1,0,1]
oMdes = pinocchio.SE3(np.eye(3), np.array([1.0, 0.0, 1.0]))

q    = pinocchio.neutral(model)  # 初始配置
eps  = 1e-4                       # 收敛阈值（误差 L2 范数）
IT_MAX = 1000                     # 最大迭代次数
DT   = 1e-1                       # 积分步长
damp = 1e-12                      # 阻尼系数（防止奇异）

i = 0
while True:
    # Step 1：正向运动学，更新 data.oMi[]
    pinocchio.forwardKinematics(model, data, q)

    # Step 2：计算局部系下的相对变换 iMd = T_i^{-1} * T_des
    iMd = data.oMi[JOINT_ID].actInv(oMdes)

    # Step 3：误差 = SE(3) 对数映射 → 6D 切向量（平移+旋转误差）
    # log(iMd) ∈ se(3) ≅ R^6，直接用作迭代误差
    err = pinocchio.log(iMd).vector

    if norm(err) < eps:
        success = True; break
    if i >= IT_MAX:
        success = False; break

    # Step 4：计算关节空间 Jacobian（在关节局部系下）
    J = pinocchio.computeJointJacobian(model, data, q, JOINT_ID)

    # Step 5：对 Jacobian 乘以 log 映射的 Jacobian 修正项 J_log6
    # 这是将 SE(3) 误差正确反向传播到关节空间所必需的链式法则
    J = -np.dot(pinocchio.Jlog6(iMd.inverse()), J)

    # Step 6：阻尼最小二乘 v = -J^T (JJ^T + λI)^{-1} e
    v = -J.T.dot(solve(J.dot(J.T) + damp * np.eye(6), err))

    # Step 7：在流形（李群）上积分，而非简单 q += v*DT
    # 对于浮动基，旋转部分是四元数，必须走李群积分
    q = pinocchio.integrate(model, q, v * DT)

    if not i % 10:
        print(f"{i}: error = {err.T}")
    i += 1

if success:
    print("Convergence achieved!")
else:
    print("Warning: algorithm did not converge")
print(f"\nresult: {q.flatten().tolist()}")
print(f"\nfinal error: {err.T}")
```

**涉及算法**：FK、几何 Jacobian、SE(3) 对数映射、$J_{\log_6}$、李群积分。

---

### 5.4 `inverse-kinematics-3d.py` — 纯位置 IK（3D 点跟踪）

```python
import numpy as np
import pinocchio
from numpy.linalg import norm, solve

model = pinocchio.buildSampleModelManipulator()
data = model.createData()

JOINT_ID = 6
oMdes = pinocchio.SE3(np.eye(3), np.array([1.0, 0.0, 1.0]))

q    = pinocchio.neutral(model)
eps  = 1e-4
IT_MAX = 1000
DT   = 1e-1
damp = 1e-12

it = 0
while True:
    pinocchio.forwardKinematics(model, data, q)
    iMd = data.oMi[JOINT_ID].actInv(oMdes)

    # 与 6D IK 的唯一区别：误差只取平移部分（R^3 而非 R^6）
    # 不关心末端姿态，只控制末端位置
    err = iMd.translation

    if norm(err) < eps:
        success = True; break
    if it >= IT_MAX:
        success = False; break

    J = pinocchio.computeJointJacobian(model, data, q, JOINT_ID)
    # 只取 Jacobian 的上 3 行（线速度部分）并取负号（局部系误差的符号约定）
    J = -J[:3, :]

    # 3×3 阻尼最小二乘
    v = -J.T.dot(solve(J.dot(J.T) + damp * np.eye(3), err))
    q = pinocchio.integrate(model, q, v * DT)

    if not it % 10:
        print(f"{it}: error = {err.T}")
    it += 1

if success:
    print("Convergence achieved!")
else:
    print("Warning: algorithm did not converge")
print(f"\nresult: {q.flatten().tolist()}")
print(f"\nfinal error: {err.T}")
```

**与 5.3 的对比**：不使用 SE(3) 对数映射，只控制末端位置；Jacobian 只取线速度行；
无需 $J_{\log_6}$ 修正（因为平移部分的误差本身是线性的）。

---

### 5.5 `inverse-dynamics.py` — RNEA 逆动力学

```python
from pathlib import Path
import numpy as np
import pinocchio as pin

pinocchio_model_dir = Path(__file__).parent.parent / "models/"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir = pinocchio_model_dir
urdf_filename = "ur5_robot.urdf"
urdf_model_path = model_path / "ur_description/urdf/" / urdf_filename

# buildModelsFromUrdf 同时返回运动学模型、碰撞模型、视觉模型
model, _, _ = pin.buildModelsFromUrdf(urdf_model_path, package_dirs=mesh_dir)
data = model.createData()

# 随机采样状态（实际使用时替换为真实状态）
q = pin.randomConfiguration(model)   # 关节角 [rad]
v = np.random.rand(model.nv, 1)       # 关节速度 [rad/s]
a = np.random.rand(model.nv, 1)       # 期望关节加速度 [rad/s²]

# RNEA：递推牛顿-欧拉，计算产生加速度 a 所需的力矩 τ
# 等价于：τ = M(q)a + C(q,v)v + g(q)
# 但用 O(n) 递推而非显式构造 M, C, g
tau = pin.rnea(model, data, q, v, a)

# 结果存在 data.tau，rnea() 也直接返回引用
print("Joint torques: " + str(tau))
```

**应用场景**：力矩前馈控制、逆动力学控制（IDC）、轨迹可行性验证。

---

### 5.6 `forward-dynamics-derivatives.py` — ABA 正向动力学 + 解析梯度

```python
import numpy as np
import pinocchio as pin

# 使用内置随机人形模型（双足，带浮动基）
model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

# 设置关节限位（样例模型默认无限位，优化器需要）
model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit =  np.ones((model.nq, 1))

q   = pin.randomConfiguration(model)         # 关节配置
v   = np.random.rand(model.nv, 1)            # 关节速度
tau = np.random.rand(model.nv, 1)            # 关节力矩

# 计算 ABA 及其对 (q, v, τ) 的解析梯度
# 内部调用三趟递推，一次性得到加速度和所有梯度
pin.computeABADerivatives(model, data, q, v, tau)

# ∂q̈/∂q：状态 q 对加速度的影响（DDP 中的 ∂f/∂q 的 q 块）
ddq_dq  = data.ddq_dq

# ∂q̈/∂v：状态 v 对加速度的影响（DDP 中的 ∂f/∂q 的 v 块）
ddq_dv  = data.ddq_dv

# ∂q̈/∂τ = M(q)^{-1}：控制输入对加速度的影响（DDP 中的 B 矩阵）
ddq_dtau = data.Minv

# 应用：DDP/iLQR 中的线性化
# f(x,u) = q̈  =>  f_x ≈ [ddq_dq, ddq_dv],  f_u = Minv
```

**应用场景**：Crocoddyl DDP 的动力学线性化、梯度下降轨迹优化。

---

### 5.7 `inverse-dynamics-derivatives.py` — RNEA 梯度

```python
import numpy as np
import pinocchio as pin

model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit =  np.ones((model.nq, 1))

q   = pin.randomConfiguration(model)
v   = np.random.rand(model.nv, 1)
a   = np.random.rand(model.nv, 1)

# 计算 RNEA 梯度：∂τ/∂q，∂τ/∂v，∂τ/∂a
pin.computeRNEADerivatives(model, data, q, v, a)

# ∂τ/∂q：配置对所需力矩的影响（重力梯度 + 科氏力梯度）
dtau_dq = data.dtau_dq

# ∂τ/∂v：速度对所需力矩的影响（科氏力 + 离心力梯度）
dtau_dv = data.dtau_dv

# ∂τ/∂a = M(q)：加速度对力矩的影响（即质量矩阵本身）
dtau_da = data.M
```

**应用场景**：逆动力学控制的线性化、运动学参数辨识的梯度计算。

---

### 5.8 `kinematics-derivatives.py` — 运动学梯度

```python
import numpy as np
import pinocchio as pin

model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit =  np.ones((model.nq, 1))

q = pin.randomConfiguration(model)
v = np.random.rand(model.nv, 1)
a = np.random.rand(model.nv, 1)

# 一次调用计算所有关节速度/加速度关于 (q,v,a) 的梯度
pin.computeForwardKinematicsDerivatives(model, data, q, v, a)

joint_name = "rleg6_joint"
joint_id = model.getJointId(joint_name)

# ∂v_i/∂q, ∂v_i/∂v（即几何 Jacobian）——世界系下
(dv_dq, dv_dv) = pin.getJointVelocityDerivatives(
    model, data, joint_id, pin.ReferenceFrame.WORLD
)
# 局部系下
(dv_dq_local, dv_dv_local) = pin.getJointVelocityDerivatives(
    model, data, joint_id, pin.ReferenceFrame.LOCAL
)

# ∂a_i/∂q, ∂a_i/∂v, ∂a_i/∂a——世界系下
(dv_dq, da_dq, da_dv, da_da) = pin.getJointAccelerationDerivatives(
    model, data, joint_id, pin.ReferenceFrame.WORLD
)
# 注意：∂v_i/∂v = ∂a_i/∂a（都等于几何 Jacobian）
```

**应用场景**：WBC 中 Frame 速度的 Jacobian、末端加速度对力矩的映射。

---

### 5.9 `simulation-contact-dynamics.py` — Talos 人形双足接触仿真

这是最重要的人形机器人 demo，演示浮动基 + 刚性接触 + 姿态控制的完整闭环仿真。

```python
import math, sys, time
from pathlib import Path
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir   = pinocchio_model_dir
urdf_filename = "talos_reduced.urdf"
srdf_filename = "talos.srdf"

# 关键：浮动基 (FreeFlyer) 使机器人基座在空间中自由运动
# 不加此参数则基座固定，无法模拟真实人形
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    model_path / "talos_data/robots" / urdf_filename,
    mesh_dir,
    pin.JointModelFreeFlyer()   # 人形必须！
)

# MeshCat 可视化器（浏览器中打开 URL 查看）
viz = MeshcatVisualizer(model, collision_model, visual_model)
viz.initViewer(open=False)
viz.loadViewerModel()

# 从 SRDF 加载参考配置（half_sitting：双腿微蹲标准站立姿）
pin.loadReferenceConfigurations(model, model_path / "talos_data/srdf" / srdf_filename)
q0 = model.referenceConfigurations["half_sitting"]
viz.display(q0)

# ---- 建立足底接触约束模型 ----
feet_name = ["left_sole_link", "right_sole_link"]
frame_ids  = [model.getFrameId(n) for n in feet_name]

contact_models = []
for frame_id in frame_ids:
    frame = model.frames[frame_id]
    # CONTACT_6D：6D 刚性接触（位置+姿态全约束，适合平足）
    contact_model = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_6D,
        model,
        frame.parentJoint,    # 接触发生在哪个关节
        frame.placement       # 接触点相对关节的位置
    )
    contact_models.append(contact_model)
    contact_datas = [cm.createData() for cm in contact_models]

num_constraints = len(feet_name)
contact_dim = 6 * num_constraints   # 每个 6D 接触提供 6 个约束

# 预分配接触约束内存
pin.initConstraintDynamics(model, data_sim, contact_models)

# ---- 仿真参数 ----
dt = 5e-3    # 5 ms 时间步（200 Hz）
# S：选择矩阵，将全身 nv 维度映射到驱动关节（浮动基 6 DOF 不可驱动）
S = np.zeros((model.nv - 6, model.nv))
S.T[6:, :] = np.eye(model.nv - 6)

# 姿态稳定控制器参数（PD 增益）
Kp_posture = 30.0
Kv_posture = 0.05 * math.sqrt(Kp_posture)  # 临界阻尼

q, v = q0.copy(), np.zeros(model.nv)
tau  = np.zeros(model.nv)

while t <= T:
    # ---- 控制律：关节空间 PD 稳定到参考姿态 ----
    pin.computeJointJacobians(model, data_control, q)

    # tau = S^T * Kp*(q_ref - q) - Kv*v  （只驱动非浮动基自由度）
    tau_control = S.T.dot(
        Kp_posture * pin.difference(model, q, q_ref) - Kv_posture * v[6:]
    )

    # ---- 带接触约束的正向动力学 ----
    # 求解 KKT 系统：
    # [M  J^T] [q̈  ] = [τ - C(q,v)v - g(q)]
    # [J  0  ] [λ  ]   [-J_dot * v          ]
    a = pin.constraintDynamics(
        model, data_sim, q, v, tau_control, contact_models, contact_datas
    )
    # a      → data_sim.ddq     (关节加速度)
    # lambda → data_sim.lambda_c (接触力)

    # ---- 半隐式欧拉积分 ----
    v += a * dt                          # 速度更新（用新加速度）
    q = pin.integrate(model, q, v * dt)  # 配置更新（流形积分！）

    viz.display(q)
    t += dt
```

**涉及算法**：浮动基建模、SRDF 参考配置、`RigidConstraintModel`、`constraintDynamics`（KKT 求解）、半隐式欧拉积分。

---

### 5.10 `talos-simulation.py` — KKT 约束姿态优化（静态）

```python
from time import sleep
import example_robot_data, numpy as np, pinocchio
from pinocchio.visualize import GepettoVisualizer

# 使用 example_robot_data 加载完整 Talos 模型
robot = example_robot_data.load("talos")
model, data = robot.model, robot.data

# 参考配置
q = robot.model.referenceConfigurations["half_sitting"]

# ---- 四肢（双足+双手）的 6D 接触约束 ----
foot_frames = ["left_sole_link", "right_sole_link",
               "gripper_left_fingertip_3_link", "gripper_right_fingertip_3_link"]
constraint_models = []
for frame_name in foot_frames:
    frame_id   = model.getFrameId(frame_name)
    joint_id   = model.frames[frame_id].parent
    contact_model = pinocchio.RigidConstraintModel(
        pinocchio.ContactType.CONTACT_6D, model,
        joint_id,
        model.frames[frame_id].placement,
        0,                    # 世界作为第二个约束体
        data.oMf[frame_id]    # 当前帧在世界系下的位姿
    )
    constraint_models.append(contact_model)

# 修改手部目标位置（抬起双臂的目标）
constraint_models[3].joint2_placement = pinocchio.SE3(
    pinocchio.rpy.rpyToMatrix(np.array([0., -np.pi/2, 0.])),
    np.array([0.6, -0.40, 1.0])
)

constraint_datas = [cm.createData() for cm in constraint_models]
kkt_constraint = pinocchio.ContactCholeskyDecomposition(model, constraint_models)
constraint_dim  = sum(cm.size() for cm in constraint_models)

# ---- squashing：迭代求解约束满足的关节配置 ----
# 在保持足底接触的同时，最小化质心高度变化
# KKT 求解方程：
#   [M   J^T] [dq ] = [rhs]
#   [J   0  ] [dy ]
def squashing(model, data, q_in):
    q = q_in.copy()
    y = np.ones(constraint_dim)   # Lagrange 乘子（接触力对偶变量）

    for k in range(N):
        pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))
        pinocchio.computeJointJacobians(model, data, q)

        com_act = data.com[0].copy()
        com_err = com_act - com_des(k)

        # 构造 Contact Cholesky 分解（稀疏 KKT 求解器）
        kkt_constraint.compute(model, data, constraint_models, constraint_datas, mu)

        # 约束残差：log6(c1Mc2) 衡量当前位姿与约束目标的 SE(3) 误差
        constraint_value = np.concatenate(
            [pinocchio.log6(cd.c1Mc2).vector for cd in constraint_datas]
        )
        J = np.vstack([
            pinocchio.getFrameJacobian(model, data, cm.joint1_id,
                                       cm.joint1_placement, cm.reference_frame)
            for cm in constraint_models
        ])

        # RHS：约束残差 + 质心误差 + 关节自由度为零
        rhs = np.concatenate([
            -constraint_value - y * mu,
            kp * mass * com_err,
            np.zeros(model.nv - 3)
        ])
        dz = kkt_constraint.solve(rhs)   # 稀疏 KKT 求解
        dy, dq = dz[:constraint_dim], dz[constraint_dim:]

        q = pinocchio.integrate(model, q, -alpha * dq)
        y -= alpha * (-dy + y)
        robot.display(q)
        sleep(0.05)
    return q
```

**涉及算法**：`ContactCholeskyDecomposition`、`computeAllTerms`、SE(3) 约束残差、质心控制。

---

### 5.11 `static-contact-dynamics.py` — 静态接触力计算

```python
# 目标：给定静止姿态 q0，计算维持静止所需的关节力矩和接触力
# 静态方程：g(q) = τ + J_c^T λ

# 策略：分块求解
# 对浮动基的 6 个方程（无驱动）：
#   g_bl = J_c_bl^T * λ  →  λ = pinv(J_c_bl^T) * g_bl
# 对驱动关节：
#   τ = g_j - J_c_j^T * λ

import numpy as np
import pinocchio as pin

# 加载 Solo-12 四足机器人（浮动基，4 个 3D 点接触）
model, _, _ = pin.buildModelsFromUrdf(urdf_model_path, mesh_dir, pin.JointModelFreeFlyer())
data = model.createData()

# 计算各项动力学量
pin.computeAllTerms(model, data, q0, np.zeros(model.nv))

# 计算各足 Jacobian（只取前 6 行对应浮动基）
Jc_feet_bl = ...   # 从 getFrameJacobian 提取

# 求接触力（伪逆）
lambda_contact = np.linalg.pinv(Jc_feet_bl.T) @ g_bl

# 求关节力矩
tau = g_j - Jc_feet_j.T @ lambda_contact
```

**涉及算法**：`computeAllTerms`、`getFrameJacobian`、静态平衡方程、伪逆。

---

### 5.12 `simulation-pendulum.py` — 倒立摆/多摆仿真

```python
# 从零手工构建多摆模型（无 URDF），演示 API 编程建模

model = pin.Model()
geom_model = pin.GeometryModel()
parent_id = 0

for k in range(N):
    # 添加绕 X 轴旋转的关节 (JointModelRX)
    joint_id = model.addJoint(parent_id, pin.JointModelRX(), joint_placement, f"joint_{k+1}")
    # 为关节绑定连杆惯量
    model.appendBodyToJoint(joint_id, body_inertia, body_placement)
    parent_id = joint_id

# 仿真主循环：ABA + 半隐式欧拉
for k in range(N_steps):
    tau_control = -damping * v          # 阻尼控制（耗散能量）
    a = pin.aba(model, data_sim, q, v, tau_control)  # ABA 正向动力学
    v += a * dt                         # 半隐式：先更新速度
    q = pin.integrate(model, q, v * dt) # 再用新速度更新配置（流形积分）
    viz.display(q)
```

**可选参数**：`--with-cart`（小车-摆系统），`-N 3`（三连摆）。  
**涉及算法**：手工建模（`addJoint` / `appendBodyToJoint`）、ABA、半隐式欧拉。

---

### 5.13 `simulation-closed-kinematic-chains.py` — 闭环运动链仿真

```python
# 演示四杆机构（closed-loop mechanism）的约束仿真
# 闭环约束通过 RigidConstraintModel 实现

# 建模：先建开链，再用约束封闭
model = pin.Model()
# ... 添加连杆 A1, B1, B2, A2 ...

# 闭环约束：要求两个关节的相对位姿满足特定条件
loop_constraint = pin.RigidConstraintModel(
    pin.ContactType.CONTACT_6D,
    model, joint1_id, placement1,
         joint2_id, placement2
)

# 仿真：带约束的正向动力学
a = pin.constraintDynamics(model, data, q, v, tau, [loop_constraint], constraint_datas)
```

**人形应用**：并联腿结构（如 Atlas 的膝关节四杆机构）即需此类约束建模。

---

### 5.14 `run-algo-in-parallel.py` — OpenMP 并行批量动力学

```python
import numpy as np
import pinocchio as pin

# 人形随机模型
model = pin.buildSampleModelHumanoid()
pool = pin.ModelPool(model)      # 内部管理多个 Data 副本，每线程一个

num_threads = pin.omp_get_max_threads()
batch_size  = 128                # 一次并行计算 128 个配置

# 准备批量输入（每列是一个配置/速度/加速度）
q   = np.column_stack([pin.randomConfiguration(model) for _ in range(batch_size)])
v   = np.zeros((model.nv, batch_size))
a   = np.zeros((model.nv, batch_size))
tau = np.zeros((model.nv, batch_size))

# 并行 RNEA：一次调用计算 128 组逆动力学（OpenMP 自动分配线程）
res_rnea = pin.rneaInParallel(num_threads, pool, q, v, a)

# 并行 ABA：一次调用计算 128 组正向动力学
res_aba  = pin.abaInParallel(num_threads, pool, q, v, tau)
```

**应用场景**：MPC 中的多轨迹并行仿真、强化学习中的批量仿真采样。

---

### 5.15 `build-reduced-model.py` — 降阶模型（锁定关节）

```python
import numpy as np
import pinocchio as pin

# 加载完整模型
model, collision_model, visual_model = pin.buildModelsFromUrdf(urdf_filename, mesh_dir)

# 要锁定（固定）的关节名称列表
jointsToLock = ["wrist_1_joint", "wrist_2_joint", "wrist_3_joint"]
jointsToLockIDs = [model.getJointId(jn) for jn in jointsToLock if model.existJointName(jn)]

# 锁定位置的参考配置
initialJointConfig = np.array([0, 0, 0, 1, 1, 1])

# Option 1：只降阶运动学模型
model_reduced = pin.buildReducedModel(model, jointsToLockIDs, initialJointConfig)

# Option 2：同时降阶视觉模型（显示用）
model_reduced, visual_model_reduced = pin.buildReducedModel(
    model, visual_model, jointsToLockIDs, initialJointConfig
)

# Option 3：同时降阶多个几何模型
model_reduced, (visual_model_reduced, collision_model_reduced) = pin.buildReducedModel(
    model,
    list_of_geom_models=[visual_model, collision_model],
    list_of_joints_to_lock=jointsToLockIDs,
    reference_configuration=initialJointConfig,
)
```

**人形应用**：将手臂锁定在某姿态，专注腿部控制研究；或将整个上身锁定，
验证纯腿部动力学。

---

### 5.16 `mimic_dynamics.py` — 关节 Mimic（传动比约束）

```python
import numpy as np
import pinocchio

# 从 URDF 直接解析 mimic 标签（如夹爪耦合关节）
model_mimic = pinocchio.buildModelFromUrdf(model_path, mimic=True)

# 或手工设置 mimic 关系：joint 10 = -1.0 × joint 9 + 0.0
model_mimic = pinocchio.transformJointIntoMimic(model_full, 9, 10, -1.0, 0.0)

# G 矩阵：将 mimic 模型的速度映射到完整模型的速度
# 对于非 mimic 关节：G[i,i]=1；对 mimic 关节：G[i,j]=scaling
G = np.zeros([model_mimic.nv, model_full.nv])
for i in range(model_full.njoints):
    scaling = getattr(model_mimic.joints[i].extract(), "scaling", None)
    G[model_mimic.joints[i].idx_v, model_full.joints[i].idx_v] = scaling if scaling else 1

# 关系：τ_mimic = G @ τ_full，M_mimic = G @ M_full @ G^T
tau_mimic = pinocchio.rnea(model_mimic, data_mimic, q_mimic, v_mimic, a_mimic)
M_mimic   = pinocchio.crba(model_mimic, data_mimic, q_mimic)
```

**人形应用**：人形机器人手部（夹爪）、踝关节并联结构中常有传动比耦合关节。

---

### 5.17 `contact-cholesky.py` — 接触 Cholesky 分解（Delassus 算子）

```python
import pinocchio as pin

# ANYmal 四足（4 个 3D 点接触）
model = pin.buildModelFromUrdf(urdf_filename)
data  = model.createData()
q0    = pin.neutral(model)

# 构造接触约束（3D 点接触）
contact_models = []
for fid in feet_frame_ids:
    frame = model.frames[fid]
    cmodel = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,
        frame.parent, frame.placement,
        pin.LOCAL_WORLD_ALIGNED,
    )
    contact_models.append(cmodel)

contact_data = [cm.createData() for cm in contact_models]

# 计算质量矩阵（CRBA）
pin.crba(model, data, q0)

# Contact Cholesky 分解：高效求解 KKT 系统
# 利用运动树稀疏性，比 dense LU 快数倍
data.contact_chol.compute(model, data, contact_models, contact_data)

# Delassus 矩阵（操作空间惯量矩阵逆）：Λ^{-1} = J M^{-1} J^T
delassus_matrix     = data.contact_chol.getInverseOperationalSpaceInertiaMatrix()
# 操作空间惯量矩阵：Λ = (J M^{-1} J^T)^{-1}
delassus_matrix_inv = data.contact_chol.getOperationalSpaceInertiaMatrix()
```

**Delassus 矩阵物理含义**：接触空间中的有效质量（等效惯量），用于求接触力响应。

---

## 6. Demo 注释：C++

### 6.1 `overview-simple.cpp` — C++ 最小示例

```cpp
#include <iostream>
#include "pinocchio/multibody/sample-models.hpp"   // 内置样例模型
#include "pinocchio/algorithm/joint-configuration.hpp"  // neutral()
#include "pinocchio/algorithm/rnea.hpp"            // rnea()

int main()
{
    pinocchio::Model model;
    // 构建内置 6-DOF 串联机械臂（等同 Python buildSampleModelManipulator）
    pinocchio::buildModels::manipulator(model);
    pinocchio::Data data(model);   // 算法数据容器

    // 零位配置（关节全零）
    Eigen::VectorXd q = pinocchio::neutral(model);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);  // 零速度
    Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);  // 零加速度

    // RNEA 逆动力学：a=0 时计算重力补偿力矩 g(q)
    const Eigen::VectorXd & tau = pinocchio::rnea(model, data, q, v, a);
    std::cout << "tau = " << tau.transpose() << std::endl;
}
```

---

### 6.2 `overview-urdf.cpp` — 从 URDF 加载并做正向运动学

```cpp
#include "pinocchio/parsers/urdf.hpp"              // URDF 解析器
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"      // forwardKinematics()
#include <iostream>

int main(int argc, char ** argv)
{
    using namespace pinocchio;

    // URDF 路径（可命令行传入，默认 UR5）
    const std::string urdf_filename = (argc <= 1)
        ? PINOCCHIO_MODEL_DIR + std::string(".../ur5_robot.urdf")
        : argv[1];

    Model model;
    // 加载 URDF → 构建运动学/动力学模型
    // 浮动基人形需改为：urdf::buildModel(urdf, JointModelFreeFlyer(), model)
    pinocchio::urdf::buildModel(urdf_filename, model);
    std::cout << "model name: " << model.name << std::endl;

    Data data(model);

    // 在关节限位内随机采样配置（四元数部分自动归一化）
    Eigen::VectorXd q = randomConfiguration(model);

    // 正向运动学：递推更新 data.oMi[i]（所有关节的 SE(3) 位姿）
    forwardKinematics(model, data, q);

    // 打印各关节在世界系下的平移坐标
    for (JointIndex id = 0; id < (JointIndex)model.njoints; ++id)
        std::cout << model.names[id] << ": "
                  << data.oMi[id].translation().transpose() << std::endl;
}
```

---

### 6.3 `overview-SE3.cpp` — SE(3) 李群运算

```cpp
#include "pinocchio/multibody/liegroup/liegroup.hpp"
using namespace pinocchio;

int main()
{
    // SE(3) 操作对象（3D 特殊欧氏群）
    // 配置向量格式：[x, y, z, qx, qy, qz, qw]（7 维，四元数表示旋转）
    typedef SpecialEuclideanOperationTpl<3, double> SE3Operation;
    SE3Operation aSE3;

    SE3Operation::ConfigVector_t pose_s, pose_g;
    // 起始位姿：位置 [1,1,1]，旋转用四元数表示
    pose_s << 1.0, 1.0, 1.0, -0.13795, 0.13795, 0.69352, 0.69352;
    // 目标位姿
    pose_g << 4.0, 3.0, 3.0, -0.46194, 0.331414, 0.800103, 0.191342;

    // 归一化四元数（确保旋转合法）
    aSE3.normalize(pose_s);
    aSE3.normalize(pose_g);

    // difference：在切空间（李代数）中计算两个 SE(3) 元素的差
    // delta_u ∈ R^6（平移 3 + 旋转 3），是 pose_g 相对 pose_s 的"速度"
    SE3Operation::TangentVector_t delta_u;
    aSE3.difference(pose_s, pose_g, delta_u);

    // integrate：从 pose_s 沿 delta_u 做李群积分，应恢复 pose_g
    SE3Operation::ConfigVector_t pose_check;
    aSE3.integrate(pose_s, delta_u, pose_check);
    // 验证：pose_check ≈ pose_g
}
```

**数学含义**：`difference(a, b)` = $\log(a^{-1} \cdot b)$；`integrate(a, u)` = $a \cdot \exp(u)$。  
**人形应用**：浮动基状态的差分（用于误差计算）和积分（用于状态更新）。

---

### 6.4 `overview-lie.cpp` — SE(2) 李群（平面运动）

```cpp
#include "pinocchio/multibody/liegroup/liegroup.hpp"
using namespace pinocchio;

int main()
{
    // SE(2)：2D 特殊欧氏群（平面移动机器人、足底在地面上的运动）
    // 配置向量：[x, y, cos(θ), sin(θ)]（4 维）
    typedef SpecialEuclideanOperationTpl<2, double> SE2Operation;
    SE2Operation aSE2;

    SE2Operation::ConfigVector_t pose_s, pose_g;
    // 起始：位置 [1,1]，朝向 45°
    pose_s << 1.0, 1.0, cos(M_PI/4.0), sin(M_PI/4.0);
    // 目标：位置 [3,-1]，朝向 -90°
    pose_g << 3.0, -1.0, cos(-M_PI/2.0), sin(-M_PI/2.0);

    // 在李代数（切空间）中计算差分
    // delta_u ∈ R^3（Δx, Δy, Δθ），代表从 pose_s 到 pose_g 的最短路径
    SE2Operation::TangentVector_t delta_u;
    delta_u.setZero();
    aSE2.difference(pose_s, pose_g, delta_u);

    // 验证：integrate(pose_s, delta_u) ≈ pose_g
    SE2Operation::ConfigVector_t pose_check;
    aSE2.integrate(pose_s, delta_u, pose_check);
}
```

**人形应用**：足步规划中步态在地面投影的差分；ZMP 轨迹在水平面的表示。

---

### 6.5 `inverse-kinematics.cpp` — C++ 6D IK

```cpp
#include "pinocchio/spatial/explog.hpp"      // log6(), Jlog6()
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"   // computeJointJacobian()
#include "pinocchio/algorithm/joint-configuration.hpp"

int main(int, char**)
{
    pinocchio::Model model;
    pinocchio::buildModels::manipulator(model);
    pinocchio::Data data(model);

    const int JOINT_ID = 6;
    // 目标位姿：单位旋转 + 平移 [1, 0, 1]
    const pinocchio::SE3 oMdes(Eigen::Matrix3d::Identity(), Eigen::Vector3d(1., 0., 1.));

    Eigen::VectorXd q = pinocchio::neutral(model);
    const double eps = 1e-4, DT = 1e-1, damp = 1e-6;

    // 6×nv Jacobian 矩阵预分配
    pinocchio::Data::Matrix6x J(6, model.nv);
    J.setZero();

    typedef Eigen::Matrix<double, 6, 1> Vector6d;
    Vector6d err;
    Eigen::VectorXd v(model.nv);

    for (int i = 0;; i++)
    {
        pinocchio::forwardKinematics(model, data, q);

        // 局部系下的相对变换：iMd = T_i^{-1} * T_des
        const pinocchio::SE3 iMd = data.oMi[JOINT_ID].actInv(oMdes);

        // 6D 误差：SE(3) 对数映射 → 切向量
        err = pinocchio::log6(iMd).toVector();
        if (err.norm() < eps) { success = true; break; }
        if (i >= IT_MAX)      { success = false; break; }

        // 计算局部系下的几何 Jacobian
        pinocchio::computeJointJacobian(model, data, q, JOINT_ID, J);

        // Jlog6：log 映射关于 SE(3) 参数的 Jacobian（6×6 矩阵）
        // 链式法则：J_eff = Jlog6(iMd^{-1}) * J_joint
        pinocchio::Data::Matrix6 Jlog;
        pinocchio::Jlog6(iMd.inverse(), Jlog);
        J = -Jlog * J;

        // 阻尼最小二乘：v = -J^T (JJ^T + λI)^{-1} e
        // JJt.ldlt().solve() 用 LDLT 分解高效求解对称正定线性系统
        pinocchio::Data::Matrix6 JJt;
        JJt.noalias() = J * J.transpose();
        JJt.diagonal().array() += damp;
        v.noalias() = -J.transpose() * JJt.ldlt().solve(err);

        // 李群积分（关键！不能用 q += v*DT）
        q = pinocchio::integrate(model, q, v * DT);
    }
    std::cout << "result: " << q.transpose() << std::endl;
}
```

---

### 6.6 `inverse-kinematics-3d.cpp` — C++ 3D 位置 IK

```cpp
// 与 6D IK 相比，只跟踪末端位置（不控制姿态）
// 误差降为 3D，Jacobian 只取线速度行

for (int i = 0;; i++)
{
    pinocchio::forwardKinematics(model, data, q);
    const pinocchio::SE3 iMd = data.oMi[JOINT_ID].actInv(oMdes);

    // 误差 = 相对变换的平移部分（R^3，非 R^6）
    err = iMd.translation();
    if (err.norm() < eps) { success = true; break; }

    pinocchio::computeJointJacobian(model, data, q, JOINT_ID, joint_jacobian);
    // 只取 Jacobian 上 3 行（线速度），取负（局部系符号约定）
    const auto J = -joint_jacobian.topRows<3>();

    // 3×3 阻尼最小二乘（比 6×6 计算量更小）
    const Eigen::Matrix3d JJt = J * J.transpose() + damp * Eigen::Matrix3d::Identity();
    v.noalias() = -J.transpose() * JJt.ldlt().solve(err);
    q = pinocchio::integrate(model, q, v * DT);
}
```

---

### 6.7 `inverse-dynamics.cpp` — C++ RNEA

```cpp
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

int main(int argc, char** argv)
{
    using namespace pinocchio;

    // 加载 URDF 模型（固定基）
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    // 随机状态
    Eigen::VectorXd q   = randomConfiguration(model);  // 随机关节角 [rad]
    Eigen::VectorXd v   = Eigen::VectorXd::Zero(model.nv); // 零速度
    Eigen::VectorXd a   = Eigen::VectorXd::Zero(model.nv); // 零加速度

    // RNEA：τ = M(q)a + C(q,v)v + g(q)
    // a=0, v=0 时结果即为重力补偿力矩 g(q)
    Eigen::VectorXd tau = pinocchio::rnea(model, data, q, v, a);
    // 结果同时存在 data.tau
    std::cout << "Joint torques: " << data.tau.transpose() << std::endl;
}
```

---

### 6.8 `inverse-dynamics-derivatives.cpp` — C++ RNEA 梯度

```cpp
#include "pinocchio/algorithm/rnea-derivatives.hpp"

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    Eigen::VectorXd q = randomConfiguration(model);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
    Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);

    // 预分配梯度矩阵（nv × nv）
    Eigen::MatrixXd dtau_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
    Eigen::MatrixXd dtau_dv = Eigen::MatrixXd::Zero(model.nv, model.nv);
    Eigen::MatrixXd dtau_da = Eigen::MatrixXd::Zero(model.nv, model.nv);

    // 一次调用同时计算 RNEA 结果和三个梯度矩阵
    // dtau_da = M(q)（质量矩阵，同时获得正向运动学结果）
    computeRNEADerivatives(model, data, q, v, a, dtau_dq, dtau_dv, dtau_da);

    // data.tau 存 RNEA 结果，dtau_dq/dv/da 存梯度
    std::cout << "Joint torque: " << data.tau.transpose() << std::endl;
}
```

---

### 6.9 `forward-dynamics-derivatives.cpp` — C++ ABA 梯度

```cpp
#include "pinocchio/algorithm/aba-derivatives.hpp"

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    Eigen::VectorXd q   = randomConfiguration(model);
    Eigen::VectorXd v   = Eigen::VectorXd::Zero(model.nv);
    Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);

    // 预分配梯度矩阵
    Eigen::MatrixXd ddq_dq   = Eigen::MatrixXd::Zero(model.nv, model.nv);
    Eigen::MatrixXd ddq_dv   = Eigen::MatrixXd::Zero(model.nv, model.nv);
    Eigen::MatrixXd ddq_dtau = Eigen::MatrixXd::Zero(model.nv, model.nv);

    // 一次调用：计算 ABA（正向动力学）及其对 (q, v, τ) 的解析梯度
    // ddq_dtau = M(q)^{-1}（质量矩阵的逆，即 DDP 中的 B 矩阵）
    computeABADerivatives(model, data, q, v, tau, ddq_dq, ddq_dv, ddq_dtau);

    std::cout << "Joint acceleration: " << data.ddq.transpose() << std::endl;
    // DDP 线性化：A = [ddq_dq, ddq_dv], B = ddq_dtau = Minv
}
```

---

### 6.10 `kinematics-derivatives.cpp` — C++ 运动学梯度

```cpp
#include "pinocchio/algorithm/kinematics-derivatives.hpp"

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    Eigen::VectorXd q = randomConfiguration(model);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
    Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);

    // 一次调用计算所有关节速度/加速度关于 (q,v,a) 的梯度
    computeForwardKinematicsDerivatives(model, data, q, v, a);

    // 取最后一个关节（末端）
    JointIndex joint_id = (JointIndex)(model.njoints - 1);
    Data::Matrix6x v_partial_dq(6, model.nv), a_partial_dq(6, model.nv),
                   a_partial_dv(6, model.nv), a_partial_da(6, model.nv);
    v_partial_dq.setZero(); a_partial_dq.setZero();
    a_partial_dv.setZero(); a_partial_da.setZero();

    // 局部系下：∂v_i/∂q, ∂a_i/∂q, ∂a_i/∂v, ∂a_i/∂a（=几何 Jacobian）
    getJointAccelerationDerivatives(
        model, data, joint_id, LOCAL,
        v_partial_dq, a_partial_dq, a_partial_dv, a_partial_da
    );

    // 世界系下：同上，但量在世界坐标系中表达
    getJointAccelerationDerivatives(
        model, data, joint_id, WORLD,
        v_partial_dq, a_partial_dq, a_partial_dv, a_partial_da
    );
    // 注意：v_partial_dv = a_partial_da（两者相等，不重复计算）
}
```

---

### 6.11 `build-reduced-model.cpp` — C++ 降阶模型

```cpp
#include "pinocchio/algorithm/model.hpp"   // buildReducedModel()

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // 方案一：指定要锁定的关节名称列表
    std::vector<std::string> joints_to_lock = {
        "elbow_joint", "wrist_3_joint", "wrist_2_joint"
    };
    std::vector<JointIndex> lock_ids;
    for (auto& name : joints_to_lock)
        if (model.existJointName(name))
            lock_ids.push_back(model.getJointId(name));

    // 锁定参考配置（从完整模型中随机采样）
    Eigen::VectorXd q_rand = randomConfiguration(model);

    // 构建降阶模型：被锁定的关节固定在 q_rand 对应的角度
    Model reduced_model = buildReducedModel(model, lock_ids, q_rand);

    // 方案二：指定要保留（不锁定）的关节
    // 先计算补集，再调用 buildReducedModel
    std::vector<std::string> joints_to_keep = {
        "shoulder_pan_joint", "shoulder_lift_joint", "wrist_1_joint"
    };
    // ... 计算 lock_ids = 全集 - joints_to_keep ...
    Model reduced_model2 = buildReducedModel(model, lock_ids, q_rand);
}
```

---

### 6.12 `interpolation-SE3.cpp` — SE(3) 插值

```cpp
#include "pinocchio/multibody/liegroup/liegroup.hpp"
using namespace pinocchio;

int main()
{
    // SE(3) 插值：在两个位姿之间做测地线（最短路径）插值
    // 类比：球面插值 SLERP 之于 SO(3)，这是 SE(3) 上的扩展
    typedef SpecialEuclideanOperationTpl<3, double> SE3Operation;
    SE3Operation aSE3;

    SE3Operation::ConfigVector_t pose_s, pose_g;
    pose_s << 1.0, 1.0, 1.0, -0.13795, 0.13795, 0.69352, 0.69352;
    pose_g << 4.0, 3.0, 3.0, -0.46194, 0.331414, 0.800103, 0.191342;
    aSE3.normalize(pose_s);
    aSE3.normalize(pose_g);

    // 在 SE(3) 流形上插值，t=0.5 即取中间位姿
    // 数学：pole_u = pose_s * exp(0.5 * log(pose_s^{-1} * pose_g))
    SE3Operation::ConfigVector_t pole_u;
    aSE3.interpolate(pose_s, pose_g, 0.5, pole_u);
    std::cout << "Interpolated: " << pole_u.transpose() << std::endl;
}
```

**人形应用**：末端轨迹在 SE(3) 上的插值（手部抓取路径平滑）、
步态规划中足底位姿序列的插值。

---

### 6.13 `geometry-models.cpp` — 碰撞/视觉几何模型

```cpp
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/geometry.hpp"   // updateGeometryPlacements()

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // 从 URDF 加载碰撞几何（用于碰撞检测，hp-fcl 格式）
    GeometryModel collision_model;
    pinocchio::urdf::buildGeom(model, urdf_filename, COLLISION, collision_model, mesh_dir);

    // 从 URDF 加载视觉几何（用于可视化渲染，mesh 格式）
    GeometryModel visual_model;
    pinocchio::urdf::buildGeom(model, urdf_filename, VISUAL, visual_model, mesh_dir);

    Data data(model);
    GeometryData collision_data(collision_model);
    GeometryData visual_data(visual_model);

    Eigen::VectorXd q = randomConfiguration(model);

    // 正向运动学更新关节位姿
    forwardKinematics(model, data, q);

    // 根据关节位姿更新几何对象在世界系下的位置（data.oMg[i]）
    updateGeometryPlacements(model, data, collision_model, collision_data);
    updateGeometryPlacements(model, data, visual_model,    visual_data);

    // data.oMg[geom_id] 即几何对象的世界系 SE(3) 变换
}
```

---

### 6.14 `collisions.cpp` — 碰撞检测（Talos + SRDF）

```cpp
#include "pinocchio/collision/collision.hpp"   // computeCollisions()
#include "pinocchio/parsers/srdf.hpp"          // removeCollisionPairs()

int main()
{
    using namespace pinocchio;

    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    GeometryModel geom_model;
    pinocchio::urdf::buildGeom(model, urdf_filename, COLLISION, geom_model, model_path);

    // 添加所有可能的碰撞对（O(n^2) 对）
    geom_model.addAllCollisionPairs();

    // 从 SRDF 移除已知的"允许自碰撞"对（如相邻连杆）
    // 减少不必要的碰撞检测计算
    pinocchio::srdf::removeCollisionPairs(model, geom_model, srdf_filename);

    GeometryData geom_data(geom_model);

    // 从 SRDF 加载 half_sitting 参考配置
    pinocchio::srdf::loadReferenceConfigurations(model, srdf_filename);
    const auto& q = model.referenceConfigurations["half_sitting"];

    // 全量碰撞检测（检测所有配对）
    computeCollisions(model, data, geom_model, geom_data, q);

    for (size_t k = 0; k < geom_model.collisionPairs.size(); ++k)
    {
        const auto& cp = geom_model.collisionPairs[k];
        const auto& cr = geom_data.collisionResults[k];
        std::cout << "pair " << cp.first << "," << cp.second
                  << ": " << (cr.isCollision() ? "COLLISION" : "ok") << std::endl;
    }

    // 早停版本：发现第一个碰撞即停止（更快）
    computeCollisions(model, data, geom_model, geom_data, q, /*stopAtFirst=*/true);

    // 单对检测（只检测第 3 对）
    const PairIndex pair_id = 2;
    Eigen::VectorXd q_neutral = neutral(model);
    updateGeometryPlacements(model, data, geom_model, geom_data, q_neutral);
    computeCollision(geom_model, geom_data, pair_id);
}
```

---

### 6.15 `multiprecision.cpp` — 多精度浮点运算

```cpp
#include "pinocchio/math/multiprecision.hpp"
#include <boost/multiprecision/cpp_dec_float.hpp>

int main(int argc, char** argv)
{
    using namespace pinocchio;

    // 标准 double 精度模型
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    // 定义 100 位十进制精度的浮点类型
    typedef boost::multiprecision::number<
        boost::multiprecision::cpp_dec_float<100>,
        boost::multiprecision::et_off> float_100;

    // 将 double 模型转换为 100 位精度模型
    // Pinocchio 的模板设计使得标量类型可以任意替换
    typedef ModelTpl<float_100> ModelMulti;
    typedef DataTpl<float_100>  DataMulti;
    ModelMulti model_multi = model.cast<float_100>();
    DataMulti  data_multi(model_multi);

    // 用 double 配置生成 100 位精度配置
    ModelMulti::ConfigVectorType q_multi  = randomConfiguration(model_multi);
    ModelMulti::TangentVectorType v_multi = ModelMulti::TangentVectorType::Random(model.nv);
    ModelMulti::TangentVectorType a_multi = ModelMulti::TangentVectorType::Random(model.nv);

    // 同时用两种精度计算 RNEA，对比结果
    rnea(model,       data,       q_multi.cast<double>(), ...);
    rnea(model_multi, data_multi, q_multi,                ...);

    // 100 位精度结果，可用于数值分析、梯度验证（有限差分基准）
    std::cout << std::setprecision(std::numeric_limits<float_100>::max_digits10)
              << data_multi.tau << std::endl;
}
```

**应用场景**：验证解析梯度的精度（与高精度数值梯度对比）、数值分析研究。

---

## 7. 人形机器人开发路径

按以下顺序逐步深入：

```
阶段 1：SE(3) 与李群基础
  → overview-SE3.cpp / overview-lie.cpp
  → 理解 difference / integrate / log / exp

阶段 2：浮动基建模 + FK
  → overview-urdf.py（加 JointModelFreeFlyer）
  → simulation-contact-dynamics.py 前半段（加载 Talos）

阶段 3：逆向动力学与梯度
  → inverse-dynamics.py / .cpp
  → inverse-dynamics-derivatives.py / .cpp（理解 ∂τ/∂q）

阶段 4：质心动量学
  → forward-dynamics-derivatives.py（buildSampleModelHumanoidRandom）
  → 重点：data.hg, data.Ag, data.com

阶段 5：接触约束动力学（核心）
  → simulation-contact-dynamics.py（完整 Talos 仿真）
  → contact-cholesky.py（Delassus 矩阵理解）
  → talos-simulation.py（KKT 约束优化）

阶段 6：轨迹优化接口
  → forward-dynamics-derivatives.py（ABA 梯度 → DDP 的 A, B 矩阵）
  → 接入 Crocoddyl：将 ABA Derivatives 作为动力学线性化来源

阶段 7：高性能应用
  → run-algo-in-parallel.py（并行批量仿真）
  → examples/codegen/（代码生成，用于嵌入式实时控制）
```

---

## 8. 与工业机器人的对比

| 维度 | 工业机器人 | 人形机器人（Pinocchio） |
|------|-----------|------------------------|
| 基座约束 | 固定（地面） | 浮动基（空间自由） |
| 配置空间 | $q \in \mathbb{R}^n$ | $q \in SE(3) \times \mathbb{R}^n$（流形） |
| 状态更新 | $q \mathrel{+}= \dot{q}\Delta t$ | `pin.integrate(model, q, v*dt)` |
| 驱动性 | 全驱动 | **欠驱动**（浮动基 6 DOF 无输入） |
| 接触 | 无/固定末端 | 多点动态接触（KKT 约束） |
| 建模格式 | URDF（固定基） | URDF + `JointModelFreeFlyer()` |
| 控制目标 | 末端轨迹跟踪 | 质心轨迹 + 接触力 + 任务层次 |
| 关键算法 | FK / IK / RNEA | 质心动量 + 接触约束动力学 + 解析梯度 |
| 主要应用框架 | ROS MoveIt | Crocoddyl, HPP, Stack-of-Tasks |

---

## 9. 关键头文件索引

| 功能 | 头文件 |
|------|--------|
| 正向运动学 | `include/pinocchio/algorithm/kinematics.hpp` |
| 运动学梯度 | `include/pinocchio/algorithm/kinematics-derivatives.hpp` |
| Jacobian 计算 | `include/pinocchio/algorithm/jacobian.hpp` |
| Frame 运动学 | `include/pinocchio/algorithm/frames.hpp` |
| Frame 梯度 | `include/pinocchio/algorithm/frames-derivatives.hpp` |
| RNEA 逆动力学 | `include/pinocchio/algorithm/rnea.hpp` |
| RNEA 梯度 | `include/pinocchio/algorithm/rnea-derivatives.hpp` |
| ABA 正向动力学 | `include/pinocchio/algorithm/aba.hpp` |
| ABA 梯度 | `include/pinocchio/algorithm/aba-derivatives.hpp` |
| 质量矩阵 CRBA | `include/pinocchio/algorithm/crba.hpp` |
| 质心动量 | `include/pinocchio/algorithm/centroidal.hpp` |
| 质心动量梯度 | `include/pinocchio/algorithm/centroidal-derivatives.hpp` |
| 质心位置 | `include/pinocchio/algorithm/center-of-mass.hpp` |
| 接触约束动力学 | `include/pinocchio/algorithm/constrained-dynamics.hpp` |
| 接触 Cholesky | `include/pinocchio/algorithm/contact-cholesky.hpp` |
| 接触信息 | `include/pinocchio/algorithm/contact-info.hpp` |
| SE(3) 类型 | `include/pinocchio/spatial/se3.hpp` |
| 空间运动（Twist） | `include/pinocchio/spatial/motion.hpp` |
| 空间力（Wrench） | `include/pinocchio/spatial/force.hpp` |
| 空间惯量 | `include/pinocchio/spatial/inertia.hpp` |
| 指数/对数映射 | `include/pinocchio/spatial/explog.hpp` |
| 李群操作 | `include/pinocchio/multibody/liegroup/liegroup.hpp` |
| 降阶模型 | `include/pinocchio/algorithm/model.hpp` |
| URDF 解析 | `include/pinocchio/parsers/urdf.hpp` |
| SRDF 解析 | `include/pinocchio/parsers/srdf.hpp` |
| MJCF 解析 | `include/pinocchio/parsers/mjcf.hpp` |
| 碰撞检测 | `include/pinocchio/collision/collision.hpp` |
| 并行算法 | `include/pinocchio/algorithm/parallel/` |
| 多精度支持 | `include/pinocchio/math/multiprecision.hpp` |

> **重要**：上表的 `.hpp` 只含**声明与文档注释**，真正的算法实现体在同名的 `.hxx` 里，
> 位于 `include/pinocchio/src/algorithm/`（例如 `kinematics.hpp` 的实现在
> `include/pinocchio/src/algorithm/kinematics.hxx`）。阅读机制见下一节 [11.1](#111-三层文件布局-hpp--hxx--cpp)。

---

## 10. 源码目录地图

> 目标：拿到仓库后，1 分钟内知道"某类东西该去哪个目录找"。

```
pinocchio/
├── include/pinocchio/            # 头文件（库的主体，几乎全是模板）
│   ├── fwd.hpp                   # 全局前置声明、宏、Options 枚举
│   ├── context.hpp               # 决定默认标量类型（double）的开关，见 11.2
│   ├── macros.hpp                # PINOCCHIO_* 宏（断言、DLL 导出、EIGEN 检查）
│   │
│   ├── spatial/                  # ★ 空间代数（数学地基，最先读）
│   │   ├── se3.hpp               #   SE3Tpl：刚体位姿 R,p + act/actInv
│   │   ├── motion.hpp            #   MotionTpl：空间速度 Twist (ω,v)
│   │   ├── force.hpp             #   ForceTpl：空间力 Wrench (τ,f)
│   │   ├── inertia.hpp           #   InertiaTpl：6×6 空间惯量
│   │   └── explog.hpp            #   log3/exp3/log6/exp6/Jlog6（李群映射）
│   │
│   ├── multibody/                # ★ 建模层：关节 / 模型 / 数据 / 访问者
│   │   ├── joint/                #   各关节类型（RX/FreeFlyer/Spherical…）
│   │   ├── joint.hpp             #   JointModel/JointData 的 variant 封装
│   │   ├── joint-motion-subspace.hpp  # 关节子空间 S（RNEA 中的 sᵢ）
│   │   ├── visitor.hpp           #   JointUnaryVisitorBase：算法遍历的引擎
│   │   ├── liegroup/             #   李群操作（integrate/difference/interpolate）
│   │   ├── sample-models.hpp     #   buildSampleModelManipulator / Humanoid
│   │   └── pool.hpp              #   ModelPool：多线程 Data 副本
│   │
│   ├── algorithm/                # ★ 算法层：声明 + 文档（实现在 src/*.hxx）
│   │   ├── kinematics.hpp        #   FK
│   │   ├── rnea.hpp / aba.hpp    #   逆/正动力学
│   │   ├── crba.hpp              #   质量矩阵
│   │   ├── jacobian.hpp          #   Jacobian
│   │   ├── centroidal.hpp        #   质心动量
│   │   ├── constrained-dynamics.hpp / contact-cholesky.hpp  # 接触
│   │   ├── model.hpp             #   ModelTpl 类本体、buildReducedModel
│   │   └── *-derivatives.hpp     #   各算法的解析梯度
│   │
│   ├── src/                      # ★★ 模板"实现体"（.hxx），与上面 .hpp 一一对应
│   │   ├── algorithm/*.hxx       #   算法真正的循环/递推在这里
│   │   ├── multibody/*.hxx       #   ModelTpl/DataTpl 的成员实现
│   │   ├── multibody/joint/*.hxx #   每种关节的 calc()/calc_aba()
│   │   └── context/              #   默认标量 (default.hxx) 等编译上下文
│   │
│   ├── parsers/                  # URDF / SRDF / MJCF / sdf 解析
│   ├── collision/                # 碰撞检测（封装 Coal / hpp-fcl）
│   ├── autodiff/                 # CppAD / CasADi 自动微分标量支持
│   ├── codegen/                  # 代码生成（编译期特化为 C 代码）
│   └── bindings/python/          # Python 绑定的头文件（Boost.Python 包装）
│
├── src/                          # ★★ .cpp：对默认标量做"显式模板实例化"
│   └── algorithm/*.cpp           #   编进 libpinocchio，见 11.1
│
├── bindings/python/              # Python 绑定的 .cpp（expose* 函数）
├── examples/                     # 本指南第 5、6 章注释的全部 demo
├── unittest/                     # 单元测试（读实现时的"用法真值"）
├── doc/                          # 官方 doxygen 文档源（a-features/ 值得读）
└── models/                       # 示例 URDF（UR5、Talos…）
```

**记忆口诀**：想看**怎么用**→`examples/`、`unittest/`；想看**接口和公式注释**→`algorithm/*.hpp`；
想看**实现细节**→`src/algorithm/*.hxx`；想看**数据结构**→`spatial/`、`src/multibody/`。

---

## 11. 源码架构三大机制

Pinocchio 是一个"模板元编程 + 零成本抽象"的库。不先理解下面三个机制，直接看实现会寸步难行。

### 11.1 三层文件布局 .hpp / .hxx / .cpp

同一个算法被拆到**三个文件**，各司其职：

| 文件 | 位置 | 内容 | 你什么时候看它 |
|------|------|------|----------------|
| `xxx.hpp` | `include/pinocchio/algorithm/` | **函数声明 + Doxygen 注释 + 公式** | 想知道"有哪些接口、参数含义、数学定义" |
| `xxx.hxx` | `include/pinocchio/src/algorithm/` | **模板函数的实现体**（真正的循环/递推） | 想知道"算法到底怎么算的" |
| `xxx.cpp` | `src/algorithm/` | **对默认标量 `double` 的显式实例化** | 几乎不用看，只是让编译加速、生成 `.so` |

三者的连接方式（以运动学为例）：

```cpp
// kinematics.hpp 末尾：把实现体拉进来
#include "pinocchio/src/algorithm/kinematics.hxx"

// kinematics.hxx 顶部：标注"这是私有实现，别直接 include 我"
// IWYU pragma: private, include "pinocchio/algorithm/kinematics.hpp"

// kinematics.cpp：对 double 版本做一次显式实例化，编进库
template ... void forwardKinematics<context::Scalar, ...>(const Model&, Data&, ...);
```

> **为什么这么拆？** 模板函数的定义必须对调用者可见（否则链接失败），所以实现放在随头文件
> 分发的 `.hxx`。但如果每个使用 Pinocchio 的项目都从头实例化 `double` 版本，编译会极慢；
> 于是库在 `.cpp` 里**预先实例化**好默认标量版本编进 `libpinocchio.so`，你的项目直接链接即可。
> **结论**：读实现永远去 `.hxx`，不是 `.cpp`。

### 11.2 Tpl 模板 + context 默认标量

库里所有核心类型的"真名"都带 `Tpl` 后缀并按 `<Scalar, Options, JointCollection>` 模板化：

```cpp
ModelTpl<Scalar, Options, JointCollectionTpl>      // 模型
DataTpl <Scalar, Options, JointCollectionTpl>      // 数据
SE3Tpl  <Scalar, Options>                          // 位姿
MotionTpl<Scalar, Options>                         // 速度
```

你平时写的 `pinocchio::Model` / `pinocchio::SE3` 只是**默认标量的 typedef 别名**。这个"默认标量是谁"
由 `include/pinocchio/context.hpp` → `src/context/default.hxx` 决定：

```cpp
#define PINOCCHIO_SCALAR_TYPE_DEFAULT double   // ← 默认就是 double
// 于是： pinocchio::Model == ModelTpl<double, 0, JointCollectionDefaultTpl>
```

**这解释了两件让新手困惑的事**：

1. 为什么函数签名总是长长一串 `template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>`
   —— 因为同一份代码要能用 `double` / `float` / `CppAD::AD<double>`（自动微分）/ 100 位高精度
   （见 [6.15 multiprecision](#615-multiprecisioncpp--多精度浮点运算)）跑，标量类型是模板参数。
2. `model.cast<NewScalar>()` 能凭空把 `double` 模型变成自动微分模型 —— 因为整棵类型树都是模板，
   换个 `Scalar` 重新实例化即可。这正是 Pinocchio 能提供**解析梯度**和**autodiff** 的架构根源。

> 读源码技巧：看签名时可以在脑子里把 `Scalar` 替换成 `double`、把 `ModelTpl<...>` 读成 `Model`，
> signature 立刻变清爽。

### 11.3 CRTP + Boost.Fusion 访问者

这是整个库最核心、也最难的机制。**问题**：`forwardKinematics` 要遍历运动树，对每个关节调用
"该关节类型特有的" `calc()`。但关节类型五花八门（旋转、移动、球形、浮动基…），运行时又存在一个
`std::vector<JointModel>` 里。怎样既保留每种关节的静态类型（零虚函数开销），又能统一遍历？

**答案**：`JointModel`/`JointData` 用 Boost.Variant 存放所有可能的关节类型，遍历时用
**Fusion 访问者（visitor）**在编译期为每种关节生成专门代码。一个算法 = 一个或多个"访问者结构体"。

以 FK 为例，`src/algorithm/kinematics.hxx` 里的实现骨架：

```cpp
namespace impl {
  // 一个"访问者"就是算法的一趟(pass)。继承 JointUnaryVisitorBase<自己>（CRTP）
  template<typename Scalar, int Options, ...>
  struct ForwardKinematicZeroStep
  : fusion::JointUnaryVisitorBase<ForwardKinematicZeroStep<...>>
  {
    // ArgsType：这趟遍历需要透传给每个关节的参数包
    typedef boost::fusion::vector<const Model&, Data&, const ConfigVectorType&> ArgsType;

    // algo()：对【每一种】具体 JointModel 类型都会实例化一份
    //   jmodel/jdata 已经是"具体关节类型"（如 JointModelRX），静态分派、无虚函数
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model, Data & data, const ConfigVectorType & q)
    {
      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      jmodel.calc(jdata.derived(), q.derived());   // ← 关节自身的运动学(每种关节不同)
      data.liMi[i] = model.jointPlacements[i] * jdata.M();   // 关节相对父的变换
      data.oMi[i]  = (parent > 0) ? data.oMi[parent] * data.liMi[i]   // 递推到世界系
                                  : data.liMi[i];
    }
  };
}

// 对外的模板函数：只是"对每个关节跑一遍这个访问者"
template<...> void forwardKinematics(const Model & model, Data & data, const ConfigVector & q)
{
  for (JointIndex i = 1; i < model.njoints; ++i)
    ForwardKinematicZeroStep<...>::run(model.joints[i], data.joints[i],
                                       typename ...::ArgsType(model, data, q));
}
```

**读懂访问者的三把钥匙**：

| 概念 | 含义 | 在哪看 |
|------|------|--------|
| `struct XxxStep : JointUnaryVisitorBase<XxxStep>` | 算法的一趟遍历（CRTP 自引用） | `visitor.hpp` |
| `static void algo(jmodel, jdata, ...)` | 对**每种关节类型**编译期特化的实体循环 | 各 `*.hxx` |
| `jmodel.calc(jdata, q)` / `calc_aba(...)` | 关节**自身**的运动学/动力学（S、变换、惯量） | `src/multibody/joint/joint-*.hxx` |

**因此，读任何一个算法实现的套路是固定的**：

1. 在 `xxx.hxx` 里找到 `namespace impl`，数一数有几个 `struct ...Step`（= 算法有几趟）；
2. 每个 `Step` 的 `algo()` 就是那一趟对单个关节做的事，对照第 4 章的公式看；
3. 遇到 `jmodel.calc()` / `jdata.S()` / `jdata.M()` 想深究，再去 `joint/joint-revolute.hxx` 等看具体关节。

RNEA（`rnea.hxx`）就有两个 Step（`Pass1` 正向、`Pass2` 反向），和 [4.3](#43-逆向动力学-rnea) 的两趟递推
严格对应；ABA（`aba.hxx`）有三个 Step，对应 [4.4](#44-正向动力学-aba) 的三趟。**公式 ↔ Step 一一对应**，
这是本指南第 4 章能直接当"实现导读"用的原因。

---

## 12. 实战：追踪一个算法从 API 到实现

以最常用的 `forwardKinematics` 为例，走通"Python 调用 → C++ 实现 → 关节细节"的完整链路。
掌握这一条，其余算法照抄即可。

**① 入口（Python）**——`examples/overview-urdf.py`：
```python
pinocchio.forwardKinematics(model, data, q)
```

**② Python→C++ 绑定**——`bindings/python/algorithm/expose-kinematics.cpp`
里用 `bp::def("forwardKinematics", ...)` 把 C++ 模板函数的 `double` 实例暴露给 Python。
（绑定机制详见第 13 节。）

**③ 声明 + 文档**——[include/pinocchio/algorithm/kinematics.hpp](include/pinocchio/algorithm/kinematics.hpp)：
读函数注释确认语义（"更新 `data.oMi[]`，$O(n)$"），拿到签名。文件末尾 `#include ".../kinematics.hxx"`。

**④ 实现体**——[include/pinocchio/src/algorithm/kinematics.hxx](include/pinocchio/src/algorithm/kinematics.hxx)：
- 对外函数体：一个 `for` 循环，对每个关节 `run` 访问者 `ForwardKinematicZeroStep`；
- 访问者 `algo()`：`jmodel.calc()` + `oMi[i] = oMi[parent] * liMi[i]`（就是 [4.1](#41-正向运动学-fk) 的递推式
  ${}^0T_i = {}^0T_{p(i)} \cdot {}^{p(i)}T_i$）。

**⑤ 关节细节**——想知道旋转关节的 `calc()` 怎么把 $q_i$ 变成变换矩阵：
[include/pinocchio/src/multibody/joint/joint-revolute.hxx](include/pinocchio/src/multibody/joint/joint-revolute.hxx)
的 `JointDataRevolute::calc()`——用 $\cos q_i,\sin q_i$ 填充 `M`（该关节的 `jdata.M()`）。

**⑥ 数据结构**——`data.oMi` 的定义在
[include/pinocchio/src/multibody/data.hxx](include/pinocchio/src/multibody/data.hxx)（`DataTpl` 成员），
`SE3` 的 `operator*`/`actInv` 在 [include/pinocchio/spatial/se3.hpp](include/pinocchio/spatial/se3.hpp)。

**通用追踪清单**（换任何算法都适用）：
```
公开函数名  →  algorithm/<name>.hpp        （看签名 + 注释 + 公式）
            →  src/algorithm/<name>.hxx    （看 impl::*Step::algo，对照第4章公式）
            →  src/multibody/joint/*.hxx   （深究单关节 calc/calc_aba）
            →  spatial/*.hpp               （深究 SE3/Motion/Force/Inertia 运算）
旁证用法    →  unittest/<name>.cpp、examples/
```

---

## 13. Python 绑定如何映射到 C++

Pinocchio 的 Python 层是 **Boost.Python** 手写绑定（不是 pybind11、不是自动生成），
所以每个 Python 符号都能在 C++ 里找到对应物。理解映射规则后，可以从 Python API 反查 C++ 实现。

**目录对应**：
```
bindings/python/<模块>/expose-<主题>.cpp      # 定义 expose 函数，调用 bp::def / bp::class_
include/pinocchio/bindings/python/<模块>/     # 绑定用的 C++ 包装头（Visitor 模式）
bindings/python/module.cpp                    # 总入口，依次调用所有 exposeXxx()
```

**映射规律**：

| Python | C++ 对应 | 说明 |
|--------|----------|------|
| `pin.forwardKinematics(...)` | `bp::def("forwardKinematics", ...)` in `expose-kinematics.cpp` | 函数：暴露的是 `double` 实例 |
| `pin.SE3` | `bp::class_<SE3>(...)` in `expose-SE3.cpp` | 类：`ModelTpl<double>` 等的别名被 `class_` 包装 |
| `model.createData()` | `ModelTpl::createData()` 成员 | 成员方法 1:1 暴露 |
| `data.oMi` | `DataTpl::oMi` 成员 | 通过 `.def_readwrite`/`add_property` 暴露 |
| `pin.ReferenceFrame.LOCAL` | `enum ReferenceFrame` | 枚举用 `bp::enum_` 暴露 |

> 找一个 Python 函数的 C++ 源头：`grep -rn "\"forwardKinematics\"" bindings/python/`，
> 命中的 `expose-*.cpp` 就是绑定点，顺着它调用的模板函数名即可跳进 `algorithm/*.hpp`。

**要点**：Python 侧看到的一切都是 C++ **默认标量（double）** 的实例；autodiff / codegen 等
非 `double` 标量的能力只在 C++ 层可用，不通过 Python 暴露。

---

## 14. 源码阅读推荐顺序

按"地基→建模→算法→扩展"的依赖顺序读，避免一上来陷进模板细节。

```
第 0 步：先读机制，别读代码
  → 本指南 §11（三大机制）+ §12（一条追踪链路）
  → 目的：建立"公式↔Step、.hpp↔.hxx、Tpl↔别名"的心智模型

第 1 步：空间代数地基（spatial/）
  → se3.hpp → motion.hpp → force.hpp → inertia.hpp → explog.hpp
  → 对照本指南 §2；这是后面一切算法的"数据类型词汇表"

第 2 步：建模层（multibody/）
  → model.hpp(ModelTpl) + src/multibody/data.hxx(DataTpl)：先看有哪些字段（对照 §3）
  → joint/joint-revolute.hxx：读懂一种最简单关节的 calc / calc_aba / S
  → visitor.hpp：读懂 JointUnaryVisitorBase::run 如何分派（§11.3）

第 3 步：运动学算法（最短、最适合练手）
  → src/algorithm/kinematics.hxx：一趟遍历，验证 §12 的追踪链
  → jacobian.hxx：几何 Jacobian 如何按列填充

第 4 步：动力学三大件（对照 §4 的公式读 Step）
  → rnea.hxx（2 趟）→ crba.hxx（复合惯量、上三角）→ aba.hxx（3 趟 + 关节化体惯量）

第 5 步：人形核心
  → centroidal.hxx（质心动量矩阵 Ag）
  → constrained-dynamics.hxx + contact-cholesky.hxx（KKT / Delassus，最硬）

第 6 步：解析梯度（读完对应正算法再读其导数）
  → rnea-derivatives.hxx → aba-derivatives.hxx → kinematics-derivatives.hxx

第 7 步：按需扩展
  → parsers/urdf（建模来源）、collision/（碰撞）、autodiff/ 与 codegen/（高级标量）
  → 想加新关节类型：doc/e-development/Implement_a_new_joint.md 是官方教程
```

**贯穿始终的两个好习惯**：
1. **对照单测读实现**——`unittest/<algo>.cpp` 给出该算法的期望输入输出，是最好的"实现真值"。
2. **对照本指南第 4 章读 `.hxx`**——每个 `impl::*Step` 都能在第 4 章找到对应公式，公式讲"为什么"，
   代码讲"怎么写"，两边对读效率最高。
