import sys
import time

import hppfcl as fcl
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# 闭环运动链（Closed Kinematic Chain）仿真
# 示例：四连杆机构（Four-Bar Linkage）
#
# 拓扑结构（开链展开 + 闭合约束）：
#   World --- [A1] --- joint1(RY) --- [B1] --- joint2(RY) --- [A2] --- joint3(RY) --- [B2]
#                                                                                       |
#                                            约束：B2末端 = World上固定点（构成闭环）
#
# 处理方法：
#   1. 将闭环结构展开为开链模型（Pinocchio 只支持树状拓扑）
#   2. 用 RigidConstraintModel(CONTACT_3D) 描述闭合约束
#   3. 先用逆几何（Inverse Geometry）找满足约束的初始配置
#   4. 再用 constraintDynamics 进行含约束的正向动力学仿真
# ============================================================

# ---- 连杆几何和惯量参数 ----
height = 0.1
width  = 0.01
radius = 0.05   # 胶囊体半径（仅用于碰撞/可视化）

# A 型连杆（较长，白色）：作为驱动臂和从动臂
mass_link_A   = 10.0
length_link_A = 1.0
shape_link_A  = fcl.Capsule(radius, length_link_A)

# B 型连杆（较短，红色）：作为连接杆
mass_link_B   = 5.0
length_link_B = 0.6
shape_link_B  = fcl.Capsule(radius, length_link_B)

# 计算各连杆的惯量（FromBox 用于动力学，形状用于碰撞检测，二者独立）
inertia_link_A = pin.Inertia.FromBox(mass_link_A, length_link_A, width, height)
# 连杆质心位于 X 轴方向上的中点
placement_center_link_A = pin.SE3.Identity()
placement_center_link_A.translation = pin.XAxis * length_link_A / 2.0
# 碰撞形状（胶囊体）的姿态：FCL 胶囊轴默认沿 Z，需旋转到沿 X
placement_shape_A = placement_center_link_A.copy()
placement_shape_A.rotation = pin.Quaternion.FromTwoVectors(
    pin.ZAxis, pin.XAxis
).matrix()

inertia_link_B = pin.Inertia.FromBox(mass_link_B, length_link_B, width, height)
placement_center_link_B = pin.SE3.Identity()
placement_center_link_B.translation = pin.XAxis * length_link_B / 2.0
placement_shape_B = placement_center_link_B.copy()
placement_shape_B.rotation = pin.Quaternion.FromTwoVectors(
    pin.ZAxis, pin.XAxis
).matrix()

# ---- 手工构建开链运动学模型 ----
# 将四连杆展开为单链：World → A1 → B1 → A2 → B2
# 闭合约束（B2末端 ↔ World上一点）稍后用 RigidConstraintModel 描述
model          = pin.Model()
collision_model = pin.GeometryModel()

RED_COLOR   = np.array([1.0, 0.0, 0.0, 1.0])
WHITE_COLOR = np.array([1.0, 1.0, 1.0, 1.0])

# A1 连杆：固定在世界坐标系（base_joint_id=0），不添加运动关节
# 只创建几何体以便可视化
base_joint_id = 0
geom_obj0 = pin.GeometryObject(
    "link_A1",
    base_joint_id,
    # 旋转使胶囊体轴与 X 轴对齐（胶囊默认轴为 Z）
    pin.SE3(pin.Quaternion.FromTwoVectors(pin.ZAxis, pin.XAxis).matrix(), np.zeros(3)),
    shape_link_A,
)
geom_obj0.meshColor = WHITE_COLOR
collision_model.addGeometryObject(geom_obj0)

# joint1：A1 末端（X=length_link_A/2 处）的绕 Y 轴旋转关节（RY）
# B1 连杆挂在 joint1 上
joint1_placement = pin.SE3.Identity()
joint1_placement.translation = pin.XAxis * length_link_A / 2.0
joint1_id = model.addJoint(
    base_joint_id, pin.JointModelRY(), joint1_placement, "link_B1"
)
model.appendBodyToJoint(joint1_id, inertia_link_B, placement_center_link_B)
geom_obj1 = pin.GeometryObject("link_B1", joint1_id, placement_shape_B, shape_link_B)
geom_obj1.meshColor = RED_COLOR
collision_model.addGeometryObject(geom_obj1)

# joint2：B1 末端（X=length_link_B 处）的绕 Y 轴旋转关节
# A2 连杆挂在 joint2 上
joint2_placement = pin.SE3.Identity()
joint2_placement.translation = pin.XAxis * length_link_B
joint2_id = model.addJoint(joint1_id, pin.JointModelRY(), joint2_placement, "link_A2")
model.appendBodyToJoint(joint2_id, inertia_link_A, placement_center_link_A)
geom_obj2 = pin.GeometryObject("link_A2", joint2_id, placement_shape_A, shape_link_A)
geom_obj2.meshColor = WHITE_COLOR
collision_model.addGeometryObject(geom_obj2)

# joint3：A2 末端的绕 Y 轴旋转关节
# B2 连杆挂在 joint3 上（B2 末端将通过约束锁定到 World）
joint3_placement = pin.SE3.Identity()
joint3_placement.translation = pin.XAxis * length_link_A
joint3_id = model.addJoint(joint2_id, pin.JointModelRY(), joint3_placement, "link_B2")
model.appendBodyToJoint(joint3_id, inertia_link_B, placement_center_link_B)
geom_obj3 = pin.GeometryObject("link_B2", joint3_id, placement_shape_B, shape_link_B)
geom_obj3.meshColor = RED_COLOR
collision_model.addGeometryObject(geom_obj3)

# 视觉模型直接复用碰撞模型（避免重复）
visual_model = collision_model

q0   = pin.neutral(model)   # 零位初始配置
data = model.createData()
pin.forwardKinematics(model, data, q0)

# ---- 定义闭环约束（CONTACT_3D：位置约束，不含姿态）----
# 约束含义：joint3（B2 根部）沿 X 轴延伸 length_link_B 的末端点
#           必须与 World 上 X = -length_link_A/2 处的点重合
# 即：B2 的末端被"钉"在世界坐标系的一个固定点上，形成闭环
constraint1_joint1_placement = pin.SE3.Identity()
constraint1_joint1_placement.translation = pin.XAxis * length_link_B

constraint1_joint2_placement = pin.SE3.Identity()
constraint1_joint2_placement.translation = -pin.XAxis * length_link_A / 2.0

# RigidConstraintModel(CONTACT_3D, model, joint1_id, joint1_placement, joint2_id, joint2_placement)
# CONTACT_3D：只约束位置（3 自由度），等效为"球铰链"接触
# joint2_id=base_joint_id=0 → 约束到世界坐标系
constraint_model = pin.RigidConstraintModel(
    pin.ContactType.CONTACT_3D,
    model,
    joint3_id,                      # 活动端：B2 末端
    constraint1_joint1_placement,
    base_joint_id,                  # 固定端：世界坐标系
    constraint1_joint2_placement,
)
constraint_data  = constraint_model.createData()
constraint_dim   = constraint_model.size()   # = 3（CONTACT_3D 的约束维度）

# ============================================================
# 阶段一：逆几何（Inverse Geometry）
# 目标：找到满足闭环约束 c(q) = 0 的初始配置 q_sol
#
# 使用增广 Lagrangian / KKT 松弛方法：
#   min_{q,λ}  ρ/2 ‖δq‖² + λᵀ c(q) + μ/2 ‖c(q)‖²
#   KKT 系统：[ρI  Jᵀ] [δq]   [Jᵀ(c+y·μ)]
#             [J   μI ] [dλ] = [-c - y·μ  ]
# ============================================================

rho = 1e-10   # 正则化系数（使 KKT 矩阵 M = ρI 正定，防止奇异）
mu  = 1e-4    # KKT 松弛参数（影响约束满足精度 vs 数值稳定性的平衡）

q   = q0.copy()
y   = np.ones(constraint_dim)   # Lagrange 乘子初始值

# data.M 在这里被借用为 KKT 质量块（正则化项），不是真实质量矩阵
data.M = np.eye(model.nv) * rho

# ContactCholeskyDecomposition：稀疏 KKT 分解器（用于求解逆几何的 KKT 系统）
kkt_constraint = pin.ContactCholeskyDecomposition(model, [constraint_model])

eps = 1e-10   # 收敛阈值
N   = 100     # 最大迭代次数

for k in range(N):
    # 计算所有关节 Jacobian（kkt_constraint.compute 需要此数据）
    pin.computeJointJacobians(model, data, q)
    # 构建 KKT 系统并分解：[ρI  Jᵀ; J  μI]
    kkt_constraint.compute(model, data, [constraint_model], [constraint_data], mu)

    # 约束残差：B2末端与 World 固定点之间的位置差（3D 向量）
    # constraint_data.c1Mc2 = joint1 接触点相对 joint2 接触点的 SE3 变换
    # .translation = 位置残差（CONTACT_3D 只关心位置）
    constraint_value = constraint_data.c1Mc2.translation

    # 接触点 Jacobian（活动端，LOCAL_WORLD_ALIGNED 参考系，只取平移部分）
    J = pin.getFrameJacobian(
        model,
        data,
        constraint_model.joint1_id,
        constraint_model.joint1_placement,
        constraint_model.reference_frame,
    )[:3, :]   # 只取前 3 行（位置 Jacobian，舍弃旋转部分）

    # 检查原始可行性（约束违反量）和对偶可行性（KKT 梯度）
    primal_feas = np.linalg.norm(constraint_value, np.inf)
    dual_feas   = np.linalg.norm(J.T.dot(constraint_value + y), np.inf)
    if primal_feas < eps and dual_feas < eps:
        print("Convergence achieved")
        break
    print("constraint_value:", np.linalg.norm(constraint_value))

    # KKT 右端向量：[−c − y·μ; 0_nv]
    # 第一块：约束残差 + Lagrange 修正项（驱动约束满足）
    # 第二块：零（无外部力矩，纯几何问题）
    rhs = np.concatenate([-constraint_value - y * mu, np.zeros(model.nv)])

    # 求解 KKT 系统：[dλ; δq] = KKT⁻¹ * rhs
    dz  = kkt_constraint.solve(rhs)
    dy  = dz[:constraint_dim]    # Lagrange 乘子增量
    dq  = dz[constraint_dim:]    # 关节角增量

    alpha = 1.0   # 步长（此处不做线搜索）
    # 在流形上积分：q ← Exp(q, -α·δq)（保证四元数归一化等约束）
    q  = pin.integrate(model, q, -alpha * dq)
    # Lagrange 乘子更新（增广 Lagrangian 对偶变量更新）
    y -= alpha * (-dy + y)

# 将角度映射到 [-π/2, π/2] 区间（消除 2π 歧义）
q_sol = (q[:] + np.pi) % np.pi - np.pi

# ============================================================
# 阶段二：含闭环约束的正向动力学仿真
# 使用 constraintDynamics（KKT 方法）：
#   [M  Jᵀ] [q̈]   [τ - C(q,v)v - g(q)]
#   [J  0 ] [λ ] = [-Kp·c(q) - Kd·Jv  ]
# 约束修正器（corrector）通过 Kp/Kd 稳定约束（Baumgarte 稳定化）
# ============================================================

q   = q_sol.copy()
v   = np.zeros(model.nv)
tau = np.zeros(model.nv)   # 无驱动力矩（自由落体 + 约束）

dt    = 5e-3     # 仿真步长 5 ms
T_sim = 10       # 仿真总时长 10 s
t     = 0

mu_sim = 1e-10   # 仿真阶段的 KKT 松弛参数（更小 → 约束更严格）

# ---- Baumgarte 稳定化参数（防止约束漂移积累）----
# 在 KKT 右端加入 PD 修正项：-Kp·c(q) - Kd·J·v
# 使约束违反量按 Kp/Kd 设定的动态特性收敛（类似 PD 控制器）
# 临界阻尼条件：Kd = 2·√Kp（无超调，快速收敛）
constraint_model.corrector.Kp[:] = 10
constraint_model.corrector.Kd[:] = 2.0 * np.sqrt(constraint_model.corrector.Kp)

# 预计算约束动力学所需的稀疏结构（内部预分配内存）
pin.initConstraintDynamics(model, data, [constraint_model])

# ProximalSettings(absolute_tol, mu, max_iter)
# 内点法（近端点法）求解 KKT 系统的精度设置
prox_settings = pin.ProximalSettings(1e-8, mu_sim, 10)

try:
    viz = MeshcatVisualizer(model, collision_model, visual_model)
    viz.initViewer(open=True)
except ImportError as error:
    print(error)
    sys.exit(0)
viz.loadViewerModel()
viz.display(q_sol)

while t <= T_sim:
    # 含闭环约束的正向动力学：求解 KKT 系统得到加速度 a
    a = pin.constraintDynamics(
        model, data, q, v, tau, [constraint_model], [constraint_data], prox_settings
    )
    # 半隐式欧拉积分：先更新速度，再在流形上积分位置
    # 半隐式比显式欧拉更稳定（速度用新加速度，位置用新速度）
    v += a * dt
    q  = pin.integrate(model, q, v * dt)
    viz.display(q)
    time.sleep(dt)
    t += dt
