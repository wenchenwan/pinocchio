import math
import sys
import time
from pathlib import Path

import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# Talos 人形机器人双足接触动力学仿真
# 演示：浮动基建模 + 刚性足底接触 + 姿态稳定控制 + KKT 求解
# 可视化：MeshCat（浏览器中打开脚本输出的 URL）
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir   = pinocchio_model_dir
urdf_filename = "talos_reduced.urdf"
urdf_model_path = model_path / "talos_data/robots" / urdf_filename
srdf_filename = "talos.srdf"
srdf_full_path = model_path / "talos_data/srdf" / srdf_filename

# 关键：传入 JointModelFreeFlyer() 使机器人基座在空间中自由运动（6 DOF 浮动基）
# 不加此参数则基座被固定在世界原点，无法模拟真实双足机器人
# 浮动基使 nq += 7（3 平移 + 4 四元数），nv += 6（3 平移速度 + 3 角速度）
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir, pin.JointModelFreeFlyer()
)

# 启动 MeshCat 可视化服务器（在浏览器中查看）
try:
    viz = MeshcatVisualizer(model, collision_model, visual_model)
    viz.initViewer(open=False)
except ImportError as err:
    print(
        "Error while initializing the viewer. "
        "It seems you should install Python meshcat"
    )
    print(err)
    sys.exit(0)

viz.loadViewerModel()

# 从 SRDF 文件加载参考配置
# "half_sitting"：Talos 的标准站立姿势（双腿微蹲，双臂自然下垂）
pin.loadReferenceConfigurations(model, srdf_full_path)
q0 = model.referenceConfigurations["half_sitting"]

# q_ref：在参考配置附近稍加扰动，作为姿态控制的目标
q_ref = pin.integrate(model, q0, 0.1 * np.random.rand(model.nv))
viz.display(q0)

# 足底 Frame 名称（Talos 的左右脚底板）
feet_name = ["left_sole_link", "right_sole_link"]
frame_ids = [model.getFrameId(frame_name) for frame_name in feet_name]

v0 = np.zeros(model.nv)
v_ref = v0.copy()

# 两份 Data：仿真用（data_sim）和控制计算用（data_control）
# 分开是为了避免控制律计算覆盖仿真中间结果
data_sim     = model.createData()
data_control = model.createData()

# ---- 构建足底接触约束模型 ----
contact_models = []
contact_datas  = []

for frame_id in frame_ids:
    frame = model.frames[frame_id]
    # CONTACT_6D：6D 刚性接触（同时约束位置和姿态，适合平足）
    # joint1_id：接触关节（足底关节）
    # joint1_placement：接触点相对该关节的 SE(3) 变换（通常为 Frame 的 placement）
    contact_model = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_6D, model, frame.parentJoint, frame.placement
    )
    contact_models.append(contact_model)
    contact_datas.append(contact_model.createData())

num_constraints = len(frame_ids)
contact_dim = 6 * num_constraints   # 每个 6D 接触贡献 6 个约束方程

# 预初始化接触约束动力学所需的内存结构
# 必须在第一次调用 constraintDynamics 之前执行
pin.initConstraintDynamics(model, data_sim, contact_models)

# ---- 仿真参数 ----
t  = 0
dt = 5e-3   # 仿真步长 5 ms（200 Hz），人形实时控制的典型频率

# 选择矩阵 S：将全身 nv 维度投影到驱动关节空间（去除浮动基的 6 个不可驱动 DOF）
# S ∈ R^{(nv-6)×nv}，S.T[6:, :] = I 表示只驱动索引 6 以后的关节
S = np.zeros((model.nv - 6, model.nv))
S.T[6:, :] = np.eye(model.nv - 6)

# PD 控制增益（关节空间姿态稳定）
Kp_posture = 30.0
# 临界阻尼系数：Kv = 2 * sqrt(Kp)（使响应不振荡不过阻尼）
Kv_posture = 0.05 * math.sqrt(Kp_posture)

q   = q0.copy()
v   = v0.copy()
tau = np.zeros(model.nv)

T = 5  # 仿真总时长（秒）

while t <= T:
    print("t:", t)
    t += dt

    tic = time.time()

    # ---- 控制律：重力补偿 + 姿态稳定 PD ----
    # Step 1：计算所有关节的几何 Jacobian（存入 data_control.J[]）
    pin.computeJointJacobians(model, data_control, q)

    # Step 2：构造接触约束 Jacobian J_constraint ∈ R^{contact_dim × nv}
    J_constraint = np.zeros((contact_dim, model.nv))
    constraint_index = 0
    for k in range(num_constraints):
        contact_model = contact_models[k]
        # getFrameJacobian：从已计算的 Jacobian 中提取指定 Frame 的 6×nv 块
        J_constraint[constraint_index : constraint_index + 6, :] = pin.getFrameJacobian(
            model,
            data_control,
            contact_model.joint1_id,
            contact_model.joint1_placement,
            contact_model.reference_frame,
        )
        constraint_index += 6

    # Step 3：构造包含驱动约束和接触约束的系数矩阵
    # A = [S; J_constraint]，b = RNEA(q,v,0)（零加速度下的动力学力矩）
    A = np.vstack((S, J_constraint))
    b = pin.rnea(model, data_control, q, v, np.zeros(model.nv))

    # Step 4：最小二乘求解 A^T * [tau_driven; lambda] = b
    # sol[:nv-6] 为驱动力矩，sol[nv-6:] 为接触力（Lagrange 乘子）
    sol = np.linalg.lstsq(A.T, b, rcond=None)[0]
    tau = np.concatenate((np.zeros((6)), sol[: model.nv - 6]))

    # Step 5：叠加关节空间 PD 稳定项（抑制漂移）
    # difference(q_ref, q) 在流形上计算配置误差（不是简单相减）
    tau[6:] += (
        -Kp_posture * (pin.difference(model, q_ref, q))[6:]
        - Kv_posture * (v - v_ref)[6:]
    )

    # ---- 带接触约束的正向动力学 ----
    # 求解 KKT 鞍点系统：
    #   [M  J^T] [q̈ ] = [τ - C(q,v)v - g(q)]
    #   [J  0  ] [λ ]   [-J̇v               ]
    # 结果：data_sim.ddq（关节加速度），data_sim.lambda_c（接触力）
    # ProximalSettings：近端正则化参数，防止约束 Jacobian 奇异时数值发散
    prox_settings = pin.ProximalSettings(1e-12, 1e-12, 10)
    a = pin.constraintDynamics(
        model, data_sim, q, v, tau, contact_models, contact_datas, prox_settings
    )
    print("a:", a.T)
    print("v:", v.T)
    print("constraint:", np.linalg.norm(J_constraint @ a))  # 约束满足残差
    print("iter:", prox_settings.iter)

    # ---- 半隐式欧拉积分 ----
    # 先用新加速度更新速度，再用新速度更新配置（比显式欧拉更稳定）
    v += a * dt
    # integrate：流形积分，保证四元数归一化等约束
    q = pin.integrate(model, q, v * dt)

    viz.display(q)
    elapsed_time = time.time() - tic

    # 等待至下一个仿真步，保持实时感
    time.sleep(max(0, dt - elapsed_time))
