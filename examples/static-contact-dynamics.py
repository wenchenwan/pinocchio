from pathlib import Path

import numpy as np
import pinocchio as pin

np.set_printoptions(linewidth=np.inf)

# ============================================================
# 静态接触力计算：给定站立配置，求各足接触力和关节力矩
# 机器人：Solo-12（四足，浮动基，12 DOF）
#
# 问题描述：
#   静止时加速度为零，运动方程退化为：
#     g(q) = τ + Jc^T * λ                              (1)
#   其中 g(q) 为重力项，τ 为关节力矩，Jc 为接触 Jacobian，λ 为接触力
#
# 求解策略：
#   将方程按浮动基（_bl, 6维）和驱动关节（_j, 12维）分块：
#     | g_bl |   |  0  |   | Jc_bl^T |   | λ |
#     | g_j  | = | τ  | + | Jc_j^T  | * | λ |        (2)
#
#   Step 1：从浮动基行求接触力（无驱动力矩，共 6 个方程求 12 个未知数 λ）
#     g_bl = Jc_bl^T * λ  →  λ = pinv(Jc_bl^T) * g_bl  (伪逆最小范数解)  (3)
#
#   Step 2：从驱动关节行求关节力矩
#     τ = g_j - Jc_j^T * λ                             (4)
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir   = pinocchio_model_dir
urdf_model_path = model_path / "solo_description/robots/solo12.urdf"

# 浮动基建模（JointModelFreeFlyer）：基座在空间自由运动
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir, pin.JointModelFreeFlyer()
)
data = model.createData()

# Solo-12 的静止站立配置
# 格式：[x, y, z, qx, qy, qz, qw（浮动基 7 维）, 12 个关节角]
q0 = np.array([
    0.0, 0.0, 0.235,          # 浮动基位置（z=0.235 m，离地高度）
    0.0, 0.0, 0.0, 1.0,       # 浮动基四元数（单位，无旋转）
    0.0,  0.8, -1.6,          # 左前腿：髋、膝、足
    0.0, -0.8,  1.6,          # 右前腿
    0.0,  0.8, -1.6,          # 左后腿
    0.0, -0.8,  1.6,          # 右后腿
])
v0 = np.zeros(model.nv)
a0 = np.zeros(model.nv)

# ---- Step 1：计算重力项 g(q) ----
# RNEA(q, v=0, a=0)：零速度零加速度时，力矩即为重力补偿项 g(q)
g_grav = pin.rnea(model, data, q0, v0, a0)
g_bl = g_grav[:6]   # 浮动基 6 维（对应平移+旋转方向的广义力）
g_j  = g_grav[6:]   # 驱动关节 12 维（每个关节的重力力矩）

# ---- Step 2：构造接触 Jacobian ----
feet_names = ["FL_FOOT", "FR_FOOT", "HL_FOOT", "HR_FOOT"]
feet_ids   = [model.getFrameId(n) for n in feet_names]
bl_id      = model.getFrameId("base_link")
ncontact   = len(feet_names)

# 计算每个足端 Frame 的 Jacobian（LOCAL 坐标系下，6×nv）
Js__feet_q = [
    np.copy(pin.computeFrameJacobian(model, data, q0, fid, pin.LOCAL))
    for fid in feet_ids
]

# 只取浮动基部分（前 6 列）的线速度行（前 3 行），得到 3×6 的 Jacobian 块
Js__feet_bl = [np.copy(J[:3, :6]) for J in Js__feet_q]

# 拼成 6×12 矩阵（转置后）：Jc_bl^T ∈ R^{6 × 3*ncontact}
Jc__feet_bl_T = np.zeros([6, 3 * ncontact])
Jc__feet_bl_T[:, :] = np.vstack(Js__feet_bl).T

# ---- Step 3：伪逆求接触力 λ（最小范数解）----
# λ = pinv(Jc_bl^T) * g_bl：在满足浮动基平衡的前提下，接触力范数最小
ls = np.linalg.pinv(Jc__feet_bl_T) @ g_bl

# 拆分为每个足端的 3D 接触力（局部坐标系）
ls__f = np.split(ls, ncontact)

# 更新 Frame 位姿，用于坐标系转换
pin.framesForwardKinematics(model, data, q0)

# 将接触力从足端局部系转换到基座（base_link）坐标系，便于验证
ls__bl = []
for l__f, foot_id in zip(ls__f, feet_ids):
    l_sp__f  = pin.Force(l__f, np.zeros(3))                        # 构造空间力（只有线力）
    l_sp__bl = data.oMf[bl_id].actInv(data.oMf[foot_id].act(l_sp__f))  # 坐标系变换
    ls__bl.append(np.copy(l_sp__bl.vector))

print("\n--- CONTACT FORCES ---")
for l__f, foot_id, name in zip(ls__bl, feet_ids, feet_names):
    print(f"Contact force at foot {name} expressed at the BL is: {l__f}")

# 验证：所有接触力之和应等于重力项（牛顿第三定律）
print(
    "Error between contact forces and gravity at base link: "
    f"{np.linalg.norm(g_bl - sum(ls__bl))}"
)

# ---- Step 4：求驱动关节力矩 ----
# 驱动关节 Jacobian 块（取前 3 行，即线速度部分，对应 3D 接触）
Js_feet_j = [np.copy(J[:3, 6:]) for J in Js__feet_q]
Jc__feet_j_T = np.zeros([12, 3 * ncontact])
Jc__feet_j_T[:, :] = np.vstack(Js_feet_j).T

# τ = g_j - Jc_j^T * λ
tau = g_j - Jc__feet_j_T @ ls

# ---- 交叉验证 1：RNEA 含外力 ----
# 将接触力作为外力传入 RNEA，验证结果是否与手动计算一致
pin.framesForwardKinematics(model, data, q0)
joint_names = ["FL_KFE", "FR_KFE", "HL_KFE", "HR_KFE"]
joint_ids   = [model.getJointId(n) for n in joint_names]

fs_ext = [pin.Force(np.zeros(6)) for _ in range(len(model.joints))]
for idx, joint in enumerate(model.joints):
    if joint.id in joint_ids:
        fext__bl = pin.Force(ls__bl[joint_ids.index(joint.id)])
        # 将接触力从基座系变换到关节局部系（RNEA 需要局部系下的外力）
        fs_ext[idx] = data.oMi[joint.id].actInv(data.oMf[bl_id].act(fext__bl))

tau_rnea = pin.rnea(model, data, q0, v0, a0, fs_ext)

print("\n--- ID: JOINT TORQUES ---")
print(f"Tau from RNEA:         {tau_rnea}")
print(f"Tau computed manually: {np.append(np.zeros(6), tau)}")
print(f"Tau error: {np.linalg.norm(np.append(np.zeros(6), tau) - tau_rnea)}")

# ---- 交叉验证 2：正向动力学 ----
# 用求得的 τ 和 λ 跑正向动力学，验证加速度是否为零
Js_feet3d_q = [np.copy(J[:3, :]) for J in Js__feet_q]
acc = pin.forwardDynamics(
    model, data, q0, v0,
    np.append(np.zeros(6), tau),   # 完整力矩向量（浮动基部分为 0）
    np.vstack(Js_feet3d_q),        # 全局接触 Jacobian（3*ncontact × nv）
    np.zeros(12),                  # 零约束漂移
)

print("\n--- FD: ACC. & CONTACT FORCES ---")
print(f"Norm of the FD acceleration: {np.linalg.norm(acc)}")   # 应接近 0
print(f"Contact forces manually: {ls}")
print(f"Contact forces FD: {data.lambda_c}")
print(f"Contact forces error: {np.linalg.norm(data.lambda_c - ls)}")
