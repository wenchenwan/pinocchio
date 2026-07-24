import numpy as np
import pinocchio as pin

# 构建内置随机人形模型（浮动基，双足）
model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

# 设置关节限位（随机人形模型默认无限位）
model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit = np.ones((model.nq, 1))

q = pin.randomConfiguration(model)     # 关节配置
v = np.random.rand(model.nv, 1)        # 关节速度
a = np.random.rand(model.nv, 1)        # 关节加速度

# 一次调用预计算所有关节速度/加速度对 (q, v, a) 的梯度
# 内部做一次完整的正向运动学递推，并保存所有中间量
# 后续的 getJointVelocityDerivatives / getJointAccelerationDerivatives 从缓存取结果
pin.computeForwardKinematicsDerivatives(model, data, q, v, a)

# 选择要查询梯度的关节（右腿末端关节）
joint_name = "rleg6_joint"
joint_id = model.getJointId(joint_name)

# ---- 空间速度对 (q, v) 的梯度 ----
# WORLD 坐标系：量在世界系中表达（绝对坐标）
# dv_dq ∈ R^{6×nq}：∂v_i/∂q（几何 Jacobian 的 q-偏导）
# dv_dv ∈ R^{6×nv}：∂v_i/∂v = 几何 Jacobian J_i（直接等于 Jacobian 本身）
(dv_dq, dv_dv) = pin.getJointVelocityDerivatives(
    model, data, joint_id, pin.ReferenceFrame.WORLD
)

# LOCAL 坐标系：量在关节自身局部系中表达
(dv_dq_local, dv_dv_local) = pin.getJointVelocityDerivatives(
    model, data, joint_id, pin.ReferenceFrame.LOCAL
)

# ---- 空间加速度对 (q, v, a) 的梯度 ----
# WORLD 坐标系：
# dv_dq ∈ R^{6×nq}：∂v_i/∂q（与上面速度导数一致）
# da_dq ∈ R^{6×nq}：∂a_i/∂q（包含科氏项和重力项贡献）
# da_dv ∈ R^{6×nv}：∂a_i/∂v（科氏加速度对速度的偏导）
# da_da ∈ R^{6×nv}：∂a_i/∂a = 几何 Jacobian J_i（加速度对关节加速度线性）
(dv_dq, da_dq, da_dv, da_da) = pin.getJointAccelerationDerivatives(
    model, data, joint_id, pin.ReferenceFrame.WORLD
)

# LOCAL 坐标系下的同一组梯度
# 注意：da_da_local 也等于局部系下的几何 Jacobian（∂v_i/∂v = ∂a_i/∂a）
(dv_dq_local, da_dq_local, da_dv_local, da_da_local) = (
    pin.getJointAccelerationDerivatives(model, data, joint_id, pin.ReferenceFrame.LOCAL)
)
