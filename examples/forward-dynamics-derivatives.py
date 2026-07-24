import numpy as np
import pinocchio as pin

##
## In this short script, we show how to compute the derivatives of the
## forward dynamics, using the algorithms explained in:
##
## Analytical Derivatives of Rigid Body Dynamics Algorithms,
## Justin Carpentier and Nicolas Mansard,
## Robotics: Science and Systems, 2018
##

# 构建内置随机人形模型（双足，带浮动基 FreeFlyer）
# 该模型包含躯干 + 双臂 + 双腿，共约 32 DOF
model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

# 设置关节限位（样例人形模型默认无限位，优化器/随机采样需要此项）
model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit = np.ones((model.nq, 1))

q   = pin.randomConfiguration(model)       # 关节配置（nq 维，含浮动基四元数）
v   = np.random.rand(model.nv, 1)          # 关节速度（nv 维）
tau = np.random.rand(model.nv, 1)          # 关节力矩（nv 维）

# 计算 ABA（关节化体算法）的正向动力学及其对 (q, v, τ) 的解析梯度
# 一次调用同时得到 q̈ 和所有梯度，无需多次调用或数值差分
pin.computeABADerivatives(model, data, q, v, tau)

# ∂q̈/∂q：配置对加速度的影响（nv×nq 矩阵）
# 在 DDP/iLQR 中对应线性化动力学的 ∂f/∂q 块（状态转移矩阵 A 的 q 列）
ddq_dq = data.ddq_dq  # Derivatives of the FD w.r.t. the joint config vector

# ∂q̈/∂v：速度对加速度的影响（nv×nv 矩阵）
# 在 DDP/iLQR 中对应状态转移矩阵 A 的 v 列
ddq_dv = data.ddq_dv  # Derivatives of the FD w.r.t. the joint velocity vector

# ∂q̈/∂τ = M(q)^{-1}：力矩对加速度的影响（nv×nv 矩阵）
# 即质量矩阵的逆，也是 DDP/iLQR 中的控制矩阵 B
# 注意：data.Minv 是稠密矩阵，由 ABA 算法的第三趟递推直接得到，无需显式求逆
ddq_dtau = data.Minv  # Derivatives of the FD w.r.t. the joint acceleration vector
