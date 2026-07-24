import numpy as np
import pinocchio as pin

##
## In this short script, we show how to compute the derivatives of the
## inverse dynamics (RNEA), using the algorithms proposed in:
##
## Analytical Derivatives of Rigid Body Dynamics Algorithms,
## Justin Carpentier and Nicolas Mansard,
## Robotics: Science and Systems, 2018
##

# 构建内置随机人形模型（浮动基，双足）
model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

# 设置关节限位（随机人形模型默认无限位）
model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit = np.ones((model.nq, 1))

q = pin.randomConfiguration(model)     # 关节配置（含浮动基四元数，nq 维）
v = np.random.rand(model.nv, 1)        # 关节速度（nv 维）
a = np.random.rand(model.nv, 1)        # 关节加速度（nv 维）

# 计算 RNEA 及其对 (q, v, a) 的解析梯度
# 等价于对 τ = M(q)a + C(q,v)v + g(q) 分别对 q, v, a 求偏导
# 解析梯度比数值差分精度高（不依赖步长）、速度快（约 10×）
pin.computeRNEADerivatives(model, data, q, v, a)

# ∂τ/∂q：配置对力矩的偏导（nv×nq 矩阵）
# 包含重力梯度 ∂g/∂q 和科氏力/离心力梯度 ∂(Cv)/∂q
# 应用：逆动力学控制的线性化、运动学参数辨识
dtau_dq = data.dtau_dq  # Derivatives of the ID w.r.t. the joint config vector

# ∂τ/∂v：速度对力矩的偏导（nv×nv 矩阵）
# 主要来自科氏力和离心力项 C(q,v)v 对 v 的导数
dtau_dv = data.dtau_dv  # Derivatives of the ID w.r.t. the joint velocity vector

# ∂τ/∂a = M(q)：加速度对力矩的偏导（nv×nv 矩阵）
# 即质量矩阵本身（因为 τ 对 a 的依赖是线性的：τ = M·a + ...）
# 同 CRBA 的计算结果，可直接用作质量矩阵
dtau_da = data.M  # Derivatives of the ID w.r.t. the joint acceleration vector
