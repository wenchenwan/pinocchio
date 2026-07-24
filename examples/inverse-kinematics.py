import numpy as np
import pinocchio
from numpy.linalg import norm, solve

# 构建内置 6-DOF 串联机械臂
model = pinocchio.buildSampleModelManipulator()
data = model.createData()

# 目标关节 ID（末端关节）
JOINT_ID = 6

# 目标位姿：旋转=单位矩阵（无旋转），平移=[1, 0, 1]
# SE3(R, p)：R 是 3×3 旋转矩阵，p 是 3D 平移向量
oMdes = pinocchio.SE3(np.eye(3), np.array([1.0, 0.0, 1.0]))

# 初始配置（零位）
q = pinocchio.neutral(model)

eps    = 1e-4   # 收敛阈值：误差 L2 范数小于此值则认为收敛
IT_MAX = 1000   # 最大迭代次数，防止不收敛时死循环
DT     = 1e-1   # 积分步长（类似梯度下降的学习率）
damp   = 1e-12  # 阻尼系数 λ，用于阻尼最小二乘防止 Jacobian 奇异

i = 0
while True:
    # 第一步：正向运动学，更新所有关节的 SE(3) 变换到 data.oMi[]
    pinocchio.forwardKinematics(model, data, q)

    # 第二步：计算末端关节局部系下的相对变换
    # iMd = T_i^{-1} * T_des，即从当前末端姿态到目标姿态的变换
    # actInv(T) 等价于 self.inverse() * T，但数值更稳定
    iMd = data.oMi[JOINT_ID].actInv(oMdes)

    # 第三步：计算 6D 误差向量（SE(3) 上的对数映射）
    # log(iMd) 将 SE(3) 元素映射到其李代数 se(3)≅R^6
    # 结果为 [平移误差(3D), 旋转误差(3D)] 的拼接，在关节局部系下表达
    err = pinocchio.log(iMd).vector

    if norm(err) < eps:
        success = True
        break
    if i >= IT_MAX:
        success = False
        break

    # 第四步：计算关节空间几何 Jacobian（在关节局部系下）
    # J ∈ R^{6×nv}，将关节速度映射到末端空间速度
    J = pinocchio.computeJointJacobian(model, data, q, JOINT_ID)

    # 第五步：修正 Jacobian（链式法则，补偿对数映射的非线性）
    # 直接用 J 作梯度下降在 SE(3) 上是错误的，需乘以 log 映射的 Jacobian J_log6
    # J_log6(iMd^{-1}) ∈ R^{6×6}，将 se(3) 中的误差正确反传到切空间
    # 前置负号：因为误差定义为 iMd（从末端到目标），梯度方向取反
    J = -np.dot(pinocchio.Jlog6(iMd.inverse()), J)

    # 第六步：阻尼最小二乘求关节速度增量
    # 最小化 ||J*v - (-err)||^2 + λ^2*||v||^2
    # 解析解：v = -J^T (J J^T + λ^2 I)^{-1} err
    # 阻尼 λ 防止 Jacobian 接近奇异时解爆炸（λ=0 即普通最小二乘）
    v = -J.T.dot(solve(J.dot(J.T) + damp * np.eye(6), err))

    # 第七步：在流形（李群）上做积分更新配置
    # 不能用普通加法 q += v*DT，因为旋转部分（四元数）不是线性空间
    # integrate(q, Δq) 等价于 q * exp(Δq)，保证四元数归一化等约束
    q = pinocchio.integrate(model, q, v * DT)

    if not i % 10:
        print(f"{i}: error = {err.T}")
    i += 1

if success:
    print("Convergence achieved!")
else:
    print(
        "\n"
        "Warning: the iterative algorithm has not reached convergence "
        "to the desired precision"
    )

print(f"\nresult: {q.flatten().tolist()}")
print(f"\nfinal error: {err.T}")
