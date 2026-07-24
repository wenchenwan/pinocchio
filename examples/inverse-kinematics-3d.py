import numpy as np
import pinocchio
from numpy.linalg import norm, solve

model = pinocchio.buildSampleModelManipulator()
data = model.createData()

JOINT_ID = 6
# 目标位姿（仅使用其平移部分，姿态被忽略）
oMdes = pinocchio.SE3(np.eye(3), np.array([1.0, 0.0, 1.0]))

q      = pinocchio.neutral(model)
eps    = 1e-4
IT_MAX = 1000
DT     = 1e-1
damp   = 1e-12  # 阻尼系数，防止近奇异时解爆炸

it = 0
while True:
    pinocchio.forwardKinematics(model, data, q)

    # 计算末端局部系下的相对变换 T_i^{-1} * T_des
    iMd = data.oMi[JOINT_ID].actInv(oMdes)

    # 与 6D IK 的关键区别：
    # 只取 iMd 的平移分量作为误差（R^3），不关心旋转误差
    # 无需 SE(3) 对数映射，平移误差本身就在线性空间中
    err = iMd.translation

    if norm(err) < eps:
        success = True
        break
    if it >= IT_MAX:
        success = False
        break

    # 计算 6×nv 几何 Jacobian（在关节局部系下）
    J = pinocchio.computeJointJacobian(model, data, q, JOINT_ID)  # in joint frame

    # 只取上 3 行（线速度部分），忽略下 3 行（角速度部分）
    # 取负号：局部系下平移误差相对于 Jacobian 的符号约定
    # 无需 Jlog6 修正（仅 6D IK 才需要，因为那里用了 log 映射）
    J = -J[:3, :]  # linear part of the Jacobian

    # 3×3 阻尼最小二乘，比 6D IK 的 6×6 系统计算量更小
    v = -J.T.dot(solve(J.dot(J.T) + damp * np.eye(3), err))

    # 流形积分（旋转部分依然需要正确的李群积分）
    q = pinocchio.integrate(model, q, v * DT)
    if not it % 10:
        print(f"{it}: error = {err.T}")
    it += 1

if success:
    print("Convergence achieved!")
else:
    print(
        "\nWarning: the iterative algorithm has not reached convergence to "
        "the desired precision"
    )

print(f"\nresult: {q.flatten().tolist()}")
print(f"\nfinal error: {err.T}")
