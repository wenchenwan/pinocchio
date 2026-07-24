from time import sleep

import numpy as np
import pinocchio
from example_robot_data import load

# ============================================================
# ANYmal 四足机器人 KKT 约束姿态优化
# 演示：四足接触约束 + Contact Cholesky KKT 求解 + 质心高度控制
# 与 talos-simulation.py 的逻辑相同，但使用四足机器人（3D 点接触）
# ============================================================

# 通过 example_robot_data 加载 ANYmal 模型（浮动基四足）
robot = load("anymal")
model = robot.model
data  = robot.data

state_name = "standing"
# 加载标准站立姿势
robot.q0 = robot.model.referenceConfigurations[state_name]

pinocchio.forwardKinematics(model, data, robot.q0)

# ANYmal 四条腿的足端 Frame 名称
lfFoot, rfFoot, lhFoot, rhFoot = "LF_FOOT", "RF_FOOT", "LH_FOOT", "RH_FOOT"
foot_frames    = [lfFoot, rfFoot, lhFoot, rhFoot]
foot_frame_ids = [robot.model.getFrameId(name) for name in foot_frames]
# 每个足端 Frame 的父关节 ID
foot_joint_ids = [
    robot.model.frames[robot.model.getFrameId(name)].parent
    for name in foot_frames
]

pinocchio.forwardKinematics(model, data, robot.q0)
pinocchio.framesForwardKinematics(model, data, robot.q0)

# ---- 四足 3D 点接触约束 ----
# 与 Talos 的 6D 接触不同，四足常用 3D 点接触（只约束位置，不约束姿态）
# joint2_id=0：约束另一端为世界（地面固定点）
# joint2_placement：当前足端在世界系下的位姿（保持足端不动）
constraint_models = []
for j, frame_id in enumerate(foot_frame_ids):
    contact_model_lf1 = pinocchio.RigidConstraintModel(
        pinocchio.ContactType.CONTACT_3D,   # 3D 点接触（仅位置约束，每足贡献 3 维）
        robot.model,
        foot_joint_ids[j],
        robot.model.frames[frame_id].placement,
        0,
        data.oMf[frame_id],
    )
    constraint_models.extend([contact_model_lf1])

robot.initViewer()
robot.loadViewerModel("pinocchio")
robot.display(robot.q0)

constraint_datas = [cm.createData() for cm in constraint_models]
q = robot.q0.copy()

# 一次性计算质量矩阵、Coriolis、重力、质心等所有动力学量
pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))

# Contact Cholesky 分解：稀疏 KKT 求解器
kkt_constraint = pinocchio.ContactCholeskyDecomposition(model, constraint_models)

# 四足 3D 接触：每足 3 维，共 12 维约束
constraint_dim = sum([cm.size() for cm in constraint_models])
N   = 100000
eps = 1e-10
mu  = 0.0   # 无正则化（约束 Jacobian 满秩时可以设为 0）

q_sol = q.copy()
robot.display(q_sol)

mass = data.mass[0]  # 总质量


def squashing(model, data, q_in):
    """
    增广拉格朗日迭代求解约束满足的关节配置
    目标：在保持四足接触的同时，使质心向目标位置移动（下蹲）
    """
    q = q_in.copy()
    y = np.ones(constraint_dim)   # Lagrange 乘子（对偶变量）初始化为 1

    N_full = 200
    com_drop_amp = 0.2   # 质心下降幅度 0.2 m
    pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))
    com_base = data.com[0].copy()  # 记录初始质心位置
    kp = 1.0
    speed = 1.0

    def com_des(k):
        # 质心目标轨迹：在 Z 方向做正弦变化（模拟下蹲/起立）
        return com_base - np.array(
            [0.0, 0.0, np.abs(com_drop_amp * np.sin(2.0 * np.pi * k * speed / N_full))]
        )

    for k in range(N):
        pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))
        pinocchio.computeJointJacobians(model, data, q)
        pinocchio.computeJointJacobians(model, data, q)

        com_act = data.com[0].copy()
        com_err = com_act - com_des(k)

        # 重新计算 Contact Cholesky 分解（q 变化后 M 和 J 也变化）
        kkt_constraint.compute(model, data, constraint_models, constraint_datas, mu)

        # 3D 点接触：约束残差为足端平移误差（不含旋转）
        constraint_value = np.concatenate(
            [cd.c1Mc2.translation for cd in constraint_datas]
        )

        # 各足端的 3D Jacobian（取 Frame Jacobian 的前 3 行）
        J = np.vstack(
            [
                pinocchio.getFrameJacobian(
                    model, data, cm.joint1_id, cm.joint1_placement, cm.reference_frame
                )[:3, :]
                for cm in constraint_models
            ]
        )

        primal_feas = np.linalg.norm(constraint_value, np.inf)
        dual_feas   = np.linalg.norm(J.T.dot(constraint_value + y), np.inf)
        print("primal_feas:", primal_feas)
        print("dual_feas:", dual_feas)
        print("constraint_value:", np.linalg.norm(constraint_value))
        print("com_error:", np.linalg.norm(com_err))

        # KKT 右端：约束残差项 + 质心误差项 + 零填充
        rhs = np.concatenate(
            [-constraint_value - y * mu, kp * mass * com_err, np.zeros(model.nv - 3)]
        )

        # 稀疏 KKT 求解：dz = [dy（对偶增量）; dq（原始增量）]
        dz = kkt_constraint.solve(rhs)
        dy = dz[:constraint_dim]
        dq = dz[constraint_dim:]

        alpha = 1.0
        q  = pinocchio.integrate(model, q, -alpha * dq)  # 流形积分更新配置
        y -= alpha * (-dy + y)                            # 对偶变量更新

        robot.display(q)
        sleep(0.05)
    return q


q_new = squashing(model, data, robot.q0)
