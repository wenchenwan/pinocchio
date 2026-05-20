from time import sleep

import example_robot_data
import numpy as np
import pinocchio
from pinocchio.visualize import GepettoVisualizer

# ============================================================
# Talos 人形机器人 KKT 约束姿态优化（静态）
# 演示：四肢多点 6D 约束 + Contact Cholesky 稀疏 KKT 求解
#       + 质心高度控制（squashing 算法）
# 可视化：GepettoViewer（需要单独启动 gepetto-gui）
# ============================================================

# 通过 example_robot_data 加载完整 Talos 人形模型（包含 URDF + SRDF）
robot = example_robot_data.load("talos")
model = robot.model
data  = robot.data

state_name = "half_sitting"
# 从 SRDF 预存的参考配置中读取标准站立姿（浮动基 + 全关节）
robot.q0 = robot.model.referenceConfigurations[state_name]

# 正向运动学：计算参考配置下各 Frame 的世界系位姿（供后续约束构建使用）
pinocchio.forwardKinematics(model, data, robot.q0)

# 四个约束末端：双足（地面接触）+ 双手（目标位置约束）
lfFoot, rfFoot, lhFoot, rhFoot = (
    "left_sole_link",
    "right_sole_link",
    "gripper_left_fingertip_3_link",
    "gripper_right_fingertip_3_link",
)
foot_frames    = [lfFoot, rfFoot, lhFoot, rhFoot]
foot_frame_ids = [robot.model.getFrameId(name) for name in foot_frames]
# 每个 Frame 对应的父关节（接触力作用点）
foot_joint_ids = [
    robot.model.frames[robot.model.getFrameId(name)].parent
    for name in foot_frames
]

# 更新所有 Frame 的世界系位姿到 data.oMf[]
pinocchio.forwardKinematics(model, data, robot.q0)
pinocchio.framesForwardKinematics(model, data, robot.q0)

# ---- 构建 6D 刚性约束模型 ----
# RigidConstraintModel(type, model, joint1_id, joint1_placement, joint2_id, joint2_placement)
# joint2_id=0 表示约束的另一端是世界（固定点）
# joint2_placement 指定约束目标位姿（当前 Frame 在世界系下的实际位姿）
constraint_models = []
for j, frame_id in enumerate(foot_frame_ids):
    contact_model_lf1 = pinocchio.RigidConstraintModel(
        pinocchio.ContactType.CONTACT_6D,
        robot.model,
        foot_joint_ids[j],
        robot.model.frames[frame_id].placement,
        0,                      # 世界坐标系作为约束另一端
        data.oMf[frame_id],    # 目标位姿 = 当前 Frame 的世界系位姿（保持不动）
    )
    constraint_models.extend([contact_model_lf1])

# 修改手部约束目标：将双手抬起到前方指定位置
# 约束索引 2/3 对应左手/右手，此处改变其 joint2_placement（目标位姿）
constraint_models[3].joint2_placement = pinocchio.SE3(
    pinocchio.rpy.rpyToMatrix(np.array([0.0, -np.pi / 2, 0.0])),  # 旋转（绕 Y 轴 -90°）
    np.array([0.6, -0.40, 1.0]),                                    # 右手目标位置
)
constraint_models[2].joint2_placement = pinocchio.SE3(
    pinocchio.rpy.rpyToMatrix(np.array([0, -np.pi / 2, 0.0])),
    np.array([0.6, 0.4, 1.0]),                                      # 左手目标位置
)

robot.setVisualizer(GepettoVisualizer())
robot.initViewer()
robot.loadViewerModel("pinocchio")
gui = robot.viewer.gui
robot.display(robot.q0)
window_id = robot.viewer.gui.getWindowID("python-pinocchio")
robot.viewer.gui.setBackgroundColor1(window_id, [1.0, 1.0, 1.0, 1.0])
robot.viewer.gui.setBackgroundColor2(window_id, [1.0, 1.0, 1.0, 1.0])
robot.viewer.gui.addFloor("hpp-gui/floor")
robot.viewer.gui.setScale("hpp-gui/floor", [0.5, 0.5, 0.5])
robot.viewer.gui.setColor("hpp-gui/floor", [0.7, 0.7, 0.7, 1.0])
robot.viewer.gui.setLightingMode("hpp-gui/floor", "OFF")
robot.display(robot.q0)

constraint_datas = [cm.createData() for cm in constraint_models]

q = robot.q0.copy()

# 计算全部动力学量（质量矩阵、质心、Jacobian 等），存入 data
pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))

# 构建 Contact Cholesky 分解器：利用运动树稀疏结构高效求解 KKT 系统
# 比朴素 LU 分解快，随约束数量线性缩放
kkt_constraint = pinocchio.ContactCholeskyDecomposition(model, constraint_models)

# 约束总维度：4 个 6D 接触各贡献 6 维 = 24 维
constraint_dim = sum([cm.size() for cm in constraint_models])

N   = 100000  # 最大迭代次数
eps = 1e-10   # 收敛阈值
mu  = 1e-8    # KKT 正则化参数（防止约束退化时数值不稳定）

q_sol = q.copy()
robot.display(q_sol)

mass = data.mass[0]  # 机器人总质量（用于质心力矩计算）


def squashing(model, data, q_in):
    """
    约束满足的关节配置求解（Squashing 算法）

    在保持所有四肢约束的同时，最小化质心高度误差。
    使用增广拉格朗日 / ADMM 迭代：
      - 原始变量 q：关节配置
      - 对偶变量 y：约束力（Lagrange 乘子）
    KKT 系统：
      [M   J^T] [dq] = [rhs_constraint + rhs_com]
      [J   μI ] [dy]
    """
    q = q_in.copy()
    y = np.ones(constraint_dim)  # 初始化对偶变量（Lagrange 乘子）

    N_full = 200
    com_drop_amp = 0.1   # 质心下降幅度（让 Talos 略微下蹲）
    pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))
    com_base = data.com[0].copy()
    kp    = 1.0   # 质心误差的比例增益
    speed = 1.0

    def com_des(k):
        # 质心目标：在竖直方向做正弦下降（模拟下蹲动作）
        return com_base - np.array(
            [0.0, 0.0, np.abs(com_drop_amp * np.sin(2.0 * np.pi * k * speed / N_full))]
        )

    for k in range(N):
        # 更新所有动力学量和 Jacobian
        pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))
        pinocchio.computeJointJacobians(model, data, q)
        pinocchio.computeJointJacobians(model, data, q)  # 调用两次确保缓存一致

        com_act = data.com[0].copy()
        com_err = com_act - com_des(k)

        # 构建 Contact Cholesky 分解（每步更新，因为 M 和 J 都随 q 变化）
        kkt_constraint.compute(model, data, constraint_models, constraint_datas, mu)

        # 约束残差：log6(c1Mc2) 衡量约束两端的 SE(3) 差异
        # c1Mc2 是约束 Frame 1 到 Frame 2 的相对变换，目标是使其为单位元
        constraint_value = np.concatenate(
            [pinocchio.log6(cd.c1Mc2) for cd in constraint_datas]
        )

        # 各约束的 Frame Jacobian（在约束参考系下）
        J = np.vstack(
            [
                pinocchio.getFrameJacobian(
                    model, data, cm.joint1_id, cm.joint1_placement, cm.reference_frame
                )
                for cm in constraint_models
            ]
        )

        # 原始可行性：约束残差的无穷范数
        primal_feas = np.linalg.norm(constraint_value, np.inf)
        print(J.shape, constraint_value.shape, y.shape)
        # 对偶可行性：J^T(constraint + y) 衡量 KKT 梯度条件
        dual_feas = np.linalg.norm(J.T.dot(constraint_value + y), np.inf)
        print("primal_feas:", primal_feas)
        print("dual_feas:", dual_feas)
        print("constraint_value:", np.linalg.norm(constraint_value))
        print("com_error:", np.linalg.norm(com_err))

        # 构造 KKT 右端向量（rhs）
        # 结构：[-约束残差 - μ*y,  质心误差力,  关节空间零（不施加额外任务）]
        rhs = np.concatenate(
            [-constraint_value - y * mu, kp * mass * com_err, np.zeros(model.nv - 3)]
        )

        # 求解 KKT 系统：dz = [dy; dq]
        dz = kkt_constraint.solve(rhs)
        dy = dz[:constraint_dim]   # 对偶变量增量
        dq = dz[constraint_dim:]   # 关节配置增量

        alpha = 1.0  # 步长（无线搜索，固定为 1）

        # 在流形上更新配置（李群积分）
        q = pinocchio.integrate(model, q, -alpha * dq)
        # 更新对偶变量（增广拉格朗日乘子更新）
        y -= alpha * (-dy + y)

        robot.display(q)
        sleep(0.05)
    return q


q_new = squashing(model, data, robot.q0)
