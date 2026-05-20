from time import sleep

import example_robot_data
import numpy as np
import pinocchio
from pinocchio.visualize import GepettoVisualizer

# ============================================================
# Cassie 双足机器人约束姿态优化
# Cassie 特点：闭环运动链（四杆机构），同时包含：
#   - URDF 内的闭链约束（robot.constraint_models 已预加载）
#   - 足底与地面的接触约束（手动添加）
# ============================================================

robot = example_robot_data.load("cassie")

# Cassie 的 URDF 包含闭链约束（四杆机构），已由 example_robot_data 解析
# robot.constraint_models：预加载的闭链 RigidConstraintModel 列表
constraint_models = robot.constraint_models

data = pinocchio.Data(robot.model)

pinocchio.forwardKinematics(robot.model, robot.data, robot.q0)

# Cassie 的足底平面关节（用于定义地面接触点）
foot_joints    = ["left-plantar-foot-joint", "right-plantar-foot-joint"]
foot_joint_ids = [robot.model.getJointId(jn) for jn in foot_joints]

# 足底前后两点的局部位置（定义足底板的接触几何）
front_placement = pinocchio.SE3(np.identity(3), np.array([-0.1,  0.11, 0.0]))
back_placement  = pinocchio.SE3(np.identity(3), np.array([ 0.03, -0.0, 0.0]))
foot_length     = np.linalg.norm(front_placement.translation - back_placement.translation)

mid_point = np.zeros(3)  # 记录所有接触点的几何中心（用于可视化）

for joint_id in foot_joint_ids:
    # 足底前点在世界系下的目标位置（当前位置投影到地面 z=0）
    joint2_placement = robot.data.oMi[joint_id] * front_placement
    joint2_placement.translation[2] = 0  # 投影到地面

    # 添加足底前点的 3D 点接触约束
    contact_model_lf1 = pinocchio.RigidConstraintModel(
        pinocchio.ContactType.CONTACT_3D,
        robot.model,
        joint_id,
        front_placement,          # 接触点在足底关节局部系下的位置
        0,                        # 约束另一端：世界（地面）
        joint2_placement,         # 地面上的目标位置
        pinocchio.ReferenceFrame.LOCAL,
    )
    mid_point += joint2_placement.translation

    # 足底后点的地面目标位置
    joint2_placement.translation[0] -= foot_length

    # 添加足底后点的 3D 点接触约束
    contact_model_lf2 = pinocchio.RigidConstraintModel(
        pinocchio.ContactType.CONTACT_3D,
        robot.model,
        joint_id,
        back_placement,
        0,
        joint2_placement,
        pinocchio.ReferenceFrame.LOCAL,
    )
    mid_point += joint2_placement.translation

    # 追加到约束列表（在原有闭链约束之后）
    constraint_models.extend([contact_model_lf1, contact_model_lf2])

mid_point /= 4.0  # 四个接触点的平均位置（几何中心）

robot.setVisualizer(GepettoVisualizer())
robot.initViewer()
robot.loadViewerModel("pinocchio")
gui = robot.viewer.gui
robot.display(robot.q0)
q0 = robot.q0.copy()

# 为所有约束创建数据容器（包括原有闭链约束 + 新增的地面接触约束）
constraint_datas = pinocchio.StdVec_RigidConstraintData()
for cm in constraint_models:
    constraint_datas.append(cm.createData())


def update_axis(q):
    """更新可视化中约束两端的坐标轴显示"""
    pinocchio.forwardKinematics(robot.model, robot.data, q)
    for j, cm in enumerate(robot.constraint_models):
        pos1 = robot.data.oMi[cm.joint1_id] * cm.joint1_placement
        pos2 = robot.data.oMi[cm.joint2_id] * cm.joint2_placement
        name1 = "hpp-gui/cm1_" + str(j)
        name2 = "hpp-gui/cm2_" + str(j)
        gui.applyConfiguration(name1, list(pinocchio.SE3ToXYZQUAT(pos1)))
        gui.applyConfiguration(name2, list(pinocchio.SE3ToXYZQUAT(pos2)))
        gui.refresh()


def check_joint(model, nj):
    """可视化单个关节的运动范围（调试用）"""
    for k in range(model.joints[nj].nv):
        for res in range(200):
            q1 = robot.q0.copy()
            theta = res * 2.0 * np.pi / 200.0
            v = np.zeros(robot.model.nv)
            v[robot.model.idx_vs[nj] + k] = theta
            q1 = pinocchio.integrate(robot.model, q1, v)
            robot.display(q1)
            update_axis(q1)
            sleep(0.005)


q = q0.copy()

pinocchio.computeAllTerms(robot.model, robot.data, q, np.zeros(robot.model.nv))

# Contact Cholesky 分解器（同时处理闭链约束和地面接触约束）
kkt_constraint = pinocchio.ContactCholeskyDecomposition(robot.model, constraint_models)
constraint_dim  = sum([cm.size() for cm in constraint_models])

N   = 1000
eps = 1e-6
mu  = 1e-4   # KKT 正则化（Cassie 闭链约束可能导致 Jacobian 退化）

q_sol = q.copy()

# --- 可视化设置 ---
window_id = robot.viewer.gui.getWindowID("python-pinocchio")
robot.viewer.gui.setBackgroundColor1(window_id, [1.0, 1.0, 1.0, 1.0])
robot.viewer.gui.setBackgroundColor2(window_id, [1.0, 1.0, 1.0, 1.0])
robot.viewer.gui.addFloor("hpp-gui/floor")
robot.viewer.gui.setScale("hpp-gui/floor", [0.5, 0.5, 0.5])
robot.viewer.gui.setColor("hpp-gui/floor", [0.7, 0.7, 0.7, 1.0])
robot.viewer.gui.setLightingMode("hpp-gui/floor", "OFF")

axis_size   = 0.08
radius      = 0.005
transparency = 0.5

# 可视化每个约束的两端坐标轴（闭链约束的两端 Frame）
for j, cm in enumerate(constraint_models):
    pos1 = robot.data.oMi[cm.joint1_id] * cm.joint1_placement
    pos2 = robot.data.oMi[cm.joint2_id] * cm.joint2_placement
    name1 = "hpp-gui/cm1_" + str(j)
    name2 = "hpp-gui/cm2_" + str(j)
    red_color = 1.0 * float(j) / float(len(constraint_models))
    robot.viewer.gui.addXYZaxis(
        name1, [red_color, 1.0, 1.0 - red_color, transparency], radius, axis_size
    )
    robot.viewer.gui.addXYZaxis(
        name2, [red_color, 1.0, 1.0 - red_color, transparency], radius, axis_size
    )
    gui.applyConfiguration(name1, list(pinocchio.SE3ToXYZQUAT(pos1)))
    gui.applyConfiguration(name2, list(pinocchio.SE3ToXYZQUAT(pos2)))
    gui.refresh()
    gui.setVisibility(name1, "OFF")
    gui.setVisibility(name2, "OFF")

mid_point_name = "hpp-gui/mid_point"
mid_point_pos  = pinocchio.SE3(np.identity(3), mid_point)
robot.viewer.gui.addXYZaxis(
    mid_point_name, [0.0, 0.0, 1.0, transparency], radius, axis_size
)
gui.applyConfiguration(mid_point_name, list(pinocchio.SE3ToXYZQUAT(mid_point_pos)))
gui.setVisibility(mid_point_name, "ALWAYS_ON_TOP")
robot.display(q_sol)

mass    = robot.data.mass[0]
com_base = robot.data.com[0].copy()
com_base[:2] = mid_point[:2]   # 质心 XY 对齐到接触点几何中心
com_drop_amp = 0.1


def com_des(k):
    # 目标质心：降低 Z 高度（下蹲），XY 保持在接触点中心上方
    return com_base - np.array([0.0, 0.0, np.abs(com_drop_amp)])


def squashing(model, data, q_in, Nin=N, epsin=eps, verbose=True):
    """
    增广拉格朗日迭代：同时满足闭链约束和地面接触约束，并控制质心位置。
    Cassie 的特殊性：constraint_models 包含两类约束：
      - 前半段（闭链）：使用 log6（6D SE(3) 残差）
      - 后半段（地面接触）：使用平移残差（3D）
    """
    q = q_in.copy()
    y = np.ones(constraint_dim)

    pinocchio.computeAllTerms(robot.model, robot.data, q, np.zeros(robot.model.nv))
    # 各方向的质心误差增益（Z 方向用更小的增益，避免过快下蹲）
    kp = np.array([1.0, 1.0, 0.1])

    for k in range(Nin):
        pinocchio.computeAllTerms(robot.model, robot.data, q, np.zeros(robot.model.nv))
        pinocchio.computeJointJacobians(robot.model, robot.data, q)
        com_act = robot.data.com[0].copy()
        com_err = com_act - com_des(k)

        kkt_constraint.compute(
            robot.model, robot.data, constraint_models, constraint_datas, mu
        )

        # 闭链约束（前段）：6D SE(3) 残差（log 映射）
        constraint_value1 = np.concatenate(
            [pinocchio.log(cd.c1Mc2) for cd in constraint_datas[:-4]]
        )
        # 地面接触约束（后段）：3D 平移残差
        constraint_value2 = np.concatenate(
            [cd.c1Mc2.translation for cd in constraint_datas[-4:]]
        )
        constraint_value = np.concatenate([constraint_value1, constraint_value2])

        # 闭链约束 Jacobian（完整 6×nv）
        J1 = np.vstack([
            pinocchio.getFrameJacobian(
                robot.model, data, cm.joint1_id, cm.joint1_placement, cm.reference_frame
            )
            for cm in constraint_models[:-4]
        ])
        # 地面接触 Jacobian（只取前 3 行：线速度）
        J2 = np.vstack([
            pinocchio.getFrameJacobian(
                robot.model, data, cm.joint1_id, cm.joint1_placement, cm.reference_frame
            )[:3, :]
            for cm in constraint_models[-4:]
        ])
        J = np.vstack([J1, J2])

        primal_feas = np.linalg.norm(constraint_value, np.inf)
        dual_feas   = np.linalg.norm(J.T.dot(constraint_value + y), np.inf)

        # 双重收敛判据：原始可行性（约束残差）+ 对偶可行性（KKT 梯度）
        if primal_feas < epsin and dual_feas < epsin:
            print("Convergence achieved")
            break
        if verbose:
            print("constraint_value:", np.linalg.norm(constraint_value))
            print("com_error:", np.linalg.norm(com_err))
            print("com_des", com_des(k))
            print("com_act", com_act)

        rhs = np.concatenate(
            [-constraint_value - y * mu, kp * com_err, np.zeros(robot.model.nv - 3)]
        )
        dz = kkt_constraint.solve(rhs)
        dy = dz[:constraint_dim]
        dq = dz[constraint_dim:]

        alpha = 1.0
        q = pinocchio.integrate(robot.model, q, -alpha * dq)
        y -= alpha * (-dy + y)

        robot.display(q)
        update_axis(q)
        # 更新质心目标位置的可视化标记
        gui.applyConfiguration(
            mid_point_name,
            list(pinocchio.SE3ToXYZQUAT(pinocchio.SE3(np.identity(3), com_des(k)))),
        )
        sleep(0.0)
    return q


# 多轮迭代直到收敛（每轮从上一轮结果重新出发）
qin = robot.q0.copy()
for k in range(1000):
    qout = squashing(robot.model, data, qin)
    # 重置浮动基旋转为单位（消除数值漂移导致的基座旋转）
    qout[3:6] = 0.0
    qout[6]   = 1.0
    qin = qout.copy()

# 最终精细化迭代（更高精度）
qout = squashing(robot.model, data, qin, Nin=100000, epsin=1e-10, verbose=False)
