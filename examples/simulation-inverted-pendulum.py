import math
import sys
import time

import hppfcl as fcl
import numpy as np
import pinocchio as pin
from pinocchio.visualize import GepettoVisualizer

# ============================================================
# 多节倒立摆仿真（N 节串联，每节绕 Y 轴旋转）
# 与 simulation-pendulum.py 的区别：
#   - 使用 GepettoVisualizer（而非 Meshcat），适合本地桌面应用
#   - 绕 Y 轴旋转（JointModelRY），在 XZ 平面内摆动
#   - 无阻尼（纯自由落体，能量守恒），适合验证积分器
#   - 无驱动力矩（τ = 0）：纯被动系统
# ============================================================

N          = 10     # 串联摆的节数（默认 10 节）
model      = pin.Model()
geom_model = pin.GeometryModel()

parent_id       = 0
joint_placement = pin.SE3.Identity()   # 第一节摆的关节在世界坐标系原点
body_mass       = 1.0
body_radius     = 0.1

# 基座球体（固定在原点，仅用于可视化）
shape0    = fcl.Sphere(body_radius)
geom0_obj = pin.GeometryObject("base", 0, pin.SE3.Identity(), shape0)
geom0_obj.meshColor = np.array([1.0, 0.1, 0.1, 1.0])   # 红色基座
geom_model.addGeometryObject(geom0_obj)

# ---- 构建 N 节串联倒立摆 ----
for k in range(N):
    joint_name = "joint_" + str(k + 1)
    # JointModelRY：绕 Y 轴旋转（1 DOF），摆在 XZ 平面内运动
    joint_id = model.addJoint(
        parent_id, pin.JointModelRY(), joint_placement, joint_name
    )

    # 每节连杆：质点在关节正上方 1.0 m（Z 方向）
    body_inertia   = pin.Inertia.FromSphere(body_mass, body_radius)
    body_placement = joint_placement.copy()
    body_placement.translation[2] = 1.0   # 摆长 = 1.0 m
    model.appendBodyToJoint(joint_id, body_inertia, body_placement)

    # 球形质点（白色）
    geom1_name = "ball_" + str(k + 1)
    shape1     = fcl.Sphere(body_radius)
    geom1_obj  = pin.GeometryObject(geom1_name, joint_id, body_placement, shape1)
    geom1_obj.meshColor = np.ones(4)
    geom_model.addGeometryObject(geom1_obj)

    # 连接杆（黑色细圆柱，从关节到质点中点）
    geom2_name      = "bar_" + str(k + 1)
    shape2          = fcl.Cylinder(body_radius / 4.0, body_placement.translation[2])
    shape2_placement = body_placement.copy()
    shape2_placement.translation[2] /= 2.0
    geom2_obj = pin.GeometryObject(geom2_name, joint_id, shape2_placement, shape2)
    geom2_obj.meshColor = np.array([0.0, 0.0, 0.0, 1.0])
    geom_model.addGeometryObject(geom2_obj)

    # 下一节的父关节为当前质心（串联结构）
    parent_id       = joint_id
    joint_placement = body_placement.copy()


visual_model = geom_model   # 视觉与碰撞共用同一几何模型
# GepettoVisualizer：基于 gepetto-viewer 的 3D 可视化（需要安装并启动 gepetto-viewer-server）
viz = GepettoVisualizer(model, geom_model, visual_model)

# Initialize the viewer.
try:
    viz.initViewer()
except ImportError as err:
    print(
        "Error while initializing the viewer. "
        "It seems you should install gepetto-viewer"
    )
    print(err)
    sys.exit(0)

try:
    viz.loadViewerModel("pinocchio")
except AttributeError as err:
    print(
        "Error while loading the viewer model. "
        "It seems you should start gepetto-viewer"
    )
    print(err)
    sys.exit(0)

# Display a robot configuration.
q0 = pin.neutral(model)
viz.display(q0)

# ---- 仿真参数 ----
dt = 0.01
T  = 5
N  = math.floor(T / dt)

# 关节角范围 [-π, +π]（用于 randomConfiguration）
model.lowerPositionLimit.fill(-math.pi)
model.upperPositionLimit.fill(+math.pi)
q = pin.randomConfiguration(model)   # 随机初始配置（倒立摆从随机姿态开始）
v = np.zeros(model.nv)               # 初始速度为零

t        = 0.0
data_sim = model.createData()

for k in range(N):
    # 无驱动力矩（τ = 0）：纯被动系统，仅受重力作用
    tau_control = np.zeros(model.nv)
    # ABA 正向动力学：q̈ = ABA(q, v, τ)
    a = pin.aba(model, data_sim, q, v, tau_control)

    # 半隐式欧拉积分
    v += a * dt
    # 注意原代码中被注释掉的 q += v*dt：
    # 对于旋转关节，直接做向量加法不正确（会破坏四元数归一化）
    # 必须使用 integrate()：在 SO(2)/SO(3) 流形上做正确的 Exp 映射
    q = pin.integrate(model, q, v * dt)

    viz.display(q)
    time.sleep(dt)
    t += dt
