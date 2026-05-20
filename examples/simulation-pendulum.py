# Parse input arguments
import argparse
import math
import sys
import time

import hppfcl as fcl
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer as Visualizer

# ============================================================
# 单/多节摆仿真（可选小车-摆系统）
# 用法：
#   python simulation-pendulum.py            → 单摆
#   python simulation-pendulum.py -N 5       → 5 节串联摆
#   python simulation-pendulum.py --with-cart → 小车-摆（Cart-Pole）
#
# 演示内容：
#   - 手工编程构建多体模型（无 URDF）
#   - ABA 正向动力学
#   - 半隐式欧拉积分
#   - Meshcat 可视化
# ============================================================

parser = argparse.ArgumentParser()
parser.add_argument(
    "--with-cart",
    help="Add a cart at the base of the pendulum to simulate a cart pole system.",
    action="store_true",
)
parser.add_argument(
    "-N", help="Number of pendulums compositing the dynamical system.", type=int
)
args = parser.parse_args()

if args.N:
    N = args.N
else:
    N = 1  # 默认单摆

model      = pin.Model()
geom_model = pin.GeometryModel()

parent_id = 0   # 当前父关节 ID（从世界/宇宙关节开始）

if args.with_cart:
    # ---- 小车（Cart）：沿 Y 轴平移的棱柱关节（JointModelPY）----
    # Cart-Pole 系统：底部小车水平滑动，顶部摆自由旋转
    cart_radius = 0.1
    cart_length = 5 * cart_radius
    cart_mass   = 2.0
    joint_name  = "joint_cart"

    # 旋转几何体使圆柱轴从 Z 轴（FCL 默认）旋转到 Y 轴（小车滑动方向）
    geometry_placement = pin.SE3.Identity()
    geometry_placement.rotation = pin.Quaternion(
        np.array([0.0, 0.0, 1.0]), np.array([0.0, 1.0, 0.0])
    ).toRotationMatrix()

    # JointModelPY：沿 Y 轴的单自由度平移关节（Prismatic joint along Y）
    joint_id = model.addJoint(
        parent_id, pin.JointModelPY(), pin.SE3.Identity(), joint_name
    )

    body_inertia   = pin.Inertia.FromCylinder(cart_mass, cart_radius, cart_length)
    body_placement = geometry_placement
    model.appendBodyToJoint(joint_id, body_inertia, body_placement)

    shape_cart = fcl.Cylinder(cart_radius, cart_length)
    geom_cart  = pin.GeometryObject("shape_cart", joint_id, geometry_placement, shape_cart)
    geom_cart.meshColor = np.array([1.0, 0.1, 0.1, 1.0])
    geom_model.addGeometryObject(geom_cart)

    parent_id = joint_id   # 摆的父关节变为小车关节
else:
    # ---- 固定基底（球体，仅用于可视化）----
    base_radius = 0.2
    shape_base  = fcl.Sphere(base_radius)
    geom_base   = pin.GeometryObject("base", 0, pin.SE3.Identity(), shape_base)
    geom_base.meshColor = np.array([1.0, 0.1, 0.1, 1.0])
    geom_model.addGeometryObject(geom_base)

# ---- 构建 N 节串联摆 ----
# 每节摆：绕 X 轴旋转关节（JointModelRX） + 质点连杆
joint_placement = pin.SE3.Identity()   # 初始关节位置（原点）
body_mass       = 1.0
body_radius     = 0.1

for k in range(N):
    joint_name = "joint_" + str(k + 1)
    # JointModelRX：绕 X 轴旋转（1 DOF），nq=1, nv=1
    joint_id = model.addJoint(
        parent_id, pin.JointModelRX(), joint_placement, joint_name
    )

    # 连杆质心在关节正上方 1.0 m（Z 轴方向）
    body_inertia            = pin.Inertia.FromSphere(body_mass, body_radius)
    body_placement          = joint_placement.copy()
    body_placement.translation[2] = 1.0   # 摆长 = 1.0 m
    model.appendBodyToJoint(joint_id, body_inertia, body_placement)

    # 连杆可视化：球体（质点）+ 细圆柱（连杆杆）
    geom1_name = "ball_" + str(k + 1)
    shape1     = fcl.Sphere(body_radius)
    geom1_obj  = pin.GeometryObject(geom1_name, joint_id, body_placement, shape1)
    geom1_obj.meshColor = np.ones(4)   # 白色球
    geom_model.addGeometryObject(geom1_obj)

    geom2_name = "bar_" + str(k + 1)
    # 细圆柱：半径 = 球半径/4，高度 = 摆长（质心 Z 坐标）
    shape2          = fcl.Cylinder(body_radius / 4.0, body_placement.translation[2])
    shape2_placement = body_placement.copy()
    shape2_placement.translation[2] /= 2.0   # 圆柱中点在关节和质心之间
    geom2_obj = pin.GeometryObject(geom2_name, joint_id, shape2_placement, shape2)
    geom2_obj.meshColor = np.array([0.0, 0.0, 0.0, 1.0])   # 黑色杆
    geom_model.addGeometryObject(geom2_obj)

    # 下一节摆的父关节 = 当前质心位置（链式串联）
    parent_id       = joint_id
    joint_placement = body_placement.copy()


visual_model = geom_model   # 共用同一几何模型用于可视化

# Initialize the viewer.
try:
    viz = Visualizer(model, geom_model, visual_model)
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
dt = 0.01     # 积分步长 10 ms
T  = 5        # 仿真总时长 5 s

N = math.floor(T / dt)

# 设置关节限位（randomConfiguration 需要有限范围）
model.lowerPositionLimit.fill(-math.pi)
model.upperPositionLimit.fill(+math.pi)

if args.with_cart:
    # 小车关节限位为 0（固定小车，只让摆自由运动）
    model.lowerPositionLimit[0] = model.upperPositionLimit[0] = 0.0

data_sim = model.createData()

t = 0.0
q            = pin.randomConfiguration(model)   # 随机初始配置
v            = np.zeros(model.nv)               # 初始速度为零
tau_control  = np.zeros(model.nv)
damping_value = 0.1   # 小阻尼系数（防止能量累积导致数值爆炸）

for k in range(N):
    tic = time.time()
    # 阻尼力矩：τ = -b·v（线性粘滞阻尼，耗散能量使仿真稳定）
    tau_control = -damping_value * v

    # ABA（Articulated Body Algorithm）：O(n) 正向动力学
    # 输入：q（配置）, v（速度）, τ（力矩）
    # 输出：q̈（加速度），内部使用铰接体惯量递推，避免显式求 M⁻¹
    a = pin.aba(model, data_sim, q, v, tau_control)

    # ---- 半隐式欧拉积分 ----
    # v(t+dt) = v(t) + a(t+dt) * dt
    #   → 速度用新加速度更新（"半隐式"的"隐"：用了当前步的加速度）
    v += a * dt
    # q(t+dt) = integrate(model, q(t), v(t+dt) * dt)
    #   → integrate() 在 Lie 群流形上积分（对旋转关节做 Exp 映射）
    #   → 等价于：q = q ⊕ (v·dt)，自动保持四元数归一化等约束
    q = pin.integrate(model, q, v * dt)

    viz.display(q)
    toc = time.time()
    ellapsed = toc - tic

    # 实时仿真：补偿计算耗时，保持显示帧率与仿真步长一致
    dt_sleep = max(0, dt - (ellapsed))
    time.sleep(dt_sleep)
    t += dt
