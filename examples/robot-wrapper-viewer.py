#
# In this short script, we show how to use RobotWrapper
# integrating different kinds of viewers
#

from pathlib import Path
from sys import argv

import pinocchio as pin
from pinocchio.robot_wrapper import RobotWrapper
from pinocchio.visualize import GepettoVisualizer, MeshcatVisualizer

# ============================================================
# RobotWrapper：高层机器人接口
# RobotWrapper 封装了 model + data + collision_model + visual_model，
# 提供比直接操作 Model/Data 更简洁的接口：
#   robot.q0         → SRDF 中的默认配置（通常是 half_sitting）
#   robot.com(q)     → 质心位置（等价于 pin.centerOfMass(model, data, q)）
#   robot.display(q) → 直接更新可视化
#   robot.buildReducedRobot() → 锁定部分关节得到降阶模型
#
# 应用场景：
#   - 快速原型开发（无需手动管理 model/data 对）
#   - 教学/调试（接口更直观）
#   - 已有完整机器人描述（URDF + SRDF）时推荐使用
# ============================================================

# If you want to visualize the robot in this example,
# you can choose which visualizer to employ
# by specifying an option from the command line:
# GepettoVisualizer: -g
# MeshcatVisualizer: -m
VISUALIZER = None
if len(argv) > 1:
    opt = argv[1]
    if opt == "-g":
        VISUALIZER = GepettoVisualizer    # 本地桌面 3D 查看器
    elif opt == "-m":
        VISUALIZER = MeshcatVisualizer    # 浏览器 3D 查看器
    else:
        raise ValueError("Unrecognized option: " + opt)

# Load the URDF model with RobotWrapper
pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path     = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir       = pinocchio_model_dir
# Talos（降阶版）：PAL Robotics 人形机器人，去掉了手部精细关节
urdf_filename   = "talos_reduced.urdf"
urdf_model_path = model_path / "talos_data/robots" / urdf_filename

# RobotWrapper.BuildFromURDF：一站式加载（解析 URDF + 创建 model/data/geom_models）
# pin.JointModelFreeFlyer()：添加浮动基（人形机器人必须）
robot = RobotWrapper.BuildFromURDF(urdf_model_path, mesh_dir, pin.JointModelFreeFlyer())

# ---- 通过 RobotWrapper 访问底层对象 ----
model = robot.model   # pin.Model
data  = robot.data    # pin.Data

# ---- 质心计算（两种等价写法）----
q0 = robot.q0   # SRDF 中 "half_sitting" 等默认配置（若没有 SRDF 则为 neutral）

# 方式 1：RobotWrapper 高层接口（内部自动调用 forwardKinematics + centerOfMass）
com = robot.com(q0)

# 方式 2：直接调用 Pinocchio 函数（等价，更底层）
# pin.centerOfMass 返回 CoM 在世界坐标系中的位置向量
com2 = pin.centerOfMass(model, data, q0)

# ---- 可视化（可选）----
if VISUALIZER:
    # setVisualizer：将查看器实例绑定到 RobotWrapper
    robot.setVisualizer(VISUALIZER())
    robot.initViewer()                         # 启动查看器服务
    robot.loadViewerModel("pinocchio")         # 加载几何模型到场景
    robot.display(q0)                          # 显示初始配置
