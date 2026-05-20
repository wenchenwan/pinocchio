# NOTE: this example needs RViz to be installed
# usage: start ROS master (roscore) and then run this test


from pathlib import Path

import pinocchio as pin
from pinocchio.visualize import RVizVisualizer

# ============================================================
# RVizVisualizer：ROS RViz 3D 可视化
# 前提：需要安装 ROS（1 或 2），并启动 roscore + RViz
# 与 Meshcat/Gepetto 的区别：
#   RViz：ROS 生态系统中的标准可视化工具，支持实时话题订阅
#         适合与 ROS 控制器/传感器集成的完整系统
#   Meshcat：独立运行，无需 ROS，更轻量
#
# RVizVisualizer 内部：
#   通过 /robot_description 话题发布 URDF 模型
#   通过 /joint_states 话题发布关节角（JointState 消息）
#   通过 TF 广播浮动基的位姿（tf_broadcaster）
#
# 演示：多机器人同场景显示（通过 rootNodeName 区分 TF namespace）
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path     = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir       = pinocchio_model_dir
# Talos（降阶版）：PAL Robotics 人形机器人，常见 ROS 机器人
urdf_filename   = "talos_reduced.urdf"
urdf_model_path = model_path / "talos_data/robots" / urdf_filename

# 浮动基加载（人形机器人）
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir, pin.JointModelFreeFlyer()
)
viz = RVizVisualizer(model, collision_model, visual_model)

# Initialize the viewer.
# initViewer：初始化 ROS 节点，发布 /robot_description 话题
viz.initViewer()
# loadViewerModel：在 RViz 中加载 robot_state_publisher 并订阅关节状态
viz.loadViewerModel("pinocchio")

# Display a robot configuration.
# 内部：发布 JointState 消息到 /joint_states 话题
q0 = pin.neutral(model)
viz.display(q0)

# ---- 第二个机器人实例（在 RViz 中使用独立的 TF namespace）----
viz2 = RVizVisualizer(model, collision_model, visual_model)
viz2.initViewer(viz.viewer)   # 共享 ROS 节点句柄
# rootNodeName="pinocchio2"：独立的 robot_description 命名空间
viz2.loadViewerModel(rootNodeName="pinocchio2")
q    = q0.copy()
q[1] = 1.0   # Y 方向偏移 1m
viz2.display(q)

input("Press enter to exit...")

# clean()：清理 RViz 中的显示对象（RViz 标记清除）
# 对于 RVizVisualizer，退出时需要显式清理，防止残留标记
viz.clean()
