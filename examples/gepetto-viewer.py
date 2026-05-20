# NOTE: this example needs gepetto-gui to be installed
# usage: launch gepetto-gui and then run this test

import sys
from pathlib import Path

import pinocchio as pin
from pinocchio.visualize import GepettoVisualizer

# ============================================================
# GepettoVisualizer：本地桌面 3D 可视化（基于 gepetto-viewer）
# 前提：需要先安装并启动 gepetto-viewer-server（桌面 GUI 程序）
# 与 Meshcat 的区别：
#   Gepetto：本地桌面程序（Qt + OpenSceneGraph），渲染质量高，支持光影
#   Meshcat：浏览器（Three.js），轻量，支持远程访问，无需安装桌面程序
#
# 演示：多机器人同场景显示（通过 rootNodeName 区分场景树节点）
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path     = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir       = pinocchio_model_dir
# Talos（降阶版）：PAL Robotics 人形机器人
urdf_filename   = "talos_reduced.urdf"
urdf_model_path = model_path / "talos_data/robots" / urdf_filename

# 浮动基加载（人形机器人必须）
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir, pin.JointModelFreeFlyer()
)
viz = GepettoVisualizer(model, collision_model, visual_model)

# Initialize the viewer.
try:
    viz.initViewer()   # 连接到已运行的 gepetto-viewer-server
except ImportError as err:
    print(
        "Error while initializing the viewer. "
        "It seems you should install gepetto-viewer"
    )
    print(err)
    sys.exit(0)

try:
    # loadViewerModel：在 Gepetto 场景树中以 "pinocchio" 为根节点加载几何体
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

# ---- 同一场景中显示第二个机器人实例 ----
# 共享同一个 viz.viewer 连接（不创建新服务器）
viz2 = GepettoVisualizer(model, collision_model, visual_model)
viz2.initViewer(viz.viewer)
# rootNodeName="pinocchio2"：第二个机器人在场景树中的节点名称（与第一个区分）
viz2.loadViewerModel(rootNodeName="pinocchio2")
q    = q0.copy()
q[1] = 1.0   # 在 Y 方向偏移 1m（两机器人并排）
viz2.display(q)
