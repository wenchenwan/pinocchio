# This examples shows how to load and move a robot in meshcat.
# Note: this feature requires Meshcat to be installed, this can be done using
# pip install --user meshcat

import sys
from pathlib import Path

import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# Meshcat 查看器：Romeo 人形机器人（DAE 网格格式）
# 与 meshcat-viewer.py 的区别：
#   - 加载 Romeo 人形机器人（而非 Solo 四足）
#   - Romeo URDF 的网格使用 .dae（COLLADA）格式（而非 .stl）
#   - Romeo URDF 中没有颜色信息，需要手工指定 color 参数
#   - 演示：多颜色实例（同一 viewer 中同时显示两个颜色不同的机器人）
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path     = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir       = pinocchio_model_dir
# Romeo：Softbank/Aldebaran 人形机器人，36 DOF（含两臂、两腿、头）
urdf_filename    = "romeo_small.urdf"
urdf_model_path  = model_path / "romeo_description/urdf" / urdf_filename

# JointModelFreeFlyer：人形机器人必须使用浮动基（nq=36+7=43, nv=36+6=42）
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir, pin.JointModelFreeFlyer()
)

# Start a new MeshCat server and client.
try:
    viz = MeshcatVisualizer(model, collision_model, visual_model)
    viz.initViewer(open=True)
except ImportError as err:
    print(
        "Error while initializing the viewer. "
        "It seems you should install Python meshcat"
    )
    print(err)
    sys.exit(0)

# Load the robot in the viewer.
# Color is needed here because the Romeo URDF doesn't contain any color, so the default
# color results in an invisible robot (alpha value set to 0).
# color=[R, G, B, A]：RGBA 颜色，用于覆盖 URDF 中缺失的颜色信息
# 黑色 + 完全不透明（alpha=1.0）
viz.loadViewerModel(color=[0.0, 0.0, 0.0, 1.0])

# ---- 手工指定 Romeo 的站立配置（half_sitting 近似）----
# Romeo 配置向量格式：
#   [0:3]   浮动基位置（x, y, z）
#   [3:7]   浮动基姿态（四元数 qx, qy, qz, qw）
#   [7:13]  左腿 6 DOF
#   [13:19] 右腿 6 DOF
#   [19]    躯干 1 DOF
#   [20:27] 左臂 7 DOF
#   [27:31] 头部 4 DOF
#   [31:38] 右臂 7 DOF
q0 = np.array(
    [
        0,
        0,
        0.840252,       # 浮动基 Z 高度（站立时约 0.84 m）
        0,
        0,
        0,
        1,              # 浮动基四元数（w=1 → 无旋转）
        0,
        0,
        -0.3490658,
        0.6981317,
        -0.3490658,
        0,              # 左腿：髋-膝-踝（约 -20°, +40°, -20°）
        0,
        0,
        -0.3490658,
        0.6981317,
        -0.3490658,
        0,              # 右腿（对称）
        0,              # 躯干
        1.5,
        0.6,
        -0.5,
        -1.05,
        -0.4,
        -0.3,
        -0.2,           # 左臂（T 形展开姿态）
        0,
        0,
        0,
        0,              # 头部
        1.5,
        -0.6,
        0.5,
        1.05,
        -0.4,
        -0.3,
        -0.2,           # 右臂（与左臂对称）
    ]
).T
viz.display(q0)

# ---- 多机器人同场景：显示半透明红色机器人 ----
# 创建第二个 visualizer，共享同一个 viewer 实例
# rootNodeName="red_robot"：在 Meshcat 场景树中用不同名称区分两个机器人
red_robot_viz = MeshcatVisualizer(model, collision_model, visual_model)
red_robot_viz.initViewer(viz.viewer)   # 共享 viewer（不新建服务器）
# color=[1, 0, 0, 0.5]：半透明红色（alpha=0.5）
red_robot_viz.loadViewerModel(rootNodeName="red_robot", color=[1.0, 0.0, 0.0, 0.5])
q    = q0.copy()
q[1] = 1.0   # 在 Y 方向偏移 1m（两机器人并排显示）
red_robot_viz.display(q)
