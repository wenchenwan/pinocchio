# This examples shows how to load several robots in panda3d_viewer.
# Note: this feature requires panda3d_viewer to be installed, this can be done using
# pip install panda3d_viewer
# ruff: noqa: E402

import sys
from pathlib import Path

# ============================================================
# Panda3D 查看器：多机器人同场景显示
# panda3d_viewer：基于 Panda3D 游戏引擎的 3D 可视化（需要 pip install panda3d_viewer）
# 与其他查看器的区别：
#   Panda3D：游戏引擎渲染，支持高质量光影效果，可交互（鼠标旋转/缩放）
#             适合展示和演示，支持多机器人同场景
#   Meshcat：浏览器内置，远程友好，更适合 Jupyter 集成
#   Gepetto：桌面 Qt 程序，最成熟但依赖重
#
# 演示：同时显示 7 个不同机器人（Talos、Romeo、ICub、Tiago、Solo8、HyQ、Hector）
# 每个机器人在 Y 方向错开排列，便于对比观察
# ============================================================

# Add path to the example-robot-data package from git submodule.
# If you have it properly installed, there is no need for this sys.path thing.
path = Path(__file__).parent.parent / "models" / "example-robot-data" / "python"
sys.path.append(str(path))
from example_robot_data.robots_loader import (
    HectorLoader,   # Hector UAV（无人机）
    HyQLoader,      # HyQ 四足机器人（IIT）
    ICubLoader,     # iCub 人形机器人（IIT）
    RomeoLoader,    # Romeo 人形机器人（SoftBank）
    Solo8Loader,    # Solo8 四足机器人（LAAS）
    TalosLoader,    # Talos 人形机器人（PAL Robotics）
    TiagoLoader,    # TIAGo 移动操作机器人（PAL Robotics）
)
from panda3d_viewer import Viewer
from pinocchio.visualize.panda3d_visualizer import Panda3dVisualizer

# 打开 Panda3D GUI 窗口（独立窗口，非浏览器）
viewer = Viewer(window_title="python-pinocchio")

# ---- 依次加载并显示各机器人 ----
# 每个 Loader 封装了对应机器人的 URDF/SRDF 路径和默认配置
loaders = (
    TalosLoader,   # 人形，36 DOF
    RomeoLoader,   # 人形，36 DOF（需要手工设置站立高度）
    ICubLoader,    # 人形，32 DOF
    TiagoLoader,   # 移动臂，12+ DOF
    Solo8Loader,   # 四足，8 DOF
    HyQLoader,     # 四足，12 DOF
    HectorLoader,  # 无人机，6 DOF
)

for i, loader in enumerate(loaders):
    # loader().robot：调用 Loader 构造函数加载机器人（返回 RobotWrapper）
    robot = loader().robot
    # 将 Panda3D 查看器绑定到 RobotWrapper
    robot.setVisualizer(Panda3dVisualizer())
    robot.initViewer(viewer=viewer)   # 共享同一个 Panda3D viewer 窗口
    # group_name：Panda3D 场景节点名称（对应机器人的 model.name）
    robot.loadViewerModel(group_name=robot.model.name)

    # 在 Y 方向按索引 i 错开各机器人（间距 3m）
    q    = robot.q0[:]
    q[1] = 3 - i   # Y 轴位置：+3m 到 -3m 依次排列
    if loader is RomeoLoader:
        q[2] = 0.87   # Romeo 需要手工设置浮动基 Z 高度（URDF 零位在地面以下）

    robot.display(q)

# 等待窗口关闭（阻塞，直到用户关闭 Panda3D 窗口）
viewer.join()
