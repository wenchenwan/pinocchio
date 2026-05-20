# This examples shows how to load and move a robot in panda3d_viewer.
# Note: this feature requires panda3d_viewer to be installed, this can be done using
# pip install panda3d_viewer
# ruff: noqa: E402


import sys
from pathlib import Path

import numpy as np

# Add path to the example-robot-data package from git submodule.
# If you have a proper install version, there is no need for this sys.path thing
path = Path(__file__).parent.parent / "models" / "example-robot-data" / "python"
sys.path.append(str(path))
from example_robot_data.robots_loader import TalosLoader
from panda3d_viewer import ViewerClosedError
from pinocchio.visualize.panda3d_visualizer import Panda3dVisualizer

# ============================================================
# Panda3D 查看器：轨迹播放（循环往复）
# 演示：
#   - Talos 人形机器人的参数化周期轨迹
#   - robot.play(traj, dt)：按轨迹逐帧播放
#   - np.flip + 循环：使轨迹来回往复（乒乓模式）
# ============================================================

# talos is a RobotWrapper object
talos = TalosLoader().robot
# Panda3D 可视化器绑定
talos.setVisualizer(Panda3dVisualizer())
talos.initViewer()
talos.loadViewerModel(group_name="talos", color=(1, 1, 1, 1))   # 白色机器人


def play_sample_trajectory():
    """
    生成并循环播放 Talos 的行走样例轨迹。
    轨迹设计：
      - 基于 talos.q0（半蹲站立配置）
      - 选取关键关节（腿部关节）施加余弦/线性函数，模拟步态
      - np.flip 翻转轨迹，实现来回往复播放（乒乓）
    """
    update_rate = 60      # 显示帧率 60 FPS
    cycle_time  = 3       # 单次轨迹时长 3 秒
    # 初始化轨迹矩阵：所有列都是 q0（静止起始）
    traj = np.repeat(talos.q0.reshape((-1, 1)), cycle_time * update_rate, axis=1)
    # beta：从 0 到 1 的线性参数（用于生成关节运动轨迹）
    beta = np.linspace(0, 1, traj.shape[1])

    # ---- 选取关节 ID 并指定轨迹 ----
    # Talos 关节索引（对应 talos_reduced.urdf 中的关节顺序）：
    #   2 → 浮动基 Z 高度（上下起伏）
    #   9,10,11 → 左腿髋/膝/踝
    #   22 → 躯干旋转（或其他）
    #   15,16,17 → 右腿髋/膝/踝
    #   30 → 另一躯干/腰关节
    # 余弦函数模拟自然的关节往复运动；线性函数模拟单调的屈伸
    traj[[2, 9, 10, 11, 22, 15, 16, 17, 30]] = (
        0.39 + 0.685 * np.cos(beta),   # Z 高度：余弦起伏（模拟重心上下移动）
        -beta,                          # 左髋：线性屈曲
        2.0 * beta,                     # 左膝：线性伸展
        -beta,                          # 左踝：线性跖屈
        0.1 + beta * 1.56,              # 躯干旋转
        -beta,                          # 右髋：线性屈曲
        2.0 * beta,                     # 右膝：线性伸展
        -beta,                          # 右踝：线性跖屈
        -0.1 - beta * 1.56,             # 另一躯干关节（与上对称）
    )

    while True:
        # robot.play(traj, dt)：
        #   traj：(N_frames, nq) 轨迹矩阵
        #   dt = 1/update_rate：每帧时长
        #   内部逐帧调用 display() 并 sleep(dt) 实现实时播放
        talos.play(traj.T, 1.0 / update_rate)
        # 翻转轨迹：来回往复播放（乒乓模式），避免轨迹不连续的跳跃
        traj = np.flip(traj, 1)


try:
    play_sample_trajectory()
except ViewerClosedError:
    # 用户关闭 Panda3D 窗口时抛出此异常，正常退出
    pass
