"""
Pose a Solo-12 robot on a surface defined through a function and displayed through an
hppfcl.HeightField.
"""

import time
from pathlib import Path

import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# HeightField 地形显示：Solo-12 站立在曲面地形上
# 演示内容：
#   - hppfcl.HeightFieldAABB：用高度图（elevation map）定义地形几何体
#   - 将地形作为 GeometryObject 添加到 Meshcat 场景
#   - addGeometryObject：动态向已初始化的 visualizer 添加几何体
#
# HeightField 应用场景：
#   - 地形碰撞检测（机器人在山地/坡道上行走）
#   - 从激光雷达点云生成高度图用于接地约束
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path     = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir       = pinocchio_model_dir
# Solo-12：12 DOF 四足机器人（每腿 3 DOF：髋内/外展、髋屈伸、膝）
urdf_filename   = "solo12.urdf"
urdf_model_path = model_path / "solo_description/robots" / urdf_filename

# 浮动基加载（四足机器人在三维空间自由运动）
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir, pin.JointModelFreeFlyer()
)

# ---- Solo-12 的特定站立配置（来自标定数据）----
# 格式：[浮动基位置(3), 浮动基四元数(4), 腿部关节角(12)] = 19 维
# 这是一个从 SRDF 或参数辨识得到的特定姿态（非零位 neutral）
q_ref = np.array(
    [
        [0.09906518],   # 浮动基 x
        [0.20099078],   # 浮动基 y
        [0.32502457],   # 浮动基 z（站立高度约 0.325 m）
        [0.19414175],   # 浮动基四元数 qx
        [-0.00524735],  # 浮动基四元数 qy
        [-0.97855773],  # 浮动基四元数 qz（接近 180° 偏转）
        [0.06860185],   # 浮动基四元数 qw
        [0.00968163],   # 腿 1 关节 1
        [0.60963582],   # 腿 1 关节 2
        [-1.61206407],  # 腿 1 关节 3
        [-0.02543309],  # 腿 2 关节 1
        [0.66709088],   # 腿 2 关节 2
        [-1.50870083],  # 腿 2 关节 3
        [0.32405118],   # 腿 3 关节 1
        [-1.15305599],  # 腿 3 关节 2
        [1.56867351],   # 腿 3 关节 3
        [-0.39097222],  # 腿 4 关节 1
        [-1.29675892],  # 腿 4 关节 2
        [1.39741073],   # 腿 4 关节 3
    ]
)

vizer = MeshcatVisualizer(model, collision_model, visual_model)
vizer.initViewer(open=True)


# ---- 定义曲面地形函数 ----
# 使用三角函数叠加模拟崎岖地形（类似沙丘/丘陵）
def ground(xy):
    return (
        np.sin(xy[0] * 3) / 5          # 低频 X 方向正弦波（振幅 0.2 m）
        + np.cos(xy[1] ** 2 * 3) / 20  # Y 方向余弦（振幅 0.05 m）
        + np.sin(xy[1] * xy[0] * 5) / 10  # XY 交叉正弦（振幅 0.1 m）
    )


def vizGround(viz, elevation_fn, space, name="ground", color=[1.0, 1.0, 0.6, 0.8]):
    """
    根据高度函数构建 HeightField 地形并添加到 Meshcat 场景。
    elevation_fn：地形高度函数 f(xy) → z，输入 (2, nx, nx) 网格坐标
    space：采样间距（越小越精细，越慢）
    """
    xg = np.arange(-2, 2, space)   # X/Y 轴采样点
    nx = xg.shape[0]
    xy_g   = np.meshgrid(xg, xg)
    xy_g   = np.stack(xy_g)        # shape (2, nx, nx)
    elev_g = np.zeros((nx, nx))
    elev_g[:, :] = elevation_fn(xy_g)   # 计算高度图矩阵

    sx = xg[-1] - xg[0]   # X 方向总长度（米）
    sy = xg[-1] - xg[0]   # Y 方向总长度（米）
    # 翻转行：HeightFieldAABB 期望 Y 轴从下到上（图像坐标系 vs. 世界坐标系）
    elev_g[:, :] = elev_g[::-1, :]
    import hppfcl

    # HeightFieldAABB：用 AABB 树加速的高度场碰撞几何体
    # 参数：sx(米), sy(米), 高度矩阵, 最小高度（z轴下界）
    heightField = hppfcl.HeightFieldAABB(sx, sy, elev_g, np.min(elev_g))
    pl  = pin.SE3.Identity()   # 地形原点在世界坐标系原点
    obj = pin.GeometryObject(name, 0, pl, heightField)
    obj.meshColor[:] = color
    # addGeometryObject：向已初始化的 visualizer 动态添加几何体
    # 等价于：向 visual_model 添加后 rebuildData()
    viz.addGeometryObject(obj)


# Load the robot in the viewer.
vizer.loadViewerModel()

# 构建并显示地形（4m×4m 范围，0.02m 采样间距，半透明淡黄色）
colorrgb = [128, 149, 255, 200]   # 淡蓝紫色（RGBA 0~255）
colorrgb = np.array(colorrgb) / 255.0   # 归一化到 [0, 1]
vizGround(vizer, ground, 0.02, color=colorrgb)

# 显示 Solo-12 站立在曲面地形上的配置
vizer.display(q_ref)
time.sleep(1.0)
