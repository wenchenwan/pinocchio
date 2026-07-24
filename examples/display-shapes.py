import sys

import numpy as np
import pinocchio as pin

# ============================================================
# 基本几何体的可视化（GepettoVisualizer）
# 演示：无关节模型 + 纯几何体（Capsule/Sphere/Box/Cylinder/Cone）
# 适用场景：
#   - 调试几何体的形状和位置（不涉及机器人动力学）
#   - 快速搭建环境障碍物或物体模型
#   - 验证几何体参数（半径、高度、颜色等）
# ============================================================

try:
    import hppfcl
except ImportError:
    print("This example requires hppfcl")
    sys.exit(0)
from pinocchio.visualize import GepettoVisualizer

# ---- 无关节的空运动学模型 ----
# 几何体不依附于任何运动关节，直接挂在世界坐标系（joint_id=0）
model = pin.Model()

# ---- 定义几何体列表 ----
geom_model = pin.GeometryModel()
geometries = [
    hppfcl.Capsule(0.1, 0.8),    # 胶囊体：半径 0.1 m，高度 0.8 m（用于杆/肢体近似）
    hppfcl.Sphere(0.5),           # 球体：半径 0.5 m
    hppfcl.Box(1, 1, 1),          # 长方体：1m × 1m × 1m
    hppfcl.Cylinder(0.1, 1.0),   # 圆柱体：半径 0.1 m，高度 1.0 m（轴沿 Z 轴）
    hppfcl.Cone(0.5, 1.0),       # 圆锥体：底半径 0.5 m，高度 1.0 m
]
for i, geom in enumerate(geometries):
    # 将每个几何体沿 X 轴等间距排列（防止重叠）
    placement = pin.SE3(np.eye(3), np.array([i, 0, 0]))
    # GeometryObject(name, parent_joint_id, parent_frame_id, placement, geometry)
    # parent_joint_id=0, parent_frame_id=0：固定在世界坐标系（不随任何关节运动）
    geom_obj  = pin.GeometryObject(f"obj{i}", 0, 0, placement, geom)
    color     = np.random.uniform(0, 1, 4)   # 随机 RGBA 颜色
    color[3]  = 1                             # 不透明度 = 1（完全不透明）
    geom_obj.meshColor = color
    geom_model.addGeometryObject(geom_obj)

# ---- GepettoVisualizer：基于 gepetto-viewer 的本地 3D 可视化 ----
# 注意：需要先启动 gepetto-viewer-server（桌面应用程序）
# 两个几何模型参数（碰撞/视觉）都用同一个 geom_model（双重用途）
viz = GepettoVisualizer(
    model=model,
    collision_model=geom_model,
    visual_model=geom_model,
)

# Initialize the viewer.
try:
    viz.initViewer()
except ImportError as error:
    print(
        "Error while initializing the viewer. "
        "It seems you should install gepetto-viewer"
    )
    print(error)
    sys.exit(0)

try:
    viz.loadViewerModel("shapes")
except AttributeError as error:
    print(
        "Error while loading the viewer model. "
        "It seems you should start gepetto-viewer"
    )
    print(error)
    sys.exit(0)

# 显示：模型无关节（nq=0），传空配置向量
viz.display(np.zeros(0))
