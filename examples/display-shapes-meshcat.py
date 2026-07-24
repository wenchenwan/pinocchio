import sys

import numpy as np
import pinocchio as pin

# ============================================================
# 基本几何体的 Meshcat 可视化
# 与 display-shapes.py 相同功能，使用 MeshcatVisualizer 替代 GepettoVisualizer
#
# MeshcatVisualizer vs GepettoVisualizer：
#   Meshcat：基于浏览器（Three.js），无需安装桌面程序，跨平台，支持远程查看
#   Gepetto：基于 Qt/OSG 桌面程序，需要本地安装 gepetto-viewer-server
#
# 注意：display-shapes-meshcat.py 中 GeometryObject 参数顺序与 display-shapes.py 略有不同：
#   Gepetto 版：GeometryObject(name, joint_id, frame_id, placement, geom)
#   Meshcat 版：GeometryObject(name, joint_id, frame_id, geom, placement)  ← geom/placement 互换
#   这是 API 版本差异，两者功能相同
# ============================================================

try:
    import hppfcl
except ImportError:
    print("This example requires hppfcl")
    sys.exit(0)
from pinocchio.visualize import MeshcatVisualizer

# 无关节的空运动学模型（几何体固定在世界坐标系）
model      = pin.Model()
geom_model = pin.GeometryModel()

# ---- 定义五种基本几何体 ----
geometries = [
    hppfcl.Capsule(0.1, 0.8),    # 胶囊体（常用于碰撞近似）
    hppfcl.Sphere(0.5),           # 球体
    hppfcl.Box(1, 1, 1),          # 长方体
    hppfcl.Cylinder(0.1, 1.0),   # 圆柱体
    hppfcl.Cone(0.5, 1.0),       # 圆锥体
]
for i, geom in enumerate(geometries):
    placement = pin.SE3(np.eye(3), np.array([i, 0, 0]))   # 沿 X 轴等间距排列
    # 注意参数顺序：(name, joint_id, frame_id, geom, placement) ← Meshcat 版本
    geom_obj  = pin.GeometryObject(f"obj{i}", 0, 0, geom, placement)
    color     = np.random.uniform(0, 1, 4)
    color[3]  = 1   # 完全不透明
    geom_obj.meshColor = color
    geom_model.addGeometryObject(geom_obj)

# ---- MeshcatVisualizer：浏览器内 3D 可视化 ----
# 碰撞模型和视觉模型都使用同一个 geom_model
viz = MeshcatVisualizer(model, geom_model, geom_model)

# Initialize the viewer.
try:
    # open=True：自动在默认浏览器中打开 Meshcat 页面（localhost:7000 或类似端口）
    viz.initViewer(open=True)
except ImportError as error:
    print(error)
    sys.exit(0)

try:
    viz.loadViewerModel("shapes")
except AttributeError as error:
    print(error)
    sys.exit(0)

# 显示：模型无关节（nq=0），传空配置向量
viz.display(np.zeros(0))
input("press enter to continue")   # 防止脚本立即退出（等待用户查看 Meshcat）
