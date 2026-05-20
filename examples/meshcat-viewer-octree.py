# This examples shows how to load and move a robot in meshcat.
# Note: this feature requires Meshcat to be installed, this can be done using
# pip install --user meshcat

import sys
import time

import hppfcl as fcl
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# OctoMap 八叉树碰撞几何体的可视化
# 应用场景：
#   - 将激光雷达/RGB-D 点云转换为体素化碰撞地图（OctoMap）
#   - 与机器人模型做碰撞检测，用于路径规划避障
#   - 在 Meshcat 中可视化 3D 点云的体素化表示
#
# 前提：hpp-fcl 需编译时开启 WITH_OCTOMAP 支持
#   通常需要安装 liboctomap-dev 并重新编译 hpp-fcl
# ============================================================

# ---- 检查 hpp-fcl 版本并确认 OctoMap 支持 ----
# hpp-fcl 3.0.0+ 通过 fcl.WITH_OCTOMAP 标志标识 OctoMap 支持
if tuple(map(int, fcl.__version__.split("."))) >= (3, 0, 0):
    with_octomap = fcl.WITH_OCTOMAP
else:
    with_octomap = False
if not with_octomap:
    print(
        "This example is skiped as HPP-FCL has not been compiled with octomap support."
    )

# 无关节运动学模型（点云/八叉树不属于任何机器人连杆）
model           = pin.Model()
collision_model = pin.GeometryModel()

# ---- 生成随机点云并转换为八叉树 ----
# np.random.rand(1000, 3)：在 [0,1]³ 立方体中随机生成 1000 个点（模拟激光雷达点云）
# 0.01：体素分辨率 0.01 m（= 1 cm，分辨率越小占用内存越多）
# makeOctree 内部：将点云插入八叉树（OcTree），合并相邻体素
octree        = fcl.makeOctree(np.random.rand(1000, 3), 0.01)
# 将八叉树包装为 GeometryObject（固定在世界坐标系原点）
octree_object = pin.GeometryObject("octree", 0, pin.SE3.Identity(), octree)
octree_object.meshColor[0] = 1.0   # 红色（R=1, G=0, B=0, A=0 → 需要注意 Alpha 可能为 0）
collision_model.addGeometryObject(octree_object)

# 视觉和碰撞模型共用（八叉树同时用于显示和碰撞检测）
visual_model = collision_model
viz = MeshcatVisualizer(model, collision_model, visual_model)

# Start a new MeshCat server and client.
try:
    viz.initViewer(open=True)
except ImportError as err:
    print(
        "Error while initializing the viewer. "
        "It seems you should install Python meshcat"
    )
    print(err)
    sys.exit(0)

viz.loadViewerModel()
time.sleep(1.0)   # 等待浏览器加载完毕（Meshcat 异步渲染）
