# This examples shows how to perform collision detection between the end-effector of a
# robot and a point cloud depicted as a Height Field
# Note: this feature requires Meshcat to be installed, this can be done using
# pip install --user meshcat

import sys
import time
from pathlib import Path

import hppfcl as fcl
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# 点云碰撞检测：将点云转换为 BVH 网格或高度场，与机械臂做碰撞检测
# 应用场景：
#   - 机械臂在激光雷达/RGB-D 点云环境中的碰撞检测
#   - 末端执行器与工作台面（HeightField）的接触检测
#   - 在线路径规划的碰撞查询
#
# 两种点云表示方式：
#   BVHModelOBBRSS：BVH（包围体层次结构）+ OBB（有向包围盒）+ RSS
#     优点：精确表示任意点云分布；缺点：构建和更新代价较高
#   HeightFieldOBBRSS：高度场（Elevation Map）+ OBB/RSS 加速
#     优点：结构简单，适合平坦地形；缺点：只能表示单值高度函数
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path     = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir       = pinocchio_model_dir
# Panda：7-DOF 机械臂，末端执行器 panda_hand 用于碰撞检测
urdf_filename   = "panda.urdf"
urdf_model_path = model_path / "panda_description/urdf" / urdf_filename

model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir
)

# ---- 生成随机点云（模拟激光雷达数据）----
num_points = 5000
points     = np.random.rand(3, num_points)   # 在 [0,1]³ 内均匀随机采样
# 点云在世界坐标系中的位置（偏移到机械臂前方偏下位置）
point_cloud_placement = pin.SE3.Identity()
point_cloud_placement.translation = np.array([0.2, 0.2, -0.5])

X = points[0, :]
Y = points[1, :]
Z = points[2, :]

# ---- 将点云转换为高度图（Elevation Map）----
# 高度图：将三维点投影到 XY 平面，记录每个格子内的最大 Z 值
nx, ny   = 20, 20
x_grid   = np.linspace(0.0, 1.0, nx)
x_half_pad = 0.5 * (x_grid[1] - x_grid[0])   # 半格偏移（防止边界效应）
x_bins   = np.digitize(X, x_grid + x_half_pad)   # 每个点的 X 格子索引
x_dim    = x_grid[-1] - x_grid[0]               # X 方向总范围（米）

y_grid   = np.linspace(0.0, 1.0, ny)
y_half_pad = 0.5 * (y_grid[1] - y_grid[0])
y_bins   = np.digitize(Y, y_grid + y_half_pad)
y_dim    = y_grid[-1] - y_grid[0]

# 线性化格子索引：point_bins[i] 是第 i 个点所在的格子编号
point_bins = y_bins * nx + x_bins
heights    = np.zeros((ny, nx))
# 对每个格子取该格子内所有点的最大 Z（堆叠高度）
np.maximum.at(heights.ravel(), point_bins, Z)

# ---- 构建 BVH 点云几何体（精确表示）----
# BVHModelOBBRSS：层次式包围盒树（OBB + RSS 双精度），支持精确碰撞查询
point_cloud = fcl.BVHModelOBBRSS()
point_cloud.beginModel(0, num_points)    # 0 个三角面片，num_points 个顶点
point_cloud.addVertices(points.T)        # 添加顶点（shape: (n, 3)）
# 注意：纯点云（无面片）的碰撞检测精度有限，通常配合凸包使用

# ---- 构建 HeightField 几何体（高效平面近似）----
# HeightFieldOBBRSS(sx, sy, heights, min_height)
# sx/sy：X/Y 方向的总尺寸（米）
# heights：(ny, nx) 高度矩阵
# min_height：高度场下界（通常为高度矩阵的最小值）
height_field = fcl.HeightFieldOBBRSS(x_dim, y_dim, heights, min(Z))
# 高度场中心偏移（使高度场原点对齐到 XY 网格中心）
height_field_placement = point_cloud_placement * pin.SE3(
    np.eye(3), 0.5 * np.array([x_grid[0] + x_grid[-1], y_grid[0] + y_grid[-1], 0.0])
)

# ---- 将点云和高度场添加到几何模型 ----
go_point_cloud = pin.GeometryObject("point_cloud", 0, point_cloud_placement, point_cloud)
go_point_cloud.meshColor = np.ones(4)   # 白色
collision_model.addGeometryObject(go_point_cloud)
visual_model.addGeometryObject(go_point_cloud)

go_height_field = pin.GeometryObject("height_field", 0, height_field_placement, height_field)
go_height_field.meshColor = np.ones(4)
# addGeometryObject 返回新添加的几何体 ID（用于后续指定碰撞对）
height_field_collision_id = collision_model.addGeometryObject(go_height_field)
visual_model.addGeometryObject(go_height_field)

# ---- 设置末端执行器与高度场之间的碰撞对 ----
# 获取 panda_hand 的碰撞几何体 ID
panda_hand_collision_id = collision_model.getGeometryId("panda_hand_0")
go_panda_hand = collision_model.geometryObjects[panda_hand_collision_id]
# buildConvexRepresentation(False)：False=不保留原始网格，节省内存
# 替换为凸包：FCL 对凸多面体的碰撞检测效率远高于 BVH 网格
go_panda_hand.geometry.buildConvexRepresentation(False)
go_panda_hand.geometry = go_panda_hand.geometry.convex   # 用凸包替换原始网格

# 添加指定碰撞对：高度场 vs panda_hand 凸包
# 不使用 addAllCollisionPairs 是因为只关心末端与地形的碰撞
collision_pair = pin.CollisionPair(height_field_collision_id, panda_hand_collision_id)
collision_model.addCollisionPair(collision_pair)

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

viz.loadViewerModel()

q0 = pin.neutral(model)
viz.display(q0)

# ---- 随机采样直到找到碰撞配置 ----
is_collision    = False
data            = model.createData()
collision_data  = collision_model.createData()
while not is_collision:
    q = pin.randomConfiguration(model)
    # computeCollisions(... , True)：True = 发现第一个碰撞立即停止（提前终止）
    # 返回值：bool，是否发生了至少一个碰撞
    is_collision = pin.computeCollisions(
        model, data, collision_model, collision_data, q, True
    )

print("Found a configuration in collision:", q)
viz.display(q)
time.sleep(1.0)
