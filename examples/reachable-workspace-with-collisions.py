import itertools
import sys
import time
from pathlib import Path

import meshcat.geometry as g
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# 带碰撞约束的可达工作空间（Collision-Aware Reachable Workspace）
# 与 reachable-workspace.py 的区别：
#   - 在采样过程中排除与障碍物碰撞的配置
#   - 使用 reachableWorkspaceWithCollisions / reachableWorkspaceWithCollisionsHull
#   - 场景中有三个胶囊体障碍物（手工放置在机械臂工作空间内）
# 应用场景：
#   - 考虑障碍物后的安全可达域分析
#   - 运动规划中的可达性地图预计算
# ============================================================


def XYZRPYtoSE3(xyzrpy):
    """将 [x, y, z, roll, pitch, yaw] 转换为 SE(3) 变换矩阵"""
    rotate = pin.utils.rotate
    # R = Rx(roll) · Ry(pitch) · Rz(yaw)（固定轴外旋顺序）
    R = rotate("x", xyzrpy[3]) @ rotate("y", xyzrpy[4]) @ rotate("z", xyzrpy[5])
    p = np.array(xyzrpy[:3])
    return pin.SE3(R, p)


pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir   = pinocchio_model_dir
# Panda：7-DOF 工业机械臂，SRDF 中定义了自碰撞排除对
urdf_path  = model_path / "panda_description/urdf/panda.urdf"
srdf_path  = model_path / "panda_description/srdf/panda.srdf"

robot, collision_model, visual_model = pin.buildModelsFromUrdf(urdf_path, mesh_dir)
data = robot.createData()

# ---- 定义障碍物（三个胶囊体）----
# XYZ-RPY 参数：[x, y, z, roll, pitch, yaw]
# 这些位置手工选取，位于 Panda 机械臂工作空间内（0.3~0.8m 半径范围）
oMobs = [
    [0.40,  0.0, 0.30, np.pi / 2, 0, 0],   # 正前方障碍物（竖置胶囊）
    [-0.08, -0.0, 0.75, np.pi / 2, 0, 0],  # 后上方障碍物
    [0.23,  -0.0, 0.04, np.pi / 2, 0, 0],  # 近前方低位障碍物
]

rad, length = 0.1, 0.4   # 胶囊体半径 0.1m，长度 0.4m
for i, xyzrpy in enumerate(oMobs):
    # CreateCapsule：创建胶囊几何对象（不绑定到特定关节）
    obs = pin.GeometryObject.CreateCapsule(rad, length)
    obs.meshColor  = np.array([1.0, 0.2, 0.2, 1.0])   # 红色
    obs.name       = f"obs{i}"
    obs.parentJoint = 0   # 固定在世界坐标系（不随任何关节运动）
    obs.placement  = XYZRPYtoSE3(xyzrpy)   # 障碍物在世界系中的位姿
    collision_model.addGeometryObject(obs)  # 加入碰撞模型（用于碰撞检测）
    visual_model.addGeometryObject(obs)     # 加入视觉模型（用于显示）

# ---- 碰撞对设置 ----
# addAllCollisionPairs：先添加所有几何体对之间的碰撞对
collision_model.addAllCollisionPairs()
# removeCollisionPairs：从 SRDF 移除不需要的自碰撞对（相邻连杆等）
pin.removeCollisionPairs(robot, collision_model, srdf_path)

# 手工添加"机器人连杆 vs 障碍物"的碰撞对
# 不添加"障碍物 vs 障碍物"（障碍物之间不需要检测互相碰撞）
nobs    = len(oMobs)
nbodies = collision_model.ngeoms - nobs   # 机器人自身几何体数量
robotBodies = range(nbodies)              # 机器人几何体 ID 范围
envBodies   = range(nbodies, nbodies + nobs)  # 障碍物 ID 范围
for a, b in itertools.product(robotBodies, envBodies):
    # 笛卡尔积：每个机器人几何体与每个障碍物配对
    collision_model.addCollisionPair(pin.CollisionPair(a, b))

# ---- 重新创建数据（修改 model 后必须重建 Data）----
# 修改 collision_model（添加几何体）后，原有 collision_data 失效
collision_data = pin.GeometryData(collision_model)
visual_data    = pin.GeometryData(visual_model)

# Start a new MeshCat server and client.
viz = MeshcatVisualizer(robot, collision_model, visual_model)
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

# ---- 工作空间计算参数 ----
q0        = (robot.upperPositionLimit.T + robot.lowerPositionLimit.T) / 2
horizon   = 0.2   # 时间范围 0.2 秒
frame     = robot.getFrameId(robot.frames[-1].name)
n_samples = 5
facet_dims = 2

# 本例默认使用非凸 CGAL 方法（convex=False，需要安装 CGAL）
convex = False

if convex:
    # reachableWorkspaceWithCollisionsHull：含碰撞过滤的凸包工作空间
    # 内部在 Monte Carlo 采样时，排除与 collision_model 中障碍物碰撞的配置
    verts, faces = pin.reachableWorkspaceWithCollisionsHull(
        robot, collision_model, q0, horizon, frame, n_samples, facet_dims
    )
    verts = verts.T

else:
    try:
        from CGAL.CGAL_Alpha_wrap_3 import *  # noqa: F403
        from CGAL.CGAL_Kernel import *  # noqa: F403
        from CGAL.CGAL_Mesh_3 import *  # noqa: F403
        from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
    except ModuleNotFoundError:
        print("To compute non convex Polytope CGAL library needs to be installed.")
        sys.exit(0)

    def vertex_to_tuple(v):
        return v.x(), v.y(), v.z()

    def halfedge_to_triangle(he):
        p1 = he
        p2 = p1.next()
        p3 = p2.next()
        return [
            vertex_to_tuple(p1.vertex().point()),
            vertex_to_tuple(p2.vertex().point()),
            vertex_to_tuple(p3.vertex().point()),
        ]

    def alpha_shape_with_cgal(coords, alpha=None):
        """用 CGAL Alpha-wrap 计算点云的非凸包络曲面（见 reachable-workspace.py 的详细注释）"""
        if alpha is None:
            bbox_diag   = np.linalg.norm(np.max(coords, 0) - np.min(coords, 0))
            alpha_value = bbox_diag / 5
        else:
            alpha_value = np.mean(alpha)
        points = [Point_3(pt[0], pt[1], pt[2]) for pt in coords]  # noqa: F405
        Q  = Polyhedron_3()
        _a = alpha_wrap_3(points, alpha_value, 0.01, Q)  # noqa: F405
        alpha_shape_vertices = np.array(
            [vertex_to_tuple(vertex.point()) for vertex in Q.vertices()]
        )
        alpha_shape_faces = np.array(
            [np.array(halfedge_to_triangle(face.halfedge())) for face in Q.facets()]
        )
        return alpha_shape_vertices, alpha_shape_faces

    # reachableWorkspaceWithCollisions：返回碰撞过滤后的可达点集
    verts = pin.reachableWorkspaceWithCollisions(
        robot, collision_model, q0, horizon, frame, n_samples, facet_dims
    )
    verts = verts.T

    alpha       = 0.1   # 比 reachable-workspace.py 的 0.2 更紧（障碍物使空间更复杂）
    verts, faces = alpha_shape_with_cgal(verts, alpha)
    verts       = faces.reshape(-1, 3)
    faces       = np.arange(len(verts)).reshape(-1, 3)

print("------------------- Display Vertex")

# ---- Meshcat 三角网格显示 ----
poly = g.TriangularMeshGeometry(vertices=verts, faces=faces)
viz.viewer["poly"].set_object(
    poly, g.MeshBasicMaterial(color=0x000000, wireframe=True, linewidth=12, opacity=0.2)
)

while True:
    viz.display(q0)
    viz.viewer["poly"].set_object(
        poly,
        g.MeshBasicMaterial(color=0x000000, wireframe=True, linewidth=2, opacity=0.2),
    )
    time.sleep(1e-2)
