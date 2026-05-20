import sys
import time
from pathlib import Path

import meshcat.geometry as g
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# 可达工作空间（Reachable Workspace）计算与可视化
# 应用场景：
#   - 末端执行器在给定时间 horizon 内可达的位置集合
#   - 用于任务规划（判断目标点是否可达）
#   - 运动能力分析（在关节限位约束下）
#
# 两种计算方法：
#   1. 凸包（Convex Hull）：reachableWorkspaceHull → 快速但保守
#      内部使用 Monte Carlo 采样 + QuickHull 算法
#   2. 非凸（Alpha Shape）：reachableWorkspace + CGAL Alpha-wrap
#      更精确但需要安装 CGAL Python 绑定
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path     = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir       = pinocchio_model_dir
# Panda（7-DOF 工业机械臂）：固定基，适合工作空间分析
urdf_filename   = "panda.urdf"
urdf_model_path = model_path / "panda_description/urdf" / urdf_filename

# buildModelsFromUrdf 返回 (model, collision_model, visual_model)
# 这里 robot 实际上是 Model（变量名误导性地取名 robot）
robot, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir
)
data = robot.createData()

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
# 起始配置：关节限位的中点（关节角范围中间值，通常接近奇异点少的中立姿态）
q0 = (robot.upperPositionLimit.T + robot.lowerPositionLimit.T) / 2

# horizon：时间范围（秒），机器人从 q0 出发，在 horizon 秒内可达的末端位置
# 等价于：速度为 1 rad/s 时各关节可转动的角度 ≈ horizon·v_max
horizon = 0.2   # 0.2 秒（小值 → 局部工作空间；大值 → 全局可达域）

# frame：末端 Frame 的 ID（这里取模型的最后一个 Frame = 末端执行器）
frame = robot.getFrameId(robot.frames[-1].name)

# n_samples：Monte Carlo 采样次数（越多越精确，越慢）
n_samples  = 5
# facet_dims：工作空间的投影维度（2 = 在某平面内，3 = 完整三维）
facet_dims = 2

# ---- 选择计算方法 ----
convex = True   # True=凸包（快速），False=非凸（精确但需要 CGAL）

if convex:
    # reachableWorkspaceHull：返回凸包的顶点和三角面片（mesh）
    # verts: (3, n_verts) 顶点坐标矩阵
    # faces: (n_faces, 3) 三角面片索引矩阵
    verts, faces = pin.reachableWorkspaceHull(
        robot, q0, horizon, frame, n_samples, facet_dims
    )
    verts = verts.T   # 转置为 (n_verts, 3)（Meshcat 期望的格式）

else:
    # 非凸工作空间：需要 CGAL Alpha-wrap 库
    # reachableWorkspace 返回采样点集，Alpha-wrap 计算其表面包络
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
        """
        用 CGAL Alpha-wrap 计算点云的非凸包络曲面。
        alpha：控制"贴合度"，越小越紧贴（但可能产生孔洞），越大越凸（接近凸包）。
        """
        if alpha is None:
            bbox_diag  = np.linalg.norm(np.max(coords, 0) - np.min(coords, 0))
            alpha_value = bbox_diag / 5   # 自动估算 alpha（包围盒对角线的 1/5）
        else:
            alpha_value = np.mean(alpha)
        points = [Point_3(pt[0], pt[1], pt[2]) for pt in coords]  # noqa: F405
        Q      = Polyhedron_3()
        _a     = alpha_wrap_3(points, alpha_value, 0.01, Q)  # noqa: F405

        alpha_shape_vertices = np.array(
            [vertex_to_tuple(vertex.point()) for vertex in Q.vertices()]
        )
        alpha_shape_faces = np.array(
            [np.array(halfedge_to_triangle(face.halfedge())) for face in Q.facets()]
        )
        return alpha_shape_vertices, alpha_shape_faces

    # reachableWorkspace：返回采样点（非凸包，只是点集）
    verts       = pin.reachableWorkspace(robot, q0, horizon, frame, n_samples, facet_dims)
    verts       = verts.T
    alpha       = 0.2
    verts, faces = alpha_shape_with_cgal(verts, alpha)
    verts       = faces.reshape(-1, 3)
    faces       = np.arange(len(verts)).reshape(-1, 3)

print("------------------- Display Vertex")

# ---- 在 Meshcat 中显示三角网格 ----
# TriangularMeshGeometry：Meshcat 原生三角网格几何体
# wireframe=True：显示线框而非实体（便于看清内部结构）
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
