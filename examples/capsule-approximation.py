"""
Copyright (c) 2020 INRIA
Inspired from Antonio El Khoury PhD:
https://tel.archives-ouvertes.fr/file/index/docid/833019/filename/thesis.pdf
Section 3.8.1 Computing minimum bounding capsules
"""

from pathlib import Path

import hppfcl
import numpy as np
import scipy.optimize as optimize

# ============================================================
# 最小包围胶囊体近似（Minimum Bounding Capsule Approximation）
# 应用场景：
#   - 将复杂 STL/OBJ 网格近似为胶囊体（Capsule），用于快速碰撞检测
#   - 人形机器人的连杆简化：用胶囊体代替精细网格，碰撞检测速度快 10~100 倍
#   - 批量处理整个 URDF 文件，生成胶囊体化的碰撞 URDF
#
# 胶囊体定义：两端点 a, b（定义轴线段）+ 半径 r
#   体积 = π·r²·|b-a| + (4/3)π·r³（圆柱 + 两个半球）
#
# 优化问题：min_{a,b,r} Volume(a,b,r)
#   subject to: max_i dist(vertex_i, segment(a,b)) ≤ r
#   → 找包围所有顶点的最小体积胶囊体
# ============================================================

EPSILON = 1e-8
# 约束松弛量：防止数值精度问题导致顶点恰好在胶囊体表面（内外判断不稳定）
CONSTRAINT_INFLATION_RATIO = 5e-3


def capsule_volume(a, b, r):
    """胶囊体体积 = 圆柱部分 + 两个半球"""
    return np.linalg.norm(b - a) * np.pi * r**2 + 4 / 3 * np.pi * r**3


def distance_points_segment(p, a, b):
    """
    计算点集 p（shape: (N, 3)）中每个点到线段 ab 的距离，返回最大值。
    用于约束函数：所有顶点到胶囊轴线的最大距离 ≤ r。
    """
    ap = p - a
    ab = b - a
    # 参数化投影：t = dot(ap, ab) / |ab|²，截断到 [0,1]（线段而非直线）
    t  = ap.dot(ab) / ab.dot(ab)
    t  = np.clip(t, 0, 1)
    # 线段上的最近点
    p_witness = a[None, :] + (b - a)[None, :] * t[:, None]
    # 各点到最近点的距离，取最大值（用于约束判断）
    dist = np.linalg.norm(p - p_witness, axis=1).max()
    return dist


def pca_approximation(vertices):
    """
    用 PCA（主成分分析）得到胶囊体的初始估计。
    思路：
      1. 找顶点分布的主方向（第一主成分）作为胶囊轴线方向
      2. 轴线端点 = 沿主方向的最小/最大投影
      3. 半径 = 在垂直于主方向的平面内的最大距离
    """
    mean     = vertices.mean(axis=0)
    vertices -= mean   # 去中心化
    # SVD 分解：U S Vᵀ，vh 的行向量为主成分方向（按方差降序）
    u, s, vh = np.linalg.svd(vertices, full_matrices=True)
    components = vh
    pca_proj   = vertices.dot(components.T)   # 各顶点在主成分坐标系中的投影
    vertices  += mean   # 恢复中心

    # 沿第一主方向（最大方差方向）的端点（加 EPSILON 保证包含所有顶点）
    a0     = mean + components[0] * (pca_proj[:, 0].min() - EPSILON)
    b0     = mean + components[0] * (pca_proj[:, 0].max() + EPSILON)
    # 在垂直于主方向的平面内的最大距离 = 初始半径估计
    radius = np.linalg.norm(pca_proj[:, 1:], axis=1).max()
    return a0, b0, radius


def capsule_approximation(vertices):
    """
    求最小包围胶囊体：以 PCA 结果为初值，用非线性约束优化精化。
    优化变量：x = [a(3), b(3), r(1)] = 7 维向量
    目标：min capsule_volume(a, b, r)
    约束：max_i dist(vertex_i, segment(a,b)) - r ≤ -inflation（所有顶点在胶囊内）
    """
    a0, b0, r0 = pca_approximation(vertices)
    constraint_inflation = CONSTRAINT_INFLATION_RATIO * r0   # 约束松弛量
    x0 = np.array(list(a0) + list(b0) + [r0])

    def constraint_cap(x):
        """约束函数：最大顶点距离 - 半径（≤ -inflation 时满足约束）"""
        return distance_points_segment(vertices, x[:3], x[3:6]) - x[6]

    def capsule_vol(x):
        return capsule_volume(x[:3], x[3:6], x[6])

    # NonlinearConstraint：-∞ ≤ constraint_cap(x) ≤ -inflation
    # 即：所有顶点到胶囊轴线的最大距离 ≤ r - inflation（严格包含，留有余量）
    constraint = optimize.NonlinearConstraint(
        constraint_cap, lb=-np.inf, ub=-constraint_inflation
    )
    res = optimize.minimize(capsule_vol, x0, constraints=constraint)

    # 验证约束满足（数值优化可能不精确）
    res_constraint = constraint_cap(res.x)
    err = (
        "The computed solution is invalid, "
        "a vertex is at a distance {:.5f} of the capsule."
    )
    assert res_constraint <= 1e-4, err.format(res_constraint)

    a, b, r = res.x[:3], res.x[3:6], res.x[6]
    return a, b, r


def approximate_mesh(filename, lMg):
    """
    加载网格文件，将所有顶点变换到连杆坐标系，计算最小包围胶囊。
    filename：STL/OBJ 网格文件路径
    lMg：网格相对连杆坐标系的 SE(3) 变换（来自 URDF <origin>）
    """
    mesh_loader = hppfcl.MeshLoader()
    mesh = mesh_loader.load(filename, np.ones(3))   # 缩放因子 = [1,1,1]
    # 将每个顶点从网格局部坐标系变换到连杆坐标系
    vertices = np.array([lMg * mesh.vertices(i) for i in range(mesh.num_vertices)])
    assert vertices.shape == (mesh.num_vertices, 3)
    a, b, r = capsule_approximation(vertices)
    return a, b, r


def parse_urdf(infile, outfile):
    """
    批量处理 URDF 文件：将所有 <collision> 中的 <mesh> 替换为 <cylinder>（胶囊近似）。
    lxml：用于解析和修改 XML；tqdm：显示处理进度条。
    注意：这里用 <cylinder> 标签代表胶囊体（URDF 标准中胶囊 = 端盖圆柱体，但 Pinocchio 可识别）
    """
    from lxml import etree

    tree = etree.parse(infile)

    def get_path(fn):
        """解析 URDF 中的 mesh filename（支持 package:// ROS 路径协议）"""
        if fn.startswith("package://"):
            relpath = fn[len("package://"):]
            import os
            for rospath in os.environ["ROS_PACKAGE_PATH"].split(":"):
                abspath = Path(rospath) / relpath
                if abspath.is_file():
                    return abspath
            raise ValueError("Could not find " + fn)
        return fn

    def get_transform(origin):
        """从 URDF <origin xyz="..." rpy="..."> 元素解析 SE(3) 变换"""
        from pinocchio import SE3, rpy
        _xyz = [float(v) for v in origin.attrib.get("xyz", "0 0 0").split(" ")]
        _rpy = [float(v) for v in origin.attrib.get("rpy", "0 0 0").split(" ")]
        return SE3(rpy.rpyToMatrix(*_rpy), np.array(_xyz))

    def set_transform(origin, a, b):
        """
        将胶囊体的两端点 a, b 写回 URDF <origin> 元素。
        胶囊体轴线方向 = (b-a)/|b-a|，从 Z 轴旋转到该方向。
        """
        from pinocchio import Quaternion, rpy
        length = np.linalg.norm(b - a)
        z      = (b - a) / length
        # 四元数旋转：将 FCL 胶囊默认轴 [0,0,1] 旋转到轴线方向 z
        R = Quaternion.FromTwoVectors(np.array([0, 0, 1]), z).matrix()
        # 胶囊体中心 = (a+b)/2
        origin.attrib["xyz"] = " ".join([str(v) for v in ((a + b) / 2)])
        origin.attrib["rpy"] = " ".join([str(v) for v in rpy.matrixToRpy(R)])

    from tqdm import tqdm

    # 遍历所有 <robot/link/collision/geometry/mesh> 节点
    for mesh in tqdm(
        tree.xpath("/robot/link/collision/geometry/mesh"), desc="Generating capsules"
    ):
        geom = mesh.getparent()    # <geometry>
        coll = geom.getparent()    # <collision>
        link = coll.getparent()    # <link>

        # 确保 <collision> 有 <origin> 子元素（若无则创建）
        if coll.find("origin") is None:
            o = etree.Element("origin")
            o.tail = geom.tail
            coll.insert(0, o)
        origin = coll.find("origin")
        lMg    = get_transform(origin)   # 网格相对连杆的变换

        meshfile = get_path(mesh.attrib["filename"])
        name     = Path(meshfile).name

        # ---- 计算最小包围胶囊 ----
        a, b, radius = approximate_mesh(meshfile, lMg)
        length       = np.linalg.norm(b - a)

        # 更新 <origin> 为胶囊体中心和轴线方向
        set_transform(origin, a, b)

        # 将 <mesh> 节点替换为 <cylinder>（用圆柱标签存储胶囊参数）
        mesh.tag = "cylinder"
        mesh.attrib.pop("filename")
        mesh.attrib["radius"] = str(radius)
        mesh.attrib["length"] = str(length)
        coll.attrib["name"]   = name

        # 在 <link> 下记录胶囊体名称（用于 Pinocchio 识别）
        if link.find("collision_checking") is None:
            link.append(etree.Element("collision_checking"))
        collision_checking = link.find("collision_checking")
        collision_checking.append(etree.Element("capsule"))
        collision_checking[-1].attrib["name"] = name

    tree.write(outfile)


if __name__ == "__main__":
    # Example for a single capsule
    # filename = "mesh.obj"
    # mesh_loader = hppfcl.MeshLoader()
    # mesh = mesh_loader.load(filename, np.ones(3))
    # vertices = mesh.vertices()
    # a, b, r = capsule_approximation(vertices)

    # Example for a whole URDF model
    # This path refers to Pinocchio source code but you can define your own directory
    # here.

    pinocchio_model_dir = Path(__file__).parent.parent / "models"
    urdf_filename = (
        pinocchio_model_dir
        + "models/example-robot-data/robots/ur_description/urdf/ur5_gripper.urdf"
    )
    parse_urdf(urdf_filename, "ur5_gripper_with_capsules.urdf")
