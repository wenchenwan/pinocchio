# This examples shows how to load and move a robot in meshcat.
# Note: this feature requires Meshcat to be installed, this can be done using
# pip install --user meshcat

import sys
from pathlib import Path

import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# Meshcat 查看器综合示例：
# 演示内容：
#   1. 加载机器人到 Meshcat（浮动基 Solo 四足机器人）
#   2. 显示帧速度（drawFrameVelocities）
#   3. 凸包（Convex Hull）提取与显示
#   4. 多机器人同场景显示（共享同一 viewer 实例）
#   5. ABA 正向动力学仿真 + 视频录制
# ============================================================

# Load the URDF model.
# Conversion with str seems to be necessary when executing this file with ipython
pinocchio_model_dir = Path(__file__).parent.parent / "models"

model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir   = pinocchio_model_dir
# Solo：四足机器人，12 DOF（每腿 3 DOF），使用浮动基
urdf_filename    = "solo.urdf"
urdf_model_path  = model_path / "solo_description/robots" / urdf_filename

# 浮动基加载（JointModelFreeFlyer）：nq=12+7=19, nv=12+6=18
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    urdf_model_path, mesh_dir, pin.JointModelFreeFlyer()
)

# ---- 启动 Meshcat 服务器并打开浏览器 ----
# initViewer(open=True)：自动启动内嵌 Meshcat 服务器并打开浏览器标签页
# 也可以单独在终端运行 "meshcat-server"，保持服务器在脚本退出后继续运行
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

# Load the robot in the viewer.
viz.loadViewerModel()

# 显示零位配置
q0 = pin.neutral(model)
viz.display(q0)
# displayVisuals(True)：启用视觉几何体渲染（默认可能只显示碰撞几何体）
viz.displayVisuals(True)

# ---- 凸包（Convex Hull）提取与显示 ----
# buildConvexRepresentation(True)：从 mesh 计算凸包（hpp-fcl 提供）
# 凸包用途：简化碰撞检测（凸体间 GJK 算法效率最高）
mesh   = visual_model.geometryObjects[0].geometry
mesh.buildConvexRepresentation(True)
convex = mesh.convex

if convex is not None:
    placement = pin.SE3.Identity()
    placement.translation[0] = 2.0   # 将凸包显示在机器人旁边（X+2m）
    geometry = pin.GeometryObject("convex", 0, placement, convex)
    geometry.meshColor = np.ones(4)   # 白色
    # ---- Phong 光照材质 ----
    # overrideMaterial=True：使用自定义材质覆盖 URDF 中的材质
    geometry.overrideMaterial = True
    geometry.meshMaterial = pin.GeometryPhongMaterial()
    geometry.meshMaterial.meshEmissionColor = np.array([1.0, 0.1, 0.1, 1.0])   # 自发光颜色（红色）
    geometry.meshMaterial.meshSpecularColor = np.array([0.1, 1.0, 0.1, 1.0])   # 镜面高光颜色（绿色）
    geometry.meshMaterial.meshShininess    = 0.8   # 高光系数（0~1，越大越亮）
    visual_model.addGeometryObject(geometry)
    # 修改 visual_model 后必须重建 visualizer 内部数据结构
    viz.rebuildData()

# ---- 多机器人同场景显示 ----
# 创建第二个 visualizer，共享同一个 Meshcat viewer 实例（viz.viewer）
# rootNodeName 区分两个机器人在场景树中的节点名称
viz2 = MeshcatVisualizer(model, collision_model, visual_model)
viz2.initViewer(viz.viewer)   # 共享 viewer（不创建新服务器）
viz2.loadViewerModel(rootNodeName="pinocchio2")
q    = q0.copy()
q[1] = 1.0   # 第二个机器人在 Y 轴上偏移 1m
viz2.display(q)

# ---- 显示帧速度箭头 ----
# standing config：Solo 站立时的关节角配置（来自 SRDF 或手工标定）
q1 = np.array(
    [0.0, 0.0, 0.235, 0.0, 0.0, 0.0, 1.0, 0.8, -1.6, 0.8, -1.6, -0.8, 1.6, -0.8, 1.6]
)
v0     = np.random.randn(model.nv) * 2   # 随机速度（用于显示速度箭头）
data   = viz.data
pin.forwardKinematics(model, data, q1, v0)   # 计算位姿+速度（含 data.v）
frame_id = model.getFrameId("HR_FOOT")       # 后右足端点 Frame ID
viz.display()
# drawFrameVelocities：在 Frame 处绘制速度箭头（线速度=蓝色，角速度=红色）
viz.drawFrameVelocities(frame_id=frame_id)

# ---- 仿真循环：ABA + 半隐式欧拉 ----
model.gravity.linear[:] = 0.0   # 关闭重力（在零重力下仿真）
dt = 0.01


def sim_loop():
    tau0   = np.zeros(model.nv)
    qs     = [q1]
    vs     = [v0]
    nsteps = 100
    for i in range(nsteps):
        q    = qs[i]
        v    = vs[i]
        # ABA 正向动力学：q̈ = ABA(q, v, τ)
        a1   = pin.aba(model, data, q, v, tau0)
        # 半隐式欧拉
        vnext = v + dt * a1
        qnext = pin.integrate(model, q, dt * vnext)
        qs.append(qnext)
        vs.append(vnext)
        viz.display(qnext)
        viz.drawFrameVelocities(frame_id=frame_id)
    return qs, vs


qs, vs = sim_loop()

fid2 = model.getFrameId("FL_FOOT")   # 前左足端点 Frame ID


def my_callback(i, *args):
    """每帧回调：刷新两个足端的速度箭头显示"""
    viz.drawFrameVelocities(frame_id)
    viz.drawFrameVelocities(fid2)


# ---- 视频录制 ----
# create_video_ctx：上下文管理器，在 with 块内的所有 display() 调用都会被录制
# viz.play(qs, dt, callback)：按仿真轨迹逐帧播放，并调用 callback 更新额外显示
with viz.create_video_ctx("../leap.mp4"):
    viz.play(qs, dt, callback=my_callback)
