import math
import sys
import time
from pathlib import Path

import hppfcl as fcl
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer as Visualizer

# ============================================================
# 模型拼接：将两个独立模型合并为一个
# 应用场景：
#   - 为机器人末端添加工具/负载
#   - 将传感器模型附加到机器人某个 Frame 上
#   - 将两个子系统（如臂 + 腕）合并
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir   = pinocchio_model_dir
urdf_model_path = model_path / "ur_description/urdf/ur5_robot.urdf"

# ---- 模型 1：从 URDF 加载 UR5 机械臂 ----
model1, collision_model1, visual_model1 = pin.buildModelsFromUrdf(
    urdf_model_path, package_dirs=mesh_dir
)

# ---- 模型 2：手工编程构建一个球形关节摆（模拟末端负载）----
model2 = pin.Model()
model2.name = "pendulum"
geom_model = pin.GeometryModel()

parent_id       = 0                      # 父关节 ID（0 = 世界/宇宙）
joint_placement = pin.SE3.Identity()     # 关节相对父关节的固定变换（单位变换）
body_mass       = 1.0                    # 连杆质量 [kg]
body_radius     = 1e-2                   # 球体半径 [m]

# 添加球形关节（3 DOF 旋转，适合末端摆）
# JointModelSpherical：3 DOF，以四元数表示，nq=4 nv=3
joint_name = "joint_spherical"
joint_id = model2.addJoint(
    parent_id, pin.JointModelSpherical(), joint_placement, joint_name
)

# 为关节绑定连杆惯量（FromSphere：均匀球体的惯量）
body_inertia = pin.Inertia.FromSphere(body_mass, body_radius)
body_placement = joint_placement.copy()
body_placement.translation[2] = 0.1   # 连杆质心在关节 Z 轴上方 0.1 m
model2.appendBodyToJoint(joint_id, body_inertia, body_placement)

# 添加球形视觉几何体
geom1_obj = pin.GeometryObject("ball", joint_id, body_placement, fcl.Sphere(body_radius))
geom1_obj.meshColor = np.ones(4)   # 白色
geom_model.addGeometryObject(geom1_obj)

# 添加杆状几何体（从关节到质心的连接杆）
shape2_placement = body_placement.copy()
shape2_placement.translation[2] /= 2.0   # 杆的中点
geom2_obj = pin.GeometryObject(
    "bar", joint_id, shape2_placement,
    fcl.Cylinder(body_radius / 4.0, body_placement.translation[2])
)
geom2_obj.meshColor = np.array([0.0, 0.0, 0.0, 1.0])  # 黑色
geom_model.addGeometryObject(geom2_obj)

visual_model2 = geom_model

# ---- 模型拼接：将 model2（摆）附加到 model1（UR5）的末端 Frame ----
# appendModel(model1, model2, geom1, geom2, frame_id, placement)
#   frame_id：model1 中作为连接点的 Frame（末端执行器 "tool0"）
#   placement：model2 相对该 Frame 的固定偏移（这里用单位变换，即直接附加到末端）
frame_id_end_effector = model1.getFrameId("tool0")
model, visual_model = pin.appendModel(
    model1, model2,
    visual_model1, visual_model2,
    frame_id_end_effector,
    pin.SE3.Identity(),   # 摆直接附加在末端，无额外偏移
)

print(
    f"Check the joints of the appended model:\n {model} \n "
    "->Notice the spherical joint at the end."
)

try:
    viz = Visualizer(model, visual_model, visual_model)
    viz.initViewer(open=True)
except ImportError as err:
    print("Error while initializing the viewer. It seems you should install Python meshcat")
    print(err)
    sys.exit(0)

viz.loadViewerModel()

# 为合并后模型设置关节限位，再随机采样配置并显示
model.lowerPositionLimit.fill(-math.pi / 2)
model.upperPositionLimit.fill(+math.pi / 2)
q = pin.randomConfiguration(model)
viz.display(q)
time.sleep(1.0)
