from pathlib import Path
from sys import argv

import pinocchio

# ============================================================
# 几何模型的加载与位姿更新
# 示例展示了 Pinocchio 几何系统的完整流程：
#   1. 加载运动学模型（Model）
#   2. 同时加载碰撞/视觉几何模型（GeometryModel）
#   3. 为三个模型创建对应的数据结构（Data/GeometryData）
#   4. 计算正向运动学后，将位姿传播到所有几何体
#
# 三模型架构：
#   Model         → 关节/连杆的运动学/动力学参数
#   collision_model → 碰撞检测用几何体（通常是简化形状：Box/Capsule）
#   visual_model   → 可视化用几何体（通常是精细网格：STL/DAE/OBJ）
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path = Path(
    # 支持命令行参数指定模型目录，默认使用内置示例
    (pinocchio_model_dir / "example-robot-data/robots") if len(argv) < 2 else argv[1]
)
mesh_dir        = pinocchio_model_dir
urdf_model_path = model_path / "ur_description/urdf/ur5_robot.urdf"

# ---- buildModelsFromUrdf：一次性加载三个模型 ----
# 返回值：(运动学模型, 碰撞几何模型, 视觉几何模型)
# 等价于分别调用：
#   model = buildModelFromUrdf(urdf_model_path)
#   collision_model = buildGeomFromUrdf(model, urdf, COLLISION, mesh_dir)
#   visual_model    = buildGeomFromUrdf(model, urdf, VISUAL,    mesh_dir)
# mesh_dir：搜索网格文件（STL/OBJ/DAE）的根目录
# 若 URDF 中的 mesh filename 是相对路径，会在 mesh_dir 下递归搜索
model, collision_model, visual_model = pinocchio.buildModelsFromUrdf(
    urdf_model_path, mesh_dir
)
print("model name: " + model.name)

# ---- createDatas：为三个模型一次性创建数据结构 ----
# 等价于：
#   data           = model.createData()
#   collision_data = pinocchio.GeometryData(collision_model)
#   visual_data    = pinocchio.GeometryData(visual_model)
# GeometryData 存储：
#   oMg[k]：第 k 个几何体在世界坐标系中的 SE(3) 位姿
#   collisionResults[k]：第 k 个碰撞对的检测结果（仅碰撞模型使用）
data, collision_data, visual_data = pinocchio.createDatas(
    model, collision_model, visual_model
)

# 采样一个随机合法配置（在关节限位内）
q = pinocchio.randomConfiguration(model)
print(f"q: {q.T}")

# ---- 正向运动学：计算所有关节在世界系中的位姿 ----
# 结果存储在 data.oMi[k]：第 k 个关节的 SE(3) 变换（origin → joint_k）
pinocchio.forwardKinematics(model, data, q)

# ---- updateGeometryPlacements：将关节位姿传播到几何体 ----
# 内部：对每个几何体 k，查找其所属关节 j：
#   geom_data.oMg[k] = data.oMi[j] * geom_model.geometryObjects[k].placement
# 其中 placement 是几何体相对所属关节的固定 SE(3) 偏移（来自 URDF）
# 必须在 forwardKinematics 之后调用（依赖 data.oMi）
pinocchio.updateGeometryPlacements(model, data, collision_model, collision_data)
pinocchio.updateGeometryPlacements(model, data, visual_model,    visual_data)

# ---- 打印各关节位置 ----
print("\nJoint placements:")
for name, oMi in zip(model.names, data.oMi):
    # oMi.translation：关节原点在世界坐标系中的位置 [x, y, z]
    print("{:<24} : {: .2f} {: .2f} {: .2f}".format(name, *oMi.translation.T.flat))

# ---- 打印碰撞几何体位置 ----
print("\nCollision object placements:")
for k, oMg in enumerate(collision_data.oMg):
    # oMg：几何体（碰撞形状）在世界坐标系中的 SE(3) 位姿
    # 碰撞体通常是简化形状（如胶囊、长方体），比视觉网格计算更快
    print("{:d} : {: .2f} {: .2f} {: .2f}".format(k, *oMg.translation.T.flat))

# ---- 打印视觉几何体位置 ----
print("\nVisual object placements:")
for k, oMg in enumerate(visual_data.oMg):
    # oMg：几何体（视觉网格）在世界坐标系中的 SE(3) 位姿
    # 视觉体通常是精细 STL/DAE 网格，用于渲染；不参与碰撞计算
    print("{:d} : {: .2f} {: .2f} {: .2f}".format(k, *oMg.translation.T.flat))
