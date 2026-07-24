from pathlib import Path

import numpy as np
import pinocchio as pin

# ============================================================
# 降阶模型：将指定关节锁定在参考配置处，构建自由度更少的模型
# 应用场景：
#   - 人形上身固定，只研究腿部动力学
#   - 手臂固定，研究躯干+腿部的平衡控制
#   - 加速仿真（减少计算量）
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir   = pinocchio_model_dir
urdf_filename = model_path / "ur_description/urdf/ur5_robot.urdf"

# 同时加载运动学模型、碰撞模型、视觉模型（降阶时可一并处理）
model, collision_model, visual_model = pin.buildModelsFromUrdf(urdf_filename, mesh_dir)

print("standard model: dim=" + str(len(model.joints)))
for jn in model.joints:
    print(jn)
print("-" * 30)

# 要锁定的关节名称列表（这些关节在降阶模型中被固定）
# 实际使用时替换为需要固定的关节名
jointsToLock = ["wrist_1_joint", "wrist_2_joint", "wrist_3_joint"]

# 将关节名转为 ID（不存在的名称会被警告跳过）
jointsToLockIDs = []
for jn in jointsToLock:
    if model.existJointName(jn):
        jointsToLockIDs.append(model.getJointId(jn))
    else:
        print("Warning: joint " + str(jn) + " does not belong to the model!")

# 锁定关节时使用的参考配置（被锁定关节固定在此角度）
# 长度必须等于完整模型的 nq
initialJointConfig = np.array([
    0, 0, 0,   # 肩部和肘部关节（保持零位）
    1, 1, 1,   # 被锁定的腕部关节（固定在 1 rad）
])

# ---- Option 1：仅降阶运动学/动力学模型（不含几何模型）----
# 最轻量，适合纯动力学计算（不需要可视化时）
model_reduced = pin.buildReducedModel(model, jointsToLockIDs, initialJointConfig)

# ---- Option 2：同时降阶视觉几何模型（用于可视化）----
# 返回降阶后的运动学模型 + 对应的视觉模型
model_reduced, visual_model_reduced = pin.buildReducedModel(
    model, visual_model, jointsToLockIDs, initialJointConfig
)

# ---- Option 3：同时降阶多个几何模型（视觉 + 碰撞）----
# 返回顺序与传入的 geom_models 列表一致
geom_models = [visual_model, collision_model]
model_reduced, geometric_models_reduced = pin.buildReducedModel(
    model,
    list_of_geom_models=geom_models,
    list_of_joints_to_lock=jointsToLockIDs,
    reference_configuration=initialJointConfig,
)
# 按顺序解包（与传入顺序一致）
visual_model_reduced, collision_model_reduced = (
    geometric_models_reduced[0],
    geometric_models_reduced[1],
)

print("joints to lock (only ids):", jointsToLockIDs)
print("reduced model: dim=" + str(len(model_reduced.joints)))
print("-" * 30)

# ---- Option 4：使用 RobotWrapper（可混用关节名和关节 ID）----
# RobotWrapper 封装了 model + collision_model + visual_model，接口更高层
# reference_configuration 可选，不传则使用零位
mixed_jointsToLockIDs = [jointsToLockIDs[0], "wrist_2_joint", "wrist_3_joint"]
robot = pin.RobotWrapper.BuildFromURDF(urdf_filename, mesh_dir)
reduced_robot = robot.buildReducedRobot(
    list_of_joints_to_lock=mixed_jointsToLockIDs,
    reference_configuration=initialJointConfig,
)

print("mixed joints to lock (names and ids):", mixed_jointsToLockIDs)
print("RobotWrapper reduced model: dim=" + str(len(reduced_robot.model.joints)))
for jn in robot.model.joints:
    print(jn)
