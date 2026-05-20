from pathlib import Path

import pinocchio as pin

# ============================================================
# Contact Cholesky 分解与 Delassus 算子
# 用于高效求解带接触约束的 KKT 鞍点系统
#
# Delassus 矩阵：Λ = J M^{-1} J^T ∈ R^{nc×nc}
#   - 接触空间的有效惯量（等效质量矩阵）
#   - 衡量单位接触力引起的接触点加速度响应
#   - 对角元素越小，接触响应越"重"（刚性更高）
#
# Contact Cholesky 比 dense LU 快，因为它利用了运动树的稀疏结构
# 复杂度：O(n·nc²) 而非 O((n+nc)³)
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
# ANYmal-B：四足机器人，固定基 URDF（足底使用 3D 点接触）
urdf_filename = (
    pinocchio_model_dir
    / "example-robot-data/robots/anymal_b_simple_description/robots/anymal.urdf"
)

# 固定基建模（ANYmal 的 URDF 本身已含浮动基关节，无需额外指定）
model = pin.buildModelFromUrdf(urdf_filename)
data  = pin.Data(model)

# 零位配置（所有关节角为 0）
q0 = pin.neutral(model)

# ---- 构造四足 3D 点接触模型 ----
feet_names    = ["LH_FOOT", "RH_FOOT", "LF_FOOT", "RF_FOOT"]
feet_frame_ids = [model.getFrameId(name) for name in feet_names]

contact_models = []
for fid in feet_frame_ids:
    frame = model.frames[fid]
    # CONTACT_3D：3D 点接触（只约束位置，不约束姿态，每个接触贡献 3 维约束）
    # LOCAL_WORLD_ALIGNED：接触 Jacobian 在原点位于足端、轴与世界系对齐的坐标系中表达
    #   优点：相比 LOCAL 更直观（Z 轴始终竖直，对应法向力方向）
    cmodel = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,
        frame.parent,        # 足端关节 ID
        frame.placement,     # 接触点相对该关节的局部位置
        pin.LOCAL_WORLD_ALIGNED,
    )
    contact_models.append(cmodel)

contact_data = [cmodel.createData() for cmodel in contact_models]

# 预初始化接触约束内存（必须在调用 constraintDynamics 或 contact_chol.compute 前执行）
pin.initConstraintDynamics(model, data, contact_models)

# 计算质量矩阵 M（CRBA），Contact Cholesky 需要 M 已计算
pin.crba(model, data, q0)

# ---- Contact Cholesky 分解 ----
# 对 KKT 系统进行稀疏 Cholesky 分解：
#   [M  J^T] = L D L^T
#   [J  0  ]
# 利用运动树的稀疏性加速分解，避免填充（fill-in）
data.contact_chol.compute(model, data, contact_models, contact_data)

# ---- 提取 Delassus 矩阵 ----
# Delassus 矩阵逆（操作空间惯量矩阵逆）：Λ^{-1} = J M^{-1} J^T ∈ R^{nc×nc}
# 物理含义：单位广义力在接触空间引起的加速度（接触空间"柔度"矩阵）
delassus_matrix = data.contact_chol.getInverseOperationalSpaceInertiaMatrix()

# 操作空间惯量矩阵：Λ = (J M^{-1} J^T)^{-1} ∈ R^{nc×nc}
# 物理含义：在接触空间中施加单位加速度所需的力（接触空间"刚度"矩阵）
# 用于求接触力：λ = Λ * (J q̈_free + J̇ v)
delassus_matrix_inv = data.contact_chol.getOperationalSpaceInertiaMatrix()
