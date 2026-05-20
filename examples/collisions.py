from pathlib import Path

import pinocchio as pin

# ============================================================
# 碰撞检测：buildGeomFromUrdf + computeCollisions
# 应用场景：
#   - 运动规划中检查是否与环境/自身发生碰撞
#   - 自碰撞检测（人形机器人手臂/腿部自交叉检查）
#   - 安全监控（实时检测机器人与障碍物的距离）
#
# Pinocchio 碰撞检测基于 hpp-fcl（高性能 FCL 分支）：
#   - 支持多种几何体：Box、Sphere、Capsule、Cylinder、Mesh 等
#   - 宽相（broadphase）：AABB 树快速排除不可能碰撞的几何对
#   - 窄相（narrowphase）：GJK/EPA 计算精确碰撞/距离信息
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path    = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir      = pinocchio_model_dir
# Romeo：Softbank/Aldebaran 的人形机器人，较复杂，适合演示自碰撞检测
urdf_filename = "romeo_small.urdf"
urdf_model_path = model_path / "romeo_description/urdf" / urdf_filename

# ---- 加载运动学模型（浮动基）----
# pin.JointModelFreeFlyer()：添加 6-DOF 浮动基（nq+=7, nv+=6）
# 人形机器人必须使用浮动基才能正确表示整体位姿
model = pin.buildModelFromUrdf(urdf_model_path, pin.JointModelFreeFlyer())

# ---- 单独加载碰撞几何模型 ----
# buildGeomFromUrdf：只加载指定类型的几何体（COLLISION 或 VISUAL）
# 与 buildModelsFromUrdf 不同，这里可以单独加载几何模型而不重新解析运动学
# GeometryType.COLLISION：解析 URDF 中 <collision> 标签下的几何体
# mesh_dir：网格文件（STL/DAE/OBJ）的搜索根目录
geom_model = pin.buildGeomFromUrdf(
    model, urdf_model_path, pin.GeometryType.COLLISION, mesh_dir
)

# ---- 添加所有可能的碰撞对 ----
# addAllCollisionPairs()：枚举所有几何体两两组合
# 对于 N 个几何体：最多 N(N-1)/2 个碰撞对（通常 100~1000 对）
# 大部分碰撞对对自碰撞检测无意义（如相邻连杆的几何体通常不会自碰撞）
geom_model.addAllCollisionPairs()
print("num collision pairs - initial:", len(geom_model.collisionPairs))

# ---- 从 SRDF 文件移除无效碰撞对 ----
# SRDF（Semantic Robot Description Format）是 MoveIt 使用的补充格式
# 其中 <disable_collisions> 标签列出了"永远不会碰撞"的几何体对（如相邻连杆）
# 移除这些对可以显著减少碰撞检测计算量（通常减少 70%~90%）
srdf_filename   = "romeo.srdf"
srdf_model_path = model_path / "romeo_description/srdf" / srdf_filename
pin.removeCollisionPairs(model, geom_model, srdf_model_path)
print(
    "num collision pairs - after removing useless collision pairs:",
    len(geom_model.collisionPairs),
)

# ---- 从 SRDF 加载参考配置（如 half_sitting）----
# SRDF 中的 <group_state> 标签存储了命名配置（如站立、坐姿、零位等）
# loadReferenceConfigurations 将这些配置加载到 model.referenceConfigurations 字典中
pin.loadReferenceConfigurations(model, srdf_model_path)

# half_sitting：Romeo 的标准半蹲站立配置（膝盖略弯，重心稳定）
# 这是碰撞检测的测试配置
q = model.referenceConfigurations["half_sitting"]

# ---- 创建数据结构 ----
data      = model.createData()
# GeometryData：碰撞检测的可变状态容器
# 存储：每个几何体在世界坐标系中的位姿、碰撞检测结果
geom_data = pin.GeometryData(geom_model)

# ---- 批量碰撞检测 ----
# computeCollisions 内部流程：
#   1. forwardKinematics(model, data, q)（更新关节位姿）
#   2. updateGeometryPlacements（将关节位姿传播到几何体）
#   3. 对所有碰撞对调用 FCL 碰撞查询
# 最后一个参数 False：不提前终止（继续检查所有对，即使已发现碰撞）
# 设为 True：一旦发现第一个碰撞立即返回（适合安全快速检查）
pin.computeCollisions(model, data, geom_model, geom_data, q, False)

# ---- 打印所有碰撞对的检测结果 ----
for k in range(len(geom_model.collisionPairs)):
    cr = geom_data.collisionResults[k]    # CollisionResult：包含碰撞/距离信息
    cp = geom_model.collisionPairs[k]     # CollisionPair：两个几何体的 ID 对
    print(
        "collision pair:",
        cp.first,       # 第一个几何体 ID
        ",",
        cp.second,      # 第二个几何体 ID
        "- collision:",
        "Yes" if cr.isCollision() else "No",
    )

# ---- 单对碰撞检测（更新几何位姿后手动检查指定对）----
# updateGeometryPlacements：仅更新几何体在世界坐标系中的位姿
# 不执行碰撞检测本身（与 computeCollisions 的区别：后者一体化）
pin.updateGeometryPlacements(model, data, geom_model, geom_data, q)
# computeCollision(geom_model, geom_data, pair_id)：只检测指定索引的碰撞对
# 适用场景：只关心特定关键对（如末端执行器与工件）
pin.computeCollision(geom_model, geom_data, 0)
