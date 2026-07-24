from pathlib import Path

import numpy as np
import pinocchio

# ============================================================
# Mimic 关节（传动比约束）动力学
# 应用场景：夹爪耦合关节、并联传动机构、差速驱动等
# 物理含义：关节 j 的运动 = scaling * 关节 i 的运动 + offset
#   例：夹爪的两根手指通过齿轮耦合，一个关节驱动两个手指同步运动
#
# 数学处理（参考 Featherstone RBDA 第 10 章 Gears）：
#   设 G ∈ R^{nv_mimic × nv_full} 为传动矩阵，则：
#     v_full = G^T * v_mimic（速度映射）
#     τ_mimic = G * τ_full  （力矩映射，力矩方向相反）
#     M_mimic = G * M_full * G^T（惯量映射）
# ============================================================

pinocchio_model_dir = Path(__file__).parent.parent / "models"
# Baxter 双臂机器人：每个手臂末端有两个耦合夹爪关节（mimic 关系）
model_path = pinocchio_model_dir / "baxter_simple.urdf"

# ---- 方法 1：从 URDF 加载，忽略 mimic 标签（默认行为）----
# 所有关节都被视为独立的，nq 更大
model_full = pinocchio.buildModelFromUrdf(model_path)

# ---- 方法 2：从 URDF 加载，解析 mimic 标签 ----
# mimic=True 时，mimic 关节被转换为受约束关节，nq 减小
# 降阶后的模型 DOF = 独立关节数（耦合关节被折叠）
model_mimic_from_urdf = pinocchio.buildModelFromUrdf(model_path, mimic=True)

print(f"{model_full.nq=}")           # 完整模型配置维度
print(f"{model_mimic_from_urdf.nq=}")  # mimic 模型配置维度（更小）

# ---- 方法 3：手工设置 mimic 关系 ----
# transformJointIntoMimic(model, joint_mimicked, joint_mimic, scaling, offset)
#   joint_mimicked：被模仿的主关节 ID
#   joint_mimic：   跟随运动的从关节 ID
#   scaling：       传动比（负号表示方向相反，如夹爪对开）
#   offset：        固定偏置（这里为 0）
model_mimic = pinocchio.transformJointIntoMimic(model_full, 9, 10, -1.0, 0.0)
model_mimic = pinocchio.transformJointIntoMimic(model_mimic, 18, 19, -1.0, 0.0)

# 验证手工设置与 URDF 解析结果一致
print(f"{(model_mimic_from_urdf == model_mimic)=}")  # True

# ---- 构造传动矩阵 G ----
# G ∈ R^{nv_mimic × nv_full}
#   对非 mimic 关节：G[i, i] = 1（恒等映射）
#   对 mimic 关节：  G[mimic_idx, i] = scaling（传动比）
G = np.zeros([model_mimic.nv, model_full.nv])
for i in range(model_full.njoints):
    mimic_scaling = getattr(model_mimic.joints[i].extract(), "scaling", None)
    if mimic_scaling:
        # mimic 关节：速度映射为 scaling 倍
        G[model_mimic.joints[i].idx_v, model_full.joints[i].idx_v] = mimic_scaling
    else:
        # 普通关节：恒等映射
        G[model_mimic.joints[i].idx_v, model_full.joints[i].idx_v] = 1
print("G = ")
print(np.array_str(G, precision=0, suppress_small=True, max_line_width=80))

# ---- 随机测试数据 ----
q_mimic = pinocchio.neutral(model_mimic)
v_mimic = np.random.random(model_mimic.nv)
a_mimic = np.random.random(model_mimic.nv)

# 利用 G 将 mimic 配置展开为完整模型配置
# 注意：offset=0 时 q_full = G^T * q_mimic 成立；offset≠0 需额外处理偏置项
q_full = G.transpose() @ q_mimic
v_full = G.transpose() @ v_mimic
a_full = G.transpose() @ a_mimic

data_mimic = model_mimic.createData()
data_full  = model_full.createData()

# ---- 验证力矩一致性：τ_mimic = G * τ_full ----
# RNEA 在两个模型上分别计算
tau_mimic = pinocchio.rnea(model_mimic, data_mimic, q_mimic, v_mimic, a_mimic)
tau_full  = pinocchio.rnea(model_full,  data_full,  q_full,  v_full,  a_full)
# 应为 True：mimic 模型的力矩 = G * 完整模型的力矩（虚功原理）
print(f"{np.allclose(tau_mimic, G @ tau_full)=}")

# ---- 验证质量矩阵一致性：M_mimic = G * M_full * G^T ----
M_mimic = pinocchio.crba(model_mimic, data_mimic, q_mimic)
M_full  = pinocchio.crba(model_full,  data_full,  q_full)
# 传动矩阵的合同变换：广义惯量的坐标变换
print(f"{np.allclose(M_mimic, G @ M_full @ G.transpose())=}")

# ---- mimic 模型的正向动力学（ABA 不直接支持 mimic，用运动方程代替）----
# 非线性项（Coriolis + 重力）：C(q,v)v + g(q)
C_mimic = pinocchio.nle(model_mimic, data_mimic, q_mimic, v_mimic)
# 从运动方程求加速度：τ = M*a + C  →  a = M^{-1}(τ - C)
a_computed = np.linalg.solve(M_mimic, tau_mimic - C_mimic)
# 验证：求得的加速度应与输入 a_mimic 一致
print(f"{np.allclose(a_mimic, a_computed)=}")
