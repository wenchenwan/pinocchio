import numpy as np
import pinocchio as pin

# ============================================================
# 并行批量运算：ModelPool + OpenMP 多线程
# 应用场景：
#   - 强化学习采样：同时对 128 个随机配置计算动力学
#   - 轨迹优化：并行计算一批轨迹节点的 RNEA/ABA
#   - Monte Carlo 仿真：批量随机配置下的力矩分布统计
#
# 核心思想：
#   model 本身是只读的（线程安全）
#   data 是可写的（每个线程需要独立 data）
#   ModelPool 内部维护 num_threads 个 data 副本，零拷贝分发
# ============================================================

# 构建标准人形模型（30+ DOF，包含浮动基）
model = pin.buildSampleModelHumanoid()
# 浮动基的前 7 维配置（位置 3D + 四元数 4D）没有上下限，需手工设置
# 否则 randomConfiguration 无法为该段生成有效随机值
model.lowerPositionLimit[:7] = -np.ones(7)
model.upperPositionLimit[:7] = +np.ones(7)

# ---- ModelPool：多线程共享模型的线程安全容器 ----
# ModelPool 内部为每个线程预分配一个独立的 Data 对象
# 线程间共享只读的 Model，各自使用私有 Data，避免数据竞争
pool = pin.ModelPool(model)

# 查询当前系统可用的最大 OpenMP 线程数
# 通常等于逻辑 CPU 核数（受 OMP_NUM_THREADS 环境变量影响）
num_threads = pin.omp_get_max_threads()
batch_size = 128   # 一批并行计算的样本数量

# ---- 批量输入：shape = (nq/nv, batch_size)，列优先存储 ----
# 每列对应一个独立的机器人配置（或速度/加速度）
q = np.empty((model.nq, batch_size))
for k in range(batch_size):
    q[:, k] = pin.randomConfiguration(model)

# v=0, a=0 的批量零向量（用于演示；实际可替换为非零批量数据）
v   = np.zeros((model.nv, batch_size))
a   = np.zeros((model.nv, batch_size))
tau = np.zeros((model.nv, batch_size))

print(f"num_threads: {num_threads}")
print(f"batch_size: {batch_size}")

# ---- 并行 RNEA（逆动力学）----
# rneaInParallel(num_threads, pool, q, v, a[, res])
# → 对每列 (q[:,k], v[:,k], a[:,k]) 并行调用 RNEA
# → 结果写入 res[:,k] = τ_k ∈ R^nv

# 形式 1：预分配输出缓冲区（避免内部分配，适合实时/高频调用）
res_rnea = np.empty((model.nv, batch_size))
pin.rneaInParallel(num_threads, pool, q, v, a, res_rnea)  # 就地写入 res_rnea

# 形式 2：内部自动分配并返回（简洁，适合原型/调试）
res_rnea2 = pin.rneaInParallel(num_threads, pool, q, v, a)

# ---- 并行 ABA（正向动力学）----
# abaInParallel(num_threads, pool, q, v, tau[, res])
# → 对每列 (q[:,k], v[:,k], tau[:,k]) 并行调用 ABA
# → 结果写入 res[:,k] = q̈_k ∈ R^nv

# 形式 1：预分配输出缓冲区
res_aba = np.empty((model.nv, batch_size))
pin.abaInParallel(num_threads, pool, q, v, tau, res_aba)  # 就地写入

# 形式 2：内部分配并返回
res_aba2 = pin.abaInParallel(num_threads, pool, q, v, tau)
