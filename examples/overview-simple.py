import pinocchio

# 构建一个内置的 6-DOF 串联机械臂模型，无需 URDF 文件
# 该模型等效于标准 6 轴工业臂的运动学/动力学结构
model = pinocchio.buildSampleModelManipulator()

# 创建算法数据容器 Data，与 model 配对
# model 存储只读的拓扑和惯量参数；data 存储每次算法调用的中间结果
# 多线程场景下，多个线程共享同一个 model，每个线程持有独立的 data
data = model.createData()

# 生成零位配置（关节角全部为 0）
# neutral() 对旋转关节返回 0，对四元数关节（如浮动基）返回单位四元数
q = pinocchio.neutral(model)

# 关节速度和加速度全为 0
v = pinocchio.utils.zero(model.nv)
a = pinocchio.utils.zero(model.nv)

# 调用 RNEA（递推牛顿-欧拉算法）计算逆动力学
# 给定 (q, v, a)，返回产生该运动所需的关节力矩 τ
# 公式：τ = M(q)a + C(q,v)v + g(q)
# 此处 v=0, a=0，故结果只反映重力补偿项 g(q)
tau = pinocchio.rnea(model, data, q, v, a)
print("tau = ", tau.T)
