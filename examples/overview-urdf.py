from pathlib import Path
from sys import argv

import pinocchio

# Pinocchio 内置模型目录，实际使用时替换为自己的路径
pinocchio_model_dir = Path(__file__).parent.parent / "models"

# URDF 路径：默认使用 UR5，也可通过命令行参数传入自定义路径
# 对于人形机器人，将此处替换为自己的人形 URDF 路径
# 注意：人形机器人还需要在 buildModelFromUrdf 中追加 JointModelFreeFlyer()
urdf_filename = (
    pinocchio_model_dir / "example-robot-data/example-robot-data/robots/ur_description/urdf/ur5_robot.urdf"
    if len(argv) < 2
    else argv[1]
)

# 从 URDF 构建运动学/动力学模型（固定基版本）
# 人形/四足浮动基版本：pinocchio.buildModelFromUrdf(urdf, pinocchio.JointModelFreeFlyer())
model = pinocchio.buildModelFromUrdf(urdf_filename)
print("model name: " + model.name)

# 创建算法数据容器，与 model 绑定
# model.nq: 配置向量维度（关节角 + 四元数自由度）
# model.nv: 速度/力矩向量维度（李代数维度，总比 nq 小 1 per SO(3) 关节）
data = model.createData()

# 在关节限位内随机采样一个合法配置
# 旋转关节: 在 [lowerPositionLimit, upperPositionLimit] 内均匀采样
# 四元数关节（FreeFlyer）: 在 SO(3) 上均匀采样并自动归一化
q = pinocchio.randomConfiguration(model)
print(f"q: {q.T}")

# 正向运动学：递推计算所有关节在世界坐标系下的 SE(3) 变换
# 结果存入 data.oMi[i]，其中 i 为关节索引（0=宇宙/固定基, 1..njoints-1=实际关节）
# 公式：oMi[i] = oMi[parent(i)] * jointPlacement[i] * jointTransform(q_i)
pinocchio.forwardKinematics(model, data, q)

# 打印每个关节在世界坐标系下的平移坐标（即关节原点的 3D 位置）
# oMi[i].translation 是 3×1 向量，oMi[i].rotation 是 3×3 旋转矩阵
for name, oMi in zip(model.names, data.oMi):
    print("{:<24} : {: .2f} {: .2f} {: .2f}".format(name, *oMi.translation.T.flat))
