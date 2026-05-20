# Copyright 2023 Inria
# SPDX-License-Identifier: BSD-2-Clause


"""
In this short script, we show how to compute inverse dynamics (RNEA), i.e. the
vector of joint torques corresponding to a given motion.
"""

from pathlib import Path

import numpy as np
import pinocchio as pin

# 模型目录，替换为你自己的 URDF 路径
pinocchio_model_dir = Path(__file__).parent.parent / "models/"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir = pinocchio_model_dir
urdf_filename = "ur5_robot.urdf"
urdf_model_path = model_path / "ur_description/urdf/" / urdf_filename

# buildModelsFromUrdf 同时返回三个模型：
#   model:           运动学/动力学模型（关节拓扑 + 惯量参数）
#   collision_model: 碰撞几何模型（用于碰撞检测，基于 hpp-fcl）
#   visual_model:    视觉几何模型（用于可视化渲染）
# 这里只用运动学模型，忽略另外两个
model, _, _ = pin.buildModelsFromUrdf(urdf_model_path, package_dirs=mesh_dir)

# 创建算法数据容器
data = model.createData()

# 随机采样关节状态（实际使用时替换为传感器读数）
q = pin.randomConfiguration(model)  # 关节角 [rad]，满足关节限位约束
v = np.random.rand(model.nv, 1)     # 关节速度 [rad/s]
a = np.random.rand(model.nv, 1)     # 期望关节加速度 [rad/s²]

# RNEA（递推牛顿-欧拉算法）：逆动力学
# 给定 (q, v, a)，计算产生该运动所需的关节力矩 τ
# 等价方程：τ = M(q)·a + C(q,v)·v + g(q)
# 算法复杂度：O(n)，远优于显式构造 M, C, g 矩阵再相乘（O(n²)~O(n³)）
# 结果同时存储在返回值和 data.tau 中
tau = pin.rnea(model, data, q, v, a)

# 打印关节力矩（单位：N·m）
print("Joint torques: " + str(tau))
