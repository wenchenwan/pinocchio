#!/usr/bin/env python
#
# Copyright (C) 2023 Inria

import numpy as np

# robot_descriptions 是一个第三方包，提供常见机器人 URDF 的快捷加载接口
# 安装：pip install robot_descriptions
try:
    from robot_descriptions.loaders.pinocchio import load_robot_description
except ModuleNotFoundError:
    print("This example loads robot descriptions from robot_descriptions.py:")
    print("\n\tpip install robot_descriptions")

print("Goal: load a legged robot from its URDF, then modify leg lengths")

print("Loading robot description from URDF...")
# 加载 Upkie 轮腿机器人（2 条腿，每条腿 3 DOF：髋、膝、轮）
robot = load_robot_description("upkie_description")
model = robot.model

# Upkie 腿长 = 关节 Y 轴 translation + 固定偏移
# 这些偏移值从 URDF 机构分析中得出
known_offsets = {"knee": 0.072, "wheel": 0.065}


def check_limb_lengths(limb_length: float) -> bool:
    """验证模型中所有腿的关节 Y 轴偏移是否与期望腿长一致"""
    print(f"Checking that limbs are {limb_length} m long... ", end="")
    for side in ("left", "right"):
        for joint in ("knee", "wheel"):
            joint_id = model.getJointId(f"{side}_{joint}")
            # jointPlacements[i]：关节 i 相对父关节的固定 SE(3) 变换
            # .translation[1]：Y 轴分量（腿长方向）
            if not np.allclose(
                model.jointPlacements[joint_id].translation[1],
                limb_length - known_offsets[joint],
            ):
                print("{side}_{joint} placement is wrong!")
                return False
    return True


# 验证原始腿长为 0.24 m
if check_limb_lengths(0.24):
    print("OK, the model is as we expect it")


def update_limb_lengths(length: float) -> None:
    """
    直接修改模型中的关节位置参数，更新腿长。
    这是 Pinocchio 模型可参数化修改的示例：
    model 加载后仍可在 Python 中直接修改关节参数，无需重新解析 URDF。
    适用场景：腿长参数辨识、设计优化、仿真中动态改变机器人结构。
    """
    for side in ("left", "right"):
        for joint in ("knee", "wheel"):
            joint_id = model.getJointId(f"{side}_{joint}")
            # 直接写入关节 Y 轴偏移（= 新腿长 - 固定结构偏移）
            model.jointPlacements[joint_id].translation[1] = (
                length - known_offsets[joint]
            )


# 将腿长更新为 0.30 m（比原始 0.24 m 长 6 cm）
update_limb_lengths(0.3)
if check_limb_lengths(0.3):
    print("OK, the update worked!")
