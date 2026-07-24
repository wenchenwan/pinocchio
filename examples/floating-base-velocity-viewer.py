from time import sleep

import hppfcl
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

# ============================================================
# 浮动基/球关节的速度分量可视化
# 演示不同根关节类型（FreeFlyer、SphericalZYX、Spherical）下：
#   - nq（配置空间维度）与 nv（速度空间维度）的差异
#   - 各速度分量对应的运动（通过 integrate 积分 Δt=0.2s 后可视化）
#
# 核心概念：
#   FreeFlyer：nq=7（位置3 + 四元数4）, nv=6（速度3 + 角速度3）
#   SphericalZYX：nq=4（四元数）, nv=3（Z-Y-X 欧拉角速度）
#   Spherical：nq=4（四元数）, nv=3（角速度，李代数表示）
#   → nq > nv 是 Lie 群约束的直接体现（配置在流形上，速度在切空间）
# ============================================================


def create_pin_cube_model(j0="freeflyer"):
    """
    创建一个单关节立方体模型，根关节类型可选。
    用于演示不同根关节的配置/速度空间维度差异。
    """
    model = pin.Model()

    if j0 == "freeflyer":
        # FreeFlyer：6-DOF 浮动基（平移 3 + 旋转 3）
        # nq=7（含四元数归一化约束），nv=6（切空间速度）
        j0 = pin.JointModelFreeFlyer()
    elif j0 == "spherical":
        # Spherical：3-DOF 球铰（纯旋转，无平移）
        # nq=4（四元数），nv=3（角速度李代数）
        j0 = pin.JointModelSpherical()
    elif j0 == "sphericalzyx":
        # SphericalZYX：3-DOF 球铰，用 Z-Y-X 欧拉角参数化
        # nq=4（四元数），nv=3（欧拉角速度）
        j0 = pin.JointModelSphericalZYX()
    else:
        raise ValueError("Unknown joint type")

    # 添加根关节（父 ID=0 → 世界坐标系）
    jointCube = model.addJoint(0, j0, pin.SE3.Identity(), "joint0")
    M = pin.SE3.Identity()
    # 绑定一个长方体惯量（0.8m × 0.4m × 0.2m，质量=1 kg）
    model.appendBodyToJoint(jointCube, pin.Inertia.FromBox(1, 0.8, 0.4, 0.2), M)
    return model


def create_pin_geometry_cube_model(model):
    """创建立方体的可视化几何模型（红色半透明长方体）"""
    # getFrameId 通过 Frame 名称获取 ID（joint0 对应的 Frame）
    jointCube  = model.getFrameId("joint0")
    geom_model = pin.GeometryModel()
    cube_shape = hppfcl.Box(0.8, 0.4, 0.2)   # 与惯量参数匹配的几何形状
    cube = pin.GeometryObject(
        "cube_shape", 0, jointCube, cube_shape, pin.SE3.Identity()
    )
    cube.meshColor = np.array([1.0, 0.1, 0.1, 0.5])   # 半透明红色
    geom_model.addGeometryObject(cube)
    return geom_model


def create_model(joint0_type="freeflyer"):
    model      = create_pin_cube_model(joint0_type)
    geom_model = create_pin_geometry_cube_model(model)
    return model, geom_model


def pin_step(model, vizer, v_index_increment, dt=0.1):
    """
    在速度向量的第 v_index_increment 维施加单位速度，
    积分 dt 秒后显示新配置，直观感受该速度分量的含义。
    """
    q = pin.neutral(model)   # 零位配置
    v = np.zeros(model.nv)
    v[v_index_increment] += 1.0   # 只在指定分量施加速度

    vizer.display(q)
    print(f"{q=}, {v=}")
    sleep(1)

    # 在 Lie 群流形上积分：q_next = Exp(q, v*dt)
    # 对于 FreeFlyer：v=[vx,vy,vz,ωx,ωy,ωz]，各分量独立测试
    # 对于 Spherical：v=[ωx,ωy,ωz]，各对应绕该轴的旋转
    q_next = pin.integrate(model, q, v * dt)
    print(f"{q_next=}")
    return q_next


if __name__ == "__main__":
    # 依次测试三种根关节类型
    list_joint0 = ["freeflyer", "sphericalzyx", "spherical"]

    for joint0_name in list_joint0:
        print(f"\njoint0_name = {joint0_name}")

        model, geom_model = create_model(joint0_name)
        # 打印 nq/nv 的差异（体现 Lie 群约束）
        print(f"{model.nq=}")   # FreeFlyer: 7, Spherical/SphericalZYX: 4
        print(f"{model.nv=}")   # FreeFlyer: 6, Spherical/SphericalZYX: 3

        q = pin.neutral(model)
        print(f"neutral q configuration: {q=}")
        # FreeFlyer 零位：[0,0,0, 0,0,0,1]（位置=0, 四元数=单位）
        # Spherical 零位：[0,0,0,1]（单位四元数）

        vizer = MeshcatVisualizer(model, geom_model, geom_model)
        vizer.initViewer(open=True, loadModel=True)
        vizer.display(q)
        sleep(2)

        # 逐个测试每个速度分量的效果
        for i in range(model.nv):
            print(f"v[{i}] += 1")
            q_next = pin_step(model, vizer, i, dt=0.2)
            vizer.display(q_next)
            sleep(2)
