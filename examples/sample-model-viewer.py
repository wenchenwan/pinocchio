from sys import argv

import pinocchio as pin
from numpy import pi
from pinocchio.visualize import GepettoVisualizer, MeshcatVisualizer, RVizVisualizer

# ============================================================
# 内置示例模型的可视化（无 URDF）
# 演示：
#   - buildSampleModelHumanoid：Pinocchio 内置的参数化人形模型
#   - buildSampleGeometryModelHumanoid：对应几何模型（简单几何体，非网格）
#   - 多种查看器切换：Gepetto / Meshcat / RViz
#   - 交互式配置切换（等待用户按回车后切换到新姿态）
# ============================================================

# GepettoVisualizer: -g   本地桌面 3D 查看器
# MeshcatVisualizer: -m   浏览器 3D 查看器
# RVizVisualizer:    -r   ROS RViz 3D 查看器（需要 ROS 环境）
VISUALIZER = None
if len(argv) > 1:
    opt = argv[1]
    if opt == "-g":
        VISUALIZER = GepettoVisualizer
    elif opt == "-m":
        VISUALIZER = MeshcatVisualizer
    elif opt == "-r":
        VISUALIZER = RVizVisualizer
    else:
        raise ValueError("Unrecognized option: " + opt)

# ---- 内置人形模型（无需 URDF）----
# buildSampleModelHumanoid：构造一个参数化人形机器人
#   - 两条腿（各 6 DOF）、两条臂（各 7 DOF）、躯干、头
#   - 浮动基（FreeFlyer），nq 约 35
model = pin.buildSampleModelHumanoid()

# buildSampleGeometryModelHumanoid：为上述模型创建简单几何体（胶囊/球体）
# 这些几何体不是真实网格，只是概念可视化（适合算法开发阶段快速验证）
visual_model    = pin.buildSampleGeometryModelHumanoid(model)
collision_model = visual_model.copy()   # 碰撞和视觉共用同一套简单几何体

q0 = pin.neutral(model)   # 零位配置（所有关节角为 0）

if VISUALIZER:
    viz = VISUALIZER(model, collision_model, visual_model)
    viz.initViewer()
    viz.loadViewerModel()
    viz.display(q0)   # 显示零位

    input("Enter to check a new configuration")

    # ---- 构造展臂站立姿态 ----
    q = q0.copy()
    # 以下关节 ID 对应 buildSampleModelHumanoid 内部定义的关节顺序
    # q[8]  → 右肩（或左肩）绕某轴旋转 90°（手臂抬起）
    # q[14] → 左肩（或右肩）绕某轴旋转 90°
    # q[23] → 右肘关节 -90°
    # q[29] → 左肘关节  90°
    # 注意：具体关节含义需对照 buildSampleModelHumanoid 的实现
    q[8]  = pi / 2
    q[14] = pi / 2
    q[23] = -pi / 2
    q[29] = pi / 2

    viz.display(q)
    input("Press enter to exit...")
