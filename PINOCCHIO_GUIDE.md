# Pinocchio 技术指南

> 面向人形机器人方向的完整参考：数学原理 · 架构设计 · 全部 Demo 注释

---

## 目录

1. [项目定位](#1-项目定位)
2. [核心数学基础](#2-核心数学基础)
   - [2.1 李群 SE(3)](#21-李群-se3)
   - [2.2 空间速度与空间力](#22-空间速度与空间力)
     - [2.2.1 空间速度（Twist）](#221-空间速度twist)
     - [2.2.2 空间力（Wrench）](#222-空间力wrench)
     - [2.2.3 为什么是对偶空间：功率不变性](#223-为什么是对偶空间功率不变性)
     - [2.2.4 坐标变换：伴随与余伴随](#224-坐标变换伴随与余伴随)
     - [2.2.5 空间加速度：最反直觉的一个](#225-空间加速度最反直觉的一个)
     - [2.2.6 运动方程中的 ω× 项](#226-运动方程中的-omegatimes-项)
     - [2.2.7 参考系：LOCAL / WORLD / LOCAL_WORLD_ALIGNED](#227-参考系local--world--local_world_aligned)
     - [2.2.8 速查对照表](#228-速查对照表)
     - [2.2.9 为什么值得接受这套形式](#229-为什么值得接受这套形式)
   - [2.3 空间惯量](#23-空间惯量)
3. [Model / Data 架构](#3-model--data-架构)
   - [3.1 Joint 与 Frame：oMi 与 oMf 的区别和联系](#31-joint-与-frameomi-与-omf-的区别和联系)
4. [核心算法原理](#4-核心算法原理)
   - [4.1 正向运动学 FK](#41-正向运动学-fk)
     - [4.1.1 Joint 层 vs Frame 层](#411-joint-层-vs-frame-层)
     - [4.1.2 Frame 位姿更新的三个 API](#412-frame-位姿更新的三个-api)
     - [4.1.3 两个必须避开的陷阱](#413-两个必须避开的陷阱)
     - [4.1.4 选用建议](#414-选用建议)
     - [4.1.5 配套的 Frame 层 API](#415-配套的-frame-层-api)
   - [4.2 逆向运动学 IK](#42-逆向运动学-ik)
     - [4.2.1 为什么误差不能写成"目标减当前"](#421-为什么误差不能写成目标减当前)
     - [4.2.2 关键推导：为什么需要 Jlog6](#422-关键推导为什么需要-jlog6)
     - [4.2.3 阻尼最小二乘的来历](#423-阻尼最小二乘的来历)
     - [4.2.4 流形积分](#424-流形积分)
     - [4.2.5 完整算法](#425-完整算法)
   - [4.3 逆向动力学 RNEA](#43-逆向动力学-rnea)
     - [4.3.1 从拉格朗日方程到递推形式](#431-从拉格朗日方程到递推形式)
     - [4.3.2 第一趟：正向传播运动学](#432-第一趟正向传播运动学根--叶)
     - [4.3.3 第二趟：反向传播力](#433-第二趟反向传播力叶--根)
     - [4.3.4 由 RNEA 派生的三个常用量](#434-由-rnea-派生的三个常用量)
   - [4.4 正向动力学 ABA](#44-正向动力学-aba)
     - [4.4.1 为什么不能直接求逆](#441-为什么不能直接求逆)
     - [4.4.2 核心概念：关节化体惯量](#442-核心概念关节化体惯量)
     - [4.4.3 递推公式的推导](#443-递推公式的推导)
     - [4.4.4 三趟结构](#444-三趟结构)
   - [4.5 质量矩阵 CRBA](#45-质量矩阵-crba)
     - [4.5.1 从动能出发](#451-从动能出发)
     - [4.5.2 稀疏性：为什么大量元素恒为零](#452-稀疏性为什么大量元素恒为零)
     - [4.5.3 复合刚体惯量与递推](#453-复合刚体惯量与递推)
     - [4.5.4 性质与后续使用](#454-性质与后续使用)
   - [4.6 质心动量学](#46-质心动量学)
     - [4.6.1 定义：把全身动量约化到质心](#461-定义把全身动量约化到质心)
     - [4.6.2 为什么它是平衡控制的基础](#462-为什么它是平衡控制的基础)
     - [4.6.3 API](#463-api)
   - [4.7 接触约束动力学](#47-接触约束动力学)
     - [4.7.1 约束的三个层次](#471-约束的三个层次位置--速度--加速度)
     - [4.7.2 从 Gauss 最小约束原理导出 KKT](#472-从-gauss-最小约束原理导出-kkt)
     - [4.7.3 Delassus 矩阵与求解](#473-delassus-矩阵与求解)
     - [4.7.4 约束类型与 Baumgarte 稳定化](#474-约束类型与-baumgarte-稳定化)
   - [4.8 解析梯度](#48-解析梯度)
     - [4.8.1 RNEA 梯度](#481-rnea-梯度)
     - [4.8.2 ABA 梯度与链式法则](#482-aba-梯度与链式法则)
     - [4.8.3 与 DDP 的状态空间矩阵的关系](#483-与-ddp-的状态空间矩阵的关系)
     - [4.8.4 运动学梯度](#484-运动学梯度)
5. [Demo 注释：Python](#5-demo-注释python)
6. [Demo 注释：C++](#6-demo-注释c)
7. [人形机器人开发路径](#7-人形机器人开发路径)
8. [与工业机器人的对比](#8-与工业机器人的对比)
9. [关键头文件索引](#9-关键头文件索引)
10. [源码目录地图](#10-源码目录地图)
11. [源码架构三大机制](#11-源码架构三大机制)
    - [11.1 三层文件布局 .hpp / .hxx / .cpp](#111-三层文件布局-hpp--hxx--cpp)
    - [11.2 Tpl 模板 + context 默认标量](#112-tpl-模板--context-默认标量)
    - [11.3 CRTP + Boost.Fusion 访问者](#113-crtp--boostfusion-访问者)
    - [11.4 第三个模板参数 JointCollectionTpl：关节类型目录](#114-第三个模板参数-jointcollectiontpl关节类型目录)
12. [实战：追踪一个算法从 API 到实现](#12-实战追踪一个算法从-api-到实现)
13. [Python 绑定如何映射到 C++](#13-python-绑定如何映射到-c)
14. [源码阅读推荐顺序](#14-源码阅读推荐顺序)

---

## 1. 项目定位

**Pinocchio** 是 INRIA/LAAS-CNRS 开发的高性能刚体动力学 C++ 模板库（v3.6.0），实现了
Roy Featherstone 体系中最先进的多体动力学算法，并提供所有主要算法的**解析梯度**。

它是以下人形/腿足机器人框架的核心计算引擎：


| 框架                      | 用途               |
| ------------------------- | ------------------ |
| **Crocoddyl**             | DDP/iLQR 轨迹优化  |
| **Humanoid Path Planner** | 人形运动规划       |
| **Stack-of-Tasks**        | 层次全身控制 (WBC) |
| **MPC in Robotics**       | 实时模型预测控制   |

**性能特点**：C++ 模板展开 + 缓存友好数据布局，单台笔记本（Intel i7 @ 2.4 GHz）可
在 1 ms 内完成百自由度机器人的动力学计算，支持 OpenMP 并行批量计算。

---

## 2. 核心数学基础

### 2.1 李群 SE(3)

刚体位姿不是普通向量，而属于**特殊欧氏群**（Special Euclidean Group）：

$$
SE(3) = \left\{ T = \begin{bmatrix} R & p \\ 0 & 1 \end{bmatrix} \;\middle|\; R \in SO(3),\; p \in \mathbb{R}^3 \right\}
$$

其中 $R$ 为旋转矩阵（满足 $R^\top R = I,\;\det R = 1$），$p$ 为平移向量。

**关键运算：**


| 运算                           | 数学含义         | Pinocchio API                         |
| ------------------------------ | ---------------- | ------------------------------------- |
| $T_1 \cdot T_2$                | 变换复合         | `T1 * T2`                             |
| $T^{-1}$                       | 逆变换           | `T.inverse()`                         |
| $T^{-1} \cdot T_{des}$         | 相对变换（误差） | `T.actInv(T_des)`                     |
| $\log(T) \in \mathfrak{se}(3)$ | 对数映射→切空间 | `pinocchio.log6(T)`                   |
| $\exp(\xi) \in SE(3)$          | 指数映射→群     | `pinocchio.exp6(xi)`                  |
| $q \oplus v\Delta t$           | 流形上积分       | `pinocchio.integrate(model, q, v*dt)` |
| $q_1 \ominus q_2$              | 流形上差分       | `pinocchio.difference(model, q1, q2)` |

**为什么需要流形积分？**
人形机器人的浮动基旋转部分以四元数表示（单位球面 $S^3$），不能直接做 $q \mathrel{+}= \dot{q}\Delta t$。
必须在李群上做指数映射才能保持旋转的约束性。这是人形机器人与固定基工业机器人
最本质的数学差异。

**状态空间维度：**

$$
n_q = n_{joints} + 7 \quad (\text{4 元数旋转 } + \text{ 3 位移})
$$

$$
n_v = n_{joints} + 6 \quad (\text{切空间维度，速度/加速度用此维度})
$$

### 2.2 空间速度与空间力

Pinocchio 采用 Featherstone **空间代数（Spatial Algebra）** 统一描述运动与力。

> **核心区别**：经典力学在"某个物体的质心 / 某个具体点"上定义速度和力；空间代数在
> **整个刚体的运动场**上定义，参考点固定为**坐标系原点**。理解这一句，下面所有公式
> 都会变得自然。

#### 2.2.1 空间速度（Twist）

刚体上任意一点 $P$ 的速度不是独立的，它由一个**速度场**决定。设刚体角速度为 $\omega$，
刚体上参考点 $O$ 的速度为 $v_O$，则任意点 $P$ 的速度为：

$$
v_P = v_O + \omega \times \overrightarrow{OP}
$$

关键观察：**$\omega$ 与参考点无关（刚体的固有属性），而 $v$ 依赖于参考点的选择。**
因此完整描述刚体运动只需 6 个数：$(v_O,\ \omega)$，其中 $v_O$ 是**坐标系原点处的速度**。

$$
\nu = \begin{bmatrix} v \\ \omega \end{bmatrix} \in \mathbb{R}^6, \qquad \nu \in \mathfrak{se}(3)
$$

> ⚠️ **排列顺序**：Pinocchio 采用**线性在前、角速度在后**（Featherstone 原书相反）。
> `Motion` 的内存布局中 `linear()` 占索引 0–2，`angular()` 占索引 3–5。
> 本文档统一使用 Pinocchio 的约定。

```cpp
pinocchio::Motion v;
v.linear()   // v ∈ R³：坐标系原点处的线速度
v.angular()  // ω ∈ R³：刚体角速度
```

**最大的陷阱：`v.linear()` 不是质心速度。**

设关节 $i$ 的局部坐标系原点为 $O_i$，连杆质心为 $C$，$c = \overrightarrow{O_iC}$：

$$
\underbrace{v_C}_{\text{质心速度}} = \underbrace{v_{O_i}}_{\texttt{v.linear()}} + \omega \times c
$$

当坐标系原点不在质心时，`v.linear()` 是"刚体延拓后，恰好经过坐标系原点的那个**虚拟
质点**的速度"，而不是任何实际物质点的速度。举例：一个绕自身固定轴自转的轮子，坐标系
原点放在轴心上，此时 $v_O = 0$ 而 $\omega \neq 0$ —— 虽然轮缘每点都在运动，但空间速度
的线性部分为零。

#### 2.2.2 空间力（Wrench）

作用在刚体上的一组力 $\{f_k\}$（作用点 $P_k$），可等效约化到坐标系原点 $O$：

$$
f = \sum_k f_k, \qquad \tau_O = \sum_k \overrightarrow{OP_k} \times f_k
$$

同样：**合力 $f$ 与约化点无关，而力矩 $\tau$ 依赖于约化点。**

$$
\phi = \begin{bmatrix} f \\ \tau_O \end{bmatrix} \in \mathbb{R}^6, \qquad \phi \in \mathfrak{se}(3)^*
$$

```cpp
pinocchio::Force f;
f.linear()   // f ∈ R³：合力（与参考点无关）
f.angular()  // τ_O ∈ R³：对坐标系原点的合力矩（随原点变化）
```

**陷阱**：`f.angular()` 不是"纯力偶"，它包含两部分贡献——真正的力偶 + 合力对原点的
力矩臂效应。换坐标系原点时 `f.angular()` 会变，`f.linear()` 不变。

#### 2.2.3 为什么是对偶空间：功率不变性

这是空间代数最本质的设计动机。**功率是标量，物理上与坐标系选择无关**：

$$
P = \phi^\top \nu = f \cdot v_O + \tau_O \cdot \omega
$$

验证其与参考点无关。换到新原点 $O'$，记 $r = \overrightarrow{OO'}$：

$$
v_{O'} = v_O + \omega\times r, \qquad \tau_{O'} = \tau_O - r\times f
$$

$$
P' = f\cdot(v_O + \omega\times r) + (\tau_O - r\times f)\cdot\omega
     = f\cdot v_O + \tau_O\cdot\omega + \underbrace{f\cdot(\omega\times r) - (r\times f)\cdot\omega}_{=\,0\ \text{（混合积轮换恒等式）}} = P
$$

因为功率必须不变，当 $\nu$ 按某规则变换时 $\phi$ 必须按其**逆转置**变换 —— 这正是
"力生活在速度的对偶空间"的准确含义，也解释了为什么 Pinocchio 中变换速度用伴随
（Adjoint），变换力用余伴随（Co-Adjoint，即伴随的转置）。

#### 2.2.4 坐标变换：伴随与余伴随

设 ${}^bM_a = (R,\ p) \in SE(3)$ 是从坐标系 $a$ 到 $b$ 的变换。

**速度变换（伴随作用 Adjoint）**：

$$
{}^b\nu = \mathrm{Ad}_{{}^bM_a}\,{}^a\nu, \qquad
\mathrm{Ad}_{M} = \begin{bmatrix} R & \hat{p}R \\ 0 & R \end{bmatrix}
$$

展开成经典形式（这两行是理解全部的关键）：

$$
{}^bv = R\,{}^av + p\times(R\,{}^a\omega), \qquad {}^b\omega = R\,{}^a\omega
$$

角速度只旋转；线速度旋转后还要加上**因参考点平移产生的牵连项** $p\times\omega$。

**力变换（余伴随作用 Co-Adjoint）**：

$$
{}^b\phi = \mathrm{Ad}_{{}^aM_b}^\top\,{}^a\phi
= \begin{bmatrix} R & 0 \\ \hat{p}R & R \end{bmatrix}{}^a\phi
$$

展开：

$$
{}^bf = R\,{}^af, \qquad {}^b\tau = R\,{}^a\tau + p\times(R\,{}^af)
$$

注意与速度变换的**对称性**：速度中 $\hat p$ 位于右上角（作用在 $\omega$ 上），力中
$\hat p$ 位于左下角（作用在 $f$ 上）—— 这正是矩阵转置的结果，是对偶性的直接体现。
物理解释也完全对应经典力学中你熟悉的**力的平移定理**：力平移后要加力矩臂效应 $p\times f$。

```cpp
Motion v_b = bMa.act(v_a);      // 速度：伴随 Ad_M
Motion v_a = bMa.actInv(v_b);   // 逆变换 Ad_{M⁻¹}
Force  f_b = bMa.act(f_a);      // 力：Pinocchio 自动使用余伴随
Force  f_a = bMa.actInv(f_b);
```

#### 2.2.5 空间加速度：最反直觉的一个

$$
a = \dot\nu = \begin{bmatrix} \dot v_O \\ \dot\omega \end{bmatrix}
$$

**`a.linear()` 不是质心的经典加速度，甚至不等于原点物质点的经典加速度**：

$$
a_{\text{classic}}(O) = \dot v_O + \omega\times v_O
$$

多出来的 $\omega\times v_O$ 项，源于空间加速度定义为空间速度的**逐分量时间导数**，而
经典加速度是物质点位置的二阶导。这个差异在 RNEA 中被系统性地吸收进递推公式，所以
日常使用不必手工处理，但读源码时必须知道。

```cpp
data.a[i]                                 // 空间加速度（李代数导数）
pinocchio::getClassicalAcceleration(...)  // 经典加速度（含 ω×v 修正）
```

#### 2.2.6 运动方程中的 $\omega\times$ 项

Newton-Euler 方程在空间代数下压缩为一行：

$$
\phi = I\,a + \nu \times^* (I\,\nu)
$$

其中 $\times^*$ 是力的叉乘算子（`Motion::cross` 的对偶）。展开即经典形式：

$$
f = m\,a_C, \qquad \tau_C = \bar I_C\,\dot\omega + \underbrace{\omega\times(\bar I_C\,\omega)}_{\text{陀螺项，来自 }\nu\times^*}
$$

$\nu \times^* (I\nu)$ 就是所有离心力 / 科氏力 / 陀螺力矩的统一来源，也是 RNEA 中
`data.f[i]` 递推的核心。

#### 2.2.7 参考系：LOCAL / WORLD / LOCAL_WORLD_ALIGNED

同一个物理量的三种表达。设关节 $i$ 的位姿 ${}^oM_i = (R,\ p)$：


| 枚举值                | 原点         | 坐标轴       | 变换                             |
| --------------------- | ------------ | ------------ | -------------------------------- |
| `LOCAL`               | 关节$i$ 原点 | 关节$i$ 的轴 | ${}^i\nu$（原生存储）            |
| `WORLD`               | **世界原点** | 世界轴       | $\mathrm{Ad}_{{}^oM_i}\,{}^i\nu$ |
| `LOCAL_WORLD_ALIGNED` | 关节$i$ 原点 | 世界轴       | $\mathrm{Ad}_{(R,\,0)}\,{}^i\nu$ |

$$
\nu^{\text{WORLD}} = \begin{bmatrix} Rv + p\times(R\omega) \\ R\omega \end{bmatrix},
\qquad
\nu^{\text{LWA}} = \begin{bmatrix} Rv \\ R\omega \end{bmatrix}
$$

关键理解：

- `WORLD` 的线速度分量含 $p\times\omega$，是"**世界原点处虚拟点**的速度"——尽管名字
  听起来最自然，它几乎不是你想要的物理量。
- `LOCAL_WORLD_ALIGNED` 的线速度分量才是**末端执行器那一点的真实速度在世界轴下的分量**。

**这就是为什么接触约束几乎总是用 `LOCAL_WORLD_ALIGNED`**（见
[`contact-cholesky.py`](examples/contact-cholesky.py)、
[`anymal-simulation.py`](examples/anymal-simulation.py) 中的 `CONTACT_3D`）：足端法向力
沿世界 $Z$ 轴，且接触点速度必须是该点的真实速度（接触约束要求它为零）。

Jacobian 同样遵循这套规则（$\nu = J(q)\,\dot q$）：
`getFrameJacobian(..., LOCAL_WORLD_ALIGNED)` 的前 3 行才是你在工业机器人中熟悉的
"几何 Jacobian 的平移部分"。

#### 2.2.8 速查对照表


| 空间代数                | 经典力学       | 关系                                               |
| ----------------------- | -------------- | -------------------------------------------------- |
| `v.angular()` $\omega$  | 角速度         | 完全相同（与参考点无关）                           |
| `v.linear()` $v_O$      | 原点处速度     | $v_C = v_O + \omega\times c$                       |
| `f.linear()` $f$        | 合力           | 完全相同（与参考点无关）                           |
| `f.angular()` $\tau_O$  | 对原点的合力矩 | $\tau_C = \tau_O - c\times f$                      |
| `a.linear()` $\dot v_O$ | ——           | $a_{\text{classic}} = \dot v_O + \omega\times v_O$ |
| $\phi^\top\nu$          | 功率           | 完全相同（标量不变量）                             |
| $\mathrm{Ad}_M$         | 速度平移公式   | $v' = Rv + p\times R\omega$                        |
| $\mathrm{Ad}_M^\top$    | 力的平移定理   | $\tau' = R\tau + p\times Rf$                       |

#### 2.2.9 为什么值得接受这套形式

从工业机器人（DH 参数 + 每个连杆单独写 Newton-Euler）转过来，会觉得空间代数增加了
心智负担。但对人形机器人它是必需的：

1. **递推公式统一**：$\nu_i = \nu_{\lambda(i)} + S_i\dot q_i$，一行覆盖旋转关节、平移
   关节、球关节、浮动基。DH 参数在球关节和浮动基上直接失效。
2. **$O(n)$ 复杂度的基础**：RNEA / ABA / CRBA 依赖这套代数的结构——力和速度的对偶让
   反向递推可以直接用 $\mathrm{Ad}^\top$，无需额外推导。
3. **质心动力学**：$h_G = A_G(q)\dot q$（CMM 矩阵）本质上就是把所有连杆的空间动量用
   $\mathrm{Ad}^\top$ 约化到质心，这是人形平衡控制（ZMP / CoM 轨迹 / Capture Point）
   的数学基础（见 [4.6](#46-质心动量学)）。
4. **接触约束**：$J\dot q = 0$ 在空间代数下自然表达 6D（面接触 `CONTACT_6D`）或
   3D（点接触 `CONTACT_3D`）约束，无需分类讨论。

### 2.3 空间惯量

空间惯量 $I \in \mathbb{R}^{6\times6}$ 连接空间速度与空间动量：$h = I\,\nu$。
在以坐标系原点为参考点（质心偏移 $c$）时，按 Pinocchio 的 $[v;\ \omega]$ 顺序：

$$
I = \begin{bmatrix} m\,\mathbb{1}_3 & -m\hat{c} \\[2pt] m\hat{c} & \bar I_C - m\hat{c}\hat{c} \end{bmatrix}
$$

其中 $m$ 为质量，$c$ 为质心相对坐标系原点的位置，$\bar I_C$ 为**绕质心**的 $3\times3$
转动惯量张量，$\hat c$ 为 $c$ 的反对称矩阵（$\hat c\,x = c\times x$）。

对应的经典表达式：

$$
\underbrace{p_{\text{lin}} = m(v_O + \omega\times c)}_{=\ m\,v_C\text{，质心线动量}},
\qquad
\underbrace{L_O = m\,c\times v_O + (\bar I_C - m\hat c\hat c)\,\omega}_{\text{对原点的角动量（含平行轴项）}}
$$

这里 $-m\hat c\hat c = m(c^\top c\,\mathbb{1} - c\,c^\top)$ 正是**平行轴定理**
（Huygens–Steiner）。当 $c = 0$（原点在质心）时，$I$ 退化为对角块形式
$\mathrm{diag}(m\mathbb{1},\ \bar I_C)$，即回到最熟悉的经典教科书形式。

```cpp
pinocchio::Inertia Y = model.inertias[i];
Y.mass()     // m：连杆质量
Y.lever()    // c：质心相对关节坐标系原点的位置
Y.inertia()  // Ī_C：绕质心的 3×3 转动惯量
```

空间力方程（Newton-Euler）：$\phi = I\,a + \nu\times^*(I\,\nu)$，详见
[2.2.6](#226-运动方程中的-omegatimes-项)。

---

## 3. Model / Data 架构

Pinocchio 强制分离不可变参数（`Model`）与可变中间结果（`Data`）：

```
Model（只读，多线程共享）
├── joints[]         关节类型列表（JointModelRX / FreeFlyer / ...）
├── inertias[]       各连杆空间惯量
├── jointPlacements[] 关节相对父连杆的固定变换
├── names[]          关节名称
├── nq / nv          配置/速度空间维度
└── referenceConfigurations  参考配置（half_sitting 等）

Data（可变，每线程独立一份）
├── oMi[]       世界系下各关节 SE(3) 变换
├── oMf[]       世界系下各 Frame 的 SE(3) 变换
├── v[]         各关节空间速度
├── a[]         各关节空间加速度
├── f[]         各关节空间力
├── tau         关节力矩（RNEA 输出）
├── ddq         关节加速度（ABA 输出）
├── M           质量矩阵（CRBA 输出）
├── Minv        质量矩阵逆（ABA 导数输出）
├── hg          质心动量
├── com[]       质心位置
├── lambda_c    接触力 Lagrange 乘子
└── ...
```

**标准用法模式：**

```python
# 固定基机器人
model = pin.buildModelFromUrdf("robot.urdf")

# 浮动基人形（必须加 FreeFlyer）
model = pin.buildModelFromUrdf("humanoid.urdf", pin.JointModelFreeFlyer())

data = model.createData()   # 每个线程各创建一份

# 算法调用：model 只读，data 存储中间结果
pin.forwardKinematics(model, data, q)
result = data.oMi[joint_id]
```

**多线程**：`ModelPool` 管理多个 `Data` 副本，支持 OpenMP 批量并行。

### 3.1 Joint 与 Frame：oMi 与 oMf 的区别和联系

阅读源码时最容易混淆的一对概念。一句话：**Joint 是会动的结构节点（有自由度），
Frame 是固定挂在关节上的具名标签（无自由度）**。

**Joint（关节）** —— 运动学树里**唯一拥有自由度（DOF）**的实体：

- 产生相对父连杆的运动；`model.joints[i]` 是关节类型，`model.parents[i]` 是父关节。
- 每个关节占用 `nq`/`nv` 维度，**FK / RNEA / ABA 等算法真正遍历的就是关节**。
- 世界位姿存于 `data.oMi[i]`（世界系 → 关节 i），数量 = `model.njoints`。

**Frame（帧）** —— 本质是一个"带标签的静态偏移 `SE3`"，本身**没有自由度**：

- 定义见 [include/pinocchio/src/multibody/frame.hxx](include/pinocchio/src/multibody/frame.hxx)（`FrameTpl`）。
  核心成员：`parentJoint`（挂在哪个关节）、`placement`（相对父关节的**固定偏移**，不变量）、
  `type`、`inertia`。
- `FrameType` 有 5 种：`JOINT` / `FIXED_JOINT` / `BODY` / `OP_FRAME`（用户自定义操作帧，如末端 TCP、
  足底）/ `SENSOR`。
- 世界位姿存于 `data.oMf[f]`，数量 = `model.nframes`（通常远多于关节数）。

> **为什么帧比关节多？** URDF 里的 `fixed` 关节没有自由度，会被"压扁"、不作为关节存在，
> 只以 `FIXED_JOINT` / `BODY` 帧的形式保留。所以典型情况下 **关节数 ≪ 帧数**。
> 当你想**用名字引用末端、工具点、足底**时（`model.getFrameId("left_sole_link")`），用的就是 Frame。

**oMi 与 oMf 的联系（本质）** —— 见 [frames.hxx:34](include/pinocchio/src/algorithm/frames.hxx#L34)：

```cpp
data.oMf[i] = data.oMi[parent] * frame.placement;
```

即 **oMf 就是在 oMi 上再乘一个静态偏移**，是 oMi 的派生量：

$$
{}^oM_f = \underbrace{{}^oM_{i(f)}}_{\texttt{oMi[parent]，随 }q\text{ 变}} \cdot \underbrace{{}^{i}M_f}_{\texttt{frame.placement，不变}}
$$


|          | `oMi`                                   | `oMf`                                               |
| -------- | --------------------------------------- | --------------------------------------------------- |
| 对象     | **关节** i（joint）                     | **帧** f（frame）                                   |
| 数量     | `model.njoints`                         | `model.nframes`（更多）                             |
| 计算来源 | FK 递推`oMi[i] = oMi[parent] * liMi[i]` | 由`oMi` 派生（乘固定偏移）                          |
| 更新函数 | `forwardKinematics(model, data, q)`     | `updateFramePlacements` / `framesForwardKinematics` |

**使用注意（有先后依赖）**：`oMf` 依赖 `oMi`，必须先做 FK 再更新帧：

```python
pin.forwardKinematics(model, data, q)      # → 填充 data.oMi
pin.updateFramePlacements(model, data)     # → 由 data.oMi 计算 data.oMf
oMf = data.oMf[model.getFrameId("left_sole_link")]

# 或一步到位（内部两件事都做）：
pin.framesForwardKinematics(model, data, q)
```

**一句话总结**：关节是"动力学计算的基本单位"（`oMi`），帧是"人类用名字引用末端 / 接触点 /
传感器的便捷标签"（`oMf`），二者由 `oMf = oMi[父关节] × 固定偏移` 相连。

---

## 4. 核心算法原理

### 4.1 正向运动学 FK

**目标**：给定关节配置 $q$，计算所有连杆在世界系中的位姿 ${}^0T_i$。

**数学原理**：沿运动链递推变换：

$$
{}^0T_i = {}^0T_{p(i)} \cdot {}^{p(i)}T_i(q_i)
$$

其中 $p(i)$ 为第 $i$ 关节的父关节，${}^{p(i)}T_i(q_i)$ 由 URDF 参数和关节类型（旋转/移动/球形）决定。

**API**：

```python
pin.forwardKinematics(model, data, q)          # 只更新 data.oMi[]
pin.forwardKinematics(model, data, q, v)       # 同时计算速度 data.v[]
pin.forwardKinematics(model, data, q, v, a)    # 同时计算加速度 data.a[]
```

**复杂度**：$O(n)$，$n$ 为关节数。

#### 4.1.1 Joint 层 vs Frame 层

Pinocchio 的运动学树中**只有关节是计算节点**，上面的递推只填充 `data.oMi[]`。但 URDF
里你真正关心的对象大多**不是关节**：末端执行器 `tool0`、足底接触点 `LF_FOOT`、相机
安装座、IMU 位置、连杆本体（Body）。这些都是 **Frame** —— 相对某个父关节的**固定**
SE(3) 偏移：

$$
\underbrace{{}^oM_f}_{\texttt{data.oMf[f]}}
= \underbrace{{}^oM_{i}}_{\texttt{data.oMi[parentJoint]}}
\cdot \underbrace{{}^{i}M_f}_{\texttt{model.frames[f].placement}}
$$

关键点：**Frame 不引入自由度**，它只是关节坐标系上的一个固定挂载点。因此
$n_{\text{frames}}$ 通常是 $n_{\text{joints}}$ 的 2 倍以上
（`buildSampleModelHumanoid()`：30 个关节 vs 70 个 Frame）。

```python
print(model.njoints, model.nframes)   # 样例人形：30 vs 70
model.frames[fid].parentJoint         # 该 Frame 挂在哪个关节上
model.frames[fid].placement           # ᶦM_f，固定偏移，与 q 无关
model.frames[fid].type                # JOINT / BODY / OP_FRAME / FIXED_JOINT / SENSOR
```

> ⚠️ 特别注意 `FIXED_JOINT` 类型。URDF 中的固定关节在建模时会被**吸收合并**（不产生
> `data.oMi` 项），只以 Frame 形式保留。所以 URDF 里很多"关节名"在 Pinocchio 中其实是
> Frame，必须用 `getFrameId` 而非 `getJointId` 查找。

#### 4.1.2 Frame 位姿更新的三个 API

```python
pin.framesForwardKinematics(model, data, q)     # FK + 刷新全部 Frame（自包含）
pin.updateFramePlacements(model, data)          # 仅刷新全部 Frame（复数，需先 FK）
pin.updateFramePlacement(model, data, frame_id) # 仅刷新单个 Frame（单数，需先 FK）
```


| 函数                                      | 内部是否调 FK       | 更新范围   | 复杂度         | 前置条件 |
| ----------------------------------------- | ------------------- | ---------- | -------------- | -------- |
| `framesForwardKinematics(model, data, q)` | ✅ 是（**仅位姿**） | 全部 Frame | $O(n_j + n_f)$ | 无       |
| `updateFramePlacements(model, data)`      | ❌ 否               | 全部 Frame | $O(n_f)$       | 需先 FK  |
| `updateFramePlacement(model, data, fid)`  | ❌ 否               | 单个 Frame | $O(1)$         | 需先 FK  |

声明在 [`algorithm/frames.hpp`](include/pinocchio/algorithm/frames.hpp)，实现在
[`src/algorithm/frames.hxx`](include/pinocchio/src/algorithm/frames.hxx)（注意路径中的
`src/` 段，这是 Pinocchio 模板库的三层布局约定，见 [11.1](#111-三层文件布局-hpp--hxx--cpp)）。

**`framesForwardKinematics` 是复合操作**（[`frames.hxx:67-68`](include/pinocchio/src/algorithm/frames.hxx#L67-L68)）：

```cpp
forwardKinematics(model, data, q);    // ← 只算位姿，不算速度/加速度
updateFramePlacements(model, data);   // ← 遍历所有 frame
```

**`updateFramePlacement` 本质上只是一次 SE(3) 乘法**（[`frames.hxx:50`](include/pinocchio/src/algorithm/frames.hxx#L50)）：

```cpp
data.oMf[frame_id] = data.oMi[frame.parentJoint] * frame.placement;
```

**`updateFramePlacements` 则是对上式的遍历**（[`frames.hxx:30-35`](include/pinocchio/src/algorithm/frames.hxx#L30-L35)），
循环从索引 1 开始 —— 索引 0 是固定的 `universe` 帧，无需更新。

#### 4.1.3 两个必须避开的陷阱

**陷阱一：`updateFramePlacement` 不调用 FK。**
它直接读取 `data.oMi[parentJoint]`，假设你已经算好了。若未先调 `forwardKinematics`，
`data.oMi` 中是上次的旧值或初始值 —— 结果**静默错误**：不报错、不抛异常，只是数值不对。
这是 Pinocchio 最常见的 bug 来源之一。

```python
pin.forwardKinematics(model, data, q)                  # 必需的前置步骤
oMf = pin.updateFramePlacement(model, data, frame_id)  # 再刷这一个
```

**陷阱二：`framesForwardKinematics` 只填位姿，不填速度。**
它内部调用的是**只带 `q` 的一元 `forwardKinematics`**，因此 `data.v[]` 和 `data.a[]`
不会被更新。随后若调用 `getFrameVelocity`，拿到的是陈旧值或零值。需要速度时必须显式
走两步：

```python
pin.forwardKinematics(model, data, q, v)   # 三参数版本，填充 data.v[]
pin.updateFramePlacements(model, data)     # 再刷 Frame 位姿
vel = pin.getFrameVelocity(model, data, fid, pin.LOCAL_WORLD_ALIGNED)
```

#### 4.1.4 选用建议

**IK 迭代循环**（每步只关心一个末端 Frame）—— 用单数版，省掉遍历全部 Frame 的开销：

```python
for i in range(IT_MAX):
    pin.forwardKinematics(model, data, q)
    pin.updateFramePlacement(model, data, tool_id)   # O(1)，只刷需要的
    iMd = data.oMf[tool_id].actInv(oMdes)
    err = pin.log6(iMd).vector
    ...
    q = pin.integrate(model, q, v * DT)
```

**多接触约束**（如四足的四个足端）—— 一次 FK + 循环刷 4 个 Frame，比全量版本更省：

```python
pin.forwardKinematics(model, data, q)
for fid in feet_frame_ids:
    pin.updateFramePlacement(model, data, fid)
```

**已调用过其它算法**（`computeJointJacobians`、`crba`、`centerOfMass` 等内部都会跑 FK）
—— 直接刷 Frame，不要重复 FK：

```python
pin.computeJointJacobians(model, data, q)    # 内部已含 forwardKinematics
pin.updateFramePlacement(model, data, fid)   # 直接刷，不重复算
```

**可视化 / 碰撞检测** —— 交给上层函数，`updateGeometryPlacements` 内部已处理，
见 [`geometry-models.py`](examples/geometry-models.py)。

#### 4.1.5 配套的 Frame 层 API

这套机制还有若干配套函数，都遵循同样的"需先 FK"约定：


| 函数                                                  | 作用                             | 前置条件                       |
| ----------------------------------------------------- | -------------------------------- | ------------------------------ |
| `getFrameJacobian(model, data, fid, rf)`              | Frame 的$6\times n_v$ Jacobian   | `computeJointJacobians`        |
| `computeFrameJacobian(model, data, q, fid, rf)`       | 一体化版本（内部含 FK）          | 无                             |
| `getFrameVelocity(model, data, fid, rf)`              | Frame 空间速度                   | `forwardKinematics(m,d,q,v)`   |
| `getFrameAcceleration(model, data, fid, rf)`          | 空间加速度（李代数导数）         | `forwardKinematics(m,d,q,v,a)` |
| `getFrameClassicalAcceleration(model, data, fid, rf)` | 经典加速度（含$\omega\times v$） | 同上                           |

参数 `rf` 为参考系枚举，取值与物理含义见 [2.2.7](#227-参考系local--world--local_world_aligned)；
末端执行器任务通常用 `LOCAL_WORLD_ALIGNED`。最后两行的区别见
[2.2.5](#225-空间加速度最反直觉的一个)。

---

### 4.2 逆向运动学 IK

**目标**：给定末端目标位姿 $T_{des}$，求关节配置 $q$。

**方法**：带阻尼最小二乘的迭代 Jacobian 法（Levenberg–Marquardt 思路）。

#### 4.2.1 为什么误差不能写成"目标减当前"

工业机器人里常把误差写成 $e = x_{des} - x$。但位姿属于 $SE(3)$ 流形，减法没有定义：
两个旋转矩阵相减不再是旋转矩阵。正确做法是用**相对变换的对数映射**：

$$
{}^iM_{des} = {}^oM_i^{-1}\cdot{}^oM_{des}, \qquad e = \log_6\!\left({}^iM_{des}\right)\in\mathbb{R}^6
$$

物理含义：$e$ 是"从当前位姿出发，沿哪个 6D 螺旋运动（旋量）走单位时间能到达目标"。
当且仅当 ${}^iM_{des} = I$ 时 $e = 0$，因此 $\|e\|\to 0$ 是正确的收敛判据。

代码对应（[`inverse-kinematics.py`](examples/inverse-kinematics.py)）：

```python
iMd = data.oMi[JOINT_ID].actInv(oMdes)   # ᶦM_des = oMᵢ⁻¹ · oM_des
err = pin.log6(iMd).vector                # e ∈ R⁶
```

#### 4.2.2 关键推导：为什么需要 Jlog6

我们要求的是 $\dfrac{\partial e}{\partial q}$，但 $e = \log_6({}^iM_{des}(q))$ 是 $q$ 的
**复合函数**，必须用链式法则拆成两段：

$$
\frac{\partial e}{\partial q}
= \underbrace{\frac{\partial \log_6(M)}{\partial M}}_{\text{Jlog6，}6\times6}
\cdot
\underbrace{\frac{\partial\, {}^iM_{des}}{\partial q}}_{\text{几何 Jacobian，}6\times n_v}
$$

**第二段**：由 ${}^iM_{des} = {}^oM_i^{-1}\,{}^oM_{des}$，只有 ${}^oM_i$ 依赖 $q$。
关节 $i$ 的局部速度为 ${}^i\nu_i = {}^iJ_i\,\dot q$，扰动 $q$ 会使 ${}^oM_i$ 右乘
$\exp(\delta)$，从而使 ${}^iM_{des}$ **左乘** $\exp(-\delta)$ —— 负号由此而来。

**第一段**：$\log_6$ 是高度非线性的（含 $\theta/\sin\theta$ 型因子）。若省略 Jlog6 而直接
用 $J$，等价于假设 $\log$ 是恒等映射 —— 在小误差时近似成立，但大角度误差下会显著拖慢
收敛甚至发散。合并两段：

$$
J_{\text{eff}} = -J_{\log_6}\!\left({}^iM_{des}^{-1}\right)\cdot {}^iJ_i(q)
$$

> **数值验证**：`Jlog6(M)` 与 $\log_6$ 的有限差分在本仓库实测吻合到 $1.7\times10^{-8}$。

```python
J = pin.computeJointJacobian(model, data, q, JOINT_ID)  # LOCAL 系
J = -np.dot(pin.Jlog6(iMd.inverse()), J)
```

#### 4.2.3 阻尼最小二乘的来历

理想情况下解 $J_{\text{eff}}\,\dot q = -e$。但 $J_{\text{eff}}$ 在奇异位形附近条件数
极大，纯最小二乘 $\dot q = -J^{+}e$ 会给出爆炸的关节速度。改为求解带正则项的问题：

$$
\min_{\dot q}\ \tfrac12\|J_{\text{eff}}\dot q + e\|^2 + \tfrac{\lambda^2}{2}\|\dot q\|^2
$$

令梯度为零：$(J^\top J + \lambda^2 I_{n_v})\dot q = -J^\top e$。利用
**推移恒等式** $(J^\top J+\lambda^2 I_{n_v})^{-1}J^\top = J^\top(JJ^\top+\lambda^2 I_6)^{-1}$
（右侧只需求逆一个 $6\times6$ 矩阵，而非 $n_v\times n_v$）：

$$
\boxed{\dot q = -J_{\text{eff}}^\top\left(J_{\text{eff}}J_{\text{eff}}^\top + \lambda^2 I_6\right)^{-1} e}
$$

对人形（$n_v\approx 36$）这个变换把求逆规模从 $36\times36$ 降到 $6\times6$。$\lambda$ 的
作用是把 $J$ 的奇异值 $\sigma$ 替换为 $\sigma/(\sigma^2+\lambda^2)$：当 $\sigma\gg\lambda$
时几乎不变，当 $\sigma\to 0$ 时增益被限制在 $1/(2\lambda)$ 而非发散。

#### 4.2.4 流形积分

$$
q_{k+1} = q_k \oplus (\dot q\,\Delta t) \equiv \texttt{pin.integrate(model, q, v*DT)}
$$

**不能写成 $q \mathrel{+}= \dot q\Delta t$** —— 浮动基的四元数分量会失去单位范数。
`integrate` 对每种关节类型分别调用其指数映射（见 [2.1](#21-李群-se3)）。

#### 4.2.5 完整算法


| 步骤        | 公式                                       | API                              |
| ----------- | ------------------------------------------ | -------------------------------- |
| ① FK       | ${}^oM_i(q)$                               | `forwardKinematics`              |
| ② 误差     | $e=\log_6({}^oM_i^{-1}{}^oM_{des})$        | `log6(oMi.actInv(oMdes))`        |
| ③ 判敛     | $\|e\|<\varepsilon$                        | 典型$\varepsilon=10^{-4}$        |
| ④ Jacobian | $J_{\text{eff}}=-J_{\log_6}\cdot J$        | `Jlog6` + `computeJointJacobian` |
| ⑤ 阻尼解   | $\dot q=-J^\top(JJ^\top+\lambda^2I)^{-1}e$ | `np.linalg.solve`                |
| ⑥ 积分     | $q\leftarrow q\oplus\dot q\Delta t$        | `integrate`                      |

**3D 位置版本的简化**：若只约束位置（[`inverse-kinematics-3d.py`](examples/inverse-kinematics-3d.py)），
误差 $e = {}^iM_{des}.\text{translation}\in\mathbb{R}^3$ 是欧氏量，**无需 Jlog6**，
且只取 $J$ 的前 3 行。这正是 2.2 节所说"位置是线性的、姿态是流形的"的直接体现。

---

### 4.3 逆向动力学 RNEA

**目标**：给定 $(q,\dot q,\ddot q)$，求所需关节力矩 $\tau$。

#### 4.3.1 从拉格朗日方程到递推形式

多体系统的运动方程为：

$$
M(q)\ddot q + \underbrace{C(q,\dot q)\dot q + g(q)}_{\text{非线性项 }n(q,\dot q)} = \tau + J_c^\top\lambda
$$

**朴素思路**是分别构造 $M$、$C$、$g$ 再相加 —— 但构造 $C$ 需要 $O(n^3)$ 的 Christoffel
符号。RNEA 的洞察是：**如果只需要 $\tau$ 这个结果，根本不必显式构造任何矩阵**，
沿运动树跑两趟递推即可，复杂度 $O(n)$。

#### 4.3.2 第一趟：正向传播运动学（根 → 叶）

对每个关节 $i$，父关节记 $\lambda(i)$，关节运动子空间记 $S_i$（旋转关节为
$[0,0,0,\ 0,0,1]^\top$ 之类的常向量）：

$$
\nu_i = {}^iX_{\lambda(i)}\,\nu_{\lambda(i)} + S_i\dot q_i
$$

**推导**：连杆 $i$ 的速度 = 父连杆速度（变换到 $i$ 系）+ 关节 $i$ 自身贡献。这就是
2.2.4 节的伴随变换 ${}^iX_{\lambda(i)} = \mathrm{Ad}_{{}^iM_{\lambda(i)}}$。

对时间求导得加速度。注意 ${}^iX_{\lambda(i)}$ **本身随时间变化**，其导数贡献
$\nu_i\times(S_i\dot q_i)$：

$$
a_i = {}^iX_{\lambda(i)}\,a_{\lambda(i)} + S_i\ddot q_i + \underbrace{\nu_i\times S_i\dot q_i}_{\text{速度积项}}
$$

最后这项就是**科氏力/离心力的几何来源**——它完全由速度的叉乘产生，不含 $\ddot q$。

> **重力的技巧**：把根节点的加速度初始化为 $a_0 = -g$（重力加速度取负），重力就自动
> 通过递推传遍全身，无需单独处理。这是 Featherstone 体系的经典手法。

#### 4.3.3 第二趟：反向传播力（叶 → 根）

对连杆 $i$ 应用 Newton–Euler 方程（见 [2.2.6](#226-运动方程中的-omegatimes-项)）：

$$
f_i^{\text{net}} = I_i a_i + \nu_i\times^* (I_i\nu_i)
$$

连杆 $i$ 受力平衡：自身惯性力 = 父关节施加的力 − 传给子关节的力 + 外力。故：

$$
f_i = I_i a_i + \nu_i\times^*(I_i\nu_i) - {}^if_i^{\text{ext}} + \sum_{j\in\text{child}(i)} {}^iX_j^*\,f_j
$$

注意子关节的力用 ${}^iX_j^* = \mathrm{Ad}^\top$ 变换（**余伴随**，见 [2.2.4](#224-坐标变换伴随与余伴随)）
—— 这正是力与速度对偶性的直接应用。

**投影到关节轴**：关节只能沿其运动子空间施力，其余分量由机械结构承担：

$$
\tau_i = S_i^\top f_i
$$

这一步是**虚功原理**：$\tau_i\dot q_i = f_i^\top(S_i\dot q_i)$ 对任意 $\dot q_i$ 成立。

#### 4.3.4 由 RNEA 派生的三个常用量

RNEA 是"万能积木"，通过特殊输入可提取动力学方程的各个部分：


| 调用                        | 得到          | 原理                           |
| --------------------------- | ------------- | ------------------------------ |
| `rnea(m,d,q,0,0)`           | $g(q)$        | 令$\dot q=\ddot q=0$，只剩重力 |
| `nle(m,d,q,v)`              | $C\dot q+g$   | 令$\ddot q=0$                  |
| `rnea(m,d,q,v,a)` $-$ `nle` | $M(q)\ddot q$ | 差分消去非线性项               |

> **数值验证**（本仓库样例人形，34 DOF）：
> `tau == M @ a + nle` ✅、`nle(v=0) == g` ✅ 均严格成立。

**复杂度**：$O(n)$，无矩阵求逆，是所有动力学算法中最快的。

---

### 4.4 正向动力学 ABA

**目标**：给定 $(q,\dot q,\tau)$，求 $\ddot q$。

#### 4.4.1 为什么不能直接求逆

形式解是 $\ddot q = M^{-1}(\tau - n(q,\dot q))$，但显式构造 $M$ 需 $O(n^2)$、
求逆需 $O(n^3)$。ABA 用**三趟 $O(n)$ 递推**得到同样结果。

#### 4.4.2 核心概念：关节化体惯量

ABA 的关键洞察是引入**关节化体惯量**（Articulated Body Inertia）$I_i^A$：

> $I_i^A$ 回答的问题是：*把以连杆 $i$ 为根的整个子树看作一个整体，当子关节都能
> **自由响应**（按各自的 $\tau_j$ 加速）时，在连杆 $i$ 上施加空间力 $f$ 会产生多大加速度？*

它满足 $f_i = I_i^A a_i + p_i^A$，其中 $p_i^A$ 是**偏置力**（速度积项与已知 $\tau$ 的贡献）。

与之对比：CRBA 中的复合刚体惯量 $I^c$ 假设子关节**全部锁死**。
"锁死"对应 $M$，"自由"对应 $M^{-1}$ —— 这是两个算法的本质分野。

#### 4.4.3 递推公式的推导

设已知子关节 $j$ 满足 $f_j = I_j^A a_j + p_j^A$。关节 $j$ 的力矩方程为：

$$
\tau_j = S_j^\top f_j = S_j^\top\left(I_j^A a_j + p_j^A\right)
$$

代入 $a_j = {}^jX_i\,a_i + S_j\ddot q_j + c_j$（$c_j$ 为速度积项），解出 $\ddot q_j$：

$$
\ddot q_j = \underbrace{\left(S_j^\top I_j^A S_j\right)^{-1}}_{\textstyle D_j^{-1}}
\left[\tau_j - S_j^\top I_j^A\left({}^jX_i a_i + c_j\right) - S_j^\top p_j^A\right]
$$

把这个 $\ddot q_j$ **回代**进 $f_j$ 的表达式，消去 $\ddot q_j$ 后整理，得到父连杆看到的
等效惯量：

$$
\boxed{I_i^A = I_i + \sum_{j\in\text{child}(i)} {}^iX_j^*\left(I_j^A - \frac{I_j^A S_j S_j^\top I_j^A}{S_j^\top I_j^A S_j}\right){}^jX_i}
$$

括号内是 **Schur 补**。它的物理意义十分直观：从子树的"锁死惯量" $I_j^A$ 中
**扣除**关节 $j$ 可自由转动的那个方向所释放的惯量。若关节 $j$ 完全刚性
（$D_j\to\infty$），修正项消失，退化为 CRBA 的复合惯量。

#### 4.4.4 三趟结构


| 趟次 | 方向     | 计算内容                                                    |
| ---- | -------- | ----------------------------------------------------------- |
| ①   | 根 → 叶 | $\nu_i$、速度积项 $c_i$（与 RNEA 第一趟同，但无 $\ddot q$） |
| ②   | 叶 → 根 | $I_i^A$、偏置力 $p_i^A$（上面的 Schur 补递推）              |
| ③   | 根 → 叶 | 用已知的$a_{\lambda(i)}$ 代入上式求 $\ddot q_i$，再求 $a_i$ |

> **数值验证**：`aba(model, data, q, v, rnea(q,v,a)) == a` 严格成立
> —— ABA 与 RNEA 互为逆运算。

**复杂度**：$O(n)$。对 $n>\!\approx 8$ 的系统快于"CRBA + Cholesky 求解"的 $O(n^3)$ 路线；
但若同时还需要 $M$ 本身（如 WBC 中的任务空间投影），则后者更划算。

---

### 4.5 质量矩阵 CRBA

**目标**：计算广义质量矩阵 $M(q)\in\mathbb{R}^{n_v\times n_v}$。

#### 4.5.1 从动能出发

$M$ 的定义来自系统动能的二次型：

$$
T = \tfrac12\dot q^\top M(q)\dot q = \sum_i \tfrac12\,\nu_i^\top I_i\,\nu_i
$$

代入 $\nu_i = J_i\dot q$（$J_i$ 为连杆 $i$ 的 $6\times n_v$ Jacobian）：

$$
M(q) = \sum_i J_i^\top I_i J_i
$$

直接按此式计算需 $O(n^2)$ 次 $6\times6$ 运算。CRBA 利用**运动树的稀疏性**降到
$O(n\,d)$（$d$ 为树深度）。

#### 4.5.2 稀疏性：为什么大量元素恒为零

关键结构性质：

$$
M_{ij}\neq 0 \iff \text{关节 }i\text{ 与 }j\text{ 在同一条从根到叶的支链上}
$$

**理由**：$M_{ij} = S_i^\top(\cdots)S_j$ 描述关节 $i$ 与 $j$ 的惯性耦合。若二者位于树的
两条不同分支（如左臂与右臂），则不存在同时被两者驱动的连杆，耦合项为零。

> **数值验证**（样例人形 34 DOF）：遍历全部 $34\times34$ 个元素，
> **非同支链却非零的元素数 = 0**，稀疏性严格成立。

#### 4.5.3 复合刚体惯量与递推

定义**复合刚体惯量** $I_i^c$ —— 把以 $i$ 为根的整个子树**锁死为单个刚体**后的空间惯量：

$$
I_i^c = I_i + \sum_{j\in\text{child}(i)} {}^iX_j^*\,I_j^c\,{}^jX_i
$$

叶到根一趟递推即可求出全部 $I_i^c$。有了它，矩阵元为：

$$
\boxed{M_{ij} = S_i^\top\,I_i^c\,{}^iX_j\,S_j,\qquad j\in\text{subtree}(i)}
$$

> ⚠️ 注意下标：用的是 **$I_i^c$（行索引 $i$ 的复合惯量）**，列 $j$ 遍历 **$i$ 的子树**。
> 源码 [`crba.hxx`](include/pinocchio/src/algorithm/crba.hxx) 中对应
> `data.M.block(idx_v(i), idx_v(i), nv(i), nvSubtree[i])` —— 行块 $i$、列跨度恰为
> $i$ 的子树，与上式一致。

**物理直觉**：$M_{ij}$ = "让关节 $j$ 单位加速时，关节 $i$ 需要提供多少力矩"。由于 $j$ 在
$i$ 的子树内，$i$ 必须带动**整个子树**，故用锁死后的复合惯量 $I_i^c$。

**复杂度**：$O(n d)$，对人形（$d\approx 6$）远快于朴素的 $O(n^2)$。

#### 4.5.4 性质与后续使用

$M$ 对称正定，故可作 $LDL^\top$ 分解。Pinocchio 的 `cholesky` 模块进一步利用上述稀疏
模式，分解复杂度亦为 $O(nd)$ 而非稠密的 $O(n^3)$：

```python
M = pin.crba(model, data, q)          # 上三角计算，返回时已对称化
pin.cholesky.decompose(model, data)   # 稀疏 LDLᵀ
pin.cholesky.solve(model, data, rhs)  # 解 M x = rhs
```

> **实现细节**：算法内部只填充上三角，Python 绑定返回前会补全为对称阵
> （实测 `crba` 返回值满足 `M == M.T`）。若直接读取 `data.M` 则需自行注意。

---

### 4.6 质心动量学

这是**人形机器人控制的核心工具**，也是 Pinocchio 区别于其他库的特色之一。

#### 4.6.1 定义：把全身动量约化到质心

系统总空间动量是各连杆动量之和。将其全部用**余伴随**变换（[2.2.4](#224-坐标变换伴随与余伴随)）
约化到**以质心 $c$ 为原点、与世界系平行**的坐标系：

$$
h_g = \sum_i {}^gX_i^*\,I_i\,\nu_i
= \begin{bmatrix} p \\ L_g \end{bmatrix} \in\mathbb{R}^6
$$

其中 $p = m\dot c$ 为总线动量，$L_g$ 为**相对质心**的总角动量。

> ⚠️ 顺序仍是**线性在前**（与 [2.2](#22-空间速度与空间力) 一致）：
> `data.hg.linear` $= p = m\dot c$，`data.hg.angular` $= L_g$。
> **数值验证**：`hg.linear == mass * vcom` 严格成立。

由于每个 $\nu_i = J_i\dot q$ 都线性依赖 $\dot q$，总动量也线性依赖 $\dot q$：

$$
\boxed{h_g = A_g(q)\,\dot q},\qquad A_g(q) = \sum_i {}^gX_i^*\,I_i\,J_i \in\mathbb{R}^{6\times n_v}
$$

$A_g$ 称为**质心动量矩阵**（CMM, Centroidal Momentum Matrix）。

> **数值验证**：`hg == Ag @ v` 严格成立；样例人形 $A_g$ 形状 $6\times34$，总质量 15.5 kg。

#### 4.6.2 为什么它是平衡控制的基础

对 $h_g$ 求时间导数，并注意**内力成对抵消**（牛顿第三定律），只剩外力：

$$
\dot h_g = \sum_k {}^gX_k^*\,\phi_k^{\text{ext}} + \begin{bmatrix} m\,\mathbf{g} \\ 0 \end{bmatrix}
$$

展开成线性/角动量两式，即人形平衡控制的出发点：

$$
\dot p = \sum_k f_k^{\text{contact}} + m\mathbf{g}
\qquad\text{（质心运动只由合外力决定）}
$$

$$
\dot L_g = \sum_k (p_k^{\text{contact}} - c)\times f_k^{\text{contact}} + \sum_k \tau_k^{\text{contact}}
$$

**三条关键推论**：

1. **欠驱动的本质**：$\dot p$ 只能通过接触力改变。关节力矩 $\tau$ 不出现在这两式中
   —— 人形在腾空相无法改变质心轨迹，这是欠驱动的数学表述。
2. **ZMP / 支撑多边形**：由于摩擦锥要求 $f_k^z\ge 0$，第二式给出对
   $\dot L_g$ 的**线性不等式约束**，等价于 ZMP 必须落在支撑多边形内。
3. **MPC 状态方程**：这两式仅含 $(c,\dot c,L_g)$ 与接触力，构成
   **低维**（通常 6~9 维）的质心动力学模型，是实时 MPC 可行的关键 —— 无需在 MPC 中
   处理全部 34 个自由度。

#### 4.6.3 API

```python
pin.computeCentroidalMomentum(model, data, q, v)
# → data.hg      : 质心动量 h_g (6D, Force 类型：linear=p, angular=L_g)
# → data.com[0]  : 质心位置 c
# → data.vcom[0] : 质心速度 ċ

pin.computeCentroidalMap(model, data, q)      # → data.Ag (6 × nv)
pin.computeCentroidalMomentumTimeVariation(model, data, q, v, a)
# → data.dhg     : ḣ_g（含 Ȧ_g v 项）
```

---

### 4.7 接触约束动力学

**目标**：人形足底接触时，在约束下同时求 $\ddot q$ 和接触力 $\lambda$。

#### 4.7.1 约束的三个层次：位置 → 速度 → 加速度

接触的物理约束写在**位置**层面：接触点不能穿透也不能滑移，即 $c(q) = 0$。
但动力学方程是关于 $\ddot q$ 的，因此需连续求导两次：

$$
c(q) = 0
\ \xrightarrow{\ \frac{d}{dt}\ }\ \underbrace{J(q)\dot q = 0}_{J \,\equiv\, \partial c/\partial q}
\ \xrightarrow{\ \frac{d}{dt}\ }\ J\ddot q + \dot J\dot q = 0
$$

**只有最后一式能直接进入动力学方程**。记 $\gamma = \dot J\dot q$（接触点的漂移加速度，
即 2.2 节的速度积项），约束为 $J\ddot q = -\gamma$。

#### 4.7.2 从 Gauss 最小约束原理导出 KKT

**Gauss 最小约束原理**指出：受约束系统的真实加速度，是在所有满足约束的加速度中，
使其与"自由加速度"的 $M$-加权距离最小者：

$$
\min_{\ddot q}\ \tfrac12\left\|\ddot q - \ddot q_{\text{free}}\right\|_{M}^2
\quad\text{s.t.}\quad J\ddot q = -\gamma
$$

其中 $\ddot q_{\text{free}} = M^{-1}(\tau - n)$ 为无约束时的加速度。构造 Lagrange 函数：

$$
\mathcal{L} = \tfrac12(\ddot q - \ddot q_{\text{free}})^\top M(\ddot q - \ddot q_{\text{free}})
- \lambda^\top(J\ddot q + \gamma)
$$

$$
\frac{\partial\mathcal{L}}{\partial\ddot q} = M\ddot q - M\ddot q_{\text{free}} - J^\top\lambda = 0
\ \Longrightarrow\ M\ddot q + n = \tau + J^\top\lambda
$$

与约束式联立，得 **KKT 鞍点系统**：

$$
\boxed{\begin{bmatrix} M & J^\top \\ J & 0\end{bmatrix}
\begin{bmatrix}\ddot q\\ -\lambda\end{bmatrix}
= \begin{bmatrix}\tau - n\\ \gamma\end{bmatrix}}
$$

**$\lambda$ 的物理意义**：Lagrange 乘子恰为**接触力**，$J^\top\lambda$ 是它映射到关节
空间的等效力矩。这不是巧合 —— 乘子的量纲由约束的量纲决定，而 $J^\top$ 正是
2.2.4 节的余伴随（力的变换）。

> **数值验证**（样例人形 + CONTACT_6D）：`constraintDynamics` 返回的 $(\ddot q,\lambda)$
> 代入 $M\ddot q + n - \tau - J^\top\lambda$ 的残差为 $1.7\times10^{-13}$。
> 注意 $J\ddot q \neq 0$（实测 0.108）—— 因为约束是 $J\ddot q = -\gamma$ 而非 $J\ddot q=0$。

#### 4.7.3 Delassus 矩阵与求解

对 KKT 系统作块消元（第一行解出 $\ddot q$ 代入第二行）：

$$
\underbrace{\left(J M^{-1} J^\top\right)}_{\textstyle \Lambda^{-1}\ \text{（Delassus 矩阵）}}\lambda
= J M^{-1}(\tau - n) + \gamma
$$

$\Lambda^{-1} = JM^{-1}J^\top\in\mathbb{R}^{n_c\times n_c}$ 称为 **Delassus 算子**，
其逆 $\Lambda$ 为**操作空间惯量矩阵**：

- $\Lambda^{-1}$：接触空间的"柔度"—— 单位接触力引起多大接触点加速度
- $\Lambda$：接触空间的"有效惯量"—— 产生单位接触加速度需多大力

**Pinocchio 的高效实现**不走上述稠密消元，而用 **Contact Cholesky 分解**直接对
$(n_v+n_c)$ 维 KKT 矩阵作稀疏 $LDL^\top$，充分利用运动树结构：


| 方法             | 复杂度                 |
| ---------------- | ---------------------- |
| 稠密 LU/LDLᵀ    | $O((n_v+n_c)^3)$       |
| Contact Cholesky | $O(n_v d + n_c^2 n_v)$ |

```python
pin.initConstraintDynamics(model, data, contact_models)   # 预分配稀疏结构
a = pin.constraintDynamics(model, data, q, v, tau, contact_models, contact_datas)
lam = data.lambda_c                                        # 接触力
# 也可直接取 Delassus：
data.contact_chol.getInverseOperationalSpaceInertiaMatrix()  # Λ⁻¹ = J M⁻¹ Jᵀ
```

完整用例见 [`contact-cholesky.py`](examples/contact-cholesky.py) 与
[`simulation-contact-dynamics.py`](examples/simulation-contact-dynamics.py)。

#### 4.7.4 约束类型与 Baumgarte 稳定化


| `ContactType` | 含义                    | 每个接触的$n_c$ |
| ------------- | ----------------------- | --------------- |
| `CONTACT_3D`  | 3D 点接触（仅位置）     | 3               |
| `CONTACT_6D`  | 6D 面接触（完整 SE(3)） | 6               |

**数值漂移问题**：我们求解的是二阶导约束 $J\ddot q = -\gamma$，积分两次后 $c(q)=0$ 只在
理论上保持。浮点误差会随时间累积，使接触点缓慢穿透或漂离。解决方法是在右端加
**PD 修正项**（Baumgarte 稳定化）：

$$
J\ddot q = -\gamma - K_p\,c(q) - K_d\,J\dot q
$$

这使约束违反量按二阶系统衰减，取 $K_d = 2\sqrt{K_p}$ 为临界阻尼（无超调）：

```python
constraint_model.corrector.Kp[:] = 10
constraint_model.corrector.Kd[:] = 2.0 * np.sqrt(constraint_model.corrector.Kp)
```

见 [`simulation-closed-kinematic-chains.py`](examples/simulation-closed-kinematic-chains.py)。

---

### 4.8 解析梯度

Pinocchio 提供所有主要算法的**解析偏导数**，这是 DDP/iLQR 等轨迹优化算法高效运行的基础。
相比数值差分，解析梯度更精确（不依赖步长选择）且通常快 10 倍以上。

#### 4.8.1 RNEA 梯度

对 $\tau = M(q)\ddot q + n(q,\dot q)$ 求偏导：

$$
\frac{\partial\tau}{\partial q},\qquad
\frac{\partial\tau}{\partial\dot q},\qquad
\frac{\partial\tau}{\partial\ddot q} = M(q)
$$

第三个等式是**精确恒等式**：$\tau$ 关于 $\ddot q$ 是线性的，系数矩阵正是 $M$。
这提供了一个免费的 $M$ 计算途径。

```python
pin.computeRNEADerivatives(model, data, q, v, a)
data.dtau_dq, data.dtau_dv, data.M    # ← 从 data 读取
```

> ⚠️ **务必从 `data` 字段读取**，而非依赖 Python 绑定的返回值顺序。
> **数值验证**：`data.dtau_dq` / `data.dtau_dv` 与有限差分吻合到 $10^{-7}$；
> `dtau_da == M` 严格成立。

#### 4.8.2 ABA 梯度与链式法则

ABA 梯度可由 RNEA 梯度推出。对恒等式 $\tau = M\ddot q + n$ 关于 $x\in\{q,\dot q\}$
求全微分，注意此时 $\ddot q$ 是 $x$ 的函数而 $\tau$ 固定：

$$
0 = \frac{\partial\tau}{\partial x}\bigg|_{\ddot q\ \text{固定}} + M\,\frac{\partial\ddot q}{\partial x}
\ \Longrightarrow\
\boxed{\frac{\partial\ddot q}{\partial x} = -M^{-1}\frac{\partial\tau}{\partial x}}
$$

关键前提：右侧的 $\partial\tau/\partial x$ 必须在 **ABA 解出的那个 $\ddot q$** 处求值。

$$
\frac{\partial\ddot q}{\partial\tau} = M(q)^{-1}
$$

> **数值验证**：在样例机械臂与人形上，
> `ddq_dq == -Minv @ dtau_dq` 与 `ddq_dv == -Minv @ dtau_dv` 均成立（$10^{-6}$ 容差）。

```python
pin.computeABADerivatives(model, data, q, v, tau)
data.ddq_dq, data.ddq_dv, data.Minv
```

> ⚠️ **`data.Minv` 只填充上三角！** 这是文档明确声明的行为
> （[`aba-derivatives.hpp:154`](include/pinocchio/algorithm/aba-derivatives.hpp#L154)），
> 实测下三角确为残留值。直接拿它做矩阵乘法会得到错误结果，须先对称化：
>
> ```python
> Minv = np.triu(data.Minv) + np.triu(data.Minv, 1).T
> ```

#### 4.8.3 与 DDP 的状态空间矩阵的关系

DDP/iLQR 需要离散状态转移的线性化 $\delta x_{k+1} = A_k\delta x_k + B_k\delta u_k$。
以 $x = (q,\dot q)$、$u=\tau$、半隐式欧拉积分为例：

$$
A_k = \begin{bmatrix}
I + \Delta t^2\,\partial\ddot q/\partial q & \Delta t\left(I + \Delta t\,\partial\ddot q/\partial\dot q\right)\\[2pt]
\Delta t\,\partial\ddot q/\partial q & I + \Delta t\,\partial\ddot q/\partial\dot q
\end{bmatrix},
\qquad
B_k = \begin{bmatrix}\Delta t^2 M^{-1}\\ \Delta t\,M^{-1}\end{bmatrix}
$$

> ⚠️ 常见误解：$\partial\ddot q/\partial q$ **本身不是** $A_k$，$M^{-1}$ 本身也不是 $B_k$。
> 它们只是构成 $A_k,B_k$ 的**分块**，还需按积分格式组装（上式对应半隐式欧拉；
> 显式欧拉的组装方式不同）。此外浮动基情形下 $q$ 的扰动须在李代数上取，
> $I$ 需替换为对应的 $\partial\,\text{integrate}/\partial q$ 项。

#### 4.8.4 运动学梯度

用于 Frame 速度/加速度对 $(q,\dot q)$ 的偏导，WBC 与 Frame 空间 MPC 必需：

$$
\frac{\partial\nu_i}{\partial q},\quad
\frac{\partial a_i}{\partial q},\quad
\frac{\partial a_i}{\partial\dot q},\quad
\frac{\partial a_i}{\partial\ddot q}
$$

两个恒等式值得记住（见 [`kinematics-derivatives.cpp`](examples/kinematics-derivatives.cpp)）：

$$
\frac{\partial\nu_i}{\partial\dot q} = \frac{\partial a_i}{\partial\ddot q} = J_i(q)
$$

即"速度对速度"与"加速度对加速度"的偏导都等于几何 Jacobian。正因二者相同，
`getJointAccelerationDerivatives` 不单独返回 $\partial\nu_i/\partial\dot q$。

```python
pin.computeForwardKinematicsDerivatives(model, data, q, v, a)
pin.getJointVelocityDerivatives(model, data, jid, pin.LOCAL)
pin.getJointAccelerationDerivatives(model, data, jid, pin.LOCAL)
```

---

## 5. Demo 注释：Python

### 5.1 `overview-simple.py` — 最小完整示例

```python
import pinocchio

# 创建内置 6-DOF 串联机械臂模型（无需 URDF）
model = pinocchio.buildSampleModelManipulator()
# 为该模型创建算法数据容器（存储中间计算结果）
data = model.createData()

# 生成零位配置 q=0（关节角全为0）
q = pinocchio.neutral(model)
# 速度与加速度置零
v = pinocchio.utils.zero(model.nv)
a = pinocchio.utils.zero(model.nv)

# 调用 RNEA：逆动力学，计算产生加速度 a 所需的关节力矩 τ
# 此处 a=0, v=0，故 τ 只反映重力补偿项 g(q)
tau = pinocchio.rnea(model, data, q, v, a)
print("tau = ", tau.T)
```

**涉及算法**：RNEA（逆动力学）。

---

### 5.2 `overview-urdf.py` — 从 URDF 加载并做正向运动学

```python
from pathlib import Path
from sys import argv
import pinocchio

pinocchio_model_dir = Path(__file__).parent.parent / "models"

# URDF 路径：默认使用 UR5，可命令行传入自定义路径
urdf_filename = (
    pinocchio_model_dir / "example-robot-data/robots/ur_description/urdf/ur5_robot.urdf"
    if len(argv) < 2
    else argv[1]
)

# 从 URDF 构建运动学/动力学模型
# 固定基机器人：不传第二个参数
# 浮动基（人形）：传入 pinocchio.JointModelFreeFlyer()
model = pinocchio.buildModelFromUrdf(urdf_filename)
print("model name: " + model.name)

# 创建算法数据容器，与 model 绑定
data = model.createData()

# 在配置空间中随机采样一个合法配置（自动满足关节限位）
q = pinocchio.randomConfiguration(model)
print(f"q: {q.T}")

# 正向运动学：递推计算所有关节在世界坐标系中的 SE(3) 变换
# 结果存入 data.oMi[i]
pinocchio.forwardKinematics(model, data, q)

# 打印每个关节的平移部分（世界系下的 3D 坐标）
for name, oMi in zip(model.names, data.oMi):
    print("{:<24} : {: .2f} {: .2f} {: .2f}".format(name, *oMi.translation.T.flat))
```

**涉及算法**：FK（正向运动学）。

---

### 5.3 `inverse-kinematics.py` — 6D 逆向运动学（SE3 误差）

```python
import numpy as np
import pinocchio
from numpy.linalg import norm, solve

model = pinocchio.buildSampleModelManipulator()
data = model.createData()

JOINT_ID = 6                                       # 末端关节 ID
# 目标位姿：旋转=单位阵（无旋转），平移=[1,0,1]
oMdes = pinocchio.SE3(np.eye(3), np.array([1.0, 0.0, 1.0]))

q    = pinocchio.neutral(model)  # 初始配置
eps  = 1e-4                       # 收敛阈值（误差 L2 范数）
IT_MAX = 1000                     # 最大迭代次数
DT   = 1e-1                       # 积分步长
damp = 1e-12                      # 阻尼系数（防止奇异）

i = 0
while True:
    # Step 1：正向运动学，更新 data.oMi[]
    pinocchio.forwardKinematics(model, data, q)

    # Step 2：计算局部系下的相对变换 iMd = T_i^{-1} * T_des
    iMd = data.oMi[JOINT_ID].actInv(oMdes)

    # Step 3：误差 = SE(3) 对数映射 → 6D 切向量（平移+旋转误差）
    # log(iMd) ∈ se(3) ≅ R^6，直接用作迭代误差
    err = pinocchio.log(iMd).vector

    if norm(err) < eps:
        success = True; break
    if i >= IT_MAX:
        success = False; break

    # Step 4：计算关节空间 Jacobian（在关节局部系下）
    J = pinocchio.computeJointJacobian(model, data, q, JOINT_ID)

    # Step 5：对 Jacobian 乘以 log 映射的 Jacobian 修正项 J_log6
    # 这是将 SE(3) 误差正确反向传播到关节空间所必需的链式法则
    J = -np.dot(pinocchio.Jlog6(iMd.inverse()), J)

    # Step 6：阻尼最小二乘 v = -J^T (JJ^T + λI)^{-1} e
    v = -J.T.dot(solve(J.dot(J.T) + damp * np.eye(6), err))

    # Step 7：在流形（李群）上积分，而非简单 q += v*DT
    # 对于浮动基，旋转部分是四元数，必须走李群积分
    q = pinocchio.integrate(model, q, v * DT)

    if not i % 10:
        print(f"{i}: error = {err.T}")
    i += 1

if success:
    print("Convergence achieved!")
else:
    print("Warning: algorithm did not converge")
print(f"\nresult: {q.flatten().tolist()}")
print(f"\nfinal error: {err.T}")
```

**涉及算法**：FK、几何 Jacobian、SE(3) 对数映射、$J_{\log_6}$、李群积分。

---

### 5.4 `inverse-kinematics-3d.py` — 纯位置 IK（3D 点跟踪）

```python
import numpy as np
import pinocchio
from numpy.linalg import norm, solve

model = pinocchio.buildSampleModelManipulator()
data = model.createData()

JOINT_ID = 6
oMdes = pinocchio.SE3(np.eye(3), np.array([1.0, 0.0, 1.0]))

q    = pinocchio.neutral(model)
eps  = 1e-4
IT_MAX = 1000
DT   = 1e-1
damp = 1e-12

it = 0
while True:
    pinocchio.forwardKinematics(model, data, q)
    iMd = data.oMi[JOINT_ID].actInv(oMdes)

    # 与 6D IK 的唯一区别：误差只取平移部分（R^3 而非 R^6）
    # 不关心末端姿态，只控制末端位置
    err = iMd.translation

    if norm(err) < eps:
        success = True; break
    if it >= IT_MAX:
        success = False; break

    J = pinocchio.computeJointJacobian(model, data, q, JOINT_ID)
    # 只取 Jacobian 的上 3 行（线速度部分）并取负号（局部系误差的符号约定）
    J = -J[:3, :]

    # 3×3 阻尼最小二乘
    v = -J.T.dot(solve(J.dot(J.T) + damp * np.eye(3), err))
    q = pinocchio.integrate(model, q, v * DT)

    if not it % 10:
        print(f"{it}: error = {err.T}")
    it += 1

if success:
    print("Convergence achieved!")
else:
    print("Warning: algorithm did not converge")
print(f"\nresult: {q.flatten().tolist()}")
print(f"\nfinal error: {err.T}")
```

**与 5.3 的对比**：不使用 SE(3) 对数映射，只控制末端位置；Jacobian 只取线速度行；
无需 $J_{\log_6}$ 修正（因为平移部分的误差本身是线性的）。

---

### 5.5 `inverse-dynamics.py` — RNEA 逆动力学

```python
from pathlib import Path
import numpy as np
import pinocchio as pin

pinocchio_model_dir = Path(__file__).parent.parent / "models/"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir = pinocchio_model_dir
urdf_filename = "ur5_robot.urdf"
urdf_model_path = model_path / "ur_description/urdf/" / urdf_filename

# buildModelsFromUrdf 同时返回运动学模型、碰撞模型、视觉模型
model, _, _ = pin.buildModelsFromUrdf(urdf_model_path, package_dirs=mesh_dir)
data = model.createData()

# 随机采样状态（实际使用时替换为真实状态）
q = pin.randomConfiguration(model)   # 关节角 [rad]
v = np.random.rand(model.nv, 1)       # 关节速度 [rad/s]
a = np.random.rand(model.nv, 1)       # 期望关节加速度 [rad/s²]

# RNEA：递推牛顿-欧拉，计算产生加速度 a 所需的力矩 τ
# 等价于：τ = M(q)a + C(q,v)v + g(q)
# 但用 O(n) 递推而非显式构造 M, C, g
tau = pin.rnea(model, data, q, v, a)

# 结果存在 data.tau，rnea() 也直接返回引用
print("Joint torques: " + str(tau))
```

**应用场景**：力矩前馈控制、逆动力学控制（IDC）、轨迹可行性验证。

---

### 5.6 `forward-dynamics-derivatives.py` — ABA 正向动力学 + 解析梯度

```python
import numpy as np
import pinocchio as pin

# 使用内置随机人形模型（双足，带浮动基）
model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

# 设置关节限位（样例模型默认无限位，优化器需要）
model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit =  np.ones((model.nq, 1))

q   = pin.randomConfiguration(model)         # 关节配置
v   = np.random.rand(model.nv, 1)            # 关节速度
tau = np.random.rand(model.nv, 1)            # 关节力矩

# 计算 ABA 及其对 (q, v, τ) 的解析梯度
# 内部调用三趟递推，一次性得到加速度和所有梯度
pin.computeABADerivatives(model, data, q, v, tau)

# ∂q̈/∂q：状态 q 对加速度的影响（DDP 中的 ∂f/∂q 的 q 块）
ddq_dq  = data.ddq_dq

# ∂q̈/∂v：状态 v 对加速度的影响（DDP 中的 ∂f/∂q 的 v 块）
ddq_dv  = data.ddq_dv

# ∂q̈/∂τ = M(q)^{-1}：控制输入对加速度的影响（DDP 中的 B 矩阵）
ddq_dtau = data.Minv

# 应用：DDP/iLQR 中的线性化
# f(x,u) = q̈  =>  f_x ≈ [ddq_dq, ddq_dv],  f_u = Minv
```

**应用场景**：Crocoddyl DDP 的动力学线性化、梯度下降轨迹优化。

---

### 5.7 `inverse-dynamics-derivatives.py` — RNEA 梯度

```python
import numpy as np
import pinocchio as pin

model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit =  np.ones((model.nq, 1))

q   = pin.randomConfiguration(model)
v   = np.random.rand(model.nv, 1)
a   = np.random.rand(model.nv, 1)

# 计算 RNEA 梯度：∂τ/∂q，∂τ/∂v，∂τ/∂a
pin.computeRNEADerivatives(model, data, q, v, a)

# ∂τ/∂q：配置对所需力矩的影响（重力梯度 + 科氏力梯度）
dtau_dq = data.dtau_dq

# ∂τ/∂v：速度对所需力矩的影响（科氏力 + 离心力梯度）
dtau_dv = data.dtau_dv

# ∂τ/∂a = M(q)：加速度对力矩的影响（即质量矩阵本身）
dtau_da = data.M
```

**应用场景**：逆动力学控制的线性化、运动学参数辨识的梯度计算。

---

### 5.8 `kinematics-derivatives.py` — 运动学梯度

```python
import numpy as np
import pinocchio as pin

model = pin.buildSampleModelHumanoidRandom()
data = model.createData()

model.lowerPositionLimit = -np.ones((model.nq, 1))
model.upperPositionLimit =  np.ones((model.nq, 1))

q = pin.randomConfiguration(model)
v = np.random.rand(model.nv, 1)
a = np.random.rand(model.nv, 1)

# 一次调用计算所有关节速度/加速度关于 (q,v,a) 的梯度
pin.computeForwardKinematicsDerivatives(model, data, q, v, a)

joint_name = "rleg6_joint"
joint_id = model.getJointId(joint_name)

# ∂v_i/∂q, ∂v_i/∂v（即几何 Jacobian）——世界系下
(dv_dq, dv_dv) = pin.getJointVelocityDerivatives(
    model, data, joint_id, pin.ReferenceFrame.WORLD
)
# 局部系下
(dv_dq_local, dv_dv_local) = pin.getJointVelocityDerivatives(
    model, data, joint_id, pin.ReferenceFrame.LOCAL
)

# ∂a_i/∂q, ∂a_i/∂v, ∂a_i/∂a——世界系下
(dv_dq, da_dq, da_dv, da_da) = pin.getJointAccelerationDerivatives(
    model, data, joint_id, pin.ReferenceFrame.WORLD
)
# 注意：∂v_i/∂v = ∂a_i/∂a（都等于几何 Jacobian）
```

**应用场景**：WBC 中 Frame 速度的 Jacobian、末端加速度对力矩的映射。

---

### 5.9 `simulation-contact-dynamics.py` — Talos 人形双足接触仿真

这是最重要的人形机器人 demo，演示浮动基 + 刚性接触 + 姿态控制的完整闭环仿真。

```python
import math, sys, time
from pathlib import Path
import numpy as np
import pinocchio as pin
from pinocchio.visualize import MeshcatVisualizer

pinocchio_model_dir = Path(__file__).parent.parent / "models"
model_path = pinocchio_model_dir / "example-robot-data/robots"
mesh_dir   = pinocchio_model_dir
urdf_filename = "talos_reduced.urdf"
srdf_filename = "talos.srdf"

# 关键：浮动基 (FreeFlyer) 使机器人基座在空间中自由运动
# 不加此参数则基座固定，无法模拟真实人形
model, collision_model, visual_model = pin.buildModelsFromUrdf(
    model_path / "talos_data/robots" / urdf_filename,
    mesh_dir,
    pin.JointModelFreeFlyer()   # 人形必须！
)

# MeshCat 可视化器（浏览器中打开 URL 查看）
viz = MeshcatVisualizer(model, collision_model, visual_model)
viz.initViewer(open=False)
viz.loadViewerModel()

# 从 SRDF 加载参考配置（half_sitting：双腿微蹲标准站立姿）
pin.loadReferenceConfigurations(model, model_path / "talos_data/srdf" / srdf_filename)
q0 = model.referenceConfigurations["half_sitting"]
viz.display(q0)

# ---- 建立足底接触约束模型 ----
feet_name = ["left_sole_link", "right_sole_link"]
frame_ids  = [model.getFrameId(n) for n in feet_name]

contact_models = []
for frame_id in frame_ids:
    frame = model.frames[frame_id]
    # CONTACT_6D：6D 刚性接触（位置+姿态全约束，适合平足）
    contact_model = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_6D,
        model,
        frame.parentJoint,    # 接触发生在哪个关节
        frame.placement       # 接触点相对关节的位置
    )
    contact_models.append(contact_model)
    contact_datas = [cm.createData() for cm in contact_models]

num_constraints = len(feet_name)
contact_dim = 6 * num_constraints   # 每个 6D 接触提供 6 个约束

# 预分配接触约束内存
pin.initConstraintDynamics(model, data_sim, contact_models)

# ---- 仿真参数 ----
dt = 5e-3    # 5 ms 时间步（200 Hz）
# S：选择矩阵，将全身 nv 维度映射到驱动关节（浮动基 6 DOF 不可驱动）
S = np.zeros((model.nv - 6, model.nv))
S.T[6:, :] = np.eye(model.nv - 6)

# 姿态稳定控制器参数（PD 增益）
Kp_posture = 30.0
Kv_posture = 0.05 * math.sqrt(Kp_posture)  # 临界阻尼

q, v = q0.copy(), np.zeros(model.nv)
tau  = np.zeros(model.nv)

while t <= T:
    # ---- 控制律：关节空间 PD 稳定到参考姿态 ----
    pin.computeJointJacobians(model, data_control, q)

    # tau = S^T * Kp*(q_ref - q) - Kv*v  （只驱动非浮动基自由度）
    tau_control = S.T.dot(
        Kp_posture * pin.difference(model, q, q_ref) - Kv_posture * v[6:]
    )

    # ---- 带接触约束的正向动力学 ----
    # 求解 KKT 系统：
    # [M  J^T] [q̈  ] = [τ - C(q,v)v - g(q)]
    # [J  0  ] [λ  ]   [-J_dot * v          ]
    a = pin.constraintDynamics(
        model, data_sim, q, v, tau_control, contact_models, contact_datas
    )
    # a      → data_sim.ddq     (关节加速度)
    # lambda → data_sim.lambda_c (接触力)

    # ---- 半隐式欧拉积分 ----
    v += a * dt                          # 速度更新（用新加速度）
    q = pin.integrate(model, q, v * dt)  # 配置更新（流形积分！）

    viz.display(q)
    t += dt
```

**涉及算法**：浮动基建模、SRDF 参考配置、`RigidConstraintModel`、`constraintDynamics`（KKT 求解）、半隐式欧拉积分。

---

### 5.10 `talos-simulation.py` — KKT 约束姿态优化（静态）

```python
from time import sleep
import example_robot_data, numpy as np, pinocchio
from pinocchio.visualize import GepettoVisualizer

# 使用 example_robot_data 加载完整 Talos 模型
robot = example_robot_data.load("talos")
model, data = robot.model, robot.data

# 参考配置
q = robot.model.referenceConfigurations["half_sitting"]

# ---- 四肢（双足+双手）的 6D 接触约束 ----
foot_frames = ["left_sole_link", "right_sole_link",
               "gripper_left_fingertip_3_link", "gripper_right_fingertip_3_link"]
constraint_models = []
for frame_name in foot_frames:
    frame_id   = model.getFrameId(frame_name)
    joint_id   = model.frames[frame_id].parent
    contact_model = pinocchio.RigidConstraintModel(
        pinocchio.ContactType.CONTACT_6D, model,
        joint_id,
        model.frames[frame_id].placement,
        0,                    # 世界作为第二个约束体
        data.oMf[frame_id]    # 当前帧在世界系下的位姿
    )
    constraint_models.append(contact_model)

# 修改手部目标位置（抬起双臂的目标）
constraint_models[3].joint2_placement = pinocchio.SE3(
    pinocchio.rpy.rpyToMatrix(np.array([0., -np.pi/2, 0.])),
    np.array([0.6, -0.40, 1.0])
)

constraint_datas = [cm.createData() for cm in constraint_models]
kkt_constraint = pinocchio.ContactCholeskyDecomposition(model, constraint_models)
constraint_dim  = sum(cm.size() for cm in constraint_models)

# ---- squashing：迭代求解约束满足的关节配置 ----
# 在保持足底接触的同时，最小化质心高度变化
# KKT 求解方程：
#   [M   J^T] [dq ] = [rhs]
#   [J   0  ] [dy ]
def squashing(model, data, q_in):
    q = q_in.copy()
    y = np.ones(constraint_dim)   # Lagrange 乘子（接触力对偶变量）

    for k in range(N):
        pinocchio.computeAllTerms(model, data, q, np.zeros(model.nv))
        pinocchio.computeJointJacobians(model, data, q)

        com_act = data.com[0].copy()
        com_err = com_act - com_des(k)

        # 构造 Contact Cholesky 分解（稀疏 KKT 求解器）
        kkt_constraint.compute(model, data, constraint_models, constraint_datas, mu)

        # 约束残差：log6(c1Mc2) 衡量当前位姿与约束目标的 SE(3) 误差
        constraint_value = np.concatenate(
            [pinocchio.log6(cd.c1Mc2).vector for cd in constraint_datas]
        )
        J = np.vstack([
            pinocchio.getFrameJacobian(model, data, cm.joint1_id,
                                       cm.joint1_placement, cm.reference_frame)
            for cm in constraint_models
        ])

        # RHS：约束残差 + 质心误差 + 关节自由度为零
        rhs = np.concatenate([
            -constraint_value - y * mu,
            kp * mass * com_err,
            np.zeros(model.nv - 3)
        ])
        dz = kkt_constraint.solve(rhs)   # 稀疏 KKT 求解
        dy, dq = dz[:constraint_dim], dz[constraint_dim:]

        q = pinocchio.integrate(model, q, -alpha * dq)
        y -= alpha * (-dy + y)
        robot.display(q)
        sleep(0.05)
    return q
```

**涉及算法**：`ContactCholeskyDecomposition`、`computeAllTerms`、SE(3) 约束残差、质心控制。

---

### 5.11 `static-contact-dynamics.py` — 静态接触力计算

```python
# 目标：给定静止姿态 q0，计算维持静止所需的关节力矩和接触力
# 静态方程：g(q) = τ + J_c^T λ

# 策略：分块求解
# 对浮动基的 6 个方程（无驱动）：
#   g_bl = J_c_bl^T * λ  →  λ = pinv(J_c_bl^T) * g_bl
# 对驱动关节：
#   τ = g_j - J_c_j^T * λ

import numpy as np
import pinocchio as pin

# 加载 Solo-12 四足机器人（浮动基，4 个 3D 点接触）
model, _, _ = pin.buildModelsFromUrdf(urdf_model_path, mesh_dir, pin.JointModelFreeFlyer())
data = model.createData()

# 计算各项动力学量
pin.computeAllTerms(model, data, q0, np.zeros(model.nv))

# 计算各足 Jacobian（只取前 6 行对应浮动基）
Jc_feet_bl = ...   # 从 getFrameJacobian 提取

# 求接触力（伪逆）
lambda_contact = np.linalg.pinv(Jc_feet_bl.T) @ g_bl

# 求关节力矩
tau = g_j - Jc_feet_j.T @ lambda_contact
```

**涉及算法**：`computeAllTerms`、`getFrameJacobian`、静态平衡方程、伪逆。

---

### 5.12 `simulation-pendulum.py` — 倒立摆/多摆仿真

```python
# 从零手工构建多摆模型（无 URDF），演示 API 编程建模

model = pin.Model()
geom_model = pin.GeometryModel()
parent_id = 0

for k in range(N):
    # 添加绕 X 轴旋转的关节 (JointModelRX)
    joint_id = model.addJoint(parent_id, pin.JointModelRX(), joint_placement, f"joint_{k+1}")
    # 为关节绑定连杆惯量
    model.appendBodyToJoint(joint_id, body_inertia, body_placement)
    parent_id = joint_id

# 仿真主循环：ABA + 半隐式欧拉
for k in range(N_steps):
    tau_control = -damping * v          # 阻尼控制（耗散能量）
    a = pin.aba(model, data_sim, q, v, tau_control)  # ABA 正向动力学
    v += a * dt                         # 半隐式：先更新速度
    q = pin.integrate(model, q, v * dt) # 再用新速度更新配置（流形积分）
    viz.display(q)
```

**可选参数**：`--with-cart`（小车-摆系统），`-N 3`（三连摆）。
**涉及算法**：手工建模（`addJoint` / `appendBodyToJoint`）、ABA、半隐式欧拉。

---

### 5.13 `simulation-closed-kinematic-chains.py` — 闭环运动链仿真

```python
# 演示四杆机构（closed-loop mechanism）的约束仿真
# 闭环约束通过 RigidConstraintModel 实现

# 建模：先建开链，再用约束封闭
model = pin.Model()
# ... 添加连杆 A1, B1, B2, A2 ...

# 闭环约束：要求两个关节的相对位姿满足特定条件
loop_constraint = pin.RigidConstraintModel(
    pin.ContactType.CONTACT_6D,
    model, joint1_id, placement1,
         joint2_id, placement2
)

# 仿真：带约束的正向动力学
a = pin.constraintDynamics(model, data, q, v, tau, [loop_constraint], constraint_datas)
```

**人形应用**：并联腿结构（如 Atlas 的膝关节四杆机构）即需此类约束建模。

---

### 5.14 `run-algo-in-parallel.py` — OpenMP 并行批量动力学

```python
import numpy as np
import pinocchio as pin

# 人形随机模型
model = pin.buildSampleModelHumanoid()
pool = pin.ModelPool(model)      # 内部管理多个 Data 副本，每线程一个

num_threads = pin.omp_get_max_threads()
batch_size  = 128                # 一次并行计算 128 个配置

# 准备批量输入（每列是一个配置/速度/加速度）
q   = np.column_stack([pin.randomConfiguration(model) for _ in range(batch_size)])
v   = np.zeros((model.nv, batch_size))
a   = np.zeros((model.nv, batch_size))
tau = np.zeros((model.nv, batch_size))

# 并行 RNEA：一次调用计算 128 组逆动力学（OpenMP 自动分配线程）
res_rnea = pin.rneaInParallel(num_threads, pool, q, v, a)

# 并行 ABA：一次调用计算 128 组正向动力学
res_aba  = pin.abaInParallel(num_threads, pool, q, v, tau)
```

**应用场景**：MPC 中的多轨迹并行仿真、强化学习中的批量仿真采样。

---

### 5.15 `build-reduced-model.py` — 降阶模型（锁定关节）

```python
import numpy as np
import pinocchio as pin

# 加载完整模型
model, collision_model, visual_model = pin.buildModelsFromUrdf(urdf_filename, mesh_dir)

# 要锁定（固定）的关节名称列表
jointsToLock = ["wrist_1_joint", "wrist_2_joint", "wrist_3_joint"]
jointsToLockIDs = [model.getJointId(jn) for jn in jointsToLock if model.existJointName(jn)]

# 锁定位置的参考配置
initialJointConfig = np.array([0, 0, 0, 1, 1, 1])

# Option 1：只降阶运动学模型
model_reduced = pin.buildReducedModel(model, jointsToLockIDs, initialJointConfig)

# Option 2：同时降阶视觉模型（显示用）
model_reduced, visual_model_reduced = pin.buildReducedModel(
    model, visual_model, jointsToLockIDs, initialJointConfig
)

# Option 3：同时降阶多个几何模型
model_reduced, (visual_model_reduced, collision_model_reduced) = pin.buildReducedModel(
    model,
    list_of_geom_models=[visual_model, collision_model],
    list_of_joints_to_lock=jointsToLockIDs,
    reference_configuration=initialJointConfig,
)
```

**人形应用**：将手臂锁定在某姿态，专注腿部控制研究；或将整个上身锁定，
验证纯腿部动力学。

---

### 5.16 `mimic_dynamics.py` — 关节 Mimic（传动比约束）

```python
import numpy as np
import pinocchio

# 从 URDF 直接解析 mimic 标签（如夹爪耦合关节）
model_mimic = pinocchio.buildModelFromUrdf(model_path, mimic=True)

# 或手工设置 mimic 关系：joint 10 = -1.0 × joint 9 + 0.0
model_mimic = pinocchio.transformJointIntoMimic(model_full, 9, 10, -1.0, 0.0)

# G 矩阵：将 mimic 模型的速度映射到完整模型的速度
# 对于非 mimic 关节：G[i,i]=1；对 mimic 关节：G[i,j]=scaling
G = np.zeros([model_mimic.nv, model_full.nv])
for i in range(model_full.njoints):
    scaling = getattr(model_mimic.joints[i].extract(), "scaling", None)
    G[model_mimic.joints[i].idx_v, model_full.joints[i].idx_v] = scaling if scaling else 1

# 关系：τ_mimic = G @ τ_full，M_mimic = G @ M_full @ G^T
tau_mimic = pinocchio.rnea(model_mimic, data_mimic, q_mimic, v_mimic, a_mimic)
M_mimic   = pinocchio.crba(model_mimic, data_mimic, q_mimic)
```

**人形应用**：人形机器人手部（夹爪）、踝关节并联结构中常有传动比耦合关节。

---

### 5.17 `contact-cholesky.py` — 接触 Cholesky 分解（Delassus 算子）

```python
import pinocchio as pin

# ANYmal 四足（4 个 3D 点接触）
model = pin.buildModelFromUrdf(urdf_filename)
data  = model.createData()
q0    = pin.neutral(model)

# 构造接触约束（3D 点接触）
contact_models = []
for fid in feet_frame_ids:
    frame = model.frames[fid]
    cmodel = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,
        frame.parent, frame.placement,
        pin.LOCAL_WORLD_ALIGNED,
    )
    contact_models.append(cmodel)

contact_data = [cm.createData() for cm in contact_models]

# 计算质量矩阵（CRBA）
pin.crba(model, data, q0)

# Contact Cholesky 分解：高效求解 KKT 系统
# 利用运动树稀疏性，比 dense LU 快数倍
data.contact_chol.compute(model, data, contact_models, contact_data)

# Delassus 矩阵（操作空间惯量矩阵逆）：Λ^{-1} = J M^{-1} J^T
delassus_matrix     = data.contact_chol.getInverseOperationalSpaceInertiaMatrix()
# 操作空间惯量矩阵：Λ = (J M^{-1} J^T)^{-1}
delassus_matrix_inv = data.contact_chol.getOperationalSpaceInertiaMatrix()
```

**Delassus 矩阵物理含义**：接触空间中的有效质量（等效惯量），用于求接触力响应。

---

## 6. Demo 注释：C++

### 6.1 `overview-simple.cpp` — C++ 最小示例

```cpp
#include <iostream>
#include "pinocchio/multibody/sample-models.hpp"   // 内置样例模型
#include "pinocchio/algorithm/joint-configuration.hpp"  // neutral()
#include "pinocchio/algorithm/rnea.hpp"            // rnea()

int main()
{
    pinocchio::Model model;
    // 构建内置 6-DOF 串联机械臂（等同 Python buildSampleModelManipulator）
    pinocchio::buildModels::manipulator(model);
    pinocchio::Data data(model);   // 算法数据容器

    // 零位配置（关节全零）
    Eigen::VectorXd q = pinocchio::neutral(model);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);  // 零速度
    Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);  // 零加速度

    // RNEA 逆动力学：a=0 时计算重力补偿力矩 g(q)
    const Eigen::VectorXd & tau = pinocchio::rnea(model, data, q, v, a);
    std::cout << "tau = " << tau.transpose() << std::endl;
}
```

---

### 6.2 `overview-urdf.cpp` — 从 URDF 加载并做正向运动学

```cpp
#include "pinocchio/parsers/urdf.hpp"              // URDF 解析器
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"      // forwardKinematics()
#include <iostream>

int main(int argc, char ** argv)
{
    using namespace pinocchio;

    // URDF 路径（可命令行传入，默认 UR5）
    const std::string urdf_filename = (argc <= 1)
        ? PINOCCHIO_MODEL_DIR + std::string(".../ur5_robot.urdf")
        : argv[1];

    Model model;
    // 加载 URDF → 构建运动学/动力学模型
    // 浮动基人形需改为：urdf::buildModel(urdf, JointModelFreeFlyer(), model)
    pinocchio::urdf::buildModel(urdf_filename, model);
    std::cout << "model name: " << model.name << std::endl;

    Data data(model);

    // 在关节限位内随机采样配置（四元数部分自动归一化）
    Eigen::VectorXd q = randomConfiguration(model);

    // 正向运动学：递推更新 data.oMi[i]（所有关节的 SE(3) 位姿）
    forwardKinematics(model, data, q);

    // 打印各关节在世界系下的平移坐标
    for (JointIndex id = 0; id < (JointIndex)model.njoints; ++id)
        std::cout << model.names[id] << ": "
                  << data.oMi[id].translation().transpose() << std::endl;
}
```

---

### 6.3 `overview-SE3.cpp` — SE(3) 李群运算

```cpp
#include "pinocchio/multibody/liegroup/liegroup.hpp"
using namespace pinocchio;

int main()
{
    // SE(3) 操作对象（3D 特殊欧氏群）
    // 配置向量格式：[x, y, z, qx, qy, qz, qw]（7 维，四元数表示旋转）
    typedef SpecialEuclideanOperationTpl<3, double> SE3Operation;
    SE3Operation aSE3;

    SE3Operation::ConfigVector_t pose_s, pose_g;
    // 起始位姿：位置 [1,1,1]，旋转用四元数表示
    pose_s << 1.0, 1.0, 1.0, -0.13795, 0.13795, 0.69352, 0.69352;
    // 目标位姿
    pose_g << 4.0, 3.0, 3.0, -0.46194, 0.331414, 0.800103, 0.191342;

    // 归一化四元数（确保旋转合法）
    aSE3.normalize(pose_s);
    aSE3.normalize(pose_g);

    // difference：在切空间（李代数）中计算两个 SE(3) 元素的差
    // delta_u ∈ R^6（平移 3 + 旋转 3），是 pose_g 相对 pose_s 的"速度"
    SE3Operation::TangentVector_t delta_u;
    aSE3.difference(pose_s, pose_g, delta_u);

    // integrate：从 pose_s 沿 delta_u 做李群积分，应恢复 pose_g
    SE3Operation::ConfigVector_t pose_check;
    aSE3.integrate(pose_s, delta_u, pose_check);
    // 验证：pose_check ≈ pose_g
}
```

**数学含义**：`difference(a, b)` = $\log(a^{-1} \cdot b)$；`integrate(a, u)` = $a \cdot \exp(u)$。
**人形应用**：浮动基状态的差分（用于误差计算）和积分（用于状态更新）。

---

### 6.4 `overview-lie.cpp` — SE(2) 李群（平面运动）

```cpp
#include "pinocchio/multibody/liegroup/liegroup.hpp"
using namespace pinocchio;

int main()
{
    // SE(2)：2D 特殊欧氏群（平面移动机器人、足底在地面上的运动）
    // 配置向量：[x, y, cos(θ), sin(θ)]（4 维）
    typedef SpecialEuclideanOperationTpl<2, double> SE2Operation;
    SE2Operation aSE2;

    SE2Operation::ConfigVector_t pose_s, pose_g;
    // 起始：位置 [1,1]，朝向 45°
    pose_s << 1.0, 1.0, cos(M_PI/4.0), sin(M_PI/4.0);
    // 目标：位置 [3,-1]，朝向 -90°
    pose_g << 3.0, -1.0, cos(-M_PI/2.0), sin(-M_PI/2.0);

    // 在李代数（切空间）中计算差分
    // delta_u ∈ R^3（Δx, Δy, Δθ），代表从 pose_s 到 pose_g 的最短路径
    SE2Operation::TangentVector_t delta_u;
    delta_u.setZero();
    aSE2.difference(pose_s, pose_g, delta_u);

    // 验证：integrate(pose_s, delta_u) ≈ pose_g
    SE2Operation::ConfigVector_t pose_check;
    aSE2.integrate(pose_s, delta_u, pose_check);
}
```

**人形应用**：足步规划中步态在地面投影的差分；ZMP 轨迹在水平面的表示。

---

### 6.5 `inverse-kinematics.cpp` — C++ 6D IK

```cpp
#include "pinocchio/spatial/explog.hpp"      // log6(), Jlog6()
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"   // computeJointJacobian()
#include "pinocchio/algorithm/joint-configuration.hpp"

int main(int, char**)
{
    pinocchio::Model model;
    pinocchio::buildModels::manipulator(model);
    pinocchio::Data data(model);

    const int JOINT_ID = 6;
    // 目标位姿：单位旋转 + 平移 [1, 0, 1]
    const pinocchio::SE3 oMdes(Eigen::Matrix3d::Identity(), Eigen::Vector3d(1., 0., 1.));

    Eigen::VectorXd q = pinocchio::neutral(model);
    const double eps = 1e-4, DT = 1e-1, damp = 1e-6;

    // 6×nv Jacobian 矩阵预分配
    pinocchio::Data::Matrix6x J(6, model.nv);
    J.setZero();

    typedef Eigen::Matrix<double, 6, 1> Vector6d;
    Vector6d err;
    Eigen::VectorXd v(model.nv);

    for (int i = 0;; i++)
    {
        pinocchio::forwardKinematics(model, data, q);

        // 局部系下的相对变换：iMd = T_i^{-1} * T_des
        const pinocchio::SE3 iMd = data.oMi[JOINT_ID].actInv(oMdes);

        // 6D 误差：SE(3) 对数映射 → 切向量
        err = pinocchio::log6(iMd).toVector();
        if (err.norm() < eps) { success = true; break; }
        if (i >= IT_MAX)      { success = false; break; }

        // 计算局部系下的几何 Jacobian
        pinocchio::computeJointJacobian(model, data, q, JOINT_ID, J);

        // Jlog6：log 映射关于 SE(3) 参数的 Jacobian（6×6 矩阵）
        // 链式法则：J_eff = Jlog6(iMd^{-1}) * J_joint
        pinocchio::Data::Matrix6 Jlog;
        pinocchio::Jlog6(iMd.inverse(), Jlog);
        J = -Jlog * J;

        // 阻尼最小二乘：v = -J^T (JJ^T + λI)^{-1} e
        // JJt.ldlt().solve() 用 LDLT 分解高效求解对称正定线性系统
        pinocchio::Data::Matrix6 JJt;
        JJt.noalias() = J * J.transpose();
        JJt.diagonal().array() += damp;
        v.noalias() = -J.transpose() * JJt.ldlt().solve(err);

        // 李群积分（关键！不能用 q += v*DT）
        q = pinocchio::integrate(model, q, v * DT);
    }
    std::cout << "result: " << q.transpose() << std::endl;
}
```

---

### 6.6 `inverse-kinematics-3d.cpp` — C++ 3D 位置 IK

```cpp
// 与 6D IK 相比，只跟踪末端位置（不控制姿态）
// 误差降为 3D，Jacobian 只取线速度行

for (int i = 0;; i++)
{
    pinocchio::forwardKinematics(model, data, q);
    const pinocchio::SE3 iMd = data.oMi[JOINT_ID].actInv(oMdes);

    // 误差 = 相对变换的平移部分（R^3，非 R^6）
    err = iMd.translation();
    if (err.norm() < eps) { success = true; break; }

    pinocchio::computeJointJacobian(model, data, q, JOINT_ID, joint_jacobian);
    // 只取 Jacobian 上 3 行（线速度），取负（局部系符号约定）
    const auto J = -joint_jacobian.topRows<3>();

    // 3×3 阻尼最小二乘（比 6×6 计算量更小）
    const Eigen::Matrix3d JJt = J * J.transpose() + damp * Eigen::Matrix3d::Identity();
    v.noalias() = -J.transpose() * JJt.ldlt().solve(err);
    q = pinocchio::integrate(model, q, v * DT);
}
```

---

### 6.7 `inverse-dynamics.cpp` — C++ RNEA

```cpp
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

int main(int argc, char** argv)
{
    using namespace pinocchio;

    // 加载 URDF 模型（固定基）
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    // 随机状态
    Eigen::VectorXd q   = randomConfiguration(model);  // 随机关节角 [rad]
    Eigen::VectorXd v   = Eigen::VectorXd::Zero(model.nv); // 零速度
    Eigen::VectorXd a   = Eigen::VectorXd::Zero(model.nv); // 零加速度

    // RNEA：τ = M(q)a + C(q,v)v + g(q)
    // a=0, v=0 时结果即为重力补偿力矩 g(q)
    Eigen::VectorXd tau = pinocchio::rnea(model, data, q, v, a);
    // 结果同时存在 data.tau
    std::cout << "Joint torques: " << data.tau.transpose() << std::endl;
}
```

---

### 6.8 `inverse-dynamics-derivatives.cpp` — C++ RNEA 梯度

```cpp
#include "pinocchio/algorithm/rnea-derivatives.hpp"

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    Eigen::VectorXd q = randomConfiguration(model);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
    Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);

    // 预分配梯度矩阵（nv × nv）
    Eigen::MatrixXd dtau_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
    Eigen::MatrixXd dtau_dv = Eigen::MatrixXd::Zero(model.nv, model.nv);
    Eigen::MatrixXd dtau_da = Eigen::MatrixXd::Zero(model.nv, model.nv);

    // 一次调用同时计算 RNEA 结果和三个梯度矩阵
    // dtau_da = M(q)（质量矩阵，同时获得正向运动学结果）
    computeRNEADerivatives(model, data, q, v, a, dtau_dq, dtau_dv, dtau_da);

    // data.tau 存 RNEA 结果，dtau_dq/dv/da 存梯度
    std::cout << "Joint torque: " << data.tau.transpose() << std::endl;
}
```

---

### 6.9 `forward-dynamics-derivatives.cpp` — C++ ABA 梯度

```cpp
#include "pinocchio/algorithm/aba-derivatives.hpp"

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    Eigen::VectorXd q   = randomConfiguration(model);
    Eigen::VectorXd v   = Eigen::VectorXd::Zero(model.nv);
    Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);

    // 预分配梯度矩阵
    Eigen::MatrixXd ddq_dq   = Eigen::MatrixXd::Zero(model.nv, model.nv);
    Eigen::MatrixXd ddq_dv   = Eigen::MatrixXd::Zero(model.nv, model.nv);
    Eigen::MatrixXd ddq_dtau = Eigen::MatrixXd::Zero(model.nv, model.nv);

    // 一次调用：计算 ABA（正向动力学）及其对 (q, v, τ) 的解析梯度
    // ddq_dtau = M(q)^{-1}（质量矩阵的逆，即 DDP 中的 B 矩阵）
    computeABADerivatives(model, data, q, v, tau, ddq_dq, ddq_dv, ddq_dtau);

    std::cout << "Joint acceleration: " << data.ddq.transpose() << std::endl;
    // DDP 线性化：A = [ddq_dq, ddq_dv], B = ddq_dtau = Minv
}
```

---

### 6.10 `kinematics-derivatives.cpp` — C++ 运动学梯度

```cpp
#include "pinocchio/algorithm/kinematics-derivatives.hpp"

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    Eigen::VectorXd q = randomConfiguration(model);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
    Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);

    // 一次调用计算所有关节速度/加速度关于 (q,v,a) 的梯度
    computeForwardKinematicsDerivatives(model, data, q, v, a);

    // 取最后一个关节（末端）
    JointIndex joint_id = (JointIndex)(model.njoints - 1);
    Data::Matrix6x v_partial_dq(6, model.nv), a_partial_dq(6, model.nv),
                   a_partial_dv(6, model.nv), a_partial_da(6, model.nv);
    v_partial_dq.setZero(); a_partial_dq.setZero();
    a_partial_dv.setZero(); a_partial_da.setZero();

    // 局部系下：∂v_i/∂q, ∂a_i/∂q, ∂a_i/∂v, ∂a_i/∂a（=几何 Jacobian）
    getJointAccelerationDerivatives(
        model, data, joint_id, LOCAL,
        v_partial_dq, a_partial_dq, a_partial_dv, a_partial_da
    );

    // 世界系下：同上，但量在世界坐标系中表达
    getJointAccelerationDerivatives(
        model, data, joint_id, WORLD,
        v_partial_dq, a_partial_dq, a_partial_dv, a_partial_da
    );
    // 注意：v_partial_dv = a_partial_da（两者相等，不重复计算）
}
```

---

### 6.11 `build-reduced-model.cpp` — C++ 降阶模型

```cpp
#include "pinocchio/algorithm/model.hpp"   // buildReducedModel()

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // 方案一：指定要锁定的关节名称列表
    std::vector<std::string> joints_to_lock = {
        "elbow_joint", "wrist_3_joint", "wrist_2_joint"
    };
    std::vector<JointIndex> lock_ids;
    for (auto& name : joints_to_lock)
        if (model.existJointName(name))
            lock_ids.push_back(model.getJointId(name));

    // 锁定参考配置（从完整模型中随机采样）
    Eigen::VectorXd q_rand = randomConfiguration(model);

    // 构建降阶模型：被锁定的关节固定在 q_rand 对应的角度
    Model reduced_model = buildReducedModel(model, lock_ids, q_rand);

    // 方案二：指定要保留（不锁定）的关节
    // 先计算补集，再调用 buildReducedModel
    std::vector<std::string> joints_to_keep = {
        "shoulder_pan_joint", "shoulder_lift_joint", "wrist_1_joint"
    };
    // ... 计算 lock_ids = 全集 - joints_to_keep ...
    Model reduced_model2 = buildReducedModel(model, lock_ids, q_rand);
}
```

---

### 6.12 `interpolation-SE3.cpp` — SE(3) 插值

```cpp
#include "pinocchio/multibody/liegroup/liegroup.hpp"
using namespace pinocchio;

int main()
{
    // SE(3) 插值：在两个位姿之间做测地线（最短路径）插值
    // 类比：球面插值 SLERP 之于 SO(3)，这是 SE(3) 上的扩展
    typedef SpecialEuclideanOperationTpl<3, double> SE3Operation;
    SE3Operation aSE3;

    SE3Operation::ConfigVector_t pose_s, pose_g;
    pose_s << 1.0, 1.0, 1.0, -0.13795, 0.13795, 0.69352, 0.69352;
    pose_g << 4.0, 3.0, 3.0, -0.46194, 0.331414, 0.800103, 0.191342;
    aSE3.normalize(pose_s);
    aSE3.normalize(pose_g);

    // 在 SE(3) 流形上插值，t=0.5 即取中间位姿
    // 数学：pole_u = pose_s * exp(0.5 * log(pose_s^{-1} * pose_g))
    SE3Operation::ConfigVector_t pole_u;
    aSE3.interpolate(pose_s, pose_g, 0.5, pole_u);
    std::cout << "Interpolated: " << pole_u.transpose() << std::endl;
}
```

**人形应用**：末端轨迹在 SE(3) 上的插值（手部抓取路径平滑）、
步态规划中足底位姿序列的插值。

---

### 6.13 `geometry-models.cpp` — 碰撞/视觉几何模型

```cpp
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/geometry.hpp"   // updateGeometryPlacements()

int main(int argc, char** argv)
{
    using namespace pinocchio;
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // 从 URDF 加载碰撞几何（用于碰撞检测，hp-fcl 格式）
    GeometryModel collision_model;
    pinocchio::urdf::buildGeom(model, urdf_filename, COLLISION, collision_model, mesh_dir);

    // 从 URDF 加载视觉几何（用于可视化渲染，mesh 格式）
    GeometryModel visual_model;
    pinocchio::urdf::buildGeom(model, urdf_filename, VISUAL, visual_model, mesh_dir);

    Data data(model);
    GeometryData collision_data(collision_model);
    GeometryData visual_data(visual_model);

    Eigen::VectorXd q = randomConfiguration(model);

    // 正向运动学更新关节位姿
    forwardKinematics(model, data, q);

    // 根据关节位姿更新几何对象在世界系下的位置（data.oMg[i]）
    updateGeometryPlacements(model, data, collision_model, collision_data);
    updateGeometryPlacements(model, data, visual_model,    visual_data);

    // data.oMg[geom_id] 即几何对象的世界系 SE(3) 变换
}
```

---

### 6.14 `collisions.cpp` — 碰撞检测（Talos + SRDF）

```cpp
#include "pinocchio/collision/collision.hpp"   // computeCollisions()
#include "pinocchio/parsers/srdf.hpp"          // removeCollisionPairs()

int main()
{
    using namespace pinocchio;

    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    GeometryModel geom_model;
    pinocchio::urdf::buildGeom(model, urdf_filename, COLLISION, geom_model, model_path);

    // 添加所有可能的碰撞对（O(n^2) 对）
    geom_model.addAllCollisionPairs();

    // 从 SRDF 移除已知的"允许自碰撞"对（如相邻连杆）
    // 减少不必要的碰撞检测计算
    pinocchio::srdf::removeCollisionPairs(model, geom_model, srdf_filename);

    GeometryData geom_data(geom_model);

    // 从 SRDF 加载 half_sitting 参考配置
    pinocchio::srdf::loadReferenceConfigurations(model, srdf_filename);
    const auto& q = model.referenceConfigurations["half_sitting"];

    // 全量碰撞检测（检测所有配对）
    computeCollisions(model, data, geom_model, geom_data, q);

    for (size_t k = 0; k < geom_model.collisionPairs.size(); ++k)
    {
        const auto& cp = geom_model.collisionPairs[k];
        const auto& cr = geom_data.collisionResults[k];
        std::cout << "pair " << cp.first << "," << cp.second
                  << ": " << (cr.isCollision() ? "COLLISION" : "ok") << std::endl;
    }

    // 早停版本：发现第一个碰撞即停止（更快）
    computeCollisions(model, data, geom_model, geom_data, q, /*stopAtFirst=*/true);

    // 单对检测（只检测第 3 对）
    const PairIndex pair_id = 2;
    Eigen::VectorXd q_neutral = neutral(model);
    updateGeometryPlacements(model, data, geom_model, geom_data, q_neutral);
    computeCollision(geom_model, geom_data, pair_id);
}
```

---

### 6.15 `multiprecision.cpp` — 多精度浮点运算

```cpp
#include "pinocchio/math/multiprecision.hpp"
#include <boost/multiprecision/cpp_dec_float.hpp>

int main(int argc, char** argv)
{
    using namespace pinocchio;

    // 标准 double 精度模型
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    Data data(model);

    // 定义 100 位十进制精度的浮点类型
    typedef boost::multiprecision::number<
        boost::multiprecision::cpp_dec_float<100>,
        boost::multiprecision::et_off> float_100;

    // 将 double 模型转换为 100 位精度模型
    // Pinocchio 的模板设计使得标量类型可以任意替换
    typedef ModelTpl<float_100> ModelMulti;
    typedef DataTpl<float_100>  DataMulti;
    ModelMulti model_multi = model.cast<float_100>();
    DataMulti  data_multi(model_multi);

    // 用 double 配置生成 100 位精度配置
    ModelMulti::ConfigVectorType q_multi  = randomConfiguration(model_multi);
    ModelMulti::TangentVectorType v_multi = ModelMulti::TangentVectorType::Random(model.nv);
    ModelMulti::TangentVectorType a_multi = ModelMulti::TangentVectorType::Random(model.nv);

    // 同时用两种精度计算 RNEA，对比结果
    rnea(model,       data,       q_multi.cast<double>(), ...);
    rnea(model_multi, data_multi, q_multi,                ...);

    // 100 位精度结果，可用于数值分析、梯度验证（有限差分基准）
    std::cout << std::setprecision(std::numeric_limits<float_100>::max_digits10)
              << data_multi.tau << std::endl;
}
```

**应用场景**：验证解析梯度的精度（与高精度数值梯度对比）、数值分析研究。

---

## 7. 人形机器人开发路径

按以下顺序逐步深入：

```
阶段 1：SE(3) 与李群基础
  → overview-SE3.cpp / overview-lie.cpp
  → 理解 difference / integrate / log / exp

阶段 2：浮动基建模 + FK
  → overview-urdf.py（加 JointModelFreeFlyer）
  → simulation-contact-dynamics.py 前半段（加载 Talos）

阶段 3：逆向动力学与梯度
  → inverse-dynamics.py / .cpp
  → inverse-dynamics-derivatives.py / .cpp（理解 ∂τ/∂q）

阶段 4：质心动量学
  → forward-dynamics-derivatives.py（buildSampleModelHumanoidRandom）
  → 重点：data.hg, data.Ag, data.com

阶段 5：接触约束动力学（核心）
  → simulation-contact-dynamics.py（完整 Talos 仿真）
  → contact-cholesky.py（Delassus 矩阵理解）
  → talos-simulation.py（KKT 约束优化）

阶段 6：轨迹优化接口
  → forward-dynamics-derivatives.py（ABA 梯度 → DDP 的 A, B 矩阵）
  → 接入 Crocoddyl：将 ABA Derivatives 作为动力学线性化来源

阶段 7：高性能应用
  → run-algo-in-parallel.py（并行批量仿真）
  → examples/codegen/（代码生成，用于嵌入式实时控制）
```

---

## 8. 与工业机器人的对比


| 维度         | 工业机器人                       | 人形机器人（Pinocchio）                   |
| ------------ | -------------------------------- | ----------------------------------------- |
| 基座约束     | 固定（地面）                     | 浮动基（空间自由）                        |
| 配置空间     | $q \in \mathbb{R}^n$             | $q \in SE(3) \times \mathbb{R}^n$（流形） |
| 状态更新     | $q \mathrel{+}= \dot{q}\Delta t$ | `pin.integrate(model, q, v*dt)`           |
| 驱动性       | 全驱动                           | **欠驱动**（浮动基 6 DOF 无输入）         |
| 接触         | 无/固定末端                      | 多点动态接触（KKT 约束）                  |
| 建模格式     | URDF（固定基）                   | URDF +`JointModelFreeFlyer()`             |
| 控制目标     | 末端轨迹跟踪                     | 质心轨迹 + 接触力 + 任务层次              |
| 关键算法     | FK / IK / RNEA                   | 质心动量 + 接触约束动力学 + 解析梯度      |
| 主要应用框架 | ROS MoveIt                       | Crocoddyl, HPP, Stack-of-Tasks            |

---

## 9. 关键头文件索引


| 功能              | 头文件                                                   |
| ----------------- | -------------------------------------------------------- |
| 正向运动学        | `include/pinocchio/algorithm/kinematics.hpp`             |
| 运动学梯度        | `include/pinocchio/algorithm/kinematics-derivatives.hpp` |
| Jacobian 计算     | `include/pinocchio/algorithm/jacobian.hpp`               |
| Frame 运动学      | `include/pinocchio/algorithm/frames.hpp`                 |
| Frame 梯度        | `include/pinocchio/algorithm/frames-derivatives.hpp`     |
| RNEA 逆动力学     | `include/pinocchio/algorithm/rnea.hpp`                   |
| RNEA 梯度         | `include/pinocchio/algorithm/rnea-derivatives.hpp`       |
| ABA 正向动力学    | `include/pinocchio/algorithm/aba.hpp`                    |
| ABA 梯度          | `include/pinocchio/algorithm/aba-derivatives.hpp`        |
| 质量矩阵 CRBA     | `include/pinocchio/algorithm/crba.hpp`                   |
| 质心动量          | `include/pinocchio/algorithm/centroidal.hpp`             |
| 质心动量梯度      | `include/pinocchio/algorithm/centroidal-derivatives.hpp` |
| 质心位置          | `include/pinocchio/algorithm/center-of-mass.hpp`         |
| 接触约束动力学    | `include/pinocchio/algorithm/constrained-dynamics.hpp`   |
| 接触 Cholesky     | `include/pinocchio/algorithm/contact-cholesky.hpp`       |
| 接触信息          | `include/pinocchio/algorithm/contact-info.hpp`           |
| SE(3) 类型        | `include/pinocchio/spatial/se3.hpp`                      |
| 空间运动（Twist） | `include/pinocchio/spatial/motion.hpp`                   |
| 空间力（Wrench）  | `include/pinocchio/spatial/force.hpp`                    |
| 空间惯量          | `include/pinocchio/spatial/inertia.hpp`                  |
| 指数/对数映射     | `include/pinocchio/spatial/explog.hpp`                   |
| 李群操作          | `include/pinocchio/multibody/liegroup/liegroup.hpp`      |
| 降阶模型          | `include/pinocchio/algorithm/model.hpp`                  |
| URDF 解析         | `include/pinocchio/parsers/urdf.hpp`                     |
| SRDF 解析         | `include/pinocchio/parsers/srdf.hpp`                     |
| MJCF 解析         | `include/pinocchio/parsers/mjcf.hpp`                     |
| 碰撞检测          | `include/pinocchio/collision/collision.hpp`              |
| 并行算法          | `include/pinocchio/algorithm/parallel/`                  |
| 多精度支持        | `include/pinocchio/math/multiprecision.hpp`              |

> **重要**：上表的 `.hpp` 只含**声明与文档注释**，真正的算法实现体在同名的 `.hxx` 里，
> 位于 `include/pinocchio/src/algorithm/`（例如 `kinematics.hpp` 的实现在
> `include/pinocchio/src/algorithm/kinematics.hxx`）。阅读机制见下一节 [11.1](#111-三层文件布局-hpp--hxx--cpp)。

---

## 10. 源码目录地图

> 目标：拿到仓库后，1 分钟内知道"某类东西该去哪个目录找"。

```
pinocchio/
├── include/pinocchio/            # 头文件（库的主体，几乎全是模板）
│   ├── fwd.hpp                   # 全局前置声明、宏、Options 枚举
│   ├── context.hpp               # 决定默认标量类型（double）的开关，见 11.2
│   ├── macros.hpp                # PINOCCHIO_* 宏（断言、DLL 导出、EIGEN 检查）
│   │
│   ├── spatial/                  # ★ 空间代数（数学地基，最先读）
│   │   ├── se3.hpp               #   SE3Tpl：刚体位姿 R,p + act/actInv
│   │   ├── motion.hpp            #   MotionTpl：空间速度 Twist (v,ω) 线性在前
│   │   ├── force.hpp             #   ForceTpl：空间力 Wrench (f,τ) 线性在前
│   │   ├── inertia.hpp           #   InertiaTpl：6×6 空间惯量
│   │   └── explog.hpp            #   log3/exp3/log6/exp6/Jlog6（李群映射）
│   │
│   ├── multibody/                # ★ 建模层：关节 / 模型 / 数据 / 访问者
│   │   ├── joint/                #   各关节类型（RX/FreeFlyer/Spherical…）
│   │   ├── joint.hpp             #   JointModel/JointData 的 variant 封装
│   │   ├── joint-motion-subspace.hpp  # 关节子空间 S（RNEA 中的 sᵢ）
│   │   ├── visitor.hpp           #   JointUnaryVisitorBase：算法遍历的引擎
│   │   ├── liegroup/             #   李群操作（integrate/difference/interpolate）
│   │   ├── sample-models.hpp     #   buildSampleModelManipulator / Humanoid
│   │   └── pool.hpp              #   ModelPool：多线程 Data 副本
│   │
│   ├── algorithm/                # ★ 算法层：声明 + 文档（实现在 src/*.hxx）
│   │   ├── kinematics.hpp        #   FK
│   │   ├── rnea.hpp / aba.hpp    #   逆/正动力学
│   │   ├── crba.hpp              #   质量矩阵
│   │   ├── jacobian.hpp          #   Jacobian
│   │   ├── centroidal.hpp        #   质心动量
│   │   ├── constrained-dynamics.hpp / contact-cholesky.hpp  # 接触
│   │   ├── model.hpp             #   ModelTpl 类本体、buildReducedModel
│   │   └── *-derivatives.hpp     #   各算法的解析梯度
│   │
│   ├── src/                      # ★★ 模板"实现体"（.hxx），与上面 .hpp 一一对应
│   │   ├── algorithm/*.hxx       #   算法真正的循环/递推在这里
│   │   ├── multibody/*.hxx       #   ModelTpl/DataTpl 的成员实现
│   │   ├── multibody/joint/*.hxx #   每种关节的 calc()/calc_aba()
│   │   └── context/              #   默认标量 (default.hxx) 等编译上下文
│   │
│   ├── parsers/                  # URDF / SRDF / MJCF / sdf 解析
│   ├── collision/                # 碰撞检测（封装 Coal / hpp-fcl）
│   ├── autodiff/                 # CppAD / CasADi 自动微分标量支持
│   ├── codegen/                  # 代码生成（编译期特化为 C 代码）
│   └── bindings/python/          # Python 绑定的头文件（Boost.Python 包装）
│
├── src/                          # ★★ .cpp：对默认标量做"显式模板实例化"
│   └── algorithm/*.cpp           #   编进 libpinocchio，见 11.1
│
├── bindings/python/              # Python 绑定的 .cpp（expose* 函数）
├── examples/                     # 本指南第 5、6 章注释的全部 demo
├── unittest/                     # 单元测试（读实现时的"用法真值"）
├── doc/                          # 官方 doxygen 文档源（a-features/ 值得读）
└── models/                       # 示例 URDF（UR5、Talos…）
```

**记忆口诀**：想看**怎么用**→`examples/`、`unittest/`；想看**接口和公式注释**→`algorithm/*.hpp`；
想看**实现细节**→`src/algorithm/*.hxx`；想看**数据结构**→`spatial/`、`src/multibody/`。

---

## 11. 源码架构三大机制

Pinocchio 是一个"模板元编程 + 零成本抽象"的库。不先理解下面三个机制，直接看实现会寸步难行。

### 11.1 三层文件布局 .hpp / .hxx / .cpp

同一个算法被拆到**三个文件**，各司其职：


| 文件      | 位置                               | 内容                                    | 你什么时候看它                         |
| --------- | ---------------------------------- | --------------------------------------- | -------------------------------------- |
| `xxx.hpp` | `include/pinocchio/algorithm/`     | **函数声明 + Doxygen 注释 + 公式**      | 想知道"有哪些接口、参数含义、数学定义" |
| `xxx.hxx` | `include/pinocchio/src/algorithm/` | **模板函数的实现体**（真正的循环/递推） | 想知道"算法到底怎么算的"               |
| `xxx.cpp` | `src/algorithm/`                   | **对默认标量 `double` 的显式实例化**    | 几乎不用看，只是让编译加速、生成`.so`  |

三者的连接方式（以运动学为例）：

```cpp
// kinematics.hpp 末尾：把实现体拉进来
#include "pinocchio/src/algorithm/kinematics.hxx"

// kinematics.hxx 顶部：标注"这是私有实现，别直接 include 我"
// IWYU pragma: private, include "pinocchio/algorithm/kinematics.hpp"

// kinematics.cpp：对 double 版本做一次显式实例化，编进库
template ... void forwardKinematics<context::Scalar, ...>(const Model&, Data&, ...);
```

> **为什么这么拆？** 模板函数的定义必须对调用者可见（否则链接失败），所以实现放在随头文件
> 分发的 `.hxx`。但如果每个使用 Pinocchio 的项目都从头实例化 `double` 版本，编译会极慢；
> 于是库在 `.cpp` 里**预先实例化**好默认标量版本编进 `libpinocchio.so`，你的项目直接链接即可。
> **结论**：读实现永远去 `.hxx`，不是 `.cpp`。

#### Demo A：用 30 行"迷你 Pinocchio"复现三层布局

下面是一个自包含的玩具工程，把 Pinocchio 的三层机制浓缩成最小可编译单元。**照着敲一遍，
就能彻底理解为什么模板库要这么拆。**

```cpp
// ---------- square.hpp（声明层，对应 algorithm/*.hpp）----------
#pragma once
template<typename Scalar> Scalar square(Scalar x);   // 只声明
#include "square.hxx"                                 // 末尾拉进实现（关键！）

// ---------- square.hxx（实现层，对应 src/algorithm/*.hxx）----------
#pragma once
template<typename Scalar> Scalar square(Scalar x) { return x * x; }  // 模板定义

// ---------- square.cpp（实例化层，对应 src/algorithm/*.cpp）----------
#include "square.hpp"
template double square<double>(double);   // 对 double 显式实例化，编进 .o/.so

// ---------- main.cpp（你的项目）----------
#include "square.hpp"
#include <iostream>
int main() {
  std::cout << square(3.0)   << "\n";   // double：直接用 square.cpp 里预实例化好的版本
  std::cout << square(3.0f)  << "\n";   // float：编译期临时实例化（因为 .hxx 可见）
}
```

**三个可自己动手验证的现象**：

1. 把 `square.hpp` 末尾的 `#include "square.hxx"` 删掉 → `square(3.0f)` **链接报错**
   `undefined reference to float square<float>(float)`。**这解释了为什么每个 `.hpp` 末尾都要
   include 对应 `.hxx`**：不这样，非默认标量就没有可见定义。
2. `square(3.0)`（double）即使删掉 `.hxx` 也能链接 —— 因为 `square.cpp` 已经预实例化了 double 版。
   **这正是 `libpinocchio.so` 的作用**：帮你把最常用的 double 版本编好，省掉重复编译。
3. 对照真实文件：`kinematics.hpp` 末尾的 `#include ".../kinematics.hxx"`、
   `kinematics.cpp` 里的 `template ... forwardKinematics<context::Scalar,...>` 就是上面②③的真身。

#### Demo B：在真实仓库里"三跳"定位任意函数的实现

```bash
# ① 声明 + 文档（看接口和公式）
grep -n "forwardKinematics" include/pinocchio/algorithm/kinematics.hpp

# ② 实现体（看真正的循环/递推）—— 注意路径里的 src/
grep -n "forwardKinematics\|struct .*Step" include/pinocchio/src/algorithm/kinematics.hxx

# ③ 实例化（确认默认标量已编进库，一般无需细看）
grep -n "forwardKinematics" src/algorithm/kinematics.cpp
```

> **记牢这条路径差异**：`include/pinocchio/algorithm/xxx.hpp`（声明）↔
> `include/pinocchio/src/algorithm/xxx.hxx`（实现）。很多人只在 `algorithm/` 里翻，找不到实现体，
> 就是因为漏了中间的 `src/`。

### 11.2 Tpl 模板 + context 默认标量

> 本节回答"实现原理"：`pinocchio::Model` 这个别名**到底是被哪几行代码、按什么顺序生成的**，
> 以及 `cast<>()` 和 autodiff 在底层是怎么落地的。读完你应能自己画出从宏到别名的完整链路。

#### 原理总览：一条从宏到别名的流水线

所有核心类型的"真名"都带 `Tpl` 后缀、按 `<Scalar, Options, JointCollection>` 模板化
（`ModelTpl` / `DataTpl` / `SE3Tpl` / `MotionTpl` …）。你平时写的 `pinocchio::Model`、`pinocchio::SE3`
只是**默认标量的 typedef 别名**。这个别名不是硬编码的 `double`，而是经过一条 **4 段预处理流水线**
在编译期"拼"出来的：

```
context.hpp                       ← 你 include 的入口
  └─(include) src/context.hxx     ← ① 定义默认宏 + 选择 context 文件
       #define PINOCCHIO_SCALAR_TYPE_DEFAULT double
       #define PINOCCHIO_CONTEXT_FILE_DEFAULT "pinocchio/src/context/default.hxx"
       #ifndef PINOCCHIO_CONTEXT_FILE                 // ← 允许外部覆盖！(关键)
         #define PINOCCHIO_CONTEXT_FILE PINOCCHIO_CONTEXT_FILE_DEFAULT
       #endif
       #include PINOCCHIO_CONTEXT_FILE
          │
          └─(include) src/context/default.hxx  ← ② 把"默认标量"接到"当前标量"
               #define PINOCCHIO_SCALAR_TYPE  PINOCCHIO_SCALAR_TYPE_DEFAULT   // = double
               #include "pinocchio/src/context/generic.hxx"
               #undef  PINOCCHIO_SCALAR_TYPE
                  │
                  └─(include) src/context/generic.hxx  ← ③ 真正产生 context 命名空间
                       namespace pinocchio::context {
                         typedef PINOCCHIO_SCALAR_TYPE Scalar;    // typedef double Scalar;
                         static constexpr int Options = 0;
                         typedef Eigen::Matrix<Scalar,Dynamic,1,Options> VectorXs;
                         typedef Eigen::Matrix<Scalar,Dynamic,Dynamic,Options> MatrixXs;
                         ... // 一批公共 Eigen 类型别名
                       }

src/multibody/fwd.hxx             ← ④ 用 context::Scalar 造出面向用户的别名
   typedef ModelTpl<context::Scalar, context::Options> Model;   // 只给 2 个参数！
   typedef DataTpl <context::Scalar, context::Options> Data;
   typedef SE3Tpl  <context::Scalar, context::Options> SE3;
```

拆开看每一段的作用：

- **① `src/context.hxx`**：定义"默认标量 = double""默认 Options = 0"，并**留了一个覆盖钩子**
  `#ifndef PINOCCHIO_CONTEXT_FILE`。如果编译时没人指定，就用 `default.hxx`。这个钩子是整个机制的灵魂
  （见下面"为什么要绕这么一圈"）。
- **② `default.hxx`**：一个薄适配层，把"默认标量宏"接到"当前标量宏"`PINOCCHIO_SCALAR_TYPE`，
  包完 `generic.hxx` 后立刻 `#undef`，保证不污染后续。
- **③ `generic.hxx`**：**真正干活的地方**。在 `namespace pinocchio::context` 里
  `typedef PINOCCHIO_SCALAR_TYPE Scalar;`，于是 `context::Scalar` 就变成了 `double`；并顺手定义
  一批以此标量为基础的公共 Eigen 类型（`VectorXs`/`MatrixXs`/…），后面 `.cpp` 实例化时都复用它们。
- **④ `src/multibody/fwd.hxx`**：`typedef ModelTpl<context::Scalar, context::Options> Model;`
  —— 注意**只传了 2 个模板参数**，第 3 个 `JointCollectionTpl` 用的是 `ModelTpl` 前向声明里的**默认实参**
  `= JointCollectionDefaultTpl`（见 `fwd.hxx:19-23`）。所以最终：

```cpp
pinocchio::Model  ==  ModelTpl<double, 0, JointCollectionDefaultTpl>
```

#### 为什么要绕这么一圈？——不是脱裤子放屁

如果只想要 `double`，大可直接 `typedef ModelTpl<double,0> Model;`。绕这一圈的**唯一目的**，是那个
覆盖钩子 `PINOCCHIO_CONTEXT_FILE`：它让人能**把整个库按另一种标量重新预编译成一个独立的 `.so`**。

看 `src/context/casadi.hxx`（构建 `libpinocchio_casadi` 时 CMake 用 `-DPINOCCHIO_CONTEXT_FILE` 指向它）：

```cpp
// casadi.hxx —— 一个"非 double"的 context
#define PINOCCHIO_TEMPLATE_INSTANTIATION_HEADER "pinocchio/autodiff/casadi.hpp"
#include <casadi/casadi.hpp>
#define PINOCCHIO_SCALAR_TYPE ::casadi::SX          // ← 当前标量换成 CasADi 符号类型
#include "pinocchio/src/context/generic.hxx"        // ← 复用同一个 generic.hxx！
#undef PINOCCHIO_SCALAR_TYPE
```

于是**同一套 `generic.hxx` + 同一套算法源码**，仅仅因为 `PINOCCHIO_SCALAR_TYPE` 变了：


| context 文件          | `context::Scalar`   | `pinocchio::Model` 变成    | 产出                  |
| --------------------- | ------------------- | -------------------------- | --------------------- |
| `default.hxx`（默认） | `double`            | `ModelTpl<double,...>`     | `libpinocchio`        |
| `casadi.hxx`          | `casadi::SX`        | `ModelTpl<casadi::SX,...>` | `libpinocchio_casadi` |
| `cppad.hxx`           | `CppAD::AD<double>` | `ModelTpl<AD<double>,...>` | `libpinocchio_cppad`  |
| `mpfr`                | 高精度浮点          | `ModelTpl<mpfr,...>`       | 高精度库              |

**这就是 context 机制的实现价值**：把"整个库默认针对哪种标量做预编译"抽成一个**可从外部一键切换的开关**，
而不必改任何算法代码。`src/**/*.cpp` 里清一色写 `context::Scalar`（而非 `double`），正是为了跟着这个开关走。

#### `.cpp` 显式实例化：也依赖 context

回顾 [§11.1](#111-三层文件布局-hpp--hxx--cpp)，`.cpp` 的作用是把默认标量版本预编译进库。它实例化的正是
`context::Scalar` 版本（`src/algorithm/kinematics.cpp`）：

```cpp
#include "pinocchio/src/context/template-instantiation.hxx"
template ... void forwardKinematics<context::Scalar, context::Options,
                                     JointCollectionDefaultTpl, ...>(const Model&, Data&, ...);
```

所以切换 context → `context::Scalar` 变 → `.cpp` 自动实例化的就是新标量版本 → 产出的库也就换了标量。
一条线全串起来了。

#### `cast<NewScalar>()` 的底层：逐字段"换标量重造"

`model.cast<float>()` 不是黑魔法，就是**新建一个 `ModelTpl<NewScalar>`，把每个字段按新标量搬过去**
（`src/multibody/model.hxx` 的 `cast()`）：

```cpp
template<typename NewScalar>
typename CastType<NewScalar, ModelTpl<...>>::type ModelTpl<...>::cast() const
{
  typedef ModelTpl<NewScalar, Options, JointCollectionTpl> ReturnType;
  ReturnType res;
  res.nq = nq;  res.nv = nv;  res.parents = parents;  res.names = names;   // 整数/拓扑：直接拷
  res.gravity = gravity.template cast<NewScalar>();                        // Eigen 量：逐个 .cast
  res.lowerPositionLimit = lowerPositionLimit.template cast<NewScalar>();
  for (size_t k = 0; k < joints.size(); ++k) {
    res.inertias[k]       = inertias[k].template cast<NewScalar>();        // 空间惯量：换标量
    res.jointPlacements[k]= jointPlacements[k].template cast<NewScalar>(); // SE3：换标量
    res.joints[k]         = joints[k].template cast<NewScalar>();          // 关节 variant：换标量
  }
  return res;                                                              // 拓扑不变，数值换类型
}
```

**关键观察**：拓扑结构（`parents`、`names`、`nq`…）与标量无关，直接整数拷贝；只有**数值容器**
（`SE3`、`Inertia`、Eigen 向量、关节参数）才 `.template cast<NewScalar>()`。因为整棵类型树都是模板，
`ModelTpl<NewScalar>` 的所有成员类型都能自动重新实例化——这就是"换个 `Scalar` 就能重造一个模型"的实现根据。

#### autodiff / 解析梯度：cast 到"会记账的标量"

把 `NewScalar` 换成 `CppAD::AD<double>` 或 `casadi::SX` 这类**会自动记录运算图的标量**后，
你照常调用 `rnea(model_ad, data_ad, q_ad, ...)`——算法一行没改，但每步 `+ - * /` 都被这些标量
记进了计算图，结束后即可对图求导。**这就是 Pinocchio "解析/自动梯度"的底层实现**：
不是为梯度单独写代码，而是**让同一套算法在"会求导的数系"上跑一遍**。所以 `autodiff/` 目录本质只是
"提供这些特殊标量 + 对应的实例化"，算法逻辑全部复用第 4 章那套。

> **读源码技巧**：看到冗长签名 `template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>`
> 时，心里把 `Scalar`→`double`、`ModelTpl<...>`→`Model` 一替换，signature 立刻清爽。反过来，
> 当你想理解 autodiff/codegen 时，再把 `Scalar` 想成 `AD<double>` / `casadi::SX` 即可，**代码是同一份**。

#### Demo C：验证"别名等价"——`Model` 就是 `ModelTpl<double,...>`

```cpp
#include "pinocchio/multibody/model.hpp"
#include <type_traits>

int main() {
  using namespace pinocchio;
  // 编译期断言：默认别名 == 指定 double 标量的模板实例
  static_assert(
    std::is_same<Model, ModelTpl<double, 0, JointCollectionDefaultTpl>>::value,
    "pinocchio::Model 就是 double 版的 ModelTpl");
  static_assert(std::is_same<SE3, SE3Tpl<double, 0>>::value, "");
  // 能编译通过，就证明了 §11.2 的等价关系
}
```

#### Demo D：同一份算法，换标量就能跑出不同能力（cast 的威力）

`model.cast<NewScalar>()` 把整棵类型树按新标量重新实例化。这就是 Pinocchio 能同时支持
**普通计算 / 高精度 / 自动微分 / 代码生成**的根源——**算法只写一遍**。

```cpp
#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/algorithm/rnea.hpp"

int main() {
  using namespace pinocchio;
  Model model;                                   // double 模型
  buildModels::manipulator(model);
  Data data(model);
  Eigen::VectorXd q = neutral(model), v = q, a = q;
  rnea(model, data, q, v, a);                    // ① 普通 double 计算

  // ② 换成 float：一行 cast，同一套 rnea 代码直接复用
  typedef ModelTpl<float> Modelf;
  Modelf model_f = model.cast<float>();
  DataTpl<float> data_f(model_f);
  rnea(model_f, data_f, q.cast<float>(), v.cast<float>(), a.cast<float>());

  // ③ 换成 CppAD::AD<double>：rnea 的输出就带上了对 (q,v,a) 的导数信息
  //    → 这正是"解析梯度/autodiff"的实现方式，见 autodiff/ 目录
  //    ModelTpl<CppAD::AD<double>> model_ad = model.cast<CppAD::AD<double>>();

  // ④ 换成 100 位高精度：见 §6.15 multiprecision.cpp
}
```

**读源码时的收获**：以后你在 `autodiff/`、`codegen/` 目录看到的"魔法"，本质都只是
"给 `ModelTpl` / `rnea` 换了个 `Scalar` 模板参数"，没有任何算法被重写。理解这一点，
这两个高级目录就不再神秘。

> **Python 侧提醒**：`pinocchio.Model` 只暴露 `double` 版；`float` / autodiff / codegen 这些
> 非默认标量的能力**只在 C++ 层可用**（见 [§13](#13-python-绑定如何映射到-c)）。

### 11.3 CRTP + Boost.Fusion 访问者

这是整个库最核心、也最难的机制。**问题**：`forwardKinematics` 要遍历运动树，对每个关节调用
"该关节类型特有的" `calc()`。但关节类型五花八门（旋转、移动、球形、浮动基…），运行时又存在一个
`std::vector<JointModel>` 里。怎样既保留每种关节的静态类型（零虚函数开销），又能统一遍历？

**答案**：`JointModel`/`JointData` 用 Boost.Variant 存放所有可能的关节类型，遍历时用
**Fusion 访问者（visitor）**在编译期为每种关节生成专门代码。一个算法 = 一个或多个"访问者结构体"。

以 FK 为例，`src/algorithm/kinematics.hxx` 里的实现骨架：

```cpp
namespace impl {
  // 一个"访问者"就是算法的一趟(pass)。继承 JointUnaryVisitorBase<自己>（CRTP）
  template<typename Scalar, int Options, ...>
  struct ForwardKinematicZeroStep
  : fusion::JointUnaryVisitorBase<ForwardKinematicZeroStep<...>>
  {
    // ArgsType：这趟遍历需要透传给每个关节的参数包
    typedef boost::fusion::vector<const Model&, Data&, const ConfigVectorType&> ArgsType;

    // algo()：对【每一种】具体 JointModel 类型都会实例化一份
    //   jmodel/jdata 已经是"具体关节类型"（如 JointModelRX），静态分派、无虚函数
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model, Data & data, const ConfigVectorType & q)
    {
      const JointIndex i = jmodel.id();
      const JointIndex parent = model.parents[i];

      jmodel.calc(jdata.derived(), q.derived());   // ← 关节自身的运动学(每种关节不同)
      data.liMi[i] = model.jointPlacements[i] * jdata.M();   // 关节相对父的变换
      data.oMi[i]  = (parent > 0) ? data.oMi[parent] * data.liMi[i]   // 递推到世界系
                                  : data.liMi[i];
    }
  };
}

// 对外的模板函数：只是"对每个关节跑一遍这个访问者"
template<...> void forwardKinematics(const Model & model, Data & data, const ConfigVector & q)
{
  for (JointIndex i = 1; i < model.njoints; ++i)
    ForwardKinematicZeroStep<...>::run(model.joints[i], data.joints[i],
                                       typename ...::ArgsType(model, data, q));
}
```

**读懂访问者的三把钥匙**：


| 概念                                              | 含义                                         | 在哪看                            |
| ------------------------------------------------- | -------------------------------------------- | --------------------------------- |
| `struct XxxStep : JointUnaryVisitorBase<XxxStep>` | 算法的一趟遍历（CRTP 自引用）                | `visitor.hpp`                     |
| `static void algo(jmodel, jdata, ...)`            | 对**每种关节类型**编译期特化的实体循环       | 各`*.hxx`                         |
| `jmodel.calc(jdata, q)` / `calc_aba(...)`         | 关节**自身**的运动学/动力学（S、变换、惯量） | `src/multibody/joint/joint-*.hxx` |

**因此，读任何一个算法实现的套路是固定的**：

1. 在 `xxx.hxx` 里找到 `namespace impl`，数一数有几个 `struct ...Step`（= 算法有几趟）；
2. 每个 `Step` 的 `algo()` 就是那一趟对单个关节做的事，对照第 4 章的公式看；
3. 遇到 `jmodel.calc()` / `jdata.S()` / `jdata.M()` 想深究，再去 `joint/joint-revolute.hxx` 等看具体关节。

RNEA（`rnea.hxx`）就有两个 Step（`Pass1` 正向、`Pass2` 反向），和 [4.3](#43-逆向动力学-rnea) 的两趟递推
严格对应；ABA（`aba.hxx`）有三个 Step，对应 [4.4](#44-正向动力学-aba) 的三趟。**公式 ↔ Step 一一对应**，
这是本指南第 4 章能直接当"实现导读"用的原因。

#### 分派链路：`run` 一次调用到底发生了什么

`ForwardKinematicZeroStep::run(model.joints[i], ...)` 里 `model.joints[i]` 是一个
`JointModelVariant`（Boost.Variant，运行时可能是 RX/RY/FreeFlyer…任意一种）。`run` 如何在
**不用虚函数**的前提下跳到正确的 `algo()`？看 [joint-unary-visitor.hxx:49](include/pinocchio/src/multibody/visitor/joint-unary-visitor.hxx#L49)：

```cpp
static ReturnType run(const JointModelTpl<...> & jmodel,      // variant
                      JointDataTpl<...> & jdata, ArgsTmp args)
{
  InternalVisitorModelAndData<...> visitor(jdata, args);
  return boost::apply_visitor(visitor, jmodel);   // ← 关键：变体分派
}
```

```
run(jmodel, jdata, args)
   └─ boost::apply_visitor(visitor, jmodel)
        └─ Boost 根据 jmodel 当前"真实类型"（如 JointModelRX），
           在【编译期已生成的】分支表里选中对应实体
              └─ 调你写的 algo<JointModelRX>(jmodel, jdata, ...)
                    └─ jmodel.calc(jdata, q)  → 进入 joint-revolute.hxx 的 calc()
```

**要点**：`apply_visitor` 对 variant 里**每一种**可能的关节类型都在编译期生成了一个 `algo` 实例，
运行时只是"选分支"，没有虚表查找。这就是"零成本抽象"——写起来像多态，跑起来像手写的 switch。

#### `jmodel.calc()` 里究竟算了什么（以旋转关节为例）

追到最底层，看 [joint-revolute.hxx](include/pinocchio/src/multibody/joint/joint-revolute.hxx) 的 `calc()`：

```cpp
template<typename ConfigVector>
void calc(JointDataDerived & data, const Eigen::MatrixBase<ConfigVector> & qs) const
{
  data.joint_q[0] = qs[idx_q()];         // 从全局 q 里取出【这个关节】的分量
  Scalar ca, sa;
  SINCOS(data.joint_q[0], &sa, &ca);     // 一次算出 sin/cos
  data.M.setValues(sa, ca);              // 填充该关节的相对变换 jdata.M()（绕轴转 q 的 SE3）
}
```

对照 [§4.1](#41-正向运动学-fk) 的 ${}^{p(i)}T_i(q_i)$：**每种关节的差异，全部封装在各自的 `calc()` 里**——
旋转关节用 $\cos/\sin$ 填旋转块，移动关节填平移块，浮动基从 q 里读 7 维位姿……而上层的
FK / RNEA 递推公式对所有关节**长得完全一样**。这就是访问者模式解耦"通用递推"与"关节特有运动学"的价值。

> `idx_q()` / `idx_v()`：每个关节在全局 `q`（维度 nq）和 `v`（维度 nv）里占的**起始下标**。
> 这是 Pinocchio 用一维大向量存所有关节配置、又能让每个关节只取自己那几维的机制。

#### Demo E：亲手写一个访问者，遍历模型打印每个关节信息

这是最能建立直觉的练习——**复刻一个 Pinocchio 内部算法的骨架**。下面代码改编自
[unittest/visitor.cpp](unittest/visitor.cpp)（库自带的可运行范例），可直接编译：

```cpp
#include "pinocchio/multibody/sample-models.hpp"
#include "pinocchio/multibody/visitor.hpp"
#include <iostream>

namespace bf = boost::fusion;

// 1) 定义一个访问者：继承 JointUnaryVisitorBase<自己>（CRTP）
struct PrintJointVisitor
: public pinocchio::fusion::JointUnaryVisitorBase<PrintJointVisitor>
{
  // 2) ArgsType：声明这趟遍历要透传给每个关节的额外参数
  typedef bf::vector<const pinocchio::Model &> ArgsType;

  // 3) algo()：会对【每种】具体关节类型各实例化一份
  template<typename JointModel>
  static void algo(const pinocchio::JointModelBase<JointModel> & jmodel,
                   const pinocchio::Model & model)
  {
    std::cout << "id=" << jmodel.id()
              << "  name="   << model.names[jmodel.id()]
              << "  type="   << jmodel.shortname()   // 如 "JointModelRX"
              << "  nq="     << jmodel.nq()           // 该关节配置维度
              << "  nv="     << jmodel.nv()           // 该关节速度维度
              << "  idx_q="  << jmodel.idx_q()        // 在全局 q 中的起始下标
              << std::endl;
  }
};

int main()
{
  pinocchio::Model model;
  pinocchio::buildModels::humanoidRandom(model);   // 带浮动基的人形

  // 4) 像 Pinocchio 内部算法一样，对每个关节 run 一遍访问者
  for (pinocchio::JointIndex i = 1; i < (pinocchio::JointIndex)model.njoints; ++i)
    PrintJointVisitor::run(model.joints[i], PrintJointVisitor::ArgsType(model));

  return 0;
}
```

**运行后你会看到**：第一个关节 `type=JointModelFreeFlyer, nq=7, nv=6, idx_q=0`（浮动基），
其余是 `JointModelRX/RY/RZ, nq=1, nv=1`。这直观印证了 [§2.1](#21-李群-se3) 的 $n_q = n_{joints}+7,\;n_v = n_{joints}+6$。

**你已经掌握了写 Pinocchio 算法的全部套路**：把上面 `algo()` 里的打印换成"计算 `oMi`"，
就是 `forwardKinematics`；换成"计算受力并反向传递"，就是 RNEA 的一个 Step。**所有内置算法都是这个模式。**

#### 为什么用 CRTP + Variant，而不是普通虚函数？


| 方案                        | 代价                                                                             |
| --------------------------- | -------------------------------------------------------------------------------- |
| 虚函数多态                  | 每个关节每次调用都有**虚表查找**，且无法内联；对每秒百万次的动力学循环是灾难     |
| CRTP + Variant（Pinocchio） | 编译期为每种关节生成专门代码，**可内联、零运行时开销**，代价是编译慢、报错信息长 |

这就是为什么 Pinocchio 能在 1ms 内算完百自由度动力学——抽象的代价全部转移到了**编译期**。
理解这个取舍，也就理解了"为什么读它的源码时模板和报错这么劝退"：这是性能换来的必然复杂度。

### 11.4 第三个模板参数 `JointCollectionTpl`：关节类型"目录"

回到那个让人困惑的模板签名：

```cpp
template<typename Scalar, int Options, template<typename, int> class JointCollectionTpl>
struct ModelTpl;
```

- 第 1 个 `Scalar`、第 2 个 `Options` 已在 [§11.2](#112-tpl-模板--context-默认标量) 讲透；
- **第 3 个 `JointCollectionTpl` 就是本节主角**——它回答了 [§11.3](#113-crtp--boostfusion-访问者)
  遗留的一个问题：**`model.joints[i]` 那个 variant，到底列了哪些关节类型？谁规定的？**

答案：由 `JointCollectionTpl` 规定。它是一份**"关节类型目录（menu）"**——一个 struct，给定
`<Scalar, Options>` 后，把库支持的所有关节类型（Model 版和 Data 版）列成 typedef，
并打包成 `apply_visitor` 要分派的那个 `boost::variant`。

> 注意签名里它是 `template<typename,int> class JointCollectionTpl`，即一个**模板的模板参数**
> （template-template parameter）：传进来的不是一个具体类型，而是一个"还差 `<Scalar,Options>`
> 才能实例化的模板"。这样 `ModelTpl` 内部可以用当前的 `Scalar` 去实例化这份目录。

#### `JointCollectionDefaultTpl`：库自带的默认目录

`JointCollectionDefaultTpl`（[joint-collection.hxx:19](include/pinocchio/src/multibody/joint/joint-collection.hxx#L19)）
就是**默认那份目录**，把所有内置关节类型枚举出来，最后压进一个 `variant`：

```cpp
template<typename _Scalar, int _Options>
struct JointCollectionDefaultTpl
{
  typedef _Scalar Scalar;
  static constexpr int Options = _Options;

  // ① 给每种内置关节起一个"用当前 Scalar 实例化好"的别名
  typedef JointModelRevoluteTpl<Scalar, Options, 0> JointModelRX;   // 绕 X 轴旋转
  typedef JointModelRevoluteTpl<Scalar, Options, 1> JointModelRY;
  typedef JointModelRevoluteTpl<Scalar, Options, 2> JointModelRZ;
  typedef JointModelFreeFlyerTpl<Scalar, Options>   JointModelFreeFlyer;  // 浮动基
  typedef JointModelSphericalTpl<Scalar, Options>   JointModelSpherical;
  typedef JointModelPrismaticTpl<Scalar, Options, 0> JointModelPX;   // 沿 X 轴平移
  // …… 还有 Planar / Translation / Helical / Universal / Composite / Mimic 等

  // ② 把它们全部装进一个 variant —— 这就是 §11.3 里 apply_visitor 分派的那个盒子！
  typedef boost::variant<
    JointModelRX, JointModelRY, JointModelRZ,
    JointModelFreeFlyer, JointModelPlanar, JointModelSpherical, /* … */
    boost::recursive_wrapper<JointModelComposite>,   // 复合关节：可嵌套，故用 recursive_wrapper
    boost::recursive_wrapper<JointModelMimic>
  > JointModelVariant;

  // ③ Data 侧同理，另有一份 JointDataVariant（每种关节的运行时数据）
  typedef boost::variant< JointDataRX, /* … */ > JointDataVariant;
};
```

而通用关节包装器 `JointModelTpl` 正是从这份目录里**取出** variant 来用
（[joint-generic.hxx:97-98](include/pinocchio/src/multibody/joint/joint-generic.hxx#L97)）：

```cpp
typedef JointCollectionTpl<Scalar, Options>       JointCollection;
typedef typename JointCollection::JointModelVariant JointModelVariant;   // ← 就是这里接上的
```

**一句话串起三节**：`JointCollectionDefaultTpl` 列出"库认识哪些关节"→ 打包成 `JointModelVariant`
→ 存进 `model.joints[i]` → [§11.3](#113-crtp--boostfusion-访问者) 的 `apply_visitor` 在这个 variant 上分派
→ 调到你 `algo()` 里对应关节类型的那份实例。**它就是整个访问者机制的"类型清单"源头。**

#### 为什么把它做成模板参数（而不是写死）？

因为这是一个**定制点（customization point）**。绝大多数用户用默认目录即可，但把它开放成模板参数后，
高级场景可以替换：

| 场景 | 做法 | 收益 |
|------|------|------|
| 只用旋转关节的机械臂 | 自定义一份只含 `JointModelRX/RY/RZ` 的精简目录 | variant 更小 → **编译更快、二进制更小、分派分支更少** |
| 需要一种库里没有的自定义关节 | 定义新关节类型并加进目录的 variant | 无需改动库源码即可扩展 |
| 不同标量的预编译库 | 目录本身是模板，`cast`/换 context 时自动跟着换 `Scalar` | 与 [§11.2](#112-tpl-模板--context-默认标量) 机制无缝配合 |

因为它是**带默认值的第 3 个参数**（`= JointCollectionDefaultTpl`，见
[fwd.hxx:19-23](include/pinocchio/src/multibody/fwd.hxx#L19)），所以：

```cpp
typedef ModelTpl<context::Scalar, context::Options> Model;   // 只写 2 个参数
//                                                            // 第 3 个自动 = JointCollectionDefaultTpl
// 于是 pinocchio::Model == ModelTpl<double, 0, JointCollectionDefaultTpl>
```

99% 的代码（包括你自己）都用这个默认值，所以平时**根本感觉不到它的存在**——但它默默决定了
"这个库到底认识哪些关节"。这也是你在 [§11.2](#112-tpl-模板--context-默认标量) 看到别名只传 2 个参数、
却等价于 3 个参数实例的原因。

> **读源码小结**：看到 `JointCollectionTpl` / `JointCollectionDefaultTpl` 时，脑子里翻译成
> "**关节类型目录**"即可。想知道"这个库支持哪些关节"，直接翻
> [joint-collection.hxx](include/pinocchio/src/multibody/joint/joint-collection.hxx) 里那个 `boost::variant<...>` 列表。

---

## 12. 实战：追踪一个算法从 API 到实现

以最常用的 `forwardKinematics` 为例，走通"Python 调用 → C++ 实现 → 关节细节"的完整链路。
掌握这一条，其余算法照抄即可。

**① 入口（Python）**——`examples/overview-urdf.py`：

```python
pinocchio.forwardKinematics(model, data, q)
```

**② Python→C++ 绑定**——`bindings/python/algorithm/expose-kinematics.cpp`
里用 `bp::def("forwardKinematics", ...)` 把 C++ 模板函数的 `double` 实例暴露给 Python。
（绑定机制详见第 13 节。）

**③ 声明 + 文档**——[include/pinocchio/algorithm/kinematics.hpp](include/pinocchio/algorithm/kinematics.hpp)：
读函数注释确认语义（"更新 `data.oMi[]`，$O(n)$"），拿到签名。文件末尾 `#include ".../kinematics.hxx"`。

**④ 实现体**——[include/pinocchio/src/algorithm/kinematics.hxx](include/pinocchio/src/algorithm/kinematics.hxx)：

- 对外函数体：一个 `for` 循环，对每个关节 `run` 访问者 `ForwardKinematicZeroStep`；
- 访问者 `algo()`：`jmodel.calc()` + `oMi[i] = oMi[parent] * liMi[i]`（就是 [4.1](#41-正向运动学-fk) 的递推式
  ${}^0T_i = {}^0T_{p(i)} \cdot {}^{p(i)}T_i$）。

**⑤ 关节细节**——想知道旋转关节的 `calc()` 怎么把 $q_i$ 变成变换矩阵：
[include/pinocchio/src/multibody/joint/joint-revolute.hxx](include/pinocchio/src/multibody/joint/joint-revolute.hxx)
的 `JointDataRevolute::calc()`——用 $\cos q_i,\sin q_i$ 填充 `M`（该关节的 `jdata.M()`）。

**⑥ 数据结构**——`data.oMi` 的定义在
[include/pinocchio/src/multibody/data.hxx](include/pinocchio/src/multibody/data.hxx)（`DataTpl` 成员），
`SE3` 的 `operator*`/`actInv` 在 [include/pinocchio/spatial/se3.hpp](include/pinocchio/spatial/se3.hpp)。

**通用追踪清单**（换任何算法都适用）：

```
公开函数名  →  algorithm/<name>.hpp        （看签名 + 注释 + 公式）
            →  src/algorithm/<name>.hxx    （看 impl::*Step::algo，对照第4章公式）
            →  src/multibody/joint/*.hxx   （深究单关节 calc/calc_aba）
            →  spatial/*.hpp               （深究 SE3/Motion/Force/Inertia 运算）
旁证用法    →  unittest/<name>.cpp、examples/
```

---

## 13. Python 绑定如何映射到 C++

Pinocchio 的 Python 层是 **Boost.Python** 手写绑定（不是 pybind11、不是自动生成），
所以每个 Python 符号都能在 C++ 里找到对应物。理解映射规则后，可以从 Python API 反查 C++ 实现。

**目录对应**：

```
bindings/python/<模块>/expose-<主题>.cpp      # 定义 expose 函数，调用 bp::def / bp::class_
include/pinocchio/bindings/python/<模块>/     # 绑定用的 C++ 包装头（Visitor 模式）
bindings/python/module.cpp                    # 总入口，依次调用所有 exposeXxx()
```

**映射规律**：


| Python                       | C++ 对应                                                       | 说明                                            |
| ---------------------------- | -------------------------------------------------------------- | ----------------------------------------------- |
| `pin.forwardKinematics(...)` | `bp::def("forwardKinematics", ...)` in `expose-kinematics.cpp` | 函数：暴露的是`double` 实例                     |
| `pin.SE3`                    | `bp::class_<SE3>(...)` in `expose-SE3.cpp`                     | 类：`ModelTpl<double>` 等的别名被 `class_` 包装 |
| `model.createData()`         | `ModelTpl::createData()` 成员                                  | 成员方法 1:1 暴露                               |
| `data.oMi`                   | `DataTpl::oMi` 成员                                            | 通过`.def_readwrite`/`add_property` 暴露        |
| `pin.ReferenceFrame.LOCAL`   | `enum ReferenceFrame`                                          | 枚举用`bp::enum_` 暴露                          |

> 找一个 Python 函数的 C++ 源头：`grep -rn "\"forwardKinematics\"" bindings/python/`，
> 命中的 `expose-*.cpp` 就是绑定点，顺着它调用的模板函数名即可跳进 `algorithm/*.hpp`。

**要点**：Python 侧看到的一切都是 C++ **默认标量（double）** 的实例；autodiff / codegen 等
非 `double` 标量的能力只在 C++ 层可用，不通过 Python 暴露。

---

## 14. 源码阅读推荐顺序

按"地基→建模→算法→扩展"的依赖顺序读，避免一上来陷进模板细节。

```
第 0 步：先读机制，别读代码
  → 本指南 §11（三大机制）+ §12（一条追踪链路）
  → 目的：建立"公式↔Step、.hpp↔.hxx、Tpl↔别名"的心智模型

第 1 步：空间代数地基（spatial/）
  → se3.hpp → motion.hpp → force.hpp → inertia.hpp → explog.hpp
  → 对照本指南 §2；这是后面一切算法的"数据类型词汇表"

第 2 步：建模层（multibody/）
  → model.hpp(ModelTpl) + src/multibody/data.hxx(DataTpl)：先看有哪些字段（对照 §3）
  → joint/joint-revolute.hxx：读懂一种最简单关节的 calc / calc_aba / S
  → visitor.hpp：读懂 JointUnaryVisitorBase::run 如何分派（§11.3）

第 3 步：运动学算法（最短、最适合练手）
  → src/algorithm/kinematics.hxx：一趟遍历，验证 §12 的追踪链
  → jacobian.hxx：几何 Jacobian 如何按列填充

第 4 步：动力学三大件（对照 §4 的公式读 Step）
  → rnea.hxx（2 趟）→ crba.hxx（复合惯量、上三角）→ aba.hxx（3 趟 + 关节化体惯量）

第 5 步：人形核心
  → centroidal.hxx（质心动量矩阵 Ag）
  → constrained-dynamics.hxx + contact-cholesky.hxx（KKT / Delassus，最硬）

第 6 步：解析梯度（读完对应正算法再读其导数）
  → rnea-derivatives.hxx → aba-derivatives.hxx → kinematics-derivatives.hxx

第 7 步：按需扩展
  → parsers/urdf（建模来源）、collision/（碰撞）、autodiff/ 与 codegen/（高级标量）
  → 想加新关节类型：doc/e-development/Implement_a_new_joint.md 是官方教程
```

**贯穿始终的两个好习惯**：

1. **对照单测读实现**——`unittest/<algo>.cpp` 给出该算法的期望输入输出，是最好的"实现真值"。
2. **对照本指南第 4 章读 `.hxx`**——每个 `impl::*Step` 都能在第 4 章找到对应公式，公式讲"为什么"，
   代码讲"怎么写"，两边对读效率最高。
