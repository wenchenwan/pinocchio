# Pinocchio 源码架构三大机制 —— 动手实验

配套 [PINOCCHIO_GUIDE.md 第 11 章](../../PINOCCHIO_GUIDE.md#11-源码架构三大机制)。

第 11 章讲的三个机制（三层文件布局、Tpl 模板、CRTP 访问者）纯读文字很抽象。
这个目录把它们做成**可编译、可运行、可破坏**的实验 —— 跑一遍比读十遍管用。

---

## 快速开始

```bash
cd doc/architecture-tutorial
make          # 编译三个实验
make run      # 依次运行
```

单独运行某个实验：`make run-01` / `make run-02` / `make run-03`

**依赖**：实验 1 只需 g++（自包含）；实验 2、3 需要已安装 Pinocchio C++ 开发包
（`pkg-config --cflags --libs pinocchio` 能输出内容即可）。

---

## 实验 1：三层文件布局 `.hpp / .hxx / .cpp`

**要回答的问题**：为什么一个算法要拆成三个文件？为什么每个 `.hpp` 末尾都要
`#include` 对应的 `.hxx`？

```bash
make run-01     # 正常情况
make break-01   # ★ 反面教材：故意删掉 include，看它怎么坏
```

四个文件与真实仓库的对应关系：

| 本实验 | 真实仓库 | 职责 |
|--------|----------|------|
| `square.hpp` | `include/pinocchio/algorithm/kinematics.hpp` | 声明 + 文档；**末尾 include .hxx** |
| `square.hxx` | `include/pinocchio/src/algorithm/kinematics.hxx` | 模板实现体（**算法真正写在这**） |
| `square.cpp` | `src/algorithm/kinematics.cpp` | 对 `double` 显式实例化，编进 `.so` |
| `main.cpp` | 你自己的项目 | 使用方 |

> ⚠️ 注意实现层路径里的 **`src/`**：`include/pinocchio/`**`src/`**`algorithm/xxx.hxx`。
> 很多人只在 `algorithm/` 目录里翻找不到实现体，就是漏了这一段。

### `make break-01` 会看到什么

注释掉 `square.hpp` 末尾的 `#include "square.hxx"` 后：

```
链接失败: `float square<float>(float)'
链接失败: `long double square<long double>(long double)'
```

而 `double` 版本**依然正常**。

**这两行输出解释了整个设计**：

- `double` 能用 → 因为 `square.cpp` 已显式实例化并编进目标文件
  （对应真实世界：`libpinocchio.so` 里预编译好了 `double` 版本，省掉你的编译时间）
- `float` / `long double` 失败 → 因为没有 `.hxx` 就看不到模板定义体，编译器无法为
  新标量生成代码
  （对应真实世界：没有这行 include，你就用不了 autodiff 标量、高精度标量）

---

## 实验 2：Tpl 模板 + context 默认标量

**要回答的问题**：`pinocchio::Model` 到底是什么？`cast<>()` 在底层做了什么？
为什么这个机制是"解析梯度"的前提？

```bash
make run-02
```

实验内容：

1. **编译期证明别名等价** —— 用 `static_assert` 断言
   `Model == ModelTpl<context::Scalar, context::Options>` 且 `context::Scalar == double`。
   能编译通过本身就是证明。
2. **同一份 `rnea()` 跑两种标量** —— `double` 版和 `long double` 版，
   **算法源码一行都没改**，只是模板参数不同。
3. **精度对比** —— 输出两者差异约 `4.4e-15`，与 `double` 的机器精度
   `2.2e-16` 同量级，说明 `double` 结果本身可靠。

**为什么这件事重要**：正因为算法是模板，把标量换成"会记账的类型"
（CppAD / CasADi 的 AD 标量）就能自动得到导数 —— 这是
[第 4.8 节解析梯度](../../PINOCCHIO_GUIDE.md#48-解析梯度)能存在的底层前提。

---

## 实验 3：CRTP + Boost.Fusion 访问者 ★ 最核心

**要回答的问题**：关节类型五花八门（旋转/移动/球形/浮动基），却存在同一个数组里，
Pinocchio 如何在**不用虚函数**的前提下对每种关节调用其特有的 `calc()`？

```bash
make run-03
```

三个递进的例子：

### ① 只读访问者：遍历打印关节信息

输出片段：

```
id= 1  root_joint      JointModelFreeFlyer   nq=7  nv=6  idx_q=0   idx_v=0
id=2   lleg1_joint     JointModelRX          nq=1  nv=1  idx_q=7   idx_v=6
```

第一行直观印证了浮动基的 `nq=7 ≠ nv=6`
（见 [§2.1](../../PINOCCHIO_GUIDE.md#21-李群-se3)）。

### ② 手写 `forwardKinematics`

复刻 `kinematics.hxx` 里 `ForwardKinematicZeroStep` 的完整骨架，核心就三行：

```cpp
jmodel.calc(jdata.derived(), q);                        // 关节自身运动学
data.liMi[i] = model.jointPlacements[i] * jdata.M();    // 相对父的变换
data.oMi[i]  = data.oMi[parent] * data.liMi[i];         // 递推到世界系
```

### ③ 与官方实现对比验证

用 SE(3) 对数映射度量差异，实测输出：

```
遍历全部 27 个关节，max‖log6(ᵒM_mine⁻¹ · ᵒM_ref)‖ = 8.354e-16
→ 完全一致 ✓
```

**做完这个实验，你就掌握了写 Pinocchio 算法的全部套路**：把 `algo()` 里的内容
换成"计算受力并反向传递"，就是 RNEA 的一个 Step。

### 读源码的固定套路

| 算法 | Step 数 | 对应章节 |
|------|---------|----------|
| `forwardKinematics` | 1 | [§4.1](../../PINOCCHIO_GUIDE.md#41-正向运动学-fk) |
| `rnea` | 2（正向 + 反向） | [§4.3](../../PINOCCHIO_GUIDE.md#43-逆向动力学-rnea) |
| `aba` | 3 | [§4.4](../../PINOCCHIO_GUIDE.md#44-正向动力学-aba) |

打开任意 `.hxx`，数一数里面有几个 `struct ...Step`，就知道算法有几趟递推
—— **公式与 Step 一一对应**。

---

## 三个机制之间的关系

```
        ┌─────────────────────────────────────────┐
        │  想要：一份算法源码，多种标量 + 多种关节  │
        └─────────────────────────────────────────┘
                          │
        ┌─────────────────┼─────────────────┐
        ▼                 ▼                 ▼
   【机制二】         【机制三】         【机制一】
   Tpl 模板          CRTP 访问者        三层文件布局
   换标量            换关节类型          解决"模板定义必须
   double/AD/…       RX/RY/FreeFlyer…    可见"与"编译要快"
        │                 │              的矛盾
        └────────┬────────┘                 │
                 ▼                          │
          零运行时开销 ◄────────────────────┘
       （抽象代价全在编译期）
```

一句话串起来：**为了让一份算法同时适配多种标量（机制二）和多种关节（机制三），
必须把实现放在头文件里；而为了不让每个下游项目都重编一遍，就有了三层布局（机制一）。**

这也解释了读 Pinocchio 源码时"模板报错又长又劝退"的根源 —— 这是把运行时开销
转移到编译期所付出的必然代价。
