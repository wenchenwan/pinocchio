# 逆动力学 RNEA 解析：`algorithm/rnea.hxx`

> **RNEA**（Recursive Newton-Euler Algorithm，递归牛顿-欧拉算法）解决**逆动力学**问题：
> 给定运动 $(q,\dot q,\ddot q)$，求需要的关节力矩 $\tau$。$O(n)$ 两趟递推。
> 本文把 [src/algorithm/rnea.hxx](../include/pinocchio/src/algorithm/rnea.hxx) 讲透：两趟递推的物理意义、重力技巧、
> 力的对偶传播、以及 `nonLinearEffects`/`computeGeneralizedGravity` 如何是同一算法的特例。
> 配套：[正向运动学解析.md](正向运动学解析.md)（速度/加速度递推）、[空间代数运算解析.md](空间代数运算解析.md)（Inertia·Motion、力的对偶变换 `act`）、[算法总览.md](算法总览.md)。

---

## 目录
- [1. 逆动力学问题](#1-逆动力学问题)
- [2. 算法骨架：牛顿-欧拉两趟递推](#2-算法骨架牛顿-欧拉两趟递推)
- [3. 前向趟 `RneaForwardStep`：传播运动、算每体的力](#3-前向趟-rneaforwardstep传播运动算每体的力)
- [4. 重力技巧：把 $g$ 塞进加速度初值](#4-重力技巧把-g-塞进加速度初值)
- [5. 后向趟 `RneaBackwardStep`：对偶传播力、投影出力矩](#5-后向趟-rneabackwardstep对偶传播力投影出力矩)
- [6. 递推如何依赖运动学树结构](#6-递推如何依赖运动学树结构)
- [7. mimic 关节的处理](#7-mimic-关节的处理)
- [8. 外力 `fext` 与电机转子惯量 `armature`](#8-外力-fext-与电机转子惯量-armature)
- [9. 三个派生量：`nle` / `g` / `staticTorque`](#9-三个派生量nle--g--statictorque)
- [10. 复杂度与和 CRBA/ABA 的关系](#10-复杂度与和-crbaaba-的关系)
- [11. 一句话总结](#11-一句话总结)

---

## 1. 逆动力学问题

**已知**关节的位置 $q$、速度 $\dot q$、期望加速度 $\ddot q$，**求**驱动它所需的关节力矩：

$$
\boxed{\,\tau=\text{RNEA}(q,\dot q,\ddot q)=M(q)\ddot q+C(q,\dot q)\dot q+g(q)\,}
$$

- 但 RNEA **不显式组装** $M$、$C$、$g$——它直接沿运动树递推出 $\tau$，$O(n)$，比"先算 $M,C,g$ 再乘"快得多。
- 物理直觉：**牛顿-欧拉**方法。先自根向叶把每个连杆的**运动**（速度、加速度）算出来，据 $f=Ia$ 得每个连杆需要的**净力**；再自叶向根把这些力**沿树累加**传回，每个关节轴上力的投影就是该关节要出的力矩。

**用途**：力矩控制（前馈力矩）、轨迹可行性检查、作为正动力学/CRBA 的子模块、算非线性项 $b$ 和重力 $g$。

---

## 2. 算法骨架：牛顿-欧拉两趟递推

```cpp
const typename Data::TangentVectorType &
rnea(model, data, q, v, a) {
  data.tau.setZero();
  data.v[0].setZero();
  data.a_gf[0] = -model.gravity;          // ★ 重力技巧，见 §4

  // Pass 1: 根 → 叶，传播运动，算每体的力 f[i]
  for (i = 1; i < njoints; ++i)
    RneaForwardStep::run(...);

  // Pass 2: 叶 → 根，对偶传播力，投影出 tau
  for (i = njoints-1; i > 0; --i)
    RneaBackwardStep::run(...);

  data.tau.array() += model.armature.array() * a.array();  // 转子惯量，见 §6
  return data.tau;
}
```

- **Pass 1 正序** `i=1..n`：因 `parents[i]<i`，父先于子算好 → 运动量自根向叶流。
- **Pass 2 逆序** `i=n..1`：子先于父算好 → 力自叶向根汇聚。
- 两趟都用**编译期访问者** `JointUnaryVisitor` 按关节类型零开销分派（见 [正向运动学解析 §4](正向运动学解析.md)）。

---

## 3. 前向趟 `RneaForwardStep`：传播运动、算每体的力

```cpp
jmodel.calc(jdata, q, v);                       // 解关节：jdata.M(), S(), v(), c()
data.liMi[i] = model.jointPlacements[i] * jdata.M();

// ① 空间速度递推（同正运动学一阶）
data.v[i] = jdata.v();
if (parent > 0) data.v[i] += data.liMi[i].actInv(data.v[parent]);

// ② 空间加速度递推（a_gf = "含重力场" 的加速度）
data.a_gf[i]  = jdata.c() + (data.v[i] ^ jdata.v());          // 偏置 + 速度耦合
data.a_gf[i] += jdata.S() * jmodel.JointMappedVelocitySelector(a);  // + S q̈
data.a_gf[i] += data.liMi[i].actInv(data.a_gf[parent]);      // + 父加速度搬进本系

// ③ 牛顿-欧拉：本体净力 f = I a + v ×* (I v)
model.inertias[i].__mult__(data.v[i],   data.h[i]);          // h = I v（本体动量）
model.inertias[i].__mult__(data.a_gf[i], data.f[i]);         // f = I a
data.f[i] += data.v[i].cross(data.h[i]);                     // + v ×* h（陀螺项）
```

### 3.1 ①② 就是运动学递推

$v_i$、$a_i$ 的递推公式与 [正向运动学解析 §5-6](正向运动学解析.md) 完全一致：$v_i={}^iX_\lambda v_\lambda+S_i\dot q_i$，$a_i={}^iX_\lambda a_\lambda+S_i\ddot q_i+\dot S_i\dot q_i+v_i\times v_{J,i}$。唯一区别是这里的加速度存进 `a_gf`（gravity field），初值含 $-g$（§4）。

### 3.2 ③ 牛顿-欧拉方程（每个刚体的力平衡）

单个刚体的空间力平衡（6D 牛顿-欧拉方程）：

$$
\boxed{\,f_i=I_i\,a_i+v_i\times^{*}(I_i v_i)\,}
$$

- $I_i a_i$：**牛顿+欧拉主项**，惯量乘空间加速度（`inertias[i] * a_gf[i]`）。因 $a_{gf}$ 含 $-g$，这一项已经把**重力惯性力**算进去了。
- $v_i\times^{*}(I_iv_i)$：**陀螺力/科氏项**。$h_i=I_iv_i$ 是本体空间动量，$v\times^{*}$ 是力的叉乘（$se(3)^*$ 上的 $\mathrm{ad}^*$，见 [空间代数 §6](空间代数运算解析.md)）。代码写成 `v.cross(h)`。
- $f_i$ 是"要让连杆 $i$ 产生加速度 $a_i$，作用在它上面的合空间力（wrench）"。此刻还**不含**子连杆通过关节传来的反作用力——那要靠 Pass 2 累加。

> `__mult__(x, y)` 是"就地惯量作用" $y\leftarrow I\cdot x$，避免临时对象。`data.h` 顺带存下来供 ABA 导数等复用。

---

## 4. 重力技巧：把 $g$ 塞进加速度初值

RNEA 处理重力用了一个经典技巧（Featherstone）：

```cpp
data.a_gf[0] = -model.gravity;   // 根（universe）的加速度初始化为 −g，而非 0
```

物理原理：**重力等效于整个系统以 $-g$ 做加速运动**（等效原理）。让虚拟的世界根带一个 $-g$ 的"假加速度"，它会通过 `actInv(a_gf[parent])` 沿树传播到每个连杆，于是每个 $a_{gf,i}$ 都自动含 $-g$ 分量。再乘惯量 $I_i a_{gf,i}$，就凭空多出 $-I_i g$——正是**重力惯性力**。

好处：**不需要单独遍历加重力项**，重力被"免费"融进加速度递推。`a_gf` 的命名（a with gravity field）就是这个意思。它不是真实加速度：真实加速度 $a_i=a_{gf,i}+g$（ABA 里 `data.oa[i]=oa_gf[i]+gravity` 就是还原）。

---

## 5. 后向趟 `RneaBackwardStep`：对偶传播力、投影出力矩

```cpp
// ① 关节力矩 = 运动子空间对本体力的投影
jmodel.JointMappedVelocitySelector(data.tau) += jdata.S().transpose() * data.f[i];

// ② 把本体力沿关节对偶传给父连杆（累加）
if (parent > 0)
  data.f[parent] += data.liMi[i].act(data.f[i]);
```

### 5.1 ① 力矩 = $S^\top f$

关节 $i$ 只能沿它的**运动子空间** $S_i$ 出力。作用在连杆 $i$ 上的总空间力 $f_i$ 在 $S_i$ 方向上的分量，就是该关节要提供的广义力矩：

$$
\boxed{\,\tau_i=S_i^\top f_i\,}
$$

这是**虚功原理**的直接结果：广义力 = 空间力对广义速度方向的对偶配对。转动关节 $S_i=(0;\hat a)$，则 $\tau_i=\hat a^\top n_i$——就是力矩沿转轴的分量，符合直觉。

> 用 `+=` 而非 `=`：为兼容 **mimic 关节**（被驱动关节把力矩累加到主动关节上），故开头 `data.tau.setZero()`。详见 §7。

### 5.2 ② 力的对偶传播（叶 → 根累加）

关键：连杆 $i$ 上的力 $f_i$ **不只由自己产生**，还要**扛住所有子连杆传来的反作用力**。逆序遍历保证：处理 $i$ 时，它的所有子节点已把力累加进了 `data.f[i]`。于是 $f_i$ 已是"$i$ 及其整个子树"的合力。把它通过关节变换搬到父连杆坐标系，累加给父：

$$
f_{\lambda}\mathrel{+}= {}^{\lambda}X_i^{*}\,f_i\qquad(\texttt{f[parent] += liMi[i].act(f[i])})
$$

注意这里用 **`act`**（不是 Pass 1 的 `actInv`）——因为力是**对偶量（wrench）**，其变换是位姿伴随的**逆转置**，恰好等于 SE3 对力的 `act`（把子系的 wrench 表达到父系）。方向也相反：Pass 1 把父的运动搬进子系用 `actInv`；Pass 2 把子的力搬进父系用 `act`。运动与力是对偶的，方向恰好互逆。见 [空间代数 §4.2/§6](空间代数运算解析.md)。

---

## 6. 递推如何依赖运动学树结构

RNEA（以及所有刚体树递推）的两趟结构，本质是**对运动树做两个方向的遍历**。理解它依赖的唯一拓扑信息，就理解了"为什么两个普通 for 循环就够了"。

### 6.1 唯一依赖：`parents[i]` + DFS 编号不变量

整个算法只读一条拓扑：

```cpp
const JointIndex parent = model.parents[i];   // 每个关节唯一的父
```

再加建树时 DFS 保证的**编号不变量**（见 [multibody §2.4](multibody子系统解析.md)）：

$$
\boxed{\;\texttt{parents[i]} < i\quad(\text{父编号一定小于子})\;}
$$

这一条是"用普通 for 循环代替显式递归"的全部理由：

- **正序** `for i=1..n`：到 $i$ 时父 $\lambda(i)$（编号更小）**必已算完** → 天然的**根→叶**遍历；
- **逆序** `for i=n..1`：到 $i$ 时所有孩子（编号更大）**必已算完** → 天然的**叶→根**遍历。

无需栈、无需递归、无需显式 `children` 列表。

### 6.2 两趟 = 两个遍历方向，读/写方向对偶

| | Pass 1（正序，根→叶） | Pass 2（逆序，叶→根） |
|---|---|---|
| 流动量 | 运动 $v,a$ **顺边下行** | 力 $f$ **逆边上汇** |
| 对父的操作 | **读** `data.v[parent]` / `a_gf[parent]` | **写** `data.f[parent] += …` |
| 边变换 | `liMi[i].actInv(·)`（父系→子系） | `liMi[i].act(·)`（子系→父系） |
| 依赖满足 | `parent<i`，父此轮早写好 | `parent<i`，父稍后才处理，先收齐孩子 |

Pass 1 里 $i$ **读父**、Pass 2 里 $i$ **写父**——数据依赖方向恰好相反，而 `parents[i]<i` 让两个方向都被循环顺序自动满足。运动/力对偶，`actInv`↔`act` 也对偶互换。

### 6.3 分支节点：多孩子的力靠 `+=` 自动合并

树的关键特征是**分支**（一父多子，如腰连着两条腿）。RNEA 靠后向那句 `+=` 无痛处理：

```
        1
       / \
      2   4      ← 关节 1 有两个孩子
      |   |
      3   5
```

逆序 `5,4,3,2,1`：`f[4]+=f[5]` → `f[1]+=f[4]`（右分支）→ `f[2]+=f[3]` → `f[1]+=f[2]`（左分支）→ 处理 1 时 `f[1]` 已收齐左右两支。`data.f[i]` 的语义是"关节 $i$ 及其**整个下游子树**的合力"，**分支合并 = `subtrees[i]` 信息被 `+=` 隐式利用**，不必显式遍历 `children[i]`。正确性来自归纳：

$$
f_i=\underbrace{I_ia_i+v_i\times^{*}I_iv_i}_{\text{本连杆}}+\sum_{c\in\text{children}(i)}{}^iX_c^{*}f_c
$$

### 6.4 复杂度 $O(n)$ 与树形无关

每条树边 $\lambda(i)\!\to\!i$ 在 Pass 1、Pass 2 各被走一次（下行 `actInv`、上行 `act`）；树有 $n-1$ 条边，每次常数运算 → 总 $O(n)$，**链式或多分支都一样**（边数恒为 $n-1$）。这正是树法相对稠密矩阵法 $O(n^3)$ 的根本优势。

---

## 7. mimic 关节的处理

**mimic（仿从）关节没有独立自由度**，其运动被另一关节（mimicked，主动关节）线性绑定：

$$
q_{\text{mimic}}=s\cdot q_{\text{mimicked}}+o,\qquad \dot q_{\text{mimic}}=s\cdot\dot q_{\text{mimicked}}
$$

$s$=`scaling`（倍率）、$o$=`offset`。典型：夹爪双指联动、带传动比耦合的关节。它在**运动树里仍是真实节点**（有 id、连杆、惯量、`parents[i]`），只是不新增 DoF。RNEA 对它的**树遍历完全照常**，特殊之处只有两点。

### 7.1 机制一：`idx_v` 指向被仿关节（故名 "Mapped"）

[joint-mimic.hxx](../include/pinocchio/src/multibody/joint/joint-mimic.hxx) 的 `setMimicIndexes` 把 mimic 的 `idx_v` 设成**被仿关节的 idx_v**（`Base::i_v = v`）。而 `JointMappedVelocitySelector` 按 `idx_v` 取段（[joint-model-base.hxx:363](../include/pinocchio/src/multibody/joint/joint-model-base.hxx#L363)）：

```cpp
return SizeDepType<NV>::segment(a, idx_v(), nvExtended());
```

普通关节 `idx_v()` 是自己的槽位；mimic 关节 `idx_v()` **指向被仿关节槽位**——这就是 "**Mapped**VelocitySelector" 的含义：把 $\dot q/\ddot q/\tau$ 映射到被仿关节的 DoF 上。（mimic 自己在展开索引 `idx_vExtended` 里另有独立位置，供雅可比用，见 [multibody §2.2 三索引系统](multibody子系统解析.md)。）

### 7.2 机制二：运动子空间带耦合倍率 $S_{\text{mimic}}=s\cdot S$

对 mimic 关节，`jdata.S()` 返回的**不是普通子空间矩阵**，而是一个 `ScaledJointMotionSubspaceTpl` 对象，内部存两样东西（[joint-mimic.hxx:222](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L222)）：

```cpp
RefJointMotionSubspace m_constraint;    // 被仿关节的原始子空间
Scalar m_scaling_factor;                // ★ 耦合倍率 s
```

调用链：`jdata.S()`（[joint-data-base.hxx:190](../include/pinocchio/src/multibody/joint/joint-data-base.hxx#L190)）→ `S_accessor()`（[joint-mimic.hxx:419](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L419)）→ 返回上面这个 `S` 成员。

**关键：倍率不是在 `calc()` 里预乘成一个数值矩阵，而是分"存 / 乘"两步、且乘是惰性的（按运算各乘一次）。**

**① 存**——构造 `JointDataMimic` 时把 $s$ 记入 `m_scaling_factor`（[joint-mimic.hxx:376](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L376)），此刻不做乘法：

```cpp
JointDataMimicTpl(const RefJointData & jdata, const Scalar & scaling, ...)
: m_jdata_mimicking(checkMimic(jdata.derived()))
, S(m_jdata_mimicking.S(), scaling)     // 打包：被仿子空间 + 倍率 s
{ ... }
```

**② 乘**——`ScaledJointMotionSubspace` 的**每个运算方法各自乘一遍** `m_scaling_factor`。RNEA 用到的两处：

| RNEA 用法 | 触发方法 | 乘倍率处 |
|---|---|---|
| Pass 1 前向 `jdata.S() * a_sel` | `__mult__(v)` | [:130](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L130) `return m_scaling_factor * jm;` |
| Pass 2 后向 `jdata.S().transpose() * f[i]` | `TransposeConst::operator*` | [:162](../include/pinocchio/src/multibody/joint/joint-mimic.hxx#L162) `return ref.m_scaling_factor * (ref.m_constraint.transpose() * f);` |

其他算法的入口同样各乘一次：`matrix_impl()`（:182，物化成 6×nv 矩阵）、`se3Action/se3ActionInverse`（:136/:142，`oMi.act(S)` 变世界系，ABA/CRBA/jacobian 用）、`Inertia * S`（:250，CRBA/ABA 力集）、`Matrix6 * S`（:275，ABA）。

**这样设计零额外存储、零预计算**：倍率随 `S` 的类型走，对 FK/RNEA/CRBA/ABA/jacobian **所有算法统一自动生效**，无需每个算法单独处理 mimic 缩放——是 Pinocchio "把语义编码进类型、靠表达式模板零开销展开"的典型（见 [multibody §4](multibody子系统解析.md)）。

### 7.3 两趟里的表现

**Pass 1 前向**：`jdata.S() * jmodel.JointMappedVelocitySelector(a)` = $s\cdot S\cdot\ddot q_{\text{mimicked}}$——mimic 连杆的 $v,a$ 照常按树节点递推，只是本关节加速度项从**被仿槽位**取、且带倍率 $s$，正好落实约束 $\dot q_{\text{mimic}}=s\dot q_{\text{mimicked}}$。

**Pass 2 后向**：力向上传播 `f[parent] += liMi.act(f[i])` 对 mimic **照常**（它是有质量的真实节点）。但力矩投影：

```cpp
jmodel.JointMappedVelocitySelector(data.tau) += jdata.S().transpose() * data.f[i];
//     ↑ 写到"被仿关节槽位"                       ↑ (s·S)ᵀ f_i
```

被映射到**被仿关节槽位**、并带倍率 $s$，即 $s\,S^\top f_{\text{mimic}}$——mimic 约束通过虚功**反射回主动关节**的力矩。

### 7.4 为什么必须 `tau.setZero()` + `+=`

[rnea.hxx:137](../include/pinocchio/src/algorithm/rnea.hxx#L137) 注释点破：被仿关节自己**和它的每个 mimic** 都往**同一个 tau 槽位**写。用赋值 `=` 会互相覆盖，必须先清零、全程 `+=` 叠加：

$$
\boxed{\;\tau_{\text{driver}}=\underbrace{S^\top f_{\text{driver}}}_{\text{主动关节自身}}+\sum_{k\in\text{mimics}}\underbrace{s_k\,S^\top f_k}_{\text{每个仿从反射回来}}\;}
$$

物理含义：主动电机的力矩 = 驱动自己那节 + 驱动所有联动仿从连杆（各按耦合比 $s_k$ 折算）。一个电机扛着整条联动链——虚功原理在耦合约束下的体现。

> 对照：[ABA](正动力学ABA解析.md) 和部分算法**不支持** mimic（`assert(MimicChecker())`），因为 mimic 破坏了每关节独立求逆的结构；RNEA 因是纯线性投影/累加，能天然容纳。[CRBA](质量矩阵CRBA解析.md) 则需专门的 `MimicStep` 补齐（被驱动关节对矩阵块的贡献）。

---

## 8. 外力 `fext` 与电机转子惯量 `armature`

### 8.1 外力版 `rnea(q,v,a,fext)`

`fext[i]` 是作用在连杆 $i$ 上的外部空间力（如接触力），在**关节 $i$ 局部系**表达。只需在 Pass 1 之后从本体力里扣掉：

```cpp
Pass1::run(...);
data.f[i] -= fext[i];     // 外力抵消一部分需要关节提供的力
```

得到的 $\tau=g(q)+C\dot q+M\ddot q-\sum J_i^\top f_i^{\text{ext}}$，即外力通过雅可比转置减轻/加重关节负担。

### 8.2 转子惯量 `armature`

```cpp
data.tau.array() += model.armature.array() * a.array();
```

`armature` 是电机转子经减速比折算的**反射惯量**（reflected inertia），在关节空间是纯对角项 $\text{diag}(\text{armature})\ddot q$，直接加到 $\tau$。它使有效质量矩阵变成 $M+\text{diag}(\text{armature})$，改善数值条件、贴合真实带减速器的关节。

---

## 9. 三个派生量：`nle` / `g` / `staticTorque`

同一套 RNEA 递推，喂不同的输入即得动力学方程的不同分块：

| 函数 | 输入设定 | 输出 | 物理量 |
|---|---|---|---|
| `rnea(q,v,a)` | 完整 | $\tau=M\ddot q+C\dot q+g$ | 完整逆动力学 |
| `nonLinearEffects(q,v)` | $\ddot q=0$ | $b=C\dot q+g$ | 非线性项（科氏+离心+重力）|
| `computeGeneralizedGravity(q)` | $\dot q=0,\ddot q=0$ | $g(q)$ | 广义重力 |
| `computeStaticTorque(q,fext)` | $\dot q=0,\ddot q=0$ + 外力 | $g-\sum J^\top f_{\text{ext}}$ | 静态平衡力矩 |

原理：RNEA 输出对 $(\ddot q,\dot q)$ 分别是线性/二次的。令 $\ddot q=0$ 抹掉 $M\ddot q$ 得 $b$；再令 $\dot q=0$ 抹掉科氏/离心（含 $c$ 和 $v\times v$ 项）只剩重力 $g$。这些专用函数内部用**简化过的单趟递推**（`NLEForwardStep` 等）省去无关计算，但数学等价。

- 由此可反推质量矩阵的一列：$M\ddot q=\text{rnea}(q,0,\ddot q)-g(q)$，取 $\ddot q=e_k$ 即得 $M$ 第 $k$ 列（不过组装整个 $M$ 用 CRBA 更快，见 [质量矩阵CRBA解析.md](质量矩阵CRBA解析.md)）。

---

## 10. 复杂度与和 CRBA/ABA 的关系

- **复杂度** $O(n)$：两趟线性遍历，每关节常数次 6D 运算。是所有动力学量里最便宜的。
- **RNEA ↔ CRBA**：CRBA 组装 $M$ 本质上是"批量 RNEA"——对每个自由度算一列。
- **RNEA ↔ ABA**：正动力学 $\ddot q=M^{-1}(\tau-b)$。其中 $b=\text{RNEA}(q,\dot q,0)=$ `nonLinearEffects`。所以"CRBA 求 $M$ + Cholesky 解 + RNEA 求 $b$"是 ABA 之外求正动力学的另一条路（$O(n^2\sim n^3)$，ABA 是 $O(n)$，见 [正动力学ABA解析.md](正动力学ABA解析.md)）。
- **RNEA ↔ 导数**：`computeRNEADerivatives` 在同一递推里附带算 $\partial\tau/\partial q,\partial\tau/\partial\dot q$，且 $\partial\tau/\partial\ddot q=M$（顺便得到质量矩阵）。见 [算法总览 §I](算法总览.md)。

---

## 11. 一句话总结

RNEA = **牛顿-欧拉两趟 $O(n)$ 递推**：

$$
\underbrace{\text{根→叶}}_{\text{传播 }v,a,\ \text{算 }f_i=I_ia_i+v_i\times^*I_iv_i}\quad\Longrightarrow\quad
\underbrace{\text{叶→根}}_{\tau_i=S_i^\top f_i,\ f_\lambda\mathrel{+}={}^\lambda X_i^* f_i}
$$

- **重力**靠"根加速度 $=-g$"的等效技巧免费融入；
- **前向用 `actInv` 搬运动、后向用 `act` 搬力**（运动/力对偶，方向互逆）；
- **树结构**：靠 `parents[i]<i` 不变量，正/逆序 for 循环即实现根↔叶两向遍历；分支靠 `+=` 自动合并；$O(n)$ 与树形无关（§6）；
- **mimic 关节**：树遍历照常，仅把 $\dot q/\ddot q/\tau$ 经 `idx_v` 映射到被仿槽位、子空间带倍率 $s$，`tau.setZero()`+`+=` 把力矩反射累加到主动关节（§7）；
- **$\tau=S^\top f$** 是虚功原理；`nle`/`g`/`staticTorque` 都是它令 $\ddot q$ 或 $\dot q$ 为零的特例。
