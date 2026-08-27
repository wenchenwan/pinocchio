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
- [6. 外力 `fext` 与电机转子惯量 `armature`](#6-外力-fext-与电机转子惯量-armature)
- [7. 三个派生量：`nle` / `g` / `staticTorque`](#7-三个派生量nle--g--statictorque)
- [8. 复杂度与和 CRBA/ABA 的关系](#8-复杂度与和-crbaaba-的关系)
- [9. 一句话总结](#9-一句话总结)

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

> 用 `+=` 而非 `=`：为兼容 **mimic 关节**（被驱动关节把力矩累加到主动关节上，见 [multibody §2 mimic](multibody子系统解析.md)），故开头 `data.tau.setZero()`。

### 5.2 ② 力的对偶传播（叶 → 根累加）

关键：连杆 $i$ 上的力 $f_i$ **不只由自己产生**，还要**扛住所有子连杆传来的反作用力**。逆序遍历保证：处理 $i$ 时，它的所有子节点已把力累加进了 `data.f[i]`。于是 $f_i$ 已是"$i$ 及其整个子树"的合力。把它通过关节变换搬到父连杆坐标系，累加给父：

$$
f_{\lambda}\mathrel{+}= {}^{\lambda}X_i^{*}\,f_i\qquad(\texttt{f[parent] += liMi[i].act(f[i])})
$$

注意这里用 **`act`**（不是 Pass 1 的 `actInv`）——因为力是**对偶量（wrench）**，其变换是位姿伴随的**逆转置**，恰好等于 SE3 对力的 `act`（把子系的 wrench 表达到父系）。方向也相反：Pass 1 把父的运动搬进子系用 `actInv`；Pass 2 把子的力搬进父系用 `act`。运动与力是对偶的，方向恰好互逆。见 [空间代数 §4.2/§6](空间代数运算解析.md)。

---

## 6. 外力 `fext` 与电机转子惯量 `armature`

### 6.1 外力版 `rnea(q,v,a,fext)`

`fext[i]` 是作用在连杆 $i$ 上的外部空间力（如接触力），在**关节 $i$ 局部系**表达。只需在 Pass 1 之后从本体力里扣掉：

```cpp
Pass1::run(...);
data.f[i] -= fext[i];     // 外力抵消一部分需要关节提供的力
```

得到的 $\tau=g(q)+C\dot q+M\ddot q-\sum J_i^\top f_i^{\text{ext}}$，即外力通过雅可比转置减轻/加重关节负担。

### 6.2 转子惯量 `armature`

```cpp
data.tau.array() += model.armature.array() * a.array();
```

`armature` 是电机转子经减速比折算的**反射惯量**（reflected inertia），在关节空间是纯对角项 $\text{diag}(\text{armature})\ddot q$，直接加到 $\tau$。它使有效质量矩阵变成 $M+\text{diag}(\text{armature})$，改善数值条件、贴合真实带减速器的关节。

---

## 7. 三个派生量：`nle` / `g` / `staticTorque`

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

## 8. 复杂度与和 CRBA/ABA 的关系

- **复杂度** $O(n)$：两趟线性遍历，每关节常数次 6D 运算。是所有动力学量里最便宜的。
- **RNEA ↔ CRBA**：CRBA 组装 $M$ 本质上是"批量 RNEA"——对每个自由度算一列。
- **RNEA ↔ ABA**：正动力学 $\ddot q=M^{-1}(\tau-b)$。其中 $b=\text{RNEA}(q,\dot q,0)=$ `nonLinearEffects`。所以"CRBA 求 $M$ + Cholesky 解 + RNEA 求 $b$"是 ABA 之外求正动力学的另一条路（$O(n^2\sim n^3)$，ABA 是 $O(n)$，见 [正动力学ABA解析.md](正动力学ABA解析.md)）。
- **RNEA ↔ 导数**：`computeRNEADerivatives` 在同一递推里附带算 $\partial\tau/\partial q,\partial\tau/\partial\dot q$，且 $\partial\tau/\partial\ddot q=M$（顺便得到质量矩阵）。见 [算法总览 §I](算法总览.md)。

---

## 9. 一句话总结

RNEA = **牛顿-欧拉两趟 $O(n)$ 递推**：

$$
\underbrace{\text{根→叶}}_{\text{传播 }v,a,\ \text{算 }f_i=I_ia_i+v_i\times^*I_iv_i}\quad\Longrightarrow\quad
\underbrace{\text{叶→根}}_{\tau_i=S_i^\top f_i,\ f_\lambda\mathrel{+}={}^\lambda X_i^* f_i}
$$

- **重力**靠"根加速度 $=-g$"的等效技巧免费融入；
- **前向用 `actInv` 搬运动、后向用 `act` 搬力**（运动/力对偶，方向互逆）；
- **$\tau=S^\top f$** 是虚功原理；`nle`/`g`/`staticTorque` 都是它令 $\ddot q$ 或 $\dot q$ 为零的特例。
