# Pinocchio 源码中的 C++ 语法与惯用技巧

> 阅读 Pinocchio 源码时反复出现的**模板 / Eigen / Boost 语法技巧**集中在这里。
> 每条：**语法是什么 → 最小示例 → 在 Pinocchio 里为什么这么用 → 典型出现处**。
>
> - 机制/算法层面的解析见 [docs/源码解析.md](源码解析.md)；
> - 空间代数运算见 [docs/空间代数运算解析.md](空间代数运算解析.md)；
> - **本文档只收"语法招式"本身**，供快速查阅。持续追加。

---

## 目录

1. [编译期尺寸守卫 `EIGEN_STATIC_ASSERT_*`](#1-编译期尺寸守卫-eigen_static_assert_)
2. [编译期成员检查 `BOOST_MPL_ASSERT`](#2-编译期成员检查-boost_mpl_assert)
3. [依赖名消歧：`typename` 与 `.template` / `::template`](#3-依赖名消歧typename-与-template--template)
4. [CRTP：`Base<Derived>` 与 `.derived()`](#4-crtpbasederived-与-derived)
5. [traits 特化：给类型挂"元数据"](#5-traits-特化给类型挂元数据)
6. [Eigen 惯用法：`MatrixBase<T>` / `Ref` / `const_cast_derived` / `noalias`](#6-eigen-惯用法matrixbaset--ref--const_cast_derived--noalias)
7. [接口双重载：in/out 参数版 + 返回值版](#7-接口双重载inout-参数版--返回值版)
8. [Tag dispatch：用空类型 `Blank` 消歧重载](#8-tag-dispatch用空类型-blank-消歧重载)
9. [模板参数技巧：模板的模板参数 + context 默认标量](#9-模板参数技巧模板的模板参数--context-默认标量)
10. [autodiff 友好分支：`if_then_else`](#10-autodiff-友好分支if_then_else)
11. [宏技巧：续行与显式实例化](#11-宏技巧续行与显式实例化)
12. [`fwd.hxx` 前向声明中心与默认模板实参"只写一次"](#12-fwdhxx-前向声明中心与默认模板实参只写一次)
13. [按输入种类分流重载 + `PlainObject` 返回类型](#13-按输入种类分流重载--plainobject-返回类型)
14. [分派 struct + 空特化 + 模板化 `run`](#14-分派-struct--空特化--模板化-run编译期选分支运行时接对象)
15. [`ref_selector`：自动选择"引用类型"](#15-ref_selector自动选择引用类型)
16. [用 traits 计算类型：运算返回类型/标量 + CRTP 基类委托](#16-用-traits-计算类型运算返回类型标量--crtp-基类委托)
17. [用户自定义转换运算符 `operator T()`](#17-用户自定义转换运算符-operator-t)

---

## 1. 编译期尺寸守卫 `EIGEN_STATIC_ASSERT_*`

**是什么**：`static_assert` + Eigen 编译期尺寸常量，在编译期检查传入的 Eigen 类型维度是否正确。

**实现**（`include/pinocchio/src/eigen-macros.hxx`）：

```cpp
#define PINOCCHIO_EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE_OR_DYNAMIC(TYPE, SIZE)  \
  static_assert(                                                                    \
    TYPE::IsVectorAtCompileTime                                                     \
    && (TYPE::SizeAtCompileTime == Eigen::Dynamic || TYPE::SizeAtCompileTime == SIZE))

#define PINOCCHIO_EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE_OR_DYNAMIC(TYPE, ROWS, COLS)  \
  static_assert(                                                                          \
    (TYPE::RowsAtCompileTime == Eigen::Dynamic || TYPE::RowsAtCompileTime == ROWS)        \
    && (TYPE::ColsAtCompileTime == Eigen::Dynamic || TYPE::ColsAtCompileTime == COLS))
```

**要点**：
- `TYPE::IsVectorAtCompileTime` / `RowsAtCompileTime` / `SizeAtCompileTime` 是 Eigen 每个类型都带的**编译期常量**。
- `Eigen::Dynamic`（哨兵值 −1）表示"编译期未知尺寸"。`== Eigen::Dynamic ||` 那一支实现"**动态尺寸放行**"。
- 单参数 `static_assert(cond)`（C++17 形式，无消息）；行尾 `\` 是宏续行；`;` 来自调用处。
- **纯编译期、零运行时开销**，类型不合就在调用行报编译错。

**代入示例**（向量宏，SIZE=3）：`Vector3d`✅、`Vector4d`❌、`VectorXd`✅（动态放行）、`Matrix3d`❌（非向量）。

**为什么用**：模板函数接受任意 Eigen 类型，用它强制"参数契约"（如 `skew` 只收 3 维向量 + 3×3 矩阵），
把维度错误在编译期拦下、给可读报错，而非在 Eigen 深处撞晦涩错误或算出错误结果。

**典型出现处**：`spatial/skew.hxx`、`spatial/explog.hxx` 等几乎每个接受 Eigen 类型的函数开头。

---

## 2. 编译期成员检查 `BOOST_MPL_ASSERT`

**是什么**：编译期断言"某类型是否在一个类型列表里"。

```cpp
// joint-generic.hxx：构造通用 JointModel 时，静态检查该具体关节确实属于当前关节目录
BOOST_MPL_ASSERT((boost::mpl::contains<typename JointModelVariant::types, JointModelDerived>));
```

- `JointModelVariant::types` 是 `boost::variant` 暴露的**编译期类型列表**（MPL 序列）；
- `boost::mpl::contains<Seq, T>` 是编译期元函数，判断"列表是否含 T"→ 产出 `true_/false_`；
- 不在列表则**编译报错**（清晰信息），而非在 variant 深处撞错。

**典型出现处**：[源码解析条目 4 特性 #7](源码解析.md)、`joint-generic.hxx` 的装箱构造函数。

---

## 3. 依赖名消歧：`typename` 与 `.template` / `::template`

**是什么**：当名字**依赖模板参数**时，编译器无法判断它是类型还是模板，需手动声明。

```cpp
typename CastType<NewScalar, ModelTpl<...>>::type   // typename：声明 ::type 是「类型」
obj.template cast<NewScalar>()                       // .template：声明 cast 是「模板成员」（对象）
Type::template algo<JointModelRX>                    // ::template：声明 algo 是「模板成员」（类型）
```

**规律**：
- `依赖类型::成员` 且成员是**类型** → 前面加 `typename`；
- `依赖对象.成员<...>` / `依赖类型::成员<...>` 且成员是**模板** → 中间加 `template`（否则 `<` 被当小于号）；
  区别只是 `.`（访问对象成员）vs `::`（访问类型成员）。

**典型出现处**：[源码解析条目 2.2/2.3](源码解析.md)（cast）、[条目 9.3](源码解析.md)（`::template algo`）。

---

## 4. CRTP：`Base<Derived>` 与 `.derived()`

**是什么**：基类以"派生类自己"为模板参数，从而在基类里静态调用派生类实现——**无虚函数的多态**。

```cpp
template<typename Derived> struct Base {
  void f() { static_cast<Derived*>(this)->f_impl(); }   // 基类回调派生实现
  Derived & derived() { return static_cast<Derived&>(*this); }
};
struct Concrete : Base<Concrete> { void f_impl() { ... } };  // 继承 Base<自己>
```

**在 Pinocchio 里**：`JointModelBase<JointModelTpl<...>>`、`SE3Base<SE3Tpl<...>>`、`MotionDense<...>` 全是 CRTP。
`.derived()` 把基类引用还原成具体类型（Eigen 也用同名手法从 `MatrixBase<T>` 还原具体表达式）。

**典型出现处**：[源码解析条目 4（joint-generic）](源码解析.md)、[PINOCCHIO_GUIDE §11.3 技巧A](../PINOCCHIO_GUIDE.md)、几乎所有 spatial 类型。

---

## 5. traits 特化：给类型挂"元数据"

**是什么**：用一个可特化的 `traits<T>` 模板，为类型 `T` 附上一组关联类型/常量，而不改动 `T` 本身。

```cpp
template<typename T> struct traits;                 // 通用声明
template<typename S,int O,...> struct traits<JointTpl<S,O,...>> {   // 针对某类型特化
  typedef S Scalar;  typedef SE3Tpl<S,O> Transformation_t;  static constexpr int NQ = ...;  // 一包元数据
};
// 用：typename traits<X>::Scalar
```

**在 Pinocchio 里**：`traits<>`（关节/spatial 类型的关联类型）、`CastType<NewScalar, C>::type`（"cast 后的类型"定制点）、
`SE3GroupAction<D>::ReturnType`（作用返回类型）都是这个套路。

**典型出现处**：[源码解析条目 2.2/4](源码解析.md)、`joint-generic.hxx`、`spatial/*-base.hxx`。

---

## 6. Eigen 惯用法：`MatrixBase<T>` / `Ref` / `const_cast_derived` / `noalias`

| 写法 | 含义 | 为什么 |
|------|------|--------|
| `const Eigen::MatrixBase<T> & x` + `x.derived()` | 收任意 Eigen 表达式的 CRTP 基类，再还原具体类型 | 一个模板接受 `VectorXd`/块/`Ref`/`Map`，零拷贝 |
| `Eigen::Ref<const context::VectorXs>` | 对某标量动态向量的**只读引用**参数 | 让 numpy 切片、块表达式也能零拷贝匹配（见显式实例化） |
| `x.const_cast_derived()` | 去掉 `const MatrixBase&` 的 const 再取具体类型 | "名义 const、实则输出"的 in/out 矩阵参数写入用 |
| `A.noalias() = B * C` | 告诉 Eigen "结果与输入不重叠"，跳过临时量 | 省一次拷贝，动力学循环里大量用 |

**典型出现处**：几乎所有 spatial/algorithm 函数；[空间代数解析各节](空间代数运算解析.md)、[源码解析条目 7.3](源码解析.md)（`const_cast_derived`）、[条目 1.2](源码解析.md)（`Eigen::Ref`）。

---

## 7. 接口双重载：in/out 参数版 + 返回值版

**是什么**：同一运算提供两个重载——一个把结果**写进给定输出参数**（无返回，省拷贝），一个**返回新对象**（写起来顺手）。

```cpp
void skew(const V & v, const M & out);   // in/out 版：结果写进 out
M    skew(const V & v) { M m; skew(v, m); return m; }   // 返回值版：内部转调 in/out 版
```

**为什么**：性能敏感处（循环内）用 in/out 版避免临时对象；一般代码用返回值版可读。**返回值版几乎总是转调 in/out 版**。

**典型出现处**：`skew`/`unSkew`/`alphaSkew`/`skewSquare`（[空间代数解析 §2](空间代数运算解析.md)）、`toActionMatrix`、exp/log 等。

---

## 8. Tag dispatch：用空类型 `Blank` 消歧重载

**是什么**：两个重载参数个数/位置会撞车时，塞一个**空标签类型**占位来区分。

```cpp
void calc(JointData & d, const ConfigVector & q);                 // (a) 只有配置 q
void calc(JointData & d, const Blank, const TangentVector & v);   // (b) 只有速度 v（Blank 占位区分）
void calc(JointData & d, const ConfigVector & q, const TangentVector & v);  // (c) 两者
```

若无 `Blank`，(a) 和 (b) 都是"data + 一个向量"，签名冲突。`Blank` 让 (b) 变成三参数，明确"这是速度那一版"。

**典型出现处**：`joint-generic.hxx` 的三个 `calc`（[源码解析条目 7.1](源码解析.md)）。

---

## 9. 模板参数技巧：模板的模板参数 + context 默认标量

**是什么**：
- **模板的模板参数**：参数本身是"还差实参才能实例化的模板"，写成 `template<typename,int> class JointCollectionTpl`；
- **context 默认标量**：核心类型都是 `XxxTpl<Scalar,...>`，`pinocchio::Model` 等是**默认标量的别名**，
  由 `context.hpp` 的一个可覆盖开关决定默认是 `double`。

```cpp
template<typename Scalar, int Options,
         template<typename,int> class JointCollectionTpl = JointCollectionDefaultTpl>  // 模板的模板参数(带默认)
struct ModelTpl;
typedef ModelTpl<context::Scalar, context::Options> Model;   // 别名（第3参用默认）
```

**为什么**：一套代码通吃 `double`/`float`/autodiff/高精度（换 `Scalar`）与不同关节目录（换 `JointCollectionTpl`）。

**典型出现处**：[源码解析条目 2/4](源码解析.md)、[PINOCCHIO_GUIDE §11.2/§11.4](../PINOCCHIO_GUIDE.md)。

---

## 10. autodiff 友好分支：`if_then_else`

**是什么**：用**函数式的三元选择** `internal::if_then_else(op, a, b, then, else)` 代替普通 `if`，
让分支在自动微分标量下也能被正确记录。

```cpp
// 普通 double：等价 (t > prec) ? sinθ/θ : 1 - θ²/6
alpha_v = internal::if_then_else(GT, t, prec, st/t, 1 - t2/6);
```

**为什么**：普通 `if` 在 AD 标量下会**丢掉未走分支的导数信息**；`if_then_else` 对 AD 标量把**两条分支都录进磁带**，
所以含分支的数学（exp/log 的小角度奇异处理等）也能自动微分。

**典型出现处**：`spatial/explog.hxx`（exp3/exp6/log 的泰勒兜底）、[空间代数解析 §7.4](空间代数运算解析.md)、[源码解析条目 3](源码解析.md)。

---

## 11. 宏技巧：续行与显式实例化

**续行 `\`**：多行宏每行末尾用 `\` 连成一条逻辑语句（见 §1 的尺寸守卫宏）。整个宏是**一条语句**，
末尾 `;` 由调用处 `MACRO(...);` 提供。

**显式模板实例化**（`.cpp` 里）：`template ... forwardKinematics<context::Scalar,...>(...);`——
开头 `template` + 无函数体 + 分号，命令编译器**为默认标量预编译一份**并导出符号，加速下游编译。
详见 [源码解析条目 1](源码解析.md)。

---

## 12. `fwd.hxx` 前向声明中心与默认模板实参"只写一次"

**是什么**：每个模块（`spatial/`、`multibody/`、`algorithm/`…）都有一个 `fwd.hxx`/`fwd.hpp`，**只声明、不实现**：
集中前向声明本模块所有类模板、定义默认标量别名、给出 traits 定制点。

```cpp
// spatial/fwd.hxx（节选）
template<typename Scalar, int Options = context::Options> class MotionTpl;   // ← 默认实参在此声明一次
template<typename _Scalar, int _Options = context::Options> struct SE3Tpl;
using SE3 = SE3Tpl<context::Scalar, context::Options>;                       // 默认标量别名
namespace internal { template<typename C,typename N,typename S> struct cast_call_normalize_method; } // 定制点
```

**核心意义：打破头文件循环依赖**。spatial 里类型互相引用（`SE3::act()` 提到 `Motion`，`Motion` 变换又提到 `SE3`）。
若直接互相 `#include` 就成循环。前向声明让各头**只需知道"类型存在"**即可用引用/指针/模板参数互相提及，
完整定义留到 `.hxx` 再 include。**看到 `fwd` 就知道"只声明不实现"。**

**关键规则：默认模板实参只能指定一次**。C++ 规定同一模板的默认实参全程只写一次，所以 `= context::Options`
**只在 fwd 头给**——后面真正定义 `MotionTpl` 时**不能再重复**（重复即编译错）。把默认实参集中在 fwd 头正是为满足此规则。

**是不是高级用法**：前向声明是基础技巧，但"默认实参只写一次""traits 定制点集中声明"是易踩坑的规范化组织手法。

**典型出现处**：`spatial/fwd.hxx`、`multibody/fwd.hxx`、各模块 `fwd.*`。

---

## 13. 按输入种类分流重载 + `PlainObject` 返回类型

**是什么**：同名函数按**输入的 Eigen 种类**（一般矩阵 vs `Map` 视图）拆成两个重载，返回类型也随之不同。

```cpp
// A：一般 Eigen 对象 —— 返回类型「跟随输入」
template<typename EigenDerived>
typename EigenDerived::PlainObject                                     // 输入求值后的具体拥有型
actOnEigenObject(const Eigen::MatrixBase<EigenDerived> & p) const
{ return (rotation() * p + translation()).eval(); }                    // R·p+t，.eval() 求值

// B：Eigen::Map 视图 —— 硬编码返回 Vector3
template<typename MapDerived>
Vector3
actOnEigenObject(const Eigen::MapBase<MapDerived> & p) const
{ return Vector3(rotation() * p + translation()); }
```

**两点区别**：
- **输入基类**：`MatrixBase`（一般对象）vs `MapBase`（`Eigen::Map`，套在外部裸内存上的视图，不拥有内存）。
  `Map` 同时 is-a `MatrixBase` 和 `MapBase`，但 `MapBase` 重载**更专一**，故传 Map 时选 B、传其它选 A。
- **返回类型**：
  - `EigenDerived::PlainObject` = 表达式**求值后拥有存储的具体类型**，**跟随输入**（保住输入的标量/选项，通用）；`.eval()` 强制求值，避免返回悬垂表达式。
  - Map 版硬编码 `Vector3`：因 Map 是视图、**不能拥有** `R·p+t` 的新结果，只能返回一个真正拥有内存的具体向量。

**为什么不合并**：只留 A → 传 Map 时 PlainObject 语义不适合当返回值；只留 B → 丢掉输入类型信息（AD 标量会被错当 double）。
故 A 走"类型跟随输入"通用路径、B 给 Map 特殊情形专门路径。

> 与 [§7 双重载](#7-接口双重载inout-参数版--返回值版) 的区别：§7 按"in/out vs 返回值"分，本条按"输入种类"分。

**补充概念 `PlainObject`**：每个 Eigen 表达式类型都带 `::PlainObject` typedef = 该表达式求值后的、拥有存储的具体矩阵类型
（`(A+B)::PlainObject` 是个真 `Matrix`）。返回 `PlainObject` = "返回和输入同形状/标量的真矩阵，而非惰性表达式"。

**典型出现处**：`spatial/se3-tpl.hxx` 的 `actOnEigenObject`。

---

## 14. 分派 struct + 空特化 + 模板化 `run`（编译期选分支、运行时接对象）

**是什么**：把"按类型选行为"做成一个可特化的 `struct`，用**不同特化**给不同类型不同实现；其中**最常见/无事可做的情形用空特化**（零开销）；`struct` 里的 `static ... run(...)` 写成**模板函数**，用实参推导接住真正要操作的对象。

**实例**：`SE3::cast<NewScalar>()` 之后要不要重新归一化旋转（`spatial/se3-tpl.hxx`）：

```cpp
// 调用点（cast 内部）：res 是【目标】SE3 —— SE3Tpl<NewScalar,Options>
internal::cast_call_normalize_method<SE3Tpl, NewScalar, Scalar>::run(res);
//                                   └源类型  └新标量 └旧标量

// 特化①：新标量 == 旧标量 —— 纯拷贝、无误差、无事可做 → 空实现（编译期消解为零代码）
template<typename Scalar, int Options>
struct cast_call_normalize_method<SE3Tpl<Scalar,Options>, Scalar, Scalar> {
  template<typename T> static void run(T &) {}          // 参数不命名 = 故意忽略
};

// 特化②：新标量 != 旧标量 —— 按精度关系决定是否 normalize
template<typename Scalar, int Options, typename NewScalar>
struct cast_call_normalize_method<SE3Tpl<Scalar,Options>, NewScalar, Scalar> {
  template<typename T> static void run(T & self) {
    if (pinocchio::cast<NewScalar>(Eigen::NumTraits<Scalar>::epsilon())
        > Eigen::NumTraits<NewScalar>::epsilon())
      self.normalize();   // 旧误差在新精度下"可见"才归一化（float→double 归一化；double→float 跳过）
  }
};
```

**三个要点**：

1. **编译期分派**：按 `<类, 新标量, 旧标量>` 三元组特化，编译器据 cast 目标类型**在编译期选中特化**——
   决策本质是类型问题（标量是否改变），放进类型系统而非运行时 `if`。
2. **空特化 = 零开销优化常见路径**：`Scalar→Scalar`（同类型 cast）是纯拷贝、不引入误差、无需归一化，故空实现。
   写成专门空特化而非"落进通用分支靠运行时判断跳过"，让最常见情形**内联后化为乌有**，一条指令都不生成。
3. **`run` 为什么是模板函数**：它要操作的对象 `res` 是 **cast 的结果（目标标量 SE3）**，与 struct 被特化的**源类型
   不是同一个类型**。写 `template<typename T> run(T&)` 让 `T` 从实参推导接住目标对象，无需在签名里重拼其类型；
   还对 `CastType` 产出的实际类型更健壮。**这样"struct 负责选分支、模板化 run 负责接对象"两件事干净分开。**

> 这是全库大量 `static ... run(...)` 分派 helper 的通用套路：**struct 定分支，`run` 接实参**。
> 与 [§5 traits 特化](#5-traits-特化给类型挂元数据)、[§10 `if_then_else`](#10-autodiff-友好分支if_then_else) 同属"把决策前移/类型化"的思路。

**典型出现处**：`spatial/se3-tpl.hxx` 的 `cast_call_normalize_method`（定制点声明在 `spatial/fwd.hxx`，见 [§12](#12-fwdhxx-前向声明中心与默认模板实参只写一次)）。

---

## 15. `ref_selector`：自动选择"引用类型"

**是什么**：包装 Eigen 内部 trait `ref_selector<D>`，自动挑选"引用 `D` 的正确类型"——泛型代码里想表达"给我一个引用 `D` 的类型"，但不知 `D` 是实体矩阵还是表达式时用它。

```cpp
// include/pinocchio/src/eigen-macros.hxx
#define PINOCCHIO_EIGEN_REF_CONST_TYPE(D) Eigen::internal::ref_selector<D>::type          // 只读引用类型
#define PINOCCHIO_EIGEN_REF_TYPE(D)       Eigen::internal::ref_selector<D>::non_const_type // 可写引用类型
```

```cpp
// 用例（motion-tpl.hxx traits）：toVector() 要零拷贝返回对内部 Vector6 的引用
typedef typename PINOCCHIO_EIGEN_REF_CONST_TYPE(Vector6) ToVectorConstReturnType;  // = const Vector6&
typedef typename PINOCCHIO_EIGEN_REF_TYPE(Vector6)       ToVectorReturnType;       // = Vector6&
```

**核心：为什么不直接写 `const D&`**——Eigen 类型分两类，引用方式应不同：

| `D` 是什么 | `ref_selector<D>::type` | 为什么 |
|-----------|--------------------------|--------|
| **实体矩阵**（`Vector6`/`Matrix<...>`，拥有存储） | `const D&`（真引用） | 有实际内存，引用零拷贝 |
| **表达式模板**（`A+B`、块，轻量代理、常是临时量） | `const D`（**按值**） | 表达式是廉价代理；持有其引用会**悬垂**，按值更安全 |

`::type` 给只读版、`::non_const_type` 给可写版。硬写 `const D&` 对表达式是错的，交给 `ref_selector` 就能自动选对。

**两个细节**：
- **为何包成宏**：`Eigen::internal::` 是 Eigen 私有命名空间；包一层，Eigen 若改内部 API 只需改这一处宏（隔离依赖）。
- **调用处要加 `typename`**：`ref_selector<D>::type` 是依赖名（见 [§3](#3-依赖名消歧typename-与-template--template)），故用处写 `typedef typename PINOCCHIO_EIGEN_REF_CONST_TYPE(...) ...`。

**典型出现处**：各 spatial 类型 traits 的 `ToVectorReturnType` 等（`motion-tpl.hxx`/`force-tpl.hxx`）。

---

## 16. 用 traits 计算类型：运算返回类型/标量 + CRTP 基类委托

**是什么**：把"某类型的返回类型/标量/拥有型是什么"这类**类型关系**集中定义在 `traits<T>` 里，运算/operator 只引用
`traits<T>::Xxx`，从而与具体类型细节解耦。有两个高频用法。

### 16.1 用 traits 定 operator 的返回类型与标量

```cpp
// motion-dense.hxx：自由函数 operator，参数收 CRTP 基类，返回类型/标量走 traits
template<typename M1, typename M2>
typename traits<M1>::MotionPlain                                   // 返回"拥有型"Motion
operator^(const MotionDense<M1> & v1, const MotionDense<M2> & v2)
{ return v1.derived().cross(v2.derived()); }

template<typename M1>
typename traits<M1>::MotionPlain
operator*(const typename traits<M1>::Scalar alpha, const MotionDense<M1> & v)  // 标量类型走 traits
{ return v * alpha; }
```

**三个价值**：
- **① 泛型**：参数 `const MotionDense<M1>&`（CRTP 基类）+ 返回 `traits<M1>::MotionPlain`，**一份模板通吃家族所有后端**
  （`MotionTpl`/`MotionRef`/表达式…），无需逐类型重载。
- **② 返回"拥有型"避免悬垂**：`traits<M1>::MotionPlain` = 该类型对应的**自有内存**类型（永远是 `MotionTpl`）。
  输入若是 `MotionRef`（视图）或表达式，结果**必须物化成拥有型**才能安全返回——`MotionPlain` 就是干这个的
  （`traits<MotionRef>::MotionPlain = MotionTpl`）。
- **③ 标量自动跟随**：`traits<M1>::Scalar` 让 `alpha` 类型跟随 Motion 的标量（double/float/AD 通用），并约束重载、避免劫持无关的 `X * 类型`。

### 16.2 返回类型 traits 对 CRTP 基类的"委托特化"

```cpp
// motion-dense.hxx：当查询目标是 CRTP 基类 MotionDense<Derived> 时，转交给具体 Derived
template<typename Derived>
struct SE3GroupAction<MotionDense<Derived>> {
  typedef typename SE3GroupAction<Derived>::ReturnType ReturnType;   // 剥掉 MotionDense<> 包装，下放给 Derived
};
template<typename Derived, typename MotionDerived>
struct MotionAlgebraAction<MotionDense<Derived>, MotionDerived> {
  typedef typename MotionAlgebraAction<Derived, MotionDerived>::ReturnType ReturnType;
};
```

- `SE3GroupAction<T>::ReturnType` = "SE3 作用在 T 上"的结果类型；`MotionAlgebraAction<T,M>` = "T 对 M 叉乘"的结果类型。
- **为什么委托**：`MotionDense` 是通用 CRTP 基类，不知道结果具体类型；真正答案由**具体后端**给（视图→`MotionPlain`）。
  故这两个偏特化**接住基类类型、转发给 `Derived`**，让"用基类类型查询"也能落到派生类的正确结果。
- 解析链示例：`SE3GroupAction<MotionDense<MotionRef>>::ReturnType` →（本特化下放）→ `SE3GroupAction<MotionRef>::ReturnType`
  →（MotionRef 特化）→ `MotionPlain`。**没有这个特化，用基类类型查询就会落到主模板而断链。**
- 本质：**"遇到 `MotionDense<Derived>` 就还原成 `Derived`"**——是运算里 `derived()`（对象层解包）在**类型层**的对应物。

> 两者共同点：**`traits` 是"类型关系的单一事实源"**。改一处 traits，所有引用它的 operator/方法自动跟着变（定制点，
> 见 [§5](#5-traits-特化给类型挂元数据)）；配合 [§4 CRTP](#4-crtpbasederived-与-derived) 用基类参数 + traits 返回类型，实现"运算逻辑与类型细节解耦"。

**典型出现处**：`spatial/motion-dense.hxx`、`force-dense.hxx` 的 operator 与 `SE3GroupAction`/`MotionAlgebraAction` 特化。

---

## 17. 用户自定义转换运算符 `operator T()`

**是什么**：一个特殊成员函数，定义"如何把**本对象**转换成**类型 T**"，让本类型能当作 T 直接使用。

```cpp
// inertia.hxx（InertiaBase）：让 Inertia 能当 6×6 矩阵用
operator Matrix6() const     // 转换运算符：目标类型 Matrix6 写在 operator 后，即返回类型（不另写）
{
  return matrix();           // 把紧凑的 (m,c,I_c) 展开成完整 6×6 空间惯量
}
```

**效果**：编译器在需要 `Matrix6` 处**自动插入转换**：

```cpp
Inertia Y = ...;
Matrix6 M = Y;                 // 隐式转换 = Y.matrix()
someFuncTakingMatrix6(Y);      // 传参自动转
```

**与转换构造函数的方向相反**：
- 转换**构造函数**（如 `MotionTpl(const MotionDense&)`）：`Target ← Source`（从别的类型**造出**本类型）；
- 转换**运算符** `operator T()`：`Source → Target`（把本类型**转成**别的类型）。二者是"进/出"一对。

**注意**：未加 `explicit` → **隐式转换**，可能悄悄发生；且每次都**物化一个 6×6**（`matrix()` 展开）有开销。
故 Pinocchio 内部高性能路径仍直接用 $(m,c,I_c)$ 分量算（见 [空间代数解析 §6](空间代数运算解析.md)），这个口子留给"确实需要显式矩阵"的场合（调试/外部线代库对接/公式验证）。可加 `explicit operator T()` 关掉隐式转换。

**典型出现处**：`spatial/inertia.hxx` 的 `operator Matrix6()`。

---

## 一句话定位

> 本文档是"招式speed本"：把 Pinocchio 里那些**第一次见会卡住**的 C++/Eigen/Boost 语法（编译期断言、依赖名消歧、
> CRTP、traits、Eigen 表达式惯用法、tag dispatch、模板的模板参数、`if_then_else`…）单独拎出来，配最小示例与出现处。
> 遇到看不懂的写法，先来这里对号入座。
