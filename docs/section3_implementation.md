---
title: "svanberg.cpp 与 athinput.figure2 逐段说明"
lang: zh-CN
fontsize: 10pt
geometry: [a4paper, margin=22mm]
mainfont: "Noto Serif CJK SC"
sansfont: "Noto Sans CJK SC"
monofont: "DejaVu Sans Mono"
colorlinks: true
header-includes:
  - |
    \usepackage{fvextra}
    \usepackage{mdframed}
    \usepackage{needspace}
    \definecolor{CodeBackground}{HTML}{EEF7FD}
    \definecolor{CodeBorder}{HTML}{D0E5F3}
    \definecolor{CodeText}{HTML}{20344B}
    \definecolor{CodeLine}{HTML}{6A8297}
    \mdfdefinestyle{codebox}{backgroundcolor=CodeBackground,
      linecolor=CodeBorder,linewidth=0.4pt,roundcorner=3pt,
      innerleftmargin=25pt,innerrightmargin=7pt,
      innertopmargin=7pt,innerbottommargin=7pt,
      skipabove=6pt,skipbelow=9pt}
    \BeforeBeginEnvironment{Verbatim}{\Needspace{7\baselineskip}\begin{mdframed}[style=codebox]}
    \AfterEndEnvironment{Verbatim}{\end{mdframed}}
    \AtBeginEnvironment{Verbatim}{\color{CodeText}}
    \renewcommand{\theFancyVerbLine}{\color{CodeLine}\scriptsize\arabic{FancyVerbLine}}
    \pretocmd{\subsection}{\Needspace{9\baselineskip}}{}{}
    \setlength{\parskip}{4pt}
    \XeTeXlinebreaklocale "zh"
    \XeTeXlinebreakskip=0pt plus 1pt
    \setlength{\emergencystretch}{3em}
---

# 一、src/svanberg.cpp

下面按原文件顺序收录代码，左侧数字是源文件行号。说明只针对这个 C++ 文件及后面的输入文件。当前输入采用球坐标、绝热能量方程和轴对称网格；代码中保留的柱坐标与其他边界分支会分别注明是否执行。

读数组时先区分：`u` 是守恒变量，保存密度、动量密度和总能量密度；`w`、`prim` 是原始变量，保存密度、速度和压强。同样的下标数字在两类数组中不代表同一种物理量。


## 1. 文件注释与头文件（第 1–36 行）

\begin{Verbatim}[numbers=left,firstnumber=1,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file disk.cpp
//! \brief Initializes stratified Keplerian accretion disk in both cylindrical and
//! spherical polar coordinates.  Initial conditions are in vertical hydrostatic eqm.

// C headers

// C++ headers
#include <algorithm>  // min
#include <cstdint>
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
\end{Verbatim}


第 1–8 行保留版权、BSD 许可和功能注释。第 6 行的 `disk.cpp` 是沿用的注释文字，实际文件名是 `svanberg.cpp`；注释不会决定编译的问题名称。第 9–12、24–25 行为空行或分类注释，不执行计算。

第 13–23 行是 C++ 标准库：

| 行 | 头文件 | 本文件中的作用 |
| --- | --- | --- |
| 13 | `algorithm` | 使用 `std::max` 限制密度下限；旁注 `min` 是遗留注释。 |
| 14 | `cstdint` | 提供 `std::uint64_t`，用于噪声散列。 |
| 15 | `cmath` | 平方根、幂、指数、对数、四舍五入及 `expm1`。 |
| 16 | `cstdlib` | 遗留包含项；当前噪声没有调用 `srand`。 |
| 17 | `cstring` | 用 `strcmp` 判断坐标系字符串。 |
| 18–19 | `fstream`、`iostream` | 保留的文件与流接口；当前没有直接用它们输出数据。 |
| 20 | `limits` | 获取浮点数最小正规格化正数，构造默认底限。 |
| 21–23 | `sstream`、`stdexcept`、`string` | 保留的字符串流、异常和字符串接口。 |

第 26–36 行接入 Athena++ 类型：`athena.hpp` 提供 `Real`、变量下标和宏；`athena_arrays.hpp` 提供数组；`bvals.hpp` 提供边界接口；`coordinates.hpp` 提供单元坐标；`eos.hpp` 提供状态方程接口；`field.hpp` 提供边界签名所需的 `FaceField`；`globals.hpp` 提供全局接口；`hydro.hpp` 访问流体数组；`mesh.hpp` 访问网格类；`orbital_advection.hpp` 访问可选轨道平流；`parameter_input.hpp` 读取输入。保留磁场类型接口不等于启用了磁场。



## 2. 局部函数、参数与 Noise（第 37–56 行）

\begin{Verbatim}[numbers=left,firstnumber=37,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
Real VelProfileCyl(const Real rad, const Real phi, const Real z);
// problem parameters which are useful to make global to this file
Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
Real dfloor;
Real Omega0, noise_amp, tau;
std::uint64_t seed;
// SplitMix64: coordinate-keyed noise independent of block decomposition.
Real Noise(Real r, Real theta, int component) {
  std::uint64_t x=seed ^ (std::uint64_t(std::llround(std::log(r)*1.e12))*0x9e3779b97f4a7c15ULL)
    ^ std::uint64_t(std::llround(theta*1.e12)) ^ (std::uint64_t(component+1)*0xbf58476d1ce4e5b9ULL);
  x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
  x=(x^(x>>27))*0x94d049bb133111ebULL; x^=x>>31;
  return 2.0*Real(x>>11)/9007199254740992.0-1.0;
}
} // namespace
\end{Verbatim}


第 38 行的匿名 `namespace` 将这些函数和变量限制在本源文件内部。第 39–42 行提前声明四个辅助函数，供后面的初始化和边界代码调用。

第 44–47 行声明共享参数：

| 变量 | 含义 |
| --- | --- |
| `gm0`, `r0`, `rho0` | $GM$、参考柱半径 $R_0$、参考中面密度 $\rho_0$。 |
| `dslope`, `pslope` | 密度和声速平方的径向幂指数，本次为 $-2$ 和 $-1$。 |
| `p0_over_r0` | 参考位置的 $P/\rho=c_0^2=0.0025$。 |
| `gamma_gas`, `dfloor` | 比热比与密度底限。 |
| `Omega0` | 参考系旋转角速度，默认 0；不是定义 $P_0$ 的开普勒频率。 |
| `noise_amp`, `tau`, `seed` | 噪声幅度、冷却时间、随机种子。 |

第 49–55 行 `Noise()` 为同一坐标和速度分量返回相同的伪随机值。第 50 行把 `log(r)` 乘 $10^{12}$ 并取整，编码径向坐标；第 51 行把极角和分量编号混入整数。`^` 是按位异或，`>>` 是右移，`ULL` 表示无符号长整型常量，均不是物理公式中的乘方。

第 52–53 行用异或、右移和固定乘数混合比特。第 54 行取高 53 位并除以 $2^{53}=9007199254740992$，先映射到 $[0,1)$，再映射到 $[-1,1)$。这些乘数用于散列，不是额外的物理参数。按坐标生成噪声，使初始化不依赖网格块的处理顺序。



## 3. 六个边界函数的声明（第 57–76 行）

\begin{Verbatim}[numbers=left,firstnumber=57,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
\end{Verbatim}


这里先声明函数，尚未设置任何边界值。`Inner/Outer` 表示内侧/外侧，`X1/X2/X3` 表示第 1/2/3 个坐标方向。球坐标下依次为 $r,\theta,\phi$。

各参数含义相同：`pmb` 是当前网格块，`pco` 是坐标对象，`prim` 是要填写的原始变量数组；`b` 为接口保留的磁场数据；`time`、`dt` 是边界调用时的时间和步长，本函数没有直接使用它们改变平衡值。`il/iu`、`jl/ju`、`kl/ku` 是相应索引范围，`ngh` 是幽灵层数。当前输入只有前四个边界函数会作为用户边界被调用，X3 使用周期条件。



## 4. InitUserMeshData：读取输入与注册边界（第 77–132 行）

\begin{Verbatim}[numbers=left,firstnumber=77,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));

  Omega0 = pin->GetOrAddReal("orbital_advection","Omega0",0.0);

  noise_amp=pin->GetOrAddReal("problem","noise_amp",0.01);
  tau=pin->GetOrAddReal("problem","tau",6.283185307179586e-10);
  seed=pin->GetOrAddInteger("problem","seed",202211);

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiskInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, DiskOuterX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiskInnerX3);
  }
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiskOuterX3);
  }

  return;
}
\end{Verbatim}


该函数在建立网格时执行，负责准备后续初始化、边界与冷却共用的参数。

- 第 87–88 行读取 `problem/GM` 和 `problem/r0`，默认分别为 0 和 1。
- 第 91 行 `GetReal` 要求输入必须提供 `rho0`；第 92 行 `GetOrAddReal` 读取 `dslope`，缺省时补入默认值 0。
- 第 95–98 行在绝热能量方程编译下读取声速平方归一化、温度幂指数以及 `hydro/gamma`。当前执行这一分支。
- 第 99–100 行是保留的等温读取分支。但末尾冷却函数直接操作能量，因此不能仅凭这个分支就认定整个文件支持直接改成等温编译。
- 第 102 行取得单精度最小正规格化正数；第 103 行据此给 `dfloor` 一个默认值。本输入已显式设为 $10^{-12}$，所以实际使用显式值。
- 第 105 行读取参考系的 `Omega0`，输入未给出，使用默认 0。
- 第 107–109 行读取 `noise_amp=0.01`、`tau=2\pi\times10^{-10}` 和 `seed=202211`。
- 第 112–129 行逐个检查六个面的边界标志。只有等于 `user` 才通过 `EnrollUserBoundaryFunction` 注册对应函数。当前 X1、X2 的四个面会注册，X3 的两个面因设置为 `periodic` 而不注册。
- 第 131 行结束函数，不创建流体初始场；初始场由下一个函数填写。



## 5. ProblemGenerator：逐单元建立初态（第 133–184 行）

\begin{Verbatim}[numbers=left,firstnumber=133,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel;
  Real x1, x2, x3;

  OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        den = DenProfileCyl(rad,phi,z);
        vel = VelProfileCyl(rad,phi,z);
        if (porb->orbital_advection_defined)
          vel -= vK(porb, x1, x2, x3);
        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM2,k,j,i) = den*vel;
          phydro->u(IM3,k,j,i) = 0.0;
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = den*vel;
        }

        // Uniform velocity perturbations in [-0.01 cs,+0.01 cs] in each component.
        Real cs=std::sqrt(PoverR(rad,phi,z));
        phydro->u(IM1,k,j,i)+=den*noise_amp*cs*Noise(x1,x2,0);
        phydro->u(IM2,k,j,i)+=den*noise_amp*cs*Noise(x1,x2,1);
        phydro->u(IM3,k,j,i)+=den*noise_amp*cs*Noise(x1,x2,2);
        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}
\end{Verbatim}


第 140–142 行声明局部坐标和临时密度、速度。`rad` 表示柱半径 $R$；球坐标运行中的 `x1` 是球半径 $r$，两者不同。第 144 行取得可选轨道速度函数的引用。

第 146–151 行遍历活动单元：`i` 对应 X1，`j` 对应 X2，`k` 对应 X3；`x1v/x2v/x3v` 读取体积平均单元中心。当前 `nx3=1`，但仍保留三重循环。

第 152 行转换到柱坐标。第 154–155 行分别计算式 (24) 的平衡密度和式 (25) 对应的旋转线速度。第 156–157 行只在启用轨道平流时减去轨道速度；当前默认未启用。

第 158 行保存密度，第 159 行把径向动量先设为 0。第 160–166 行按坐标系选择旋转分量：柱坐标时方位动量在 `IM2`；球坐标时在 `IM3`。当前执行第 164–165 行，得到 $m_\theta=0$、$m_\phi=\rho v_\phi$。

## 6. GetCylCoord：把网格坐标转换为 R、phi、z（第 185–201 行）

\begin{Verbatim}[numbers=left,firstnumber=185,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(k);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}
\end{Verbatim}


第 190 行的 `&rad`、`&phi`、`&z` 是引用输出：函数改变这些变量后，调用者能直接读到结果。

第 191–194 行用于柱坐标运行，直接读取 $R,\phi,z$。第 195–198 行用于当前球坐标运行：

$$
R=|r\sin\theta|,\qquad \phi=x_3,\qquad z=r\cos\theta.
$$

`abs` 保证柱半径非负；`cos(theta)` 保留上下半盘的高度符号。这里转换的是坐标，没有转换速度，也没有改变流体数组。



## 7. DenProfileCyl：式 (24) 的密度（第 202–214 行）

\begin{Verbatim}[numbers=left,firstnumber=202,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! computes density in cylindrical coordinates

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*std::pow(rad/r0,dslope);
  Real dentem = denmid*std::exp(gm0/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = dentem;
  return std::max(den,dfloor);
}
\end{Verbatim}


该函数逐项实现

$$
\rho=\rho_0(R/R_0)^{-q_D}
\exp\!\left[\frac{GM}{c_s^2(R)}\left(\frac1r-\frac1R\right)\right].
$$

- 第 207–208 行声明临时密度，并先赋参考声速平方。
- 第 209 行在当前绝热编译下调用 `PoverR()`，换成当前位置的 $c_s^2(R)$。
- 第 210 行是中面密度 $\rho_0(R/R_0)^{-q_D}$，因为 `dslope=-q_D=-2`。
- 第 211 行中 `sqrt(SQR(rad)+SQR(z))` 就是球半径 $r$；`gm0/p_over_r` 就是 $GM/c_s^2$；指数内部的差是 $1/r-1/R$。
- 第 212 行复制临时结果。第 213 行返回密度与 `dfloor` 的较大值；低于数值底限时才改变解析密度。

`phi` 保留在接口中，但这个轴对称密度解不依赖它。



## 8. PoverR 声速（第 215–223 行）

\begin{Verbatim}[numbers=left,firstnumber=215,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! computes pressure/density in cylindrical coordinates

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, pslope);
  return poverr;
}
\end{Verbatim}


第 220 行声明返回值。第 221 行计算

$$
\frac P\rho=c_s^2(R)
=c_0^2(R/R_0)^{-q_T}.
$$

本输入给出 `p0_over_r0=0.0025`、`r0=1`、`pslope=-1`，所以返回 $0.0025/R$。第 222 行返回该值。它不是声速；需要声速时再开平方。`phi` 和 `z` 没进入表达式，因此目标温度在给定柱半径上不随高度变化。



## 9. VelProfileCyl：式 (25) 对应的旋转线速度（第 224–235 行）

\begin{Verbatim}[numbers=left,firstnumber=224,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! computes rotational velocity in cylindrical coordinates

Real VelProfileCyl(const Real rad, const Real phi, const Real z) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+pslope)*p_over_r/(gm0/rad) + (1.0+pslope)
             - pslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(gm0/rad)*std::sqrt(vel) - rad*Omega0;
  return vel;
}
} // namespace
\end{Verbatim}


第 229 行读取局部 $c_s^2$。第 230–231 行把式 (25) 的无量纲括号存入临时变量 `vel`：

$$
B=-(q_D+q_T)h^2+(1-q_T)+q_T R/r.
$$

对应关系为 `dslope=-q_D`、`pslope=-q_T`，以及

$$
\frac{\mathrm{p\_over\_r}}{\mathrm{gm0}/\mathrm{rad}}
=\frac{c_s^2}{GM/R}=h^2.
$$

因此第 230 行的第一项是 $-(q_D+q_T)h^2$，第二项是 $1-q_T$；第 231 行是 $q_TR/r$。

第 232 行再计算 $\sqrt{GM/R}\sqrt B-R\Omega_{\rm frame}$。这里同名变量 `vel` 从无量纲括号变为速度，返回值为 $v_\phi=R\Omega$，不是角速度。当前 `Omega0=0`，故不作旋转参考系修正。取本输入后得到 $v_\phi=\sqrt{GM/R}\sqrt{-3h^2+R/r}$。高度依赖的 $R/r$ 给出了垂直剪切。



## 10. DiskInnerX1：径向内边界（第 236–281 行）

\begin{Verbatim}[numbers=left,firstnumber=236,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = vel;
          prim(IM3,k,j,il-i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = 0.0;
          prim(IM3,k,j,il-i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  }
}
\end{Verbatim}


这个函数填充内侧径向面的幽灵单元。当前输入把这个面设置为 `user`，因此会注册并调用该函数。

函数开头声明柱坐标临时量与旋转速度，并取得可选轨道速度函数。第 246–263 行是柱坐标分支；第 263–280 行是球坐标分支。当前球坐标编译选择后者；两个分支的主要区别是旋转速度存入哪个原始变量分量。

循环让 `i` 从 1 遍历到 `ngh`，定位 `il-i`，即由活动区边缘向外数的每一层幽灵单元。其余两个方向 `k`、`j` 遍历传入的横向索引范围。该索引表示幽灵位置，不是在重新初始化整个活动区域。

每个幽灵单元依次执行：

1. `GetCylCoord()` 取得该幽灵位置的 $R,\phi,z$，不会简单复制相邻活动单元的值。
2. `DenProfileCyl()` 按式 (24) 重新计算平衡密度。
3. `VelProfileCyl()` 按式 (25) 计算旋转线速度；若轨道平流启用，再减去对应轨道速度。当前这项附加分支不执行。
4. 柱坐标分支第 255/256/257 行依次设置 $(v_R,v_\phi,v_z)=(0,v_\phi,0)$。
5. 球坐标分支第 272/273/274 行依次设置 $(v_r,v_\theta,v_\phi)=(0,0,v_\phi)$。
6. 第 276 行在绝热能量编译下设置 $P=\rho c_s^2$，之后各层循环结束，函数返回。

这里 `prim` 是原始变量数组，速度无需乘密度；`prim(IEN)` 使用第 4 号压强槽位，因为 `IEN=IPR=4`。这与 `u(IEN)` 代表总能量的含义不同。函数没有重新加入随机噪声，也没有使用 `time`、`dt` 改变平衡解。



## 11. DiskOuterX1：径向外边界（第 282–327 行）

\begin{Verbatim}[numbers=left,firstnumber=282,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = vel;
          prim(IM3,k,j,iu+i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = 0.0;
          prim(IM3,k,j,iu+i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  }
}
\end{Verbatim}


这个函数填充外侧径向面的幽灵单元。当前输入把这个面设置为 `user`，因此会注册并调用该函数。

函数开头声明柱坐标临时量与旋转速度，并取得可选轨道速度函数。第 292–309 行是柱坐标分支；第 309–326 行是球坐标分支。当前球坐标编译选择后者；两个分支的主要区别是旋转速度存入哪个原始变量分量。

循环让 `i` 从 1 遍历到 `ngh`，定位 `iu+i`，即由活动区边缘向外数的每一层幽灵单元。其余两个方向 `k`、`j` 遍历传入的横向索引范围。该索引表示幽灵位置，不是在重新初始化整个活动区域。

每个幽灵单元依次执行：

1. `GetCylCoord()` 取得该幽灵位置的 $R,\phi,z$，不会简单复制相邻活动单元的值。
2. `DenProfileCyl()` 按式 (24) 重新计算平衡密度。
3. `VelProfileCyl()` 按式 (25) 计算旋转线速度；若轨道平流启用，再减去对应轨道速度。当前这项附加分支不执行。
4. 柱坐标分支第 301/302/303 行依次设置 $(v_R,v_\phi,v_z)=(0,v_\phi,0)$。
5. 球坐标分支第 318/319/320 行依次设置 $(v_r,v_\theta,v_\phi)=(0,0,v_\phi)$。
6. 第 322 行在绝热能量编译下设置 $P=\rho c_s^2$，之后各层循环结束，函数返回。

这里 `prim` 是原始变量数组，速度无需乘密度；`prim(IEN)` 使用第 4 号压强槽位，因为 `IEN=IPR=4`。这与 `u(IEN)` 代表总能量的含义不同。函数没有重新加入随机噪声，也没有使用 `time`、`dt` 改变平衡解。



## 12. DiskInnerX2：极角内边界（第 328–373 行）

\begin{Verbatim}[numbers=left,firstnumber=328,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = vel;
          prim(IM3,k,jl-j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = 0.0;
          prim(IM3,k,jl-j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  }
}
\end{Verbatim}


这个函数填充内侧极角面的幽灵单元。当前输入把这个面设置为 `user`，因此会注册并调用该函数。

函数开头声明柱坐标临时量与旋转速度，并取得可选轨道速度函数。第 338–355 行是柱坐标分支；第 355–372 行是球坐标分支。当前球坐标编译选择后者；两个分支的主要区别是旋转速度存入哪个原始变量分量。

循环让 `j` 从 1 遍历到 `ngh`，定位 `jl-j`，即由活动区边缘向外数的每一层幽灵单元。其余两个方向 `k`、`i` 遍历传入的横向索引范围。该索引表示幽灵位置，不是在重新初始化整个活动区域。

每个幽灵单元依次执行：

1. `GetCylCoord()` 取得该幽灵位置的 $R,\phi,z$，不会简单复制相邻活动单元的值。
2. `DenProfileCyl()` 按式 (24) 重新计算平衡密度。
3. `VelProfileCyl()` 按式 (25) 计算旋转线速度；若轨道平流启用，再减去对应轨道速度。当前这项附加分支不执行。
4. 柱坐标分支第 347/348/349 行依次设置 $(v_R,v_\phi,v_z)=(0,v_\phi,0)$。
5. 球坐标分支第 364/365/366 行依次设置 $(v_r,v_\theta,v_\phi)=(0,0,v_\phi)$。
6. 第 368 行在绝热能量编译下设置 $P=\rho c_s^2$，之后各层循环结束，函数返回。

这里 `prim` 是原始变量数组，速度无需乘密度；`prim(IEN)` 使用第 4 号压强槽位，因为 `IEN=IPR=4`。这与 `u(IEN)` 代表总能量的含义不同。函数没有重新加入随机噪声，也没有使用 `time`、`dt` 改变平衡解。



## 13. DiskOuterX2：极角外边界（第 374–419 行）

\begin{Verbatim}[numbers=left,firstnumber=374,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = vel;
          prim(IM3,k,ju+j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = 0.0;
          prim(IM3,k,ju+j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  }
}
\end{Verbatim}


这个函数填充外侧极角面的幽灵单元。当前输入把这个面设置为 `user`，因此会注册并调用该函数。

函数开头声明柱坐标临时量与旋转速度，并取得可选轨道速度函数。第 384–401 行是柱坐标分支；第 401–418 行是球坐标分支。当前球坐标编译选择后者；两个分支的主要区别是旋转速度存入哪个原始变量分量。

循环让 `j` 从 1 遍历到 `ngh`，定位 `ju+j`，即由活动区边缘向外数的每一层幽灵单元。其余两个方向 `k`、`i` 遍历传入的横向索引范围。该索引表示幽灵位置，不是在重新初始化整个活动区域。

每个幽灵单元依次执行：

1. `GetCylCoord()` 取得该幽灵位置的 $R,\phi,z$，不会简单复制相邻活动单元的值。
2. `DenProfileCyl()` 按式 (24) 重新计算平衡密度。
3. `VelProfileCyl()` 按式 (25) 计算旋转线速度；若轨道平流启用，再减去对应轨道速度。当前这项附加分支不执行。
4. 柱坐标分支第 393/394/395 行依次设置 $(v_R,v_\phi,v_z)=(0,v_\phi,0)$。
5. 球坐标分支第 410/411/412 行依次设置 $(v_r,v_\theta,v_\phi)=(0,0,v_\phi)$。
6. 第 414 行在绝热能量编译下设置 $P=\rho c_s^2$，之后各层循环结束，函数返回。

这里 `prim` 是原始变量数组，速度无需乘密度；`prim(IEN)` 使用第 4 号压强槽位，因为 `IEN=IPR=4`。这与 `u(IEN)` 代表总能量的含义不同。函数没有重新加入随机噪声，也没有使用 `time`、`dt` 改变平衡解。



## 14. DiskInnerX3：方位角内边界（第 420–465 行）

\begin{Verbatim}[numbers=left,firstnumber=420,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = vel;
          prim(IM3,kl-k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = 0.0;
          prim(IM3,kl-k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  }
}
\end{Verbatim}


这个函数填充内侧方位角面的幽灵单元。当前输入将 X3 设置为 `periodic`，所以这个用户边界函数不会被注册或调用。函数仍保留在文件中，只有把对应面显式设为 `user` 才使用。

函数开头声明柱坐标临时量与旋转速度，并取得可选轨道速度函数。第 430–447 行是柱坐标分支；第 447–464 行是球坐标分支。当前球坐标编译选择后者；两个分支的主要区别是旋转速度存入哪个原始变量分量。

循环让 `k` 从 1 遍历到 `ngh`，定位 `kl-k`，即由活动区边缘向外数的每一层幽灵单元。其余两个方向 `j`、`i` 遍历传入的横向索引范围。该索引表示幽灵位置，不是在重新初始化整个活动区域。

每个幽灵单元依次执行：

1. `GetCylCoord()` 取得该幽灵位置的 $R,\phi,z$，不会简单复制相邻活动单元的值。
2. `DenProfileCyl()` 按式 (24) 重新计算平衡密度。
3. `VelProfileCyl()` 按式 (25) 计算旋转线速度；若轨道平流启用，再减去对应轨道速度。当前这项附加分支不执行。
4. 柱坐标分支第 439/440/441 行依次设置 $(v_R,v_\phi,v_z)=(0,v_\phi,0)$。
5. 球坐标分支第 456/457/458 行依次设置 $(v_r,v_\theta,v_\phi)=(0,0,v_\phi)$。
6. 第 460 行在绝热能量编译下设置 $P=\rho c_s^2$，之后各层循环结束，函数返回。

这里 `prim` 是原始变量数组，速度无需乘密度；`prim(IEN)` 使用第 4 号压强槽位，因为 `IEN=IPR=4`。这与 `u(IEN)` 代表总能量的含义不同。函数没有重新加入随机噪声，也没有使用 `time`、`dt` 改变平衡解。



## 15. DiskOuterX3：方位角外边界（第 466–511 行）

\begin{Verbatim}[numbers=left,firstnumber=466,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = vel;
          prim(IM3,ku+k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = 0.0;
          prim(IM3,ku+k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  }
}
\end{Verbatim}


这个函数填充外侧方位角面的幽灵单元。当前输入将 X3 设置为 `periodic`，所以这个用户边界函数不会被注册或调用。函数仍保留在文件中，只有把对应面显式设为 `user` 才使用。

函数开头声明柱坐标临时量与旋转速度，并取得可选轨道速度函数。第 476–493 行是柱坐标分支；第 493–510 行是球坐标分支。当前球坐标编译选择后者；两个分支的主要区别是旋转速度存入哪个原始变量分量。

循环让 `k` 从 1 遍历到 `ngh`，定位 `ku+k`，即由活动区边缘向外数的每一层幽灵单元。其余两个方向 `j`、`i` 遍历传入的横向索引范围。该索引表示幽灵位置，不是在重新初始化整个活动区域。

每个幽灵单元依次执行：

1. `GetCylCoord()` 取得该幽灵位置的 $R,\phi,z$，不会简单复制相邻活动单元的值。
2. `DenProfileCyl()` 按式 (24) 重新计算平衡密度。
3. `VelProfileCyl()` 按式 (25) 计算旋转线速度；若轨道平流启用，再减去对应轨道速度。当前这项附加分支不执行。
4. 柱坐标分支第 485/486/487 行依次设置 $(v_R,v_\phi,v_z)=(0,v_\phi,0)$。
5. 球坐标分支第 502/503/504 行依次设置 $(v_r,v_\theta,v_\phi)=(0,0,v_\phi)$。
6. 第 506 行在绝热能量编译下设置 $P=\rho c_s^2$，之后各层循环结束，函数返回。

这里 `prim` 是原始变量数组，速度无需乘密度；`prim(IEN)` 使用第 4 号压强槽位，因为 `IEN=IPR=4`。这与 `u(IEN)` 代表总能量的含义不同。函数没有重新加入随机噪声，也没有使用 `time`、`dt` 改变平衡解。



## 16. UserWorkInLoop：式 (27) 的逐步热松弛（第 512–528 行）

\begin{Verbatim}[numbers=left,firstnumber=512,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]

// Eq. (27), exact exponential relaxation once per completed timestep.
// Update both U and W because Athena computes the next timestep from W.
void MeshBlock::UserWorkInLoop() {
  const Real f=-std::expm1(-pmy_mesh->dt/tau);
  for(int k=ks;k<=ke;++k) for(int j=js;j<=je;++j) for(int i=is;i<=ie;++i) {
    const Real rho=phydro->u(IDN,k,j,i);
    const Real kinetic=0.5*(SQR(phydro->u(IM1,k,j,i))+
      SQR(phydro->u(IM2,k,j,i))+SQR(phydro->u(IM3,k,j,i)))/rho;
    Real R,phi,z; GetCylCoord(pcoord,R,phi,z,i,j,k);
    const Real target=rho*PoverR(R,phi,z)/(gamma_gas-1.0);
    const Real internal=phydro->u(IEN,k,j,i)-kinetic;
    const Real cooled=internal+f*(target-internal);
    phydro->u(IEN,k,j,i)=kinetic+cooled;
    phydro->w(IPR,k,j,i)=(gamma_gas-1.0)*cooled;
  }
}
\end{Verbatim}


该函数在每个完整流体时间步之后执行。第 513–514 行说明作用：实施冷却，同时保持两类热变量一致。

第 516 行计算 $f=1-\exp(-\Delta t/\tau)$。`expm1(x)` 返回 $e^x-1$，所以前面的负号得到论文式 (27) 的正松弛系数。`pmy_mesh->dt` 是当前完整步长，不是代码自行设定的额外时间。

第 517 行遍历活动单元；第 518 行从守恒数组读取密度。第 519–520 行使用三个动量密度计算动能密度；代入 `u(IM*)` 保存的 $\rho v$，可写成：

$$
K=\frac12\rho\left(v_r^2+v_\theta^2+v_\phi^2\right).
$$

第 521 行转换单元坐标。第 522 行用当前密度计算目标内能密度：

$$
\epsilon_0=\frac{\rho c_s^2(R)}{\gamma-1}.
$$

第 523 行由 $E-K$ 取得当前内能。第 524 行进行热松弛：

$$
\epsilon_{\rm new}=\epsilon+
(\epsilon_0-\epsilon)[1-\exp(-\Delta t/\tau)].
$$

在这个子步骤中 $\rho$ 不变且 $\epsilon\propto\rho T$，因此等价于论文对温度的更新。第 525 行把原动能加回总能量，第 526 行更新原始压强 $P=(\gamma-1)\epsilon_{\rm new}$；第 527–528 行结束循环与函数。

这个函数不改变密度或动量，也不会把总能量直接覆盖成内能。对于当前极小的 `tau`，通常 $f$ 在数值上接近 1，温度迅速回到局部目标。同步更新 `w(IPR)` 是为了后续读取原始变量的时间步估计使用更新后的压强。



\newpage

# 二、inputs/athinput.figure2

下面按输入文件的原顺序解释每一项。尖括号包围的块标题（例如 `<mesh>`）开始一个参数块，直到下一个块标题；同名参数在不同块中的含义可以不同，例如 `<mesh>` 的 `nx1` 是总网格数，`<meshblock>` 的 `nx1` 是单块网格数。下面所有时间值均为代码单位，$GM=R_0=1$ 时 $P_0=2\pi$。


## 1. comment 和 job（第 1–8 行）

\begin{Verbatim}[numbers=left,firstnumber=1,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<comment>
problem = Svanberg et al. 2022, Section 3, Figure 2


<job>
problem_id = vsi


\end{Verbatim}


- 第 1 行 `<comment>`：开始描述块。
- 第 2 行 `problem=...`：问题说明文字，不参与动力学计算，也不会自动载入论文参数。
- 第 5 行 `<job>`：开始任务命名块。
- 第 6 行 `problem_id=vsi`：输出文件的基础名称为 `vsi`，不改变物理初态。



## 2. time：步长控制与终止时间（第 9–17 行）

\begin{Verbatim}[numbers=left,firstnumber=9,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<time>
cfl_number = 0.3
nlim = -1
tlim = 3141.592653589793
integrator = vl2
xorder = 2
ncycle_out = 100


\end{Verbatim}


- 第 9 行 `<time>`：开始时间控制块。
- 第 10 行 `cfl_number=0.3`：控制流体 CFL 步长的安全系数，不表示固定时间步长为 0.3。
- 第 11 行 `nlim=-1`：不设置最大步数限制。
- 第 12 行 `tlim=3141.592653589793`：$500\times2\pi=500P_0$，即 Figure 2 的目标时刻。论文整段 $3000P_0$ 对应 $18849.55592153876$，不是这里当前的值。
- 第 13 行 `integrator=vl2`：采用两阶段、二阶时间精度的 van Leer 时间积分。
- 第 14 行 `xorder=2`：采用二阶空间重构，在这里对应 PLM。
- 第 15 行 `ncycle_out=100`：每 100 个流体步在标准输出打印一次循环信息；不是每 100 步保存一次完整快照。



## 3. mesh：整个计算域（第 18–37 行）

\begin{Verbatim}[numbers=left,firstnumber=18,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<mesh>
nx1 = 704
x1min = 1.0
x1max = 30.0
x1rat = 1.005
ix1_bc = user
ox1_bc = user
nx2 = 448
x2min = 1.3207963267948966
x2max = 1.8207963267948966
ix2_bc = user
ox2_bc = user
nx3 = 1
x3min = 0.0
x3max = 6.283185307179586
ix3_bc = periodic
ox3_bc = periodic
num_threads = 4


\end{Verbatim}


- 第 18 行 `<mesh>`：整个网格的参数块。
- 第 19 行 `nx1=704`：径向活动单元总数为 704。
- 第 20 行 `x1min=1.0`：球半径内边界 $r=1$。
- 第 21 行 `x1max=30.0`：球半径外边界 $r=30$。
- 第 22 行 `x1rat=1.005`：相邻径向单元宽度之比为 1.005；端点保持为 1 和 30。对应单元面可写为 $r_i=1+29(q^i-1)/(q^{704}-1)$，$q=1.005$。
- 第 23 行 `ix1_bc=user`：径向内边界调用 `DiskInnerX1()`。
- 第 24 行 `ox1_bc=user`：径向外边界调用 `DiskOuterX1()`。
- 第 25 行 `nx2=448`：极角方向活动单元总数为 448。
- 第 26 行 `x2min=1.3207963267948966`：极角内边界 $\pi/2-5h$，其中 $h=0.05$。
- 第 27 行 `x2max=1.8207963267948966`：极角外边界 $\pi/2+5h$。这些角度使用弧度。
- 第 28 行 `ix2_bc=user`：极角内边界调用 `DiskInnerX2()`。
- 第 29 行 `ox2_bc=user`：极角外边界调用 `DiskOuterX2()`。
- 第 30 行 `nx3=1`：只有一个方位角单元，不解析方位角变化；仍然演化方位速度与动量。
- 第 31 行 `x3min=0.0`：方位角坐标起点为 0。
- 第 32 行 `x3max=6.283185307179586`：方位角坐标终点为 $2\pi$。
- 第 33 行 `ix3_bc=periodic`：方位角内边界为周期条件，不调用 `DiskInnerX3()`。
- 第 34 行 `ox3_bc=periodic`：方位角外边界为周期条件，不调用 `DiskOuterX3()`。
- 第 35 行 `num_threads=4`：设置网格处理使用的 OpenMP 线程数为 4；需配合启用 OpenMP 的可执行程序。

输入未写 `x2rat`，因此采用默认值 1，极角均匀划分。极角跨度 $0.5=10h$，448 格在中面附近对应约 44.8 格/H，与论文文字中的约 12 格/H 不一致；这里实际执行的是明确列出的范围和 448 格。也没有另外设置“10 格/H”的参数来覆盖径向网格。



## 4. meshblock：单个网格块（第 38–43 行）

\begin{Verbatim}[numbers=left,firstnumber=38,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<meshblock>
nx1 = 176
nx2 = 448
nx3 = 1


\end{Verbatim}


- 第 38 行 `<meshblock>`：开始单块网格配置。
- 第 39 行 `nx1=176`：每块有 176 个径向活动单元。
- 第 40 行 `nx2=448`：每块覆盖全部 448 个极角活动单元。
- 第 41 行 `nx3=1`：每块只有一个方位角活动单元。

总块数为 $(704/176)(448/448)(1/1)=4$。网格分块用于计算组织，没有改变总分辨率；幽灵单元不计入上述活动单元数。



## 5. hydro：状态方程与底限（第 44–49 行）

\begin{Verbatim}[numbers=left,firstnumber=44,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<hydro>
gamma = 1.4
dfloor = 1.e-12
pfloor = 1.e-16


\end{Verbatim}


- 第 44 行 `<hydro>`：开始流体参数块。
- 第 45 行 `gamma=1.4`：$\gamma=7/5$，用于 $P=(\gamma-1)\epsilon$。`InitUserMeshData()` 读取为 `gamma_gas`，初始化能量和冷却都使用它。
- 第 46 行 `dfloor=1.e-12`：密度数值底限。`DenProfileCyl()` 返回值也通过该底限截断。
- 第 47 行 `pfloor=1.e-16`：压强数值底限，由流体状态处理使用；本文件没有把它读成另一个独立全局变量。

底限用于数值保护，不是额外的密度或温度幂律。



## 6. problem：盘模型、噪声与冷却（第 50–61 行）

\begin{Verbatim}[numbers=left,firstnumber=50,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<problem>
GM = 1.0
r0 = 1.0
rho0 = 1.0
p0_over_r0 = 0.0025
pslope = -1.0
dslope = -2.0
noise_amp = 0.01
seed = 202211
tau = 6.283185307179586e-10


\end{Verbatim}


- 第 50 行 `<problem>`：开始问题生成器参数块。
- 第 51 行 `GM=1.0`：中心引力参数。既用于本文件的解析平衡，也启用内置点质量引力；本文件没有再手写一次重复引力。
- 第 52 行 `r0=1.0`：柱半径参考位置 $R_0=1$。
- 第 53 行 `rho0=1.0`：$R=R_0$ 的中面密度归一化为 1。
- 第 54 行 `p0_over_r0=0.0025`：目标 $P/\rho$ 在参考半径的值，即 $c_0^2$；对应 $c_0=0.05$，不是 $c_0=0.0025$。
- 第 55 行 `pslope=-1.0`：$c_s^2\propto R^{-1}$，对应论文 $q_T=1$。
- 第 56 行 `dslope=-2.0`：中面密度 $\rho_{\rm mid}\propto R^{-2}$，对应论文 $q_D=2$。
- 第 57 行 `noise_amp=0.01`：每个速度分量噪声的最大绝对幅度为局部声速的 1%，由第 170–172 行乘入动量增量。
- 第 58 行 `seed=202211`：散列噪声的种子，决定可重复的初始扰动。
- 第 59 行 `tau=6.283185307179586e-10`：$10^{-10}P_0$ 换算后的代码冷却时间，不是 $10^{-10}$ 代码时间。

这些参数给出 $h^2=c_s^2/(GM/R)=0.0025$，故 $h=0.05$。输入没有单独的 `h` 字段；它已由声速、引力与径向幂律组合确定。



## 7. output1：完整流场快照（第 62–67 行）

\begin{Verbatim}[numbers=left,firstnumber=62,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<output1>
file_type = vtk
variable = prim
dt = 6.283185307179586


\end{Verbatim}


- 第 62 行 `<output1>`：第一组输出配置。
- 第 63 行 `file_type=vtk`：以 VTK 格式保存空间场。
- 第 64 行 `variable=prim`：保存原始流体变量，包括密度、压强与速度；不是只保存二维图片。
- 第 65 行 `dt=6.283185307179586`：输出间隔 $2\pi=1P_0$。

这里的 `dt` 是该输出块的采样间隔，不是 `<time>` 控制的流体积分步长。常规输出在首个达到或跨过计划时刻的步写出，所以准确输出时间需读取文件头。



## 8. output2：检查点（第 68–72 行）

\begin{Verbatim}[numbers=left,firstnumber=68,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<output2>
file_type = rst
dt = 62.83185307179586


\end{Verbatim}


- 第 68 行 `<output2>`：第二组输出配置。
- 第 69 行 `file_type=rst`：保存可供重启的检查点，而非可视化快照。
- 第 70 行 `dt=62.83185307179586`：间隔 $20\pi=10P_0$。

这与 VTK 的用途不同：检查点保存继续积分需要的状态。



## 9. output3：history（第 73–75 行）

\begin{Verbatim}[numbers=left,firstnumber=73,numbersep=6pt,fontsize=\footnotesize,breaklines=true,breakanywhere=true]
<output3>
file_type = hst
dt = 0.6283185307179586
\end{Verbatim}


- 第 73 行 `<output3>`：第三组输出配置。
- 第 74 行 `file_type=hst`：保存 history 时间序列，例如积分量。
- 第 75 行 `dt=0.6283185307179586`：间隔 $0.2\pi=0.1P_0$。

history 不保存每个空间单元的完整流场，不能代替 VTK 快照。
