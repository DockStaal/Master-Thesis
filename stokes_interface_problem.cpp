/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "finiteElement.hpp"
#include "baseProblem.hpp"
#include "paraview.hpp"
#include "../tool.hpp"
#include "../num/matlab.hpp"
#include "../external/matplotlibcpp.h"

// #define MACRO
// #define NON_SYMMETRIC

#define BDM1_ELEMENTS
// #define RT0_ELEMENTS
// #define RT1_ELEMENTS

#define STATIC_DROP
// #define DISCONTINUOU_CONTINOUS_PROBLEM
// #define CONTINUOUS_PROBLEM
// #define DYNAMIC_DROP
// #define DISCONTINUOUS_PROBLEM
// #define VORTEX

std::string macro_string = "nm";

#ifdef NON_SYMMETRIC
    const std::string non_symmetric_string = "ns";
#else
    const std::string non_symmetric_string = "s";
#endif

#ifdef BDM1_ELEMENTS
    const std::string elements = "bdm1";
#endif

#ifdef RT0_ELEMENTS
    const std::string elements = "rt0";
#endif

#ifdef RT1_ELEMENTS
    const std::string elements = "rt1";
#endif


#ifdef STATIC_DROP

std::string problemName = "sd";

// Options mesh and domain
int nx = 30; // Number of elements in x direction
int ny = 30; // Number of elements in y direction

double bottom_left_x = -1; // Bottom left corner x coordinate
double bottom_left_y = -1; // Bottom left corner y coordinate

double width_domain = 2; // Width of the domain
double height_domain = 2; // Height of the domain

double Rad = 0.6; // Radius of the drop
// Paramers in problem
const R sigma = 5;
double interface_curvature = 1 / Rad;

double mu1 = 1;
double mu2 = 1;

double kappa1 = 1. / 2.;
double kappa2 = 1. / 2.;


// External force function in the momentum conservation equation
//
// @param P:         2D coordinate
// @param component: First or second component in the momentum conservation
//                   equation. i.e. f_1(x, y) or f_2(x, y).
// @param domain:    First or second domain i.e. \Omega_1 or \Omega_2.
//
// @return:          Value of the external force function f at coordinate P
//                   in component 1 or 2.
R rhs(double *P, int component, int domain) {
   // R mu=1;
   R x = P[0];
   R y = P[1];
   if (component == 0)
      return 0;
   else if (component == 1)
      return 0;
   else
      return 0;
}

// Dirichlet boundary condition for the velocity.
//
// @param P:         2D coordinate
// @param component: First or second component of the velocity.
//
// @return:          Value of the Dirichlet boundary condition at coordinate P
//                   in component 1 or 2.
static R u_exact(double* P, int component, int domain) {
    R x = P[0];
    R y = P[1];
    return 0;
}

static R p_exact(double* P, int component, int domain) {
    R x = P[0];
    R y = P[1];
    return (domain == 0)? 0 : sigma / (Rad);
}

// Levelset function that describes the interface \Gamma.
//
// @param P:         2D coordinate
//
// @return:          Value of the levelset function at coordinate P.
//                   The levelset describes a circle centered at (0, 0) with
//                   radius 0.5;
R Gamma(double* P, int i) {
    R2 shift(0., 0.);
   return sqrt((P[0] - shift.x) * (P[0] - shift.x) +
               (P[1] - shift.y) * (P[1] - shift.y)) - Rad;
}

#endif

#ifdef CONTINUOUS_PROBLEM

std::string problemName = "cp";

// Options mesh and domain
int nx = 100; // Number of elements in x direction
int ny = 100; // Number of elements in y direction

double bottom_left_x = 0; // Bottom left corner x coordinate
double bottom_left_y = 0; // Bottom left corner y coordinate

double width_domain = 1; // Width of the domain
double height_domain = 1; // Height of the domain

// Paramers in problem
const R sigma = 0;
double interface_curvature = 0;

double kappa1 = 0.5;
double kappa2 = 0.5;

double mu1 = 1;
double mu2 = 1;

// External force function in the momentum conservation equation
//
// @param P:         2D coordinate
// @param component: First or second component in the momentum conservation
//                   equation. i.e. f_1(x, y) or f_2(x, y).
// @param domain:    First or second domain i.e. \Omega_1 or \Omega_2.
//
// @return:          Value of the external force function f at coordinate P
//                   in component 1 or 2.
R rhs(double *P, int component, int domain) {
   // R mu=1;
   R x = P[0];
   R y = P[1];
   if (component == 0)
      return 0;
   else if (component == 1)
      return 0;
   else
      return 0;
}

// Dirichlet boundary condition for the velocity.
//
// @param P:         2D coordinate
// @param component: First or second component of the velocity.
//
// @return:          Value of the Dirichlet boundary condition at coordinate P
//                   in component 1 or 2.
R u_exact(double* P, int component, int domain) {
    R x = P[0];
    R y = P[1];
    if (component == 0)
      return 20 * x * y * y * y;
    else
        return 5 * x * x * x * x - 5 * y * y * y * y;
}

R p_exact(double* P, int component, int domain) {
    R x = P[0];
    R y = P[1];
    return 60 * x * x * y - 20 * y * y * y - 5;
}

// Levelset function that describes the interface \Gamma.
//
// @param P:         2D coordinate
//
// @return:          Value of the levelset function at coordinate P.
//                   The levelset describes a circle centered at (0.5, 0.5) with
//                   radius 0.3
R Gamma(double* P, int i) {
    R2 shift(0.5, 0.5);
    return sqrt((P[0] - shift.x) * (P[0] - shift.x) +
               (P[1] - shift.y) * (P[1] - shift.y)) - 0.3;
}

#endif

#ifdef DISCONTINUOU_CONTINOUS_PROBLEM

std::string problemName = "Discontinuous-continuous problem";

// Options mesh and domain
int nx = 100; // Number of elements in x direction
int ny = 100; // Number of elements in y direction

double bottom_left_x = 0; // Bottom left corner x coordinate
double bottom_left_y = 0; // Bottom left corner y coordinate

double width_domain = 1; // Width of the domain
double height_domain = 1; // Height of the domain

// Paramers in problem
const R sigma = 0.00001 * 10;
double interface_curvature = 1. / 0.3;

double kappa1 = 0.5;
double kappa2 = 0.5;

double mu1 = 1;
double mu2 = 1;

// External force function in the momentum conservation equation
//
// @param P:         2D coordinate
// @param component: First or second component in the momentum conservation
//                   equation. i.e. f_1(x, y) or f_2(x, y).
// @param domain:    First or second domain i.e. \Omega_1 or \Omega_2.
//
// @return:          Value of the external force function f at coordinate P
//                   in component 1 or 2.
R rhs(double *P, int component, int domain) {
   // R mu=1;
   R x = P[0];
   R y = P[1];
   if (component == 0)
      return 0;
   else if (component == 1)
      return 0;
   else
      return 0;
}

// Dirichlet boundary condition for the velocity.
//
// @param P:         2D coordinate
// @param component: First or second component of the velocity.
//
// @return:          Value of the Dirichlet boundary condition at coordinate P
//                   in component 1 or 2.
R u_exact(double* P, int component, int domain) {
    R x = P[0];
    R y = P[1];
    if (component == 0)
      return 20 * x * y * y * y;
    else
        return 5 * x * x * x * x - 5 * y * y * y * y;
}

R p_exact(double* P, int component, int domain) {
    R x = P[0];
    R y = P[1];
    if (domain == 0)
        return 60 * x * x * y - 20 * y * y * y;
    else
        return 60 * x * x * y - 20 * y * y * y + sigma * interface_curvature;
}

// Levelset function that describes the interface \Gamma.
//
// @param P:         2D coordinate
//
// @return:          Value of the levelset function at coordinate P.
//                   The levelset describes a circle centered at (0.5, 0.5) with
//                   radius 0.3
R Gamma(double* P, int i) {
    R2 shift(0.5, 0.5);
    return sqrt((P[0] - shift.x) * (P[0] - shift.x) +
               (P[1] - shift.y) * (P[1] - shift.y)) - 0.3;
}

#endif

#ifdef DYNAMIC_DROP

std::string problemName = "dd";

int nx = 30; // Number of elements in x direction
int ny = 30; // Number of elements in y direction

double bottom_left_x = -1.; // Bottom left corner x coordinate
double bottom_left_y = -1.; // Bottom left corner y coordinate

double width_domain = 2; // Width of the domain
double height_domain = 2; // Height of the domain

// Paramers in problem
const R sigma = 10.;
double interface_curvature = 3./ 2.;

double kappa1 = 0.5;
double kappa2 = 0.5;

double mu1 = 100.;
double mu2 = 1.;

R2 shift(0, 0);

R Gamma(double* P, int i) {
   return sqrt((P[0] - shift.x) * (P[0] - shift.x) +
               (P[1] - shift.y) * (P[1] - shift.y)) -
          2. / 3.;
}

static R falpha1(const R r) {
   // const R MU1 = 1e-1;
   // const R MU2 = 1e-3;
   const R MU1 = mu1;
   const R MU2 = mu2;
   const R r0  = 2. / 3;
   return 1. / MU1 + (1. / MU2 - 1. / MU1) * exp(r * r - r0 * r0);
}
static R falpha2(const R r) {
   R MU2 = mu2; // 1e-3;
   return 1. / MU2;
}

static R falpha(const R r) {
   const R r0 = 2. / 3;
   return (r < r0) ? falpha2(r) : falpha1(r);
}

static R rhs(double* P, int i) {
   const R r2 = Norme2_2(R2(P[0], P[1]));
   const R s  = exp(-r2);
   R2 R(4 * s * (r2 - 2) * P[1] + 3 * P[0] * P[0], -4 * s * (r2 - 2) * P[0]);
   return (i < 2) ? R[i] : 0;
}

static R2 fun_velocity1(double* P) {
   R r = Norme2(R2(P[0], P[1]));
   R2 R(-P[1], P[0]);
   R = falpha1(r) * exp(-r * r) * R;
   return R;
}
static R2 fun_velocity2(double* P) {
   R r = Norme2(R2(P[0], P[1]));
   R2 R(-P[1], P[0]);
   R = falpha2(r) * exp(-r * r) * R;
   return R;
}
static R fun_pressure1(double* P) { return pow(P[0], 3); }
static R fun_pressure2(double* P) {
   // R sigma = 1;//24.5;//700;
   return pow(P[0], 3) + sigma * 3. / 2.;
}

static R u_exact(double* P, int ci, int domain) {
   if (domain == 0)
      return fun_velocity1(P)[ci];
   else
      return fun_velocity2(P)[ci];
}
static R p_exact(double* P, int ci, int domain) {
   if (domain == 0)
      return fun_pressure1(P);
   else
      return fun_pressure2(P);
}

R2 fparam(double t) {
   return R2(2. / 3 * cos(t + 1. / 3), 2. / 3 * sin(t + 1. / 3));
}

#endif

#ifdef VORTEX

std::string problemName = "v";

int nx = 60; // Number of elements in x direction
int ny = 60; // Number of elements in y direction

double bottom_left_x = -1.; // Bottom left corner x coordinate
double bottom_left_y = -1.; // Bottom left corner y coordinate

double width_domain = 2.; // Width of the domain
double height_domain = 2.; // Height of the domain

// Paramers in problem
const R sigma = 2;
double interface_curvature = 3. / 2.;

double mu1 = 1.;
double mu2 = 1.;

double kappa1 = 1. / 2.;
double kappa2 = 1. / 2.;

double Re = 10;

R2 shift(0, 0);

R Gamma(double* P, int i) {
   return sqrt((P[0] - shift.x) * (P[0] - shift.x) +
               (P[1] - shift.y) * (P[1] - shift.y)) -
          2. / 3.;
}

 static R rhs(double* P, int component) {
    return (component == 0)? Re * P[0] : Re * P[1];
 }

 static R u_exact(double* P, int component, int domain) {
     return (component == 0)? -P[1] : P[0];
 }

 static R p_exact(double* P, int component, int domain) {
     R x = P[0];
     R y = P[1];
     if (domain == 1)
        return Re * ((x * x + y * y)/2 - 1./3.) + sigma * interface_curvature;
     else
        return Re * ((x * x + y * y)/2 - 1./3.);
 }

#endif

#ifdef DISCONTINUOUS_PROBLEM

std::string problemName = "dp";

int nx = 120; // Number of elements in x direction
int ny = 60; // Number of elements in y direction

double bottom_left_x = 0.; // Bottom left corner x coordinate
double bottom_left_y = -1.; // Bottom left corner y coordinate

double width_domain = 4.; // Width of the domain
double height_domain = 2.; // Height of the domain

// Paramers in problem
const R sigma = 10.;
double interface_curvature = 1;

double mu1 = 1.;
double mu2 = 100.;

double kappa1 = 1. / 2.;
double kappa2 = 1. / 2.;


R2 shift(0, 0);

R Gamma(double* P, int i) {
   return P[1];
}

 static R rhs(double* P, int component) {
    return (component == 0)? 2 * P[0] : 4 * P[0];
 }

 static R u_exact(double* P, int component, int domain) {
     R x = P[0];
     R y = P[1];
     if (component == 0 and domain == 0)
       return x * x * y / (2*mu1);
     else if (component == 0 and domain == 1)
         return x * x * y / (2*mu2);
     else if (component == 1 and domain == 0)
         return -x * y * y / (2*mu1);
     else if (component == 1 and domain == 1)
         return -x * y * y / (2 * mu2);
     else
         return 0;
 }
 static R p_exact(double* P, int component, int domain) {
     R x = P[0];
     R y = P[1];
     if (domain == 1)
        return 2 * x * y + x * x + 10 * x;
     else
        return 2 * x * y + x * x;
 }
#endif

R divergence(double *, int i, int j) {
   return 0;
}

template <typename T>
std::string to_string(const T a_value, const int n)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

std::vector<double> method(
        double jump_penalty,
        double boundary_penalty,
        double tangential_penalty,
        double ghost_penalty,
        double delta,
        std::string paraview_filename,
        std::string export_filename)
{
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    Mesh mesh(nx
            , ny
            , bottom_left_x
            , bottom_left_y
            , width_domain
            , height_domain);

    const R h = mesh[0].lenEdge(0); //element size
    CutFEMParameter mu(mu1, mu2);
    CutFEMParameter mu_inv(1. / (2. * mu1), 1. / (2 * mu2));
    int thread_count = 1; 

    jump_penalty = jump_penalty / h;
    boundary_penalty =  boundary_penalty / h;
    tangential_penalty = tangential_penalty / h;

    Normal n;
    Tangent t;

    // Create the cut mesh
    Space Qh(mesh, DataFE<Mesh2>::P1); // Used to interpolate the levelset function
    Fun_h Gamma_interpolation(Qh, Gamma); // Interpolation of the levelset function
    InterfaceLevelSet<Mesh> interface(mesh, Gamma_interpolation); // Levelset function
    ActiveMesh<Mesh> active_mesh(mesh, interface); //Create cut mesh

    // No defualt constructor unfortunately
    #ifdef BDM1_ELEMENTS
        Space Lh(mesh, DataFE<Mesh2>::BDM1); // Space for the velocity
        Space Wh(mesh, DataFE<Mesh>::P0); // Space for the pressure
    #endif

    #ifdef RT0_ELEMENTS
        Space Lh(mesh, DataFE<Mesh2>::RT0); // Space for the velocity
        Space Wh(mesh, DataFE<Mesh>::P0); // Space for the pressure
    #endif

    #ifdef RT1_ELEMENTS
        Space Lh(mesh, DataFE<Mesh2>::RT1); // Space for the velocity
        Space Wh(mesh, DataFE<Mesh>::P1dc); // Space for the pressure
    #endif

    // Create the cut finite element space
    CutSpace Vh(active_mesh, Lh); // Cut velocity space
    CutSpace Ph(active_mesh, Wh); // Cut pressure space

    Fun_h f_interpolated(Vh, rhs); // interpolates fun_rhs
    Fun_h bc_interpolated(Vh, u_exact); // interpolates fun_bc
    Fun_h p_interpolated(Ph, p_exact); // interpolates fun_p

    // Create placeholders for the finite elements and the test functions
    FunTest u(Vh, 2), p(Ph, 1), v(Vh, 2), q(Ph, 1);

    ProblemOption optionProblem;
    optionProblem.solver_name_  = "mumps";
    optionProblem.clear_matrix_ = true;

    // Create the stokes problem object
    CutFEM<Mesh2> stokes(Vh,thread_count, optionProblem);
    stokes.add(Ph);

    // Add macro element ghost penalty
    MacroElement<Mesh> macro(active_mesh, delta);

    stokes.addBilinear(
          + contractProduct(Eps(u), 2 * mu * Eps(v))
          - innerProduct(p, div(v))
          + innerProduct(div(u), q)
          , active_mesh);

    #ifdef NON_SYMMETRIC
    stokes.addBilinear(
          - innerProduct(average(2 * mu * Eps(u) * n, kappa1, kappa2), jump(v))
          + innerProduct(jump(u), average(2 * mu * Eps(v) * n, kappa1, kappa2))
          + innerProduct(average(p, kappa1, kappa2), jump(v * n))
          + innerProduct(jump_penalty * jump(u), jump(v))
          , interface);

    stokes.addBilinear(
          + innerProduct(boundary_penalty * u, v) // Weak enforcement for u \cdot t = g \cdot t on the boundary
          - innerProduct(2 * mu * Eps(u) * n, v) // natural
          + innerProduct(u, 2 * mu * Eps(v) * n) // natural
          + innerProduct(p, v * n) // natural
          , active_mesh
          , INTEGRAL_BOUNDARY);

    stokes.addBilinear(
          - innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t))
          + innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2))
          + innerProduct(tangential_penalty * (jump(u * t)), jump(v * t))
          , active_mesh
          , INTEGRAL_INNER_EDGE_2D);
    #else
    stokes.addBilinear(
          - innerProduct(average(2 * mu * Eps(u) * n, kappa1, kappa2), jump(v))
          - innerProduct(jump(u), average(2 * mu * Eps(v) * n, kappa1, kappa2))
          + innerProduct(average(p, kappa1, kappa2), jump(v * n))
          + innerProduct(jump_penalty * jump(u), jump(v))
          , interface);

    stokes.addBilinear(
          + innerProduct(boundary_penalty * u, v) // Weak enforcement for u \cdot t = g \cdot t on the boundary
          - innerProduct(2 * mu * Eps(u) * n, v) // natural
          - innerProduct(u, 2 * mu * Eps(v) * n) // natural
          + innerProduct(p, v * n) // natural
          , active_mesh
          , INTEGRAL_BOUNDARY);
    stokes.addBilinear(
          - innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t))
          - innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2))
          + innerProduct(tangential_penalty * (jump(u * t)), jump(v * t))
          , active_mesh
          , INTEGRAL_INNER_EDGE_2D);
    #endif

    stokes.addLinear(
          + innerProduct(f_interpolated.exprList(), v)
          , active_mesh);

    stokes.addLinear(
          + innerProduct(interface_curvature * sigma, average(v * n, kappa2, kappa1))
          , interface);

    #ifdef NON_SYMMETRIC
    stokes.addLinear(
          + innerProduct(bc_interpolated.exprList(), boundary_penalty * v)
          + innerProduct(bc_interpolated.exprList(), 2 * mu * Eps(v) * n)
          , active_mesh
          , INTEGRAL_BOUNDARY);
    #else
    stokes.addLinear(
          + innerProduct(bc_interpolated.exprList(), boundary_penalty * v)
          - innerProduct(bc_interpolated.exprList(), 2 * mu * Eps(v) * n)
          , active_mesh
          , INTEGRAL_BOUNDARY);
    #endif

    #if (defined(MACRO) && defined(RT1_ELEMENTS))
    FunTest grad2un = grad(grad(u) * n) * n;
    stokes.addFaceStabilization(
          + innerProduct(ghost_penalty * h * jump(u), jump(v))
          + innerProduct(ghost_penalty * pow(h, 3) * jump(grad(u) * n), jump(grad(v) * n))
          + innerProduct(ghost_penalty * pow(h, 5) * jump(grad2un), jump(grad2un))
          - innerProduct(ghost_penalty * h * jump(p), jump(div(v)))
          - innerProduct(ghost_penalty * pow(h,3) * jump(grad(p)), jump(grad(div(v))))
          + innerProduct(ghost_penalty * h * jump(div(u)), jump(q))
          + innerProduct(ghost_penalty * pow(h,3) * jump(grad(div(u))), jump(grad(q)))
          , active_mesh
          , macro);
    #elif (defined(MACRO))
    stokes.addFaceStabilization(
          + innerProduct(ghost_penalty * h * jump(u), jump(v))
          + innerProduct(ghost_penalty * pow(h, 3) * jump(grad(u) * n), jump(grad(v) * n))
          - innerProduct(ghost_penalty * h * jump(p), jump(div(v)))
          + innerProduct(ghost_penalty * h * jump(div(u)), jump(q))
          , active_mesh
          , macro);
    #elif defined(RT1_ELEMENTS)
    FunTest grad2un = grad(grad(u) * n) * n;
    stokes.addFaceStabilization(
          + innerProduct(ghost_penalty * h * jump(u), jump(v))
          + innerProduct(ghost_penalty * pow(h, 3) * jump(grad(u) * n), jump(grad(v) * n))
          + innerProduct(ghost_penalty * pow(h, 5) * jump(grad2un), jump(grad2un))
          - innerProduct(ghost_penalty * h * jump(p), jump(div(v)))
          - innerProduct(ghost_penalty * pow(h,3) * jump(grad(p)), jump(grad(div(v))))
          + innerProduct(ghost_penalty * h * jump(div(u)), jump(q))
          + innerProduct(ghost_penalty * pow(h,3) * jump(grad(div(u))), jump(grad(q)))
          , active_mesh);
    #else
    stokes.addFaceStabilization(
          + innerProduct(ghost_penalty * h * jump(u), jump(v))
          + innerProduct(ghost_penalty * pow(h, 3) * jump(grad(u) * n), jump(grad(v) * n))
          - innerProduct(ghost_penalty * h * jump(p), jump(div(v)))
          + innerProduct(ghost_penalty * h * jump(div(u)), jump(q))
          , active_mesh);
    #endif
    // Note we only set the normal component, because of the degrees of freedom of BDM1.
    // stokes.setDirichlet(bc_interpolated, active_mesh);

    // We add the Lagrangian
    double mean_p = integral(active_mesh, p_interpolated, 0);

    CutFEM<Mesh2> lagr(Vh); lagr.add(Ph);
    Rn zero_vec = lagr.rhs_;
    lagr.addLinear(innerProduct(1., p), active_mesh);

    // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;
    std::vector<double> lag_row(lagr.rhs_.begin(), lagr.rhs_.end());
    lagr.rhs_ = 0.;

    lagr.addLinear(innerProduct(1, v * n), active_mesh, INTEGRAL_BOUNDARY);
    // lagr.addLinear(innerProduct(1., q), active_mesh);

    // stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_ ,mean_p);
    stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhsSpan(), mean_p);


    export_filename = (export_filename + non_symmetric_string + "_" + macro_string + "_"  + to_string(h, 3) + ".dat");
    matlab::Export(stokes.mat_[0], export_filename);

    stokes.solve(stokes.mat_[0], stokes.rhs_);

    // Where in the solution of Ax = b we have u and p
    Rn_ data_u = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
    Rn_ data_p = stokes.rhs_(SubArray(Ph.get_nb_dof(), Vh.get_nb_dof()));

    // Now we create the finite element solution.
    Fun_h uh(Vh, data_u);
    Fun_h ph(Ph, data_p);

    // Paraview things
    paraview_filename = (paraview_filename + non_symmetric_string + "_" + macro_string + "_" + to_string(h, 3) + ".vtk");
    Paraview<Mesh> paraview(active_mesh, paraview_filename);

    paraview.add(uh, "uh", 0, 2);
    paraview.add(ph, "ph", 0, 1);

    // Error calculations
    Fun_h u_error(Vh, u_exact);
    Fun_h p_error(Ph, p_exact);
    Fun_h u_exact_interpolated(Vh, u_exact);
    Fun_h p_exact_interpolated(Ph, p_exact);

    auto dxu = dx(uh.expr(0));
    auto dyu = dy(uh.expr(1));

    u_error.v -= uh.v;
    p_error.v -= ph.v;

    paraview.add(u_exact_interpolated, "u_exact", 0, 2);
    paraview.add(p_exact_interpolated, "p_exact", 0, 1);
    paraview.add(u_error, "u_error", 0, 2);
    paraview.add(p_error, "p_error", 0, 1);
    paraview.add(fabs(dxu + dyu), "div_u");

    double error_u = L2normCut(uh, u_exact, 0, 2);
    double error_p = L2normCut(ph, p_exact, 0, 1);
    double error_div = L2normCut(fabs(dxu + dyu), divergence, active_mesh);

    return std::vector<double> {h, error_u, error_p, error_div};
}


int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);

    std::string filename_paraview = "../../data/paraview/stokes/" + elements + "/" + problemName + "/" ;
    std::string filename_error_data = "../../data/error/stokes/"+ elements + "/" + problemName + "/" ;
    std::string filename_export = "../../data/matrices/stokes/"+ elements + "/" + problemName + "/" ;

    double jump_penalty = 5000.; //Is Niche term, so divide by h
    double boundary_penalty = 20.; //Is Niche term, so divide by h
    double tangential_penalty = 7.; //Is Niche term, so divide by h
    double ghost_penalty = 1.;
    double delta = 0.25;

    #ifdef MACRO 
        macro_string = to_string(delta, 2);
    #endif

    nx = ny = 10;
    int it = 5;

    std::vector<std::vector<double>> results;

    for (int i = 0; i < it; i++) {
        results.push_back(method(jump_penalty, boundary_penalty, tangential_penalty, ghost_penalty, delta, filename_paraview, filename_export));
        // nx += 10;
        // ny += 10;
        nx *= 2;
        ny *= 2;
    }

    std::cout << '\n';
    std::cout << std::string(25, ' ') << "The problem that is being solved is " + problemName << std::endl;
    std::cout << std::string(35, ' ') << "The elements used are " + elements << std::endl;
    std::cout << std::string(37, ' ') << "The delta used is " << delta << std::endl;
    std::cout << std::endl;

    std::cout << std::setw(15) << "h" << " "
              << std::setw(15) << "error_u" << " "
              << std::setw(15) << "order_u" << " "
              << std::setw(15) << "error_p" << " "
              << std::setw(15) << "order_p" << " "
              << std::setw(15) << "error_div" << " "
              << '\n' << '\n';

    for (int i = 0; i < results.size(); i++) {
        for (int j = 0; j < results[i].size(); j++) {
            std::cout << std::setw(15) << results[i][j] << " ";
            if (j != 0 && j != 3) {
                if (i == 0) {
                    std::cout << std::setw(15) << "----  " << " ";
                }
                else {
                    std::cout << std::setw(15) << std::log(results[i - 1][j] / results[i][j]) / std::log(results[i-1][0]/results[i][0]) << " ";
                }
            }
        }
            std::cout << std::endl;
    }

    // Drawing
    std::vector<double> h;
    std::vector<double> error_u;
    std::vector<double> error_p;
    std::vector<double> error_div;
    for (auto result : results) {
        h.push_back(result[0]);
        error_u.push_back(result[1]);
        error_p.push_back(result[2]);
        error_div.push_back(result[3]);
    }

    std::ofstream file;
    file.open(filename_error_data + non_symmetric_string + "_" + macro_string + "_" + to_string(h[0], 3) + "_" + to_string(h.back(), 3) + "_" + std::to_string(it) + ".csv");
    file << "#problem," << problemName << std::endl;
    file << "#elements," << elements << std::endl;
    file << "#delta," << delta << std::endl;
    file << "#jump_penalty," << jump_penalty << std::endl;
    file << "#boundary_penalty," << boundary_penalty << std::endl;
    file << "#tangential_penalty," << tangential_penalty << std::endl;
    file << "#ghost_penalty," << ghost_penalty << std::endl;
    file << "h,error_u,error_p,error_div" << std::endl;
    for (int i = 0; i < h.size(); i++) {
        file << h[i] << "," << error_u[i] << "," << error_p[i] << "," << error_div[i] << std::endl;
    }
}
