#pragma once


double affineMap(const double& r, const double& xL, const double& xR);
double mapAffine(const double& x, const double& xL, const double& xR);


// Gaussian.cpp
template <size_t N> // order
struct storedGL;

template<size_t N>
storedGL<N> GaussLegendre();

// GaussianLobatto.cpp
template <size_t N> // order
struct storedGLobatto;

//double L(double r, const double* gLnodes, const double* uCell);


//main.cpp
class Grid;
class Solution;

//Lagrangian.cpp
template <size_t Np>
class LagrangianPoly;


//#include "Gaussian.cpp"
#include "GaussianLobatto.cpp"
#include "Lagrangian.cpp"
