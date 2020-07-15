#include<iostream>
#include<cmath>
#include<fstream>

#include"random64.h"

// Geometry constants
const int N = 5;
const int square = 16;
const int Lx = square, Ly = square;
const int L2 = Lx*Ly;

// Time constants
const int t_eq = (int)(600*std::pow(std::sqrt(L2)/8.0, 2.1));
const int t_corr = (int)(5*10.8*std::pow(std::sqrt(L2)/8.0, 2.125));

// Thermodynamic constants
const double kB = 1.0;
const double rho = 0.3;
