#include<iostream>
#include<cmath>
#include<fstream>

#include"random64.h"

// Geometry constants
const int D = 2;

const int N = 80;
const double L = 16;
const double L2 = L*L;

const double dr = 1;

// Time constants
const int t_eq = (int)(600*std::pow(std::sqrt(L2)/8.0, 2.1));
const int t_corr = (int)(5*10.8*std::pow(std::sqrt(L2)/8.0, 2.125));

// Thermodynamic constants
const double kB = 1.0;
