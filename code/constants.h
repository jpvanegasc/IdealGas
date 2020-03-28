#include<iostream>
#include<cmath>
#include<fstream>

// Geometry costants
const double square = 10.0;
const double Lx = square, Ly = square;

// PEFRL
const double Zi = 0.1786178958448091e0;
const double Lambda = 0.2123418310626054*(-1);
const double Xi = 0.06626458266981849*(-1);

const double coef1 = (1 - 2*Lambda)/2;
const double coef2 = (1 - 2*(Xi+Zi));

// Implementation
const int N = 1;