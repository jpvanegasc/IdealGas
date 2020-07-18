/**
 * Main  module. This is not for creating animations, its sole purpose is saving & processing data
 * For anmations please use animate.cpp
 */
#include"metropolis2D.h"

int main(void){
    CRandom ran64(1);
    Metropolis Ising;

    double T = 0.01;
    double beta = 1.0/(kB*T);

    Ising.initialize(ran64);

    for(int t=0; t<t_eq; t++){
        Ising.metropolis_step(beta, ran64);
    }

    return 0;
}