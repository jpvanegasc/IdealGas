/**
 * Animation module.
 * For calculating and saving data please use main.cpp
 */
#include"metropolis2D.h"
#include"animation.hpp"

int main(int argc, char *argv[]){
    CRandom ran64(1);
    Metropolis Ising;

    double T = 0.01;
    double beta = 1.0/(kB*T);

    Ising.initialize(ran64);
    start_animation(argc);

    for(int t=0; t<t_eq/2; t++){
        begin_frame(argc);
        Ising.print();
        end_frame(argc);

        Ising.metropolis_step(beta, ran64);
    }

    return 0;
}