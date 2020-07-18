/**
 * Animation module.
 * For calculating and saving data please use main.cpp
 */
#include"metropolis2D.h"
#include"animation.hpp"

int main(int argc, char *argv[]){
    CRandom ran64(1);
    Metropolis Gas;

    double T = 1000;
    double beta = 1.0/(kB*T);

    Gas.initialize(ran64);
    start_animation(argc);


    for(int t=0; t<t_eq; t++)
        for(int mcs=0; mcs<N; mcs++)
            Gas.metropolis_step(beta, ran64);

    for(int t=0; t<100; t++){
        begin_frame(argc);
        Gas.print();
        end_frame(argc);

        Gas.metropolis_step(beta, ran64);
    }

    return 0;
}