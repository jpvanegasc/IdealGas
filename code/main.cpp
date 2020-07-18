/**
 * Main  module. This is not for creating animations, its sole purpose is saving & processing data
 * For anmations please use animate.cpp
 */
#include"metropolis2D.h"

int main(void){
    CRandom ran64(1);
    Metropolis Ising;
    std::ofstream File("data.dat");

    double T = 10;
    double beta = 1.0/(kB*T);
    double t_max = 1e4;

    Ising.initialize(ran64);

    for(int t=0; t<t_eq; t++)
        for(int mcs=0; mcs<N; mcs++)
            Ising.metropolis_step(beta, ran64);

    for(int t=0; t<t_max; t++){
        for(int mcs=0; mcs<N; mcs++)
            Ising.metropolis_step(beta, ran64);
        File << t << '\t' << Ising.get_mean_r() << '\n';
    }

    File.close();

    return 0;
}