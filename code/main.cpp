/**
 * Main  module. This is not for creating animations, its sole purpose is saving & processing data
 * For animations please use animate.cpp
 */
#include"metropolis2D.h"

int main(void){
    CRandom ran64(1);
    Metropolis Gas;
    std::ofstream File("data.dat");

    double T = 10;
    double beta = 1.0/(kB*T);
    double t_max = 1e4;

    Gas.initialize(ran64);

    for(int t=0; t<t_eq; t++)
        for(int mcs=0; mcs<N; mcs++)
            Gas.metropolis_step(beta, ran64);

    for(int t=0; t<t_max; t++){
        for(int mcs=0; mcs<N; mcs++)
            Gas.metropolis_step(beta, ran64);
        File << t << '\t' << Gas.get_mean_r() << '\n';
    }

    File.close();

    return 0;
}