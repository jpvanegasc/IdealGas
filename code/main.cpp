/**
 * Basic Monte Carlo implementation using Metropolis 
 */
#include"metropolis2D.h"
#include"animation.hpp"

int main(void){
    CRandom ran64(1);
    Metropolis Ising;

    double T = 0.01;
    double beta = 1.0/(kB*T);

    Ising.initialize(ran64);
    int N1_old = Ising.get_N(0), N2_old = Ising.get_N(1);

    for(int t=0; t<t_eq; t++){

        Ising.metropolis_translation(beta, ran64, 0);
        Ising.metropolis_translation(beta, ran64, 1);

        Ising.metropolis_transfer(0, 1, beta, ran64);
        Ising.metropolis_transfer(1, 0, beta, ran64);
    }

    // Process data & print
    int N1_new = Ising.get_N(0), N2_new = Ising.get_N(1);
    double rho1_O = ((double)N1_old)/L2, rho2_O = ((double)N2_old)/L2;
    double rho1_n = ((double)N1_new)/L2, rho2_n = ((double)N2_new)/L2;

    std::cout << "old: " << rho1_O << ", " << rho2_O << '\n'
        << "new: " << rho1_n << ", " << rho2_n << std::endl;

    std::cout << t_eq << std::endl;

    return 0;
}