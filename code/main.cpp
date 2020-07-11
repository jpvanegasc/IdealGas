/**
 * Basic Monte Carlo implementation using Metropolis 
 */
#include"metropolis2D.h"

int main(void){
    CRandom ran64(1);
    Metropolis Ising;
    std::ofstream File("data.dat");

    double E_prom, E2_prom, M_prom, M2_prom, M4_prom;
    int t_data = 1e4;

    for(double T=0.1; T<10; T += 0.1){
        double beta = 1.0/(kB*T);
        E_prom = E2_prom = M_prom = M2_prom = M4_prom = 0;
        Ising.initialize_down();

        // Get to equilibrium
        for(int t=0; t<t_eq; t++){
            for(int mcs=0; mcs<L2; mcs++)
                Ising.metropolis_step(beta, ran64);
        }

        // Get & save data, increment time
        for(int t=0; t<t_data; t++){
            double E = Ising.get_E(), M = Ising.get_M();
            E_prom += E; E2_prom += E*E; M_prom += M; M2_prom += M*M; M4_prom += M*M*M*M;

            for(int mcs=0; mcs<L2; mcs++)
                Ising.metropolis_step(beta, ran64);
        }

        // Process data & print
        E_prom /= t_data; E2_prom /= t_data;
        M_prom /= t_data; M2_prom /= t_data; M4_prom /= t_data;
        double Xs = beta*(M2_prom - M_prom*M_prom);
        double Cv = beta*beta*kB*(E2_prom - E_prom*E_prom);
        double UBinder = 1.0 -(1.0/3.0)*(M4_prom/(M2_prom*M2_prom));

        File << 
            T << '\t' << 
            E_prom << '\t' << 
            Cv << '\t' << 
            M_prom << '\t' << 
            Xs << '\t' << 
            UBinder << std::endl;
    }

    File.close();

    return 0;
}