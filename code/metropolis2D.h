/**
 * Monte Carlo Ising 2D by Metropolis algorithm
 */
#include"constants.h"

class Metropolis{
    private:
        int m[Lx][Ly]; int E, M;
    public:
        void initialize_down(void);
        void metropolis_step(double beta, CRandom &ran);
        double get_E(void);
        double get_M(void);
};

inline double Metropolis::get_E(void){return (double)E;}
inline double Metropolis::get_M(void){return (double)std::abs(M);}

void Metropolis::initialize_down(void){
    for(int i=0; i<Lx; i++)
        for(int j=0; j<Ly; j++)
            m[i][j] = -1;
        M = -L2; E = -2*L2;
}

void Metropolis::metropolis_step(double beta, CRandom &ran){
    int n = (int)L2*ran.r(); int i = n%Lx, j = n/Ly; // Just one random number to save computing time
    int dE = 2*m[i][j]*(m[(i+1)%Lx][j] + m[(Lx+i-1)%Lx][j] + m[i][(j+1)%Ly] + m[i][(Ly+j-1)%Ly]);

    if(dE <= 0){
        m[i][j] *= -1; E += dE; M += 2*m[i][j];
        return;
    }
    else if(ran.r() < std::exp(-beta*dE)){
        m[i][j] *= -1; E += dE; M += 2*m[i][j];
        return;
    }
}
