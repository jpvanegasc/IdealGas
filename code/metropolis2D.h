/**
 * Monte Carlo Ising 2D by Metropolis algorithm
 */
#include"constants.h"

class Metropolis{
    private:
        int b[Lx][Ly][2];
    public:
        void initialize(CRandom &ran);
        double calculate_energy(int ix, int iy, int box);
        void metropolis_step(double beta, CRandom &ran, int box);
        void print(void);
};

void Metropolis::initialize(CRandom &ran){
    int n = N;
    while(n < 0){
        int r = (int)L2*ran.r(); int i = r%Lx, j = r/Ly;
        if(b[i][j][0] == 0){
            b[i][j][0] == 1; n--;
        }
    }
    n = N;
    while(n < 0){
        int r = (int)L2*ran.r(); int i = r%Lx, j = r/Ly;
        if(b[i][j][1] == 0){
            b[i][j][1] == 1; n--;
        }
    }
}

double Metropolis::calculate_energy(int ix, int iy, int box){
    double E = 0;
    for(int i=0; i<Lx; i++)
        for(int j=0; j<Ly; j++){
            if(i == ix && j == iy) continue;
            if(b[i][j][box] == 0) continue;
            else if(b[i][j][box] == 1){
                double r = std::sqrt((i-ix)*(i-ix) + (j-iy)*(j-iy));
                E += 4*(std::pow(r, -12.0) - std::pow(r, -6.0));
            }
        }

    return E;
}

void Metropolis::metropolis_step(double beta, CRandom &ran, int box){
    bool flag = true;
    while(flag){
        double r = ran.r(); // Just one random number to save computing time
        int n = (int)L2*r; int i = n%Lx, j = n/Ly;
        if(b[i][j][box] == 0) continue;
        else flag = false;

        double dE = calculate_energy(i, j, box);

        if(dE <= 0){
            b[i][j][box] *= -1;
            return;
        }
        else if(ran.r() < std::exp(-beta*dE)){
            b[i][j][box] *= -1;
            return;
        }
    }
}
