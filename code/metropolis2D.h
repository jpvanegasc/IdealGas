/**
 * Monte Carlo Ising 2D by Metropolis algorithm
 */
#include"constants.h"

class Metropolis{
    private:
        double r[2*N][2];

        double calculate_energy(int n_particle, double x, double y);
    public:
        void initialize(CRandom &ran);
        void metropolis_translation(double beta, CRandom &ran, int box);
        void metropolis_transfer(int from, int destination, double beta, CRandom &ran);
        int get_N(int box);
        void print(void);
};

double Metropolis::calculate_energy(int n_particle, double x, double y){
    double E = 0;

    #pragma omp parallel for shared(E, r) reduction(+: E)
    for(int n=0; n<2*N; n++){
        if(n == n_particle) continue;

        double ix = r[n][0];

        bool flag = ((x<=Lx)&&(ix<=Lx)) || ((x>Lx)&&(ix>Lx));
        if(flag == true){

            double iy = r[n][1];
            double r = std::sqrt((ix-x)*(ix-x) + (iy-y)*(iy-y));
            E += 4*(std::pow(r, -12.0) - std::pow(r, -6.0));
        }
    }

    return E;
}


void Metropolis::initialize(CRandom &ran){
    for(int n=0; n<N; n++){
        r[n][0] = Lx*ran.r(); r[n][1] = Ly*ran.r();
    }
    for(int n=N; n<2*N; n++){
        r[n][0] = Lx*ran.r() + Lx; r[n][1] = Ly*ran.r();
    }
}

void Metropolis::metropolis_translation(double beta, CRandom &ran, int box){
    int n = (int)N*ran.r() + box*N; double x = r[n][0], y = r[n][1];
    double drx = dr*(2*ran.r()-1), dry = dr*(2*ran.r()-1);

    // Periodic boundaries
    double x_new = x + drx, y_new = y + dry;

    if(x_new < 0) x_new += Lx;
    else if(x_new > (Lx + box*Lx)) x_new -= Lx;

    if(y_new < 0) y_new += Ly;
    if(y_new > Ly) y_new -= Ly;

    // Metropolis criteria
    double dE = calculate_energy(n, x_new, y_new) - calculate_energy(n, x, y);

    if(dE <= 0){
        r[n][0] = x_new; r[n][1] = y_new;
        return;
    }
    else if(ran.r() < std::exp(-beta*dE)){
        r[n][0] = x_new; r[n][1] = y_new;
        return;
    }
}

void Metropolis::metropolis_transfer(int from, int destination, double beta, CRandom &ran){
    int n = (int)N*ran.r() + from*N; double x = r[n][0], y = r[n][1];
    double x_new = Lx*ran.r() + destination*Lx, y_new = Ly*ran.r();

    // Metropolis criteria
    double dE = calculate_energy(n, x_new, y_new) - calculate_energy(n, x, y);

    if(dE <= 0){
        r[n][0] = x_new; r[n][1] = y_new;
        return;
    }
    else if(ran.r() < std::exp(-beta*dE)){
        r[n][0] = x_new; r[n][1] = y_new;
        return;
    }
}

int Metropolis::get_N(int box){
    int N_particles = 0;
    if(box == 0){
        #pragma omp parallel for shared(N, r) reduction(+, N)
        for(int n=0; n<2*N; n++){
            if(r[n][0] <= Lx) N_particles++;
        }
    }
    else if(box == 1){
        #pragma omp parallel for shared(N, r) reduction(+, N)
        for(int n=0; n<2*N; n++){
            if(r[n][0] > Lx) N_particles++;
        }
    }

    return N_particles;
}

void Metropolis::print(void){
    for(int n=0; n<2*N; n++)
        std::cout << " , " << r[n][0] << "+0.1*cos(t)," << r[n][1] << "+0.1*sin(t)";
}
