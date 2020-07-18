/**
 * Monte Carlo Ising 2D by Metropolis algorithm
 */
#include"constants.h"

class Metropolis{
    private:
        double r[N][2];

        double calculate_energy(int n_particle, double x, double y);
    public:
        void initialize(CRandom &ran);
        void metropolis_step(double beta, CRandom &ran);
        void print(void);
};

/**
 * Calculates the energy for a given particle on a given postion using the Lennard-Jones interaction
 * potential.
 * 
 * @param n_particle: particle id that is being calculated
 * @params x, y: x and y positions of particle that is being calculated
 * 
 * @return Energy for a given particle
 */
double Metropolis::calculate_energy(int n_particle, double x, double y){
    double E = 0;

    #pragma omp parallel for shared(E, r) reduction(+: E)
    for(int n=0; n<N; n++){
        if(n == n_particle) continue;

        double ix = r[n][0], iy = r[n][1];
        double r = std::sqrt((ix-x)*(ix-x) + (iy-y)*(iy-y));
        E += 4*(std::pow(r, -12.0) - std::pow(r, -6.0));
    }

    return E;
}

/**
 * Initialize both boxes with random positions
 */
void Metropolis::initialize(CRandom &ran){
    for(int n=0; n<N; n++){
        r[n][0] = Lx*ran.r(); r[n][1] = Ly*ran.r();
    }
}

/**
 * Moves a single particle within its box using the Metropolis criteria, with periodic boundaries.
 * -> Not very good yet. Selection from within box is not fully correct
 * 
 * @param beta: Thermodynamic beta
 * @param box: id of box frmom which the particle is being selected
 */
void Metropolis::metropolis_step(double beta, CRandom &ran){
    int n = (int)N*ran.r(); double x = r[n][0], y = r[n][1];
    double drx = dr*(2*ran.r()-1), dry = dr*(2*ran.r()-1);

    // Periodic boundaries
    double x_new = x + drx, y_new = y + dry;

    if(x_new < 0) x_new += Lx;
    else if(x_new > Lx) x_new -= Lx;

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

void Metropolis::print(void){
    for(int n=0; n<N; n++)
        std::cout << " , " << r[n][0] << "+0.1*cos(t)," << r[n][1] << "+0.1*sin(t)";
}
