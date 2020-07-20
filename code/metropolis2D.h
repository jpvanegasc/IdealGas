/**
 * 2D Monte Carlo by Metropolis algorithm
 */
#include"constants.h"

class Metropolis{
    private:
        double r[N][D] = {0};

        double calculate_energy(int &n_particle, double &x0, double &y0);
    public:
        void initialize(CRandom &ran);
        void metropolis_step(double beta, CRandom &ran);
        double get_mean_r(void);
        void print(void);
};

/**
 * Calculates the energy for a given particle on a given postion using the Lennard-Jones interaction
 * potential. This function is the most called of all the module, so it's heavily optimized. Sorry 
 * if it's not very readable.
 * 
 * @param n_particle: particle id that is being calculated
 * @param x: x position of particle that is being calculated
 * @param y: y position of particle that is being calculated
 * 
 * @return Energy for a given particle
 */
double Metropolis::calculate_energy(int &n_particle, double &x0, double &y0){
    double E = 0, dx = 0, dy = 0, dr2 = 0; int n = 0;

    for(n=0; n<n_particle; n++){
        dx = r[n][0] - x0; dy = r[n][1] - y0;
        dx -= L*std::fabs(dx/L);
        dy -= L*std::fabs(dy/L);

        dr2 = dx*dx + dy*dy;
        dr2 *= dr2*dr2;
        dr2 = 1/dr2;

        E += 4*dr2*(dr2-1);
    }
    for(n=n_particle+1; n<N; n++){

        dx = r[n][0] - x0; dy = r[n][1] - y0;
        dx -= L*std::fabs(dx/L);
        dy -= L*std::fabs(dy/L);

        dr2 = dx*dx + dy*dy;
        dr2 *= dr2*dr2;
        dr2 = 1/dr2;

        E += 4*dr2*(dr2-1);
    }

    return E;
}

/* Initialize with random positions */
void Metropolis::initialize(CRandom &ran){
    for(int n=0; n<N; n++){
        r[n][0] = L*ran.r(); r[n][1] = L*ran.r();
    }
}

/**
 * Moves a single particle using the Metropolis criteria, with periodic boundaries.
 * 
 * @param beta: Thermodynamic beta
 */
void Metropolis::metropolis_step(double beta, CRandom &ran){
    int n = (int)N*ran.r(); double x = r[n][0], y = r[n][1];
    double drx = dr*(2*ran.r()-1), dry = dr*(2*ran.r()-1);

    // Periodic boundaries
    double x_new = x + drx, y_new = y + dry;

    if(x_new < 0) x_new += L;
    else if(x_new > L) x_new -= L;

    if(y_new < 0) y_new += L;
    else if(y_new > L) y_new -= L;

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

/**
 * Calculates the mean minimum distance between all particles of the system in a given state
 * 
 * @return minimum mean distance
 */
double Metropolis::get_mean_r(void){
    double mean_r = 0, dx = 0, dy = 0, dr2 = 0;

    for(int n1=0; n1<N; n1++)
        for(int n2=n1+1; n2<N; n2++){
            dx = r[n1][0] - r[n2][0]; dy = r[n1][1] - r[n2][1];
            dx -= L*std::fabs(dx/L);
            dy -= L*std::fabs(dy/L);

            dr2 = dx*dx + dy*dy;

            mean_r += std::sqrt(dr2);
        }

    mean_r /= N_connections;

    return mean_r;
}

void Metropolis::print(void){
    for(int n=0; n<N; n++)
        std::cout << " , " << r[n][0] << "+0.1*cos(t)," << r[n][1] << "+0.1*sin(t)";
}
