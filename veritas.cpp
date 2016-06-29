#include "veritas.hpp"
#include "Settings.hpp"
#include "SolverManager.hpp"
#include "EMSolver.hpp"


void initialConditions(Input &grid, Particles &particles, Output &output) {

    grid.minEfficiency = 0.6, grid.dx = 0.5, grid.k = 0.01; //double
    grid.refinementCriteria = 1e-8, grid.cfl = 0.5, grid.sizeWeight = 0.0; //double
    grid.preLength = 0, grid.postLength = 0; //double
    grid.nx = 76, grid.r = 2, grid.Lfinest = 5, grid.loadBalanceLength = 30; //unsigned int
    grid.tempEM = {0.0}; //std::vector<double>

    particles.mass = {9.10938291e-31,9.10938291e-31*1836}; //std::vector<double>
    particles.charge = {-1.60217657e-19,1.60217657e-19}; //std::vector<double>
    particles.misc = {{0.0,0.01},{0.0,0.01}}; //std::vector<std::vector<double>>
    particles.np = {250,50}; //std::vector<unsigned int>
    particles.dp = {0.1,0.1}; //std::vector<double>
    particles.pmin = {0.1,0.1}; //std::vector<double>

    output.rectangleData = true; //bool
    output.charge = true; //bool
    output.energy = true; //bool
    output.potential = true; //bool
    output.EFieldLongitudinal = true; //bool
    output.EFieldTransverse = true; //bool
    output.BFieldTransverse = true; //bool
    output.AFieldSquared = true; //bool
    output.precision = 15; //int

}


void Settings::settingsOverride() {
    // Note that n_x and n_p must be multiples of the refinement ration.

    pmin [0]= -40 * m[0] * cs;
    pmin [1]= -400 * m[0] * cs;

    dp[0] = 80 * m[0] * cs / (p_size[0] - 1);
    dp[1] = 800 * m[0] * cs / (p_size[1] - 1);

    double lambda = 1e-6;
    dx = 10 * lambda / (x_size);

    sizeWeight = refinementCriteria * 10000;

    time = 0.0;

    quadratureDepth = 2;

    double T = lambda / cs;
    double omega = DPI / T;
    double Ec = m[0] * cs / std::fabs(q[0]);
    double Nc = omega * omega * m[0] * eps0 / (q[0] * q[0]);
    double a0 = 1.0;
    double n = 2.0;

    temp[0].resize(2, 0.0);
    temp[0][0] = n * Nc;
    temp[0][1] = 5e-4 * std::pow(m[0] * cs, 2.0);

    temp[1].resize(2, 0.0);
    temp[1][0] = n * Nc;
    temp[1][1] = 5e-4 * std::pow(m[0] * cs, 2.0)*1836.0;

    tempEM.resize(2, 0.0);
    tempEM[0] = lambda;
    tempEM[1] = a0 * Ec / std::sqrt(2.0);

    //end overrides.
}

bool Settings::RefinementOverride(double x, double p, int depth, int particleType) {
    return (x > 2.7e-6) && (x < 7.3e-6) && (std::fabs(p*p/(2)) < (temp[particleType][1]));
}

double Settings::GetBY(double x, double t) {
    double T = tempEM[0] / cs;
    double k = DPI / tempEM[0];
    double omega = cs * k;

    // sin^2 ramp
    if (t < (4 * T)) {
        return -tempEM[1] * (-k * std::pow(std::sin((DPI / 16.0) * ((t - x / cs) / T)), 2.0) * std::cos(omega * t - k * x) - 2 * DPI / (16.0 * T * cs) * std::pow(std::sin((DPI / 16.0) * (t - x / cs) / T), 1.0) * std::cos((DPI / 16.0) * (t - x / cs) / T) * std::sin(omega * t - k * x));
    }

    return tempEM[1] * (k * std::cos(omega * t - k * x));
}

double Settings::GetBZ(double x, double t) {
    double T = tempEM[0] / cs;
    double k = DPI / tempEM[0];
    double omega = cs * k;

    // sin^2 ramp
    if (t < (4 * T)) {
        return tempEM[1] * (k * std::pow(std::sin((DPI / 16.0) * ((t - x / cs) / T)), 2.0) * std::sin(omega * t - k * x) - 2 * DPI / (16.0 * T * cs) * std::pow(std::sin((DPI / 16.0) * (t - x / cs) / T), 1.0) * std::cos((DPI / 16.0) * (t - x / cs) / T) * std::cos(omega * t - k * x));
    }

    return tempEM[1] * (k * std::sin(omega * t - k * x));
}

double Settings::InitialDistribution(double x, double p, int particleType) {

    double n0=0.0;

    if ((x > 3.0e-6)&&(x < 7.0e-6)) {
        n0=1.0;
    }

    return n0*temp[particleType][0]  * std::exp(-(p*p) / (2.0 * temp[particleType][1])) / std::sqrt(DPI * temp[particleType][1]);

}

bool LOUD = false, NOISY = false;
int main() {
    omp_set_nested(1);

    auto start = std::chrono::steady_clock::now();

    Input grid;
    Particles particles;
    Output output;

    initialConditions(grid, particles, output);

    Settings settings(grid, particles, output);
    SolverManager SM(settings);

    double T = settings.tempEM[0] / cs;
    double dt = T / 200, t = dt, T_stop = 8.5*T;
    double tLast = 0.0;
    int counter = 0;

    while (t < T_stop) {

        double dt_adaptive = std::min(SM.CalculateDt(settings.cfl), dt);

        if (t > 3 * T) {
            SM.Advance(dt_adaptive);
        } else {
            dt_adaptive=dt;
            SM.AdvanceFields(dt_adaptive);
        }

        if ((counter > 200)&&(t > 3 * T)) {
            SM.reGrid(t);
            counter = 0;
        } else {
            counter++;
        }

        tLast += dt_adaptive;

        if (tLast > T / 10) {
            SM.fileOutput(t);
            tLast = 0.0;
        }

        t += dt_adaptive;
    }
    auto end = std::chrono::steady_clock::now();
    std::cout << "Test completed in " << std::fixed << std::setprecision(2) << std::chrono::duration<double, std::ratio<60, 1> >(end - start).count() << " mins (" << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s)." << std::endl;

    return(EXIT_SUCCESS);
}
