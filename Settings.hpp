#ifndef __settings_hpp__
#define __settings_hpp__

struct Input {
    double minEfficiency = 0.6, dx = 0.5, k = 0.01;
    double refinementCriteria = 1e-8, cfl = 0.9, sizeWeight = 0.0;
    double preLength = 2, postLength = 2;
    unsigned int nx = 76*4, r = 2, Lfinest = 5;
    double plasma_xl_bound = 3.0e-6, plasma_xr_bound = 7.0e-6;
    std::vector<double> tempEM = {0.0};
};

struct Particles {
    std::vector<double> mass = {9.10938291e-31};
    std::vector<double> charge = {-1.60217657e-19};
    std::vector<std::vector<double>> misc = {{0.0,0.01}};
    std::vector<unsigned int> np = {50};
    std::vector<double> dp = {0.1};
    std::vector<double> pmin = {0.1};
};

struct Output {
    bool time = true;
    bool rectangleData = true;
    bool charge = true;
    bool energy = true;
    bool potential = true;
    bool EFieldLongitudinal = true;
    bool EFieldTransverse = true;
    bool BFieldTransverse = true;
    bool AFieldSquared = true;
    int precision = 15;
};

class Settings {
    EMFieldSolver *EMSolver;
public:
    Output output;
    double dx, time, minEfficiency,preLength,postLength,refinementCriteria,cfl,sizeWeight,plasma_xl_bound,plasma_xr_bound;
    std::vector<double> m, q, m_inv, dp, fMax, tempEM, pmin;
    std::vector<std::vector<double>> temp;
    unsigned int x_size_finest, x_size, refinementRatio;
    int maxDepth,quadratureDepth;
    std::vector<unsigned int> p_size, p_size_finest;
    Settings(const Input &grid, const Particles &particles, const Output &out);
    void title(std::string output, char spacer);
    double InitialDistribution(double x, double p, int particleType);
    double GetDp(int level, int particleType);
    double GetDx(int level);
    int GetXSize(int level);
    int GetPSize(int level, int particleType);
    double GetMass(int i);
    double GetCharge(int i);
    double GetVectorPotentialY(double x, double t);
    double GetVectorPotentialZ(double x, double t);
    void UpdateTime(int step, double dt);
    double GetEY(double x, double t);
    double GetEZ(double x, double t);
    double GetBY(double x, double t);
    double GetBZ(double x, double t);
    void DetermineMaximum();
    double GetfMax(int i);
    void settingsOverride();
    bool RefinementOverride(double x,double p,int depth,int particleType);
};

#endif /* __settings_hpp__ */
