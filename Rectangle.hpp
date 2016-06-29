#ifndef __Rectangle_hpp__
#define __Rectangle_hpp__

#include "Settings.hpp"

class Rectangle {
    std::vector<double> errorWeights,interpolationMatrix,interpolatCoefsREL,interpolatCoefsREF,fp,fx,ep,ex,FxL,FxH,FpL,FpH,FxLS,FpLS,FxDS,FpDS,Rp,Rm,Cx,Cp;
    std::vector<std::shared_ptr<Rectangle>> finer_level,finer_level_x,finer_level_p;
    std::vector<std::weak_ptr<Rectangle>> exterior_boundary_xm, exterior_boundary_xp, exterior_boundary_pm, exterior_boundary_pp;
    std::vector<bool> is_nested,is_interpolated;
    std::vector<bool> is_interrior_level_boundary_x, is_interrior_level_boundary_p;
    std::vector<bool> is_exterrior_boundary_same_level_xm, is_exterrior_boundary_same_level_xp, is_exterrior_boundary_same_level_pm, is_exterrior_boundary_same_level_pp;
    std::vector<double> chargeRNested,currentRNested,chargeRNoNested,currentRNoNested;
    int refinementRatio, depth, particleType;
    double dp, dx, m_inv;
    Settings &settings;
    EMFieldSolver *EMSolver;
    bool up, down, left, right;
public:
    int n_x, n_p, x_pos, p_pos;
    double relativeToBottom;
    std::vector<double> f, chargeR, energyR, currentR;
    Rectangle(int n_x, int n_p, int x_pos, int p_pos, int depth, Settings &settings, const std::shared_ptr<Rectangle> &bc, bool up, bool down, bool left, bool right, int particleType);
    Rectangle();
    ~Rectangle();
    void UpdateInterriorPoints(int val=1);
    void UpdateSameLevelBoundaries(int val=1);
    void UpdateDifferentLevelBoundaries(int val=1);
    inline int Index(int i, int j, int state);
    inline int IndexNS(int i, int j);
    inline int GetFinestIndex(int i);
    inline double Momentum(double i);
    inline double CalculateFluxToCoarseP(int i, int j);
    inline double CalculateFluxToCoarseX(int i, int j);
    inline double CalculateFluxToCoarsePL(int i, int j);
    inline double CalculateFluxToCoarseXL(int i, int j);
    double GetValueFromFinerLevel(int i, int j,int val=1);
    virtual double GetValueFromSameLevel(int i, int j,int val=1);
    void CalculateRhoAndJ();
    void CalculateEnergy();
    void InitializeDistribution();
    void SetFieldSolver(EMFieldSolver *solver);
    std::vector<double> GetWenoValueFromCoarseLevel(int i, int j, int d,int val=1);
    double RGKGetFluxP(int i, int j);
    double RGKGetFluxX(int i, int j);
    double GetEfield(int i);
    void CalculateConnectivitySame(const std::shared_ptr<Rectangle> &rectangle);
    void CalculateConnectivityFromFiner(const std::shared_ptr<Rectangle> &current, const std::shared_ptr<Rectangle> &rectangle);
    double GetWenoEdgeValue(double f1,double f2,double f3,double f4,bool right);
    void GetDataFromCoarseLevelRectangle(const std::shared_ptr<Rectangle> &rectangle);
    void GetDataFromSameLevelRectangle(const std::shared_ptr<Rectangle> &rectangle);
    inline double Gamma(double p_x,double a_squared);
    void getError(std::vector<coords> &flaggedCells, int particleType);
    bool OnLine(int x,int y1,int y2);
    inline double ErrorEstimate(int i,int j);
    void GetDataFromCoarseNewLevelRectangle(const std::shared_ptr<Rectangle> &rectangle);
    std::vector<double> GetInterpolantsREF(double f1, double f2, double f3, double f4, double f5);
    std::vector<double> GetInterpolantsREL(double f1, double f2, double f3, double f4, double f5);
    void UpdateCornerPoints(int val=1);
    void SetCFromDifferentLevel(int i, int j,Rectangle *rectangle,int t);
    void FCTTimeStep(double timestep,int step,int subStep);
    double RGKGetFluxPL(int i, int j);
    double RGKGetFluxXL(int i, int j);
    inline int Index6(int i, int j, int state);
    inline int Index3(int i, int j, int state);
    void CalculateDifferentBoundaryC();
    double GetMagneticForce(int i);
    double GetCellAverageASquared(int i);
    inline double RGKGetFluxPLN(int i,int j);
    inline double RGKGetFluxXLN(int i,int j);
    inline double RGKGetFluxPN(int i,int j);
    inline double RGKGetFluxXN(int i,int j);
    void CalculateSameBoundaryC();
    void SetCFromSameLevel(int i, int j,Rectangle *rectangle,int t);
    void UpdateCFromSameLevel(int i, int j,Rectangle *rectangle,int t);
    void UpdateSameBoundaryC();
};

inline int Rectangle::GetFinestIndex(int i) {
    return relativeToBottom * (i + x_pos);
}

inline double Rectangle::Momentum(double i) {
    return settings.pmin[particleType] + dp * (i + p_pos);
}

inline int Rectangle::Index(int i, int j, int state) {
    return 8*((n_p + 4) * (i + 2) + 2 + j) + state;
}

inline int Rectangle::Index3(int i, int j, int state) {
    return 3*((n_p + 4) * (i + 2) + 2 + j) + state;
}

inline int Rectangle::Index6(int i, int j, int state) {
    return 6*((n_p + 4) * (i + 2) + 2 + j) + state;
}

inline int Rectangle::IndexNS(int i, int j) {
    return (n_p + 4) * (i + 2) + 2 + j;
}

inline double Rectangle::ErrorEstimate(int i,int j) {
    return errorWeights[0]*std::fabs(f[Index3(i+1,j,1)]-f[Index3(i-1,j,1)])+errorWeights[1]*std::fabs(f[Index3(i,j+1,1)]-f[Index3(i,j-1,1)])+errorWeights[2]*std::fabs(f[Index3(i+1,j,1)]-2*f[Index3(i,j,1)]+f[Index3(i-1,j,1)])+errorWeights[3]*std::fabs(f[Index3(i,j+1,1)]-2*f[Index3(i,j,1)]+f[Index3(i,j-1,1)])+errorWeights[4]*std::fabs(f[Index3(i,j,1)]);
}

inline double Rectangle::CalculateFluxToCoarseP(int i,int j) {
    i=i*refinementRatio-x_pos;
    j=j*refinementRatio-p_pos;
    double temp = 0.0;

    for (int k=0; k<refinementRatio; k++) {
        temp+=RGKGetFluxP(i+k,j);
    }
    temp*=(1.0/(refinementRatio*refinementRatio));
    return temp;
}

inline double Rectangle::CalculateFluxToCoarseX(int i,int j) {
    i=i*refinementRatio-x_pos;
    j=j*refinementRatio-p_pos;
    double temp = 0.0;

    for (int k=0; k<refinementRatio; k++) {
        temp+=RGKGetFluxX(i,j+k);
    }
    temp*=(1.0/(refinementRatio*refinementRatio));
    return temp;
}

inline double Rectangle::CalculateFluxToCoarsePL(int i,int j) {
    i=i*refinementRatio-x_pos;
    j=j*refinementRatio-p_pos;
    double temp = 0.0;

    for (int k=0; k<refinementRatio; k++) {
        temp+=RGKGetFluxPL(i+k,j);
    }
    temp*=(1.0/(refinementRatio*refinementRatio));
    return temp;
}

inline double Rectangle::CalculateFluxToCoarseXL(int i,int j) {
    i=i*refinementRatio-x_pos;
    j=j*refinementRatio-p_pos;
    double temp = 0.0;

    for (int k=0; k<refinementRatio; k++) {
        temp+=RGKGetFluxXL(i,j+k);
    }
    temp*=(1.0/(refinementRatio*refinementRatio));
    return temp;
}

inline double Rectangle::Gamma(double p_x,double a_squared) {
    return std::sqrt(1+((p_x*p_x)+a_squared)*((m_inv*c_inv)*(m_inv*c_inv)));
}

inline double Rectangle::RGKGetFluxPLN(int i,int j){
    i=i*refinementRatio-x_pos;
    j=j*refinementRatio-p_pos;

    return RGKGetFluxPL(i,j);
}

inline double Rectangle::RGKGetFluxPN(int i,int j){
    i=i*refinementRatio-x_pos;
    j=j*refinementRatio-p_pos;

    return RGKGetFluxP(i,j);
}

inline double Rectangle::RGKGetFluxXLN(int i,int j){
    i=i*refinementRatio-x_pos;
    j=j*refinementRatio-p_pos;

    return RGKGetFluxXL(i,j);

}

inline double Rectangle::RGKGetFluxXN(int i,int j){
    i=i*refinementRatio-x_pos;
    j=j*refinementRatio-p_pos;

    return RGKGetFluxX(i,j);

}

#endif /* __Rectangle_hpp__ */
