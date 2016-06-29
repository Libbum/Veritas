#ifndef __EMSolver_hpp__
#define __EMSolver_hpp__

#if !defined(USINGMKL)
extern "C" int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern "C" int dgetrs_(char *trans, int *n, int *m, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
#endif

class EMFieldSolver {
    Settings &settings;
    #if defined(USINGMKL)
    lapack_int *IPIV;
    #else
    int *IPIV;
    #endif
    int n_prepad, n_postpad;
    unsigned int x_size;
    std::vector<std::shared_ptr<Mesh>> meshes;
    std::vector<std::unique_ptr<std::ofstream>> energyStreams;
    std::ofstream *chargeStream, *ELongStream, *ETransStream, *potentialStream, *BStream, *AsqStream, *timeStream;
    std::vector<double> charge, J, By, Bz, Ey, Ez, Ay, Az, a_squared, neutralizationCharge;
    std::vector<std::vector<double>> charges, energies;
    double *PHI, *pM, fieldCoef, Ex0;
public:
    EMFieldSolver(Settings &settings, const std::vector<std::shared_ptr<Mesh> > &meshes);
    ~EMFieldSolver();
    void AssembleRhoAndJ();
    void AssembleEnergy();
    double GetASquared(int i);
    double GetEfield(int i);
    void DumpCharge();
    void DumpEnergy();
    void DumpEFieldLongitudinal();
    void DumpEFieldTransverse();
    void DumpPotential();
    void DumpBFieldTransverse();
    void DumpAsqField();
    void DumpTime(double time);
    void UpdatePotential();
    void RGKStep(int step, double timestep);
    void RGKCalculateRHS(int step);
    void RGKUpdateIntermediateSolution(int step, double timestep);
    inline int Index(int i, int step);
    void InterpolateToFaces();
    void EnforceChargeNeutralization();
    __attribute__((vector, nothrow)) double GetCellAverageASquared(int i);
    double EstimateCFLBound();
    double GetMagneticForce(int i);
};

inline int EMFieldSolver::Index(int i, int step) {
    return step * (x_size + n_prepad + n_postpad) + i;
}

#endif /* __EMSolver_hpp__ */
