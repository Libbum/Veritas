#ifndef __solvermanager_hpp__
#define __solvermanager_hpp__

class SolverManager {
    EMFieldSolver *EMSolver;
    Settings &settings;
    std::vector<std::shared_ptr<Mesh>> meshes;
public:
    SolverManager(Settings &settings);
    ~SolverManager();
    void Advance(double timeStep);
    void AdvanceFields(double timeStep);
    void reGrid(double time);
    std::string centeredOutput(std::string const& original, int targetSize);
    void screenOutput(const std::shared_ptr<Mesh> &mesh);
    void fileOutput(double time);
    double CalculateDt(double cfl);
};

#endif /* __solvermanager_hpp__ */
