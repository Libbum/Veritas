#ifndef __solvermanager_hpp__
#define __solvermanager_hpp__
#include "EMSolver.hpp"

class SolverManager {
    std::shared_ptr<EMFieldSolver> EMSolver;
    Settings &settings;
    std::vector<std::shared_ptr<Mesh>> meshes;
public:
    SolverManager(Settings &settings);
    void Advance(double timeStep);
    void AdvanceFields(double timeStep);
    void reGrid(double t);
    std::string centeredOutput(std::string const& original, int targetSize);
    void screenOutput(const std::shared_ptr<Mesh> &mesh);
    void OutputRectangles(double t);
    void fileOutput(double t);
    double CalculateDt(double cfl);
};

#endif /* __solvermanager_hpp__ */
