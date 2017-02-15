#ifndef __Mesh_hpp__
#define __Mesh_hpp__

class Mesh {
    std::vector<level> hierarchy;
    Settings &settings;
    std::shared_ptr<EMFieldSolver> EMSolver;
    std::shared_ptr<Rectangle> bc;
public:
    int particleType;
    std::vector<std::unique_ptr<Level>> levels;
    Mesh(int particleType, Settings &settings);
    void PushData(int val=1);
    void Advance(double timeStep, int step);
    void InterpolateRhoAndJToFinestMesh(std::vector<double> &charge, std::vector<double> &J);
    void InterpolateEnergyToFinestMesh(std::vector<double> &energy);
    void SetUpVlasovPoisson(int i);
    void SetFieldSolver(const std::shared_ptr<EMFieldSolver> &solver);
    void updateHierarchy(bool init=false);
    bool static choose_second(const coords &lhs, const coords &rhs);
    bool static choose_first(const coords &lhs, const coords &rhs);
    int sgn(int &x);
    void printRectangle(rect &r);
    void printRectangle(rect &r, int level);
    int countCells(const rect &span);
    void isFlaggedInside(std::vector<coords> &flagged, std::vector<bool> &isInside, rect &r);
    void getSigs(rect &rectangle, std::vector<coords> &flagged, std::vector<int> &sigX, std::vector<int> &sigP);
    std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int> > computeSignatures(rect &rectangle, std::vector<coords> &flagged);
    coords identifyInflection(rect &rectangle, std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int> > &signatures);
    std::tuple<bool, int, int> hasHole(std::vector<int> &sig);
    level splitRectangle(rect &rectangle, std::vector<coords> &flagged, const double &minEfficiency);
    void getError(const int &lvl, bool init, std::vector<coords> &flaggedCells);
    void interpRectanglesUp(level &identified, const int &lvl);
    void mergeDownFlaggedData(const int &lvl, const rect &r, std::vector<coords> &foundCells);
    void getExtrema(rect &extrema, const std::vector<coords> &flaggedCells);
    void outputRectangleData(double tidx);
    inline double getNe(double xp,double tempp0,double plasma_xl_bound,double plasma_xr_bound);
    void promoteHierarchyToMesh(bool init);
    void InterMeshDataTransfer(const std::vector<std::unique_ptr<Level>> &levels_n);
    void PushBoundaryC();
};

#pragma omp declare simd simdlen(8)
inline double Mesh::getNe(double xp,double tempp0,double plasma_xl_bound,double plasma_xr_bound) {
    double ne=0.0;
    if ((xp > plasma_xl_bound)&&(xp < plasma_xr_bound)) {
        ne =  tempp0;
    }
    return ne;
}

#endif /* __Mesh_hpp__ */
