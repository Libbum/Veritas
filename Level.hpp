#ifndef __Level_hpp__
#define __Level_hpp__

class Level {
    std::vector<double> chargeL, currentL, energyL;
    int x_size, p_size, depth, particleType;
    Settings &settings;
    EMFieldSolver *EMSolver;
public:
    std::vector<std::shared_ptr<Rectangle>> rectangles;
    Level(int particleType, int depth, Settings &settings);
    void PushData(int updateType,int val=1);
    void CollectRhoAndJ();
    void CollectEnergy();
    void InterpolateRhoAndJToFinestMesh(std::vector<double> &charge, std::vector<double> &J);
    void InterpolateEnergyToFinestMesh(std::vector<double> &energy);
    void SetUpVlasovPoisson(int i, bool option = false);
    void SetFieldSolver(EMFieldSolver *solver);
    void GetDataFromSameLevel(const std::unique_ptr<Level> &level);
    void GetDataFromCoarserLevel(const std::unique_ptr<Level> &level);
    void GetDataFromCoarseNewLevel(const std::unique_ptr<Level> &level);
    void FCTTimeStep(double timestep,int step,int subStep);
};

#endif /* __Level_hpp__ */
