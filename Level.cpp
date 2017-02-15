#include "veritas.hpp"
#include "Settings.hpp"
#include "Level.hpp"
#include "Rectangle.hpp"

Level::Level(int particleType, int depth, Settings &settings) : x_size(settings.GetXSize(settings.maxDepth - depth)), p_size(settings.GetPSize(settings.maxDepth - depth, particleType)), depth(depth), particleType(particleType), settings(settings), EMSolver(NULL) {
    chargeL.resize(settings.x_size_finest, 0.0);
    energyL.resize(settings.p_size_finest[particleType], 0.0);
    currentL.resize(settings.x_size_finest, 0.0);
}

void Level::FCTTimeStep(double timestep,int step,int subStep) {
    //pragma omp parallel for schedule(guided)
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles[i]->FCTTimeStep(timestep,step,subStep);
    }
}

void Level::InterpolateRhoAndJToFinestMesh(std::vector<double> &charge, std::vector<double> &J) {
    CollectRhoAndJ();

    unsigned int iloop = settings.x_size_finest;

    #pragma ivdep
    for (unsigned int i = 0; i < iloop; i++) {
        charge[i] += chargeL[i];
        J[i] += currentL[i];
    }
}

void Level::InterpolateEnergyToFinestMesh(std::vector<double> &energy) {
    CollectEnergy();

    unsigned int iloop = settings.p_size_finest[particleType];

    #pragma ivdep
    for (unsigned int i = 0; i < iloop; i++) {
        energy[i] += energyL[i];
    }
}

void Level::CollectRhoAndJ() {
    unsigned int iloop = rectangles.size();
    for (unsigned int i = 0; i < iloop; i++) {
        rectangles[i]->CalculateRhoAndJ();
    }

    std::fill(chargeL.begin(), chargeL.end(), 0.0);
    std::fill(currentL.begin(), currentL.end(), 0.0);

    //no omp
    for (unsigned int i = 0; i < iloop; i++) {
        int shift = rectangles[i]->x_pos;
        unsigned int jloop = rectangles[i]->chargeR.size();
        #pragma ivdep
        for (unsigned int j = 0; j < jloop; j++) {
            int rtb = (int)rectangles[i]->relativeToBottom;
            chargeL[shift * rtb + j] += rectangles[i]->chargeR[j];
            currentL[shift * rtb + j] += rectangles[i]->currentR[j];
        }
    }
}

void Level::CollectEnergy() {
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles[i]->CalculateEnergy();
    }

    std::fill(energyL.begin(), energyL.end(), 0.0);

    for (auto &r: rectangles) {
        int shift = r->p_pos;
        int rtb = r->relativeToBottom;
        for (unsigned int i = 0; i < r->energyR.size(); i++) {
            energyL[shift * rtb + i] += r->energyR[i];
        }
    }
}

void Level::SetFieldSolver(const std::shared_ptr<EMFieldSolver> &solver) {
    EMSolver = solver;

    for (auto &r: rectangles) {
        r->SetFieldSolver(EMSolver);
    }
}

void Level::PushData(int updateType,int val) {
    switch (updateType) {
        case 0:
            for (auto &r: rectangles) {
                r->UpdateInterriorPoints(val);
            }
            break;
        case 1:
            for (auto &r: rectangles) {
                r->UpdateSameLevelBoundaries(val);
            }
            break;
        case 2:
            for (auto &r: rectangles) {
                r->UpdateDifferentLevelBoundaries(val);
            }
            break;
        case 3:
            for (auto &r: rectangles) {
                r->UpdateCornerPoints(val);
            }
            break;
        case 4:
            for (auto &r: rectangles) {
                r->CalculateSameBoundaryC();
            }
            break;
        case 5:
            for (auto &r: rectangles) {
                r->CalculateDifferentBoundaryC();
            }
            break;
        case 6:
            for (auto &r: rectangles) {
                r->UpdateSameBoundaryC();
            }
            break;
    }
}

void Level::GetDataFromSameLevel(const std::unique_ptr<Level> &level) {
    for (unsigned int i = 0; i < level->rectangles.size(); i++) {
        for (unsigned int j = 0; j < rectangles.size(); j++) {
            level->rectangles[i]->GetDataFromSameLevelRectangle(rectangles[j]);
        }
    }
}

void Level::GetDataFromCoarserLevel(const std::unique_ptr<Level> &level) {
    for (unsigned int i = 0; i < level->rectangles.size(); i++) {
        for (unsigned int j = 0; j < rectangles.size(); j++) {
            level->rectangles[i]->GetDataFromCoarseLevelRectangle(rectangles[j]);
        }
    }
}

void Level::GetDataFromCoarseNewLevel(const std::unique_ptr<Level> &level) {
    for (unsigned int i = 0; i < level->rectangles.size(); i++) {
        for (unsigned int j = 0; j < rectangles.size(); j++) {
            level->rectangles[i]->GetDataFromCoarseNewLevelRectangle(rectangles[j]);
        }
    }
}
