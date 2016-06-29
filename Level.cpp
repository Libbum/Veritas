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
#pragma omp parallel for schedule(dynamic,1) num_threads(settings.numThreads)
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles[i]->FCTTimeStep(timestep,step,subStep);
    }
}

void Level::InterpolateRhoAndJToFinestMesh(std::vector<double> &charge, std::vector<double> &J) {
    CollectRhoAndJ();

    for (unsigned int i = 0; i < settings.x_size_finest; i++) {
        charge[i] += chargeL[i];
        J[i] += currentL[i];
    }
}

void Level::InterpolateEnergyToFinestMesh(std::vector<double> &energy) {
    CollectEnergy();

    for (unsigned int i = 0; i < settings.p_size_finest[particleType]; i++) {
        energy[i] += energyL[i];
    }
}

void Level::CollectRhoAndJ() {
    #pragma omp parallel for schedule(dynamic,1) num_threads(settings.numThreads)
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles[i]->CalculateRhoAndJ();
    }

    std::fill(chargeL.begin(), chargeL.end(), 0.0);
    std::fill(currentL.begin(), currentL.end(), 0.0);

    for (auto &r: rectangles) {
        int shift = r->x_pos;

        for (unsigned int j = 0; j < r->chargeR.size(); j++) {
            chargeL.at(shift * r->relativeToBottom + j) += r->chargeR.at(j);
            currentL.at(shift * r->relativeToBottom + j) += r->currentR.at(j);
        }
    }
}

void Level::CollectEnergy() {
    #pragma omp parallel for schedule(dynamic,1) num_threads(settings.numThreads)
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles[i]->CalculateEnergy();
    }

    std::fill(energyL.begin(), energyL.end(), 0.0);

    for (auto &r: rectangles) {
        int shift = r->p_pos;
        int rtb = r->relativeToBottom;
        for (unsigned int i = 0; i < r->energyR.size(); i++) {
            energyL.at(shift * rtb + i) += r->energyR.at(i);
        }
    }
}

void Level::SetFieldSolver(EMFieldSolver *solver) {
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
