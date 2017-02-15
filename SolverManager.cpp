#include "veritas.hpp"
#include "SolverManager.hpp"
#include "Settings.hpp"
#include "EMSolver.hpp"
#include "Mesh.hpp"
#include "Level.hpp"

SolverManager::SolverManager(Settings &settings) : settings(settings) {
    std::cout << std::setfill('-') << std::setw(3) << '-' << " Time (T) " << std::setw(4) << '-' << '|' << std::setw(3) << '-' << " Particle " << std::setw(3) << '-' << '|' << std::setw(4) << '-' << " Levels " << std::setw(4) << '-' << '|' << std::setw(8) << '-' << " Rectangles " << std::setw(9) << '-' << std::endl;
    int numSpecies = settings.p_size.size();

    for (int i = 0; i < numSpecies; i++) {
        meshes.push_back(std::make_shared<Mesh>(i, settings));
    }

    EMSolver = std::make_shared<EMFieldSolver>(settings, meshes);
    std::cout << std::setfill(' ') << std::setw(4) << ' ' << std::fixed << 0.0 << std::setw(5) << ' '; //NOTE: Time is hardcoded here

    for (auto & mesh: meshes) {
        mesh->SetFieldSolver(EMSolver);
        mesh->PushData();
        screenOutput(mesh);
    }

    EMSolver->EnforceChargeNeutralization();
}

void SolverManager::Advance(double timeStep) {
    for (int i = 0; i < 6; i++) {  //RK Steps
        EMSolver->AssembleRhoAndJ();
        EMSolver->UpdatePotential();

        for (unsigned int j = 0; j < settings.q.size(); j++) {
            meshes[j]->Advance(timeStep, i);
        }
        settings.UpdateTime(i, timeStep);
        EMSolver->RGKStep(i, timeStep);
    }
}

void SolverManager::AdvanceFields(double timeStep) {
    for (int i = 0; i < 6; i++) {  //RK Steps
        settings.UpdateTime(i, timeStep);
        EMSolver->RGKStep(i, timeStep);
    }
}

void SolverManager::reGrid(double t) {
    double T = settings.tempEM[0] / cs;

    std::cout << std::setfill(' ') << std::setw(4) << ' ' << std::fixed << std::setprecision(5) << t/T << std::setw(5) << ' ';

    for (auto & mesh: meshes) {
        mesh->updateHierarchy();
        screenOutput(mesh);
    }
}

std::string SolverManager::centeredOutput(std::string const& original, int targetSize) {
    int padding = targetSize - original.size();

    return padding > 0 ? std::string(padding / 2, ' ') + original + std::string(padding / 2, ' ') : original;
}

void SolverManager::OutputRectangles(double t) {
    if (settings.output.rectangleData) {
        for (auto & mesh: meshes) {
            mesh->outputRectangleData(t);
        }
    }
}

void SolverManager::screenOutput(const std::shared_ptr<Mesh> &mesh) {
    std::stringstream lvlDetails;
    std::stringstream pt;
    std::stringstream rectNums;

    if (mesh->particleType > 0) {
        std::cout << std::setw(17) << ' ';
    }

    pt << mesh->particleType;
    std::cout << std::setw(16) << centeredOutput(pt.str(), 15) << std::setw(1) << ' ';

    int full = 0;

    for (auto& lvl: mesh->levels) {
        if (!lvl->rectangles.empty()) {
            full++;
            rectNums << lvl->rectangles.size() << ", ";
        }
    }

    std::string rectNumss = rectNums.str();
    rectNumss.erase(rectNumss.size() - 2, 2);

    lvlDetails << full << " (" << mesh->levels.size() << ")";
    std::cout << std::setw(16) << centeredOutput(lvlDetails.str(), 15) << std::setw(1) << ' ';
    std::cout << std::setw(30) << centeredOutput(rectNumss, 29) << std::endl;
}

void SolverManager::fileOutput(double t) {
    //Write all file contents
   //rectangleData is handled separately
   if (settings.output.energy) {
       EMSolver->AssembleEnergy();
   }

   #pragma omp parallel sections
   {
       #pragma omp section
       {
           if (settings.output.charge) {
               EMSolver->DumpCharge();
           }
       }
       #pragma omp section
       {
           if (settings.output.energy) {
               EMSolver->DumpEnergy();
           }
       }
       #pragma omp section
       {
           if (settings.output.potential) {
               EMSolver->DumpPotential();
           }
       }
       #pragma omp section
       {
           if (settings.output.EFieldLongitudinal) {
               EMSolver->DumpEFieldLongitudinal();
           }
       }
       #pragma omp section
       {
           if (settings.output.EFieldTransverse) {
               EMSolver->DumpEFieldTransverse();
           }
       }
       #pragma omp section
       {
           if (settings.output.BFieldTransverse) {
               EMSolver->DumpBFieldTransverse();
           }
       }
       #pragma omp section
       {
           if (settings.output.AFieldSquared) {
               EMSolver->DumpAsqField();
           }
       }
       #pragma omp section
       {
           if (settings.output.time) {
               EMSolver->DumpTime(t);
           }
       }
   }
}

double SolverManager::CalculateDt(double cfl) {
    return cfl * EMSolver->EstimateCFLBound();
}
