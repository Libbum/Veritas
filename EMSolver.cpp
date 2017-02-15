#include "veritas.hpp"
#include "EMSolver.hpp"
#include "Settings.hpp"
#include "Mesh.hpp"

EMFieldSolver::EMFieldSolver(Settings &settings, const std::vector<std::shared_ptr<Mesh> > &meshes) : settings(settings), x_size(settings.x_size_finest), meshes(meshes), chargeStream(NULL), ELongStream(NULL), ETransStream(NULL), potentialStream(NULL), BStream(NULL), AsqStream(NULL), timeStream(NULL), Ex0(0.0) {
    n_prepad = std::max(std::floor(settings.preLength / settings.GetDx(0)), 2.0);
    n_postpad = std::max(std::floor(settings.postLength / settings.GetDx(0)), 2.0);

    charge.resize(x_size, 0);
    J.resize(x_size, 0);
    neutralizationCharge.resize(x_size, 0.0);

    charges.resize(settings.q.size(),charge);
    for (unsigned int i=0; i<settings.p_size.size();i++) {
        energies.push_back({0.0});
        energies[i].resize(settings.p_size_finest[i],0.0);
    }

    By.resize((x_size + n_prepad + n_postpad) * 8, 0.0);
    Bz.resize((x_size + n_prepad + n_postpad) * 8, 0.0);
    Ey.resize((x_size + n_prepad + n_postpad) * 8, 0.0);
    Ez.resize((x_size + n_prepad + n_postpad) * 8, 0.0);
    Ay.resize((x_size + n_prepad + n_postpad) * 8, 0.0);
    Az.resize((x_size + n_prepad + n_postpad) * 8, 0.0);
    a_squared.resize(x_size + 1, 0.0);

    double w1(-16.0 / 30.0 / 0.4), w2(1.0 / 30.0 / 0.4), w3(1 / 0.4);

    PHI = new double[x_size];
    pM = new double[x_size * x_size];

    for (unsigned int i = 0; i < x_size; i++) {
        PHI[i] = 0.0;
    }

    for (unsigned int i = 0; i < x_size * x_size; i++) {
        pM[i] = 0.0;
    }

    for (unsigned int i = 2; i < x_size - 2; i++) {
        pM[i * x_size + i] = w3;
        pM[i * x_size + i + 1] = w1;
        pM[i * x_size + i - 1] = w1;
        pM[i * x_size + i + 2] = w2;
        pM[i * x_size + i - 2] = w2;
    }

    pM[0] = 1.0;

    pM[x_size + 1] = w3;
    pM[x_size + 2] = w1;
    pM[x_size + 0] = w1;
    pM[x_size + 3] = w2;
    pM[2 * x_size - 1] = w2;

    pM[(x_size - 2) * x_size + x_size - 2] = w3;
    pM[(x_size - 2) * x_size + x_size - 1] = w1;
    pM[(x_size - 2) * x_size + x_size - 3] = w1;
    pM[(x_size - 2) * x_size + 0] = w2;
    pM[(x_size - 2) * x_size + x_size - 4] = w2;

    pM[(x_size - 1) * x_size + (x_size - 1)] = w3;
    pM[(x_size - 1) * x_size + 0] = w1;
    pM[(x_size - 1) * x_size + (x_size - 2)] = w1;
    pM[(x_size - 1) * x_size + 1] = w2;
    pM[(x_size - 1) * x_size + (x_size - 3)] = w2;

    #if defined(USINGMKL)
    lapack_int INFO = 3;
    lapack_int LDA = x_size;
    lapack_int N = x_size;
    IPIV = new lapack_int[x_size];

    INFO = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, pM, LDA, IPIV);

    if (INFO != 0) std::cerr << "WARNING: Illegal parameter or singular result from dgetrf" << std::endl;

    #else
    int INFO = 3;
    int LDA = x_size;
    int N = x_size;
    IPIV = new int[x_size];

    dgetrf_(&N, &N, pM, &LDA, IPIV, &INFO);
    #endif

    fieldCoef = (1.0 / (12 * settings.GetDx(0)));
}

EMFieldSolver::~EMFieldSolver() {
    delete[] PHI;
    delete[] pM;
    delete[] IPIV;
    delete chargeStream;
    delete ELongStream;
    delete ETransStream;
    delete potentialStream;
    delete BStream;
    delete AsqStream;
    delete timeStream;
}

void EMFieldSolver::AssembleRhoAndJ() {
    std::fill(charge.begin(), charge.end(), 0.0);
    std::fill(J.begin(), J.end(), 0.0);

    for (unsigned int i=0; i<meshes.size();i++) {
        std::fill(charges[i].begin(), charges[i].end(),0.0);
    }

    for (unsigned int i=0; i<meshes.size();i++) {
        meshes[i]->InterpolateRhoAndJToFinestMesh(charges[i], J);
    }

    for (unsigned int i=0; i<meshes.size();i++) {
        #pragma ivdep
        for (unsigned int j=0; j<x_size; j++) {
            charge[j]+=charges[i][j];
        }
    }
}

void EMFieldSolver::AssembleEnergy() {
    for (unsigned int i=0; i<meshes.size();i++) {
        std::fill(energies[i].begin(), energies[i].end(),0.0);
    }
    for (unsigned int i=0; i<meshes.size();i++) {
        meshes[i]->InterpolateEnergyToFinestMesh(energies[i]);
    }
}

double EMFieldSolver::GetASquared(int i) {
    return a_squared[std::min(std::max(i, 0), (int)x_size)];
}

double EMFieldSolver::GetEfield(int i) {
    int ip1 = i + 1;
    int im1 = i - 1;
    int ip2 = i + 2;
    int im2 = i - 2;

    ip1 = (ip1 > -1) ? ip1 : (ip1 + x_size);
    im1 = (im1 > -1) ? im1 : (im1 + x_size);
    ip2 = (ip2 > -1) ? ip2 : (ip2 + x_size);
    im2 = (im2 > -1) ? im2 : (im2 + x_size);

    ip1 = (ip1 < (int)x_size) ? ip1 : (ip1 - x_size);
    im1 = (im1 < (int)x_size) ? im1 : (im1 - x_size);
    ip2 = (ip2 < (int)x_size) ? ip2 : (ip2 - x_size);
    im2 = (im2 < (int)x_size) ? im2 : (im2 - x_size);

    return -fieldCoef * (8 * (PHI[ip1] - PHI[im1]) - PHI[ip2] + PHI[im2]) + Ex0;
}

void EMFieldSolver::UpdatePotential() {
    double te = eps0_inv;

    for (unsigned int i = 0; i < x_size; i++) {
        PHI[i] = te * (charge[i] + neutralizationCharge[i]);
    }

    double w = std::pow(settings.GetDx(0), 2.0);

    for (unsigned int i = 0; i < x_size; i++) {
        PHI[i] *= w;
    }

    char TRANS = 'N';
    #if defined(USINGMKL)
    lapack_int INFO = 3;
    lapack_int LDA = x_size;
    lapack_int LDB = x_size;
    lapack_int N = x_size;
    lapack_int NRHS = 1;

    INFO = LAPACKE_dgetrs(LAPACK_COL_MAJOR, TRANS, N, NRHS, pM, LDA, IPIV, PHI, LDB);

    if (INFO != 0) std::cerr << "WARNING: Illegal parameter or singular result from dgetrs" << std::endl;

    #else
    int INFO = 3;
    int LDA = x_size;
    int LDB = x_size;
    int N = x_size;
    int NRHS = 1;

    dgetrs_(&TRANS, &N, &NRHS, pM, &LDA, IPIV, PHI, &LDB, &INFO);
    #endif

    Ex0 += -(GetEfield(-1) + GetEfield(0)) * 0.5;
}

void EMFieldSolver::RGKStep(int step, double timestep) {

    RGKCalculateRHS(step);

    RGKUpdateIntermediateSolution(step, timestep);

    InterpolateToFaces();

}

void EMFieldSolver::RGKUpdateIntermediateSolution(int step, double timestep) {

    unsigned int iloop = x_size + n_prepad + n_postpad;

    if (step == -1) {
    } else if (step == 0) {
        double a0 = 0.5 * timestep;
        #pragma ivdep
        for (unsigned int i = 0; i < iloop; i++) {
            int idx0 = Index(i, 0);
            int idx1 = Index(i, 1);
            int idx2 = Index(i, 2);
            By[idx1] = By[idx0] + a0 * By[idx2];
            Bz[idx1] = Bz[idx0] + a0 * Bz[idx2];
            Ey[idx1] = Ey[idx0] + a0 * Ey[idx2];
            Ez[idx1] = Ez[idx0] + a0 * Ez[idx2];
            Ay[idx1] = Ay[idx0] + a0 * Ay[idx2];
            Az[idx1] = Az[idx0] + a0 * Az[idx2];
        }
    } else if (step == 1) {
        double a0 = 0.221776 * timestep;
        double a1 = 0.110224 * timestep;
        #pragma ivdep
        for (unsigned int i = 0; i < iloop; i++) {
            int idx0 = Index(i, 0);
            int idx1 = Index(i, 1);
            int idx2 = Index(i, 2);
            int idx3 = Index(i, 3);
            By[idx1] = By[idx0] + a0 * By[idx2] + a1 * By[idx3];
            Bz[idx1] = Bz[idx0] + a0 * Bz[idx2] + a1 * Bz[idx3];
            Ey[idx1] = Ey[idx0] + a0 * Ey[idx2] + a1 * Ey[idx3];
            Ez[idx1] = Ez[idx0] + a0 * Ez[idx2] + a1 * Ez[idx3];
            Ay[idx1] = Ay[idx0] + a0 * Ay[idx2] + a1 * Ay[idx3];
            Az[idx1] = Az[idx0] + a0 * Az[idx2] + a1 * Az[idx3];
        }
    } else if (step == 2) {
        double a0 = -0.04884659515311857 * timestep;
        double a1 = -0.17772065232640102 * timestep;
        double a2 = 0.8465672474795197 * timestep;
        #pragma ivdep
        for (unsigned int i = 0; i < iloop; i++) {
            int idx0 = Index(i, 0);
            int idx1 = Index(i, 1);
            int idx2 = Index(i, 2);
            int idx3 = Index(i, 3);
            int idx4 = Index(i, 4);
            By[idx1] = By[idx0] + a0 * By[idx2] + a1 * By[idx3] + a2 * By[idx4];
            Bz[idx1] = Bz[idx0] + a0 * Bz[idx2] + a1 * Bz[idx3] + a2 * Bz[idx4];
            Ey[idx1] = Ey[idx0] + a0 * Ey[idx2] + a1 * Ey[idx3] + a2 * Ey[idx4];
            Ez[idx1] = Ez[idx0] + a0 * Ez[idx2] + a1 * Ez[idx3] + a2 * Ez[idx4];
            Ay[idx1] = Ay[idx0] + a0 * Ay[idx2] + a1 * Ay[idx3] + a2 * Ay[idx4];
            Az[idx1] = Az[idx0] + a0 * Az[idx2] + a1 * Az[idx3] + a2 * Az[idx4];
        }
    } else if (step == 3) {
        double a0 = -0.15541685842491548 * timestep;
        double a1 = -0.3567050098221991 * timestep;
        double a2 = 1.0587258798684427 * timestep;
        double a3 = 0.30339598837867193 * timestep;
        #pragma ivdep
        for (unsigned int i = 0; i < iloop; i++) {
            int idx0 = Index(i, 0);
            int idx1 = Index(i, 1);
            int idx2 = Index(i, 2);
            int idx3 = Index(i, 3);
            int idx4 = Index(i, 4);
            int idx5 = Index(i, 5);
            By[idx1] = By[idx0] + a0 * By[idx2] + a1 * By[idx3] + a2 * By[idx4] + a3 * By[idx5];
            Bz[idx1] = Bz[idx0] + a0 * Bz[idx2] + a1 * Bz[idx3] + a2 * Bz[idx4] + a3 * Bz[idx5];
            Ey[idx1] = Ey[idx0] + a0 * Ey[idx2] + a1 * Ey[idx3] + a2 * Ey[idx4] + a3 * Ey[idx5];
            Ez[idx1] = Ez[idx0] + a0 * Ez[idx2] + a1 * Ez[idx3] + a2 * Ez[idx4] + a3 * Ez[idx5];
            Ay[idx1] = Ay[idx0] + a0 * Ay[idx2] + a1 * Ay[idx3] + a2 * Ay[idx4] + a3 * Ay[idx5];
            Az[idx1] = Az[idx0] + a0 * Az[idx2] + a1 * Az[idx3] + a2 * Az[idx4] + a3 * Az[idx5];
        }
    } else if (step == 4) {
        double a0 = 0.2014243506726763 * timestep;
        double a1 = 0.008742057842904185 * timestep;
        double a2 = 0.15993995707168115 * timestep;
        double a3 = 0.4038290605220775 * timestep;
        double a4 = 0.22606457389066084 * timestep;
        #pragma ivdep
        for (unsigned int i = 0; i < iloop; i++) {
            int idx0 = Index(i, 0);
            int idx1 = Index(i, 1);
            int idx2 = Index(i, 2);
            int idx3 = Index(i, 3);
            int idx4 = Index(i, 4);
            int idx5 = Index(i, 5);
            int idx6 = Index(i, 6);
            By[idx1] = By[idx0] + a0 * By[idx2] + a1 * By[idx3] + a2 * By[idx4] + a3 * By[idx5] + a4 * By[idx6];
            Bz[idx1] = Bz[idx0] + a0 * Bz[idx2] + a1 * Bz[idx3] + a2 * Bz[idx4] + a3 * Bz[idx5] + a4 * Bz[idx6];
            Ey[idx1] = Ey[idx0] + a0 * Ey[idx2] + a1 * Ey[idx3] + a2 * Ey[idx4] + a3 * Ey[idx5] + a4 * Ey[idx6];
            Ez[idx1] = Ez[idx0] + a0 * Ez[idx2] + a1 * Ez[idx3] + a2 * Ez[idx4] + a3 * Ez[idx5] + a4 * Ez[idx6];
            Ay[idx1] = Ay[idx0] + a0 * Ay[idx2] + a1 * Ay[idx3] + a2 * Ay[idx4] + a3 * Ay[idx5] + a4 * Ay[idx6];
            Az[idx1] = Az[idx0] + a0 * Az[idx2] + a1 * Az[idx3] + a2 * Az[idx4] + a3 * Az[idx5] + a4 * Az[idx6];
        }
    } else if (step==5) {
        double a0 = 0.2014243506726763 * timestep;
        double a1 = 0.008742057842904185 * timestep;
        double a2 = 0.15993995707168115 * timestep;
        double a3 = 0.4038290605220775 * timestep;
        double a4 = 0.22606457389066084 * timestep;
        double a5 = 0;

        double b0 = 0.15791629516167136 * timestep - a0;
        double b1 = 0.0 * timestep - a1;
        double b2 = 0.18675894052400077 * timestep - a2;
        double b3 = 0.6805652953093346 * timestep - a3;
        double b4 = -0.27524053099500667 * timestep - a4;
        double b5 = 0.25 * timestep - a5;
        #pragma ivdep
        for (unsigned int i = 0; i < iloop; i++) {
            int idx0 = Index(i, 0);
            int idx1 = Index(i, 1);
            int idx2 = Index(i, 2);
            int idx3 = Index(i, 3);
            int idx4 = Index(i, 4);
            int idx5 = Index(i, 5);
            int idx6 = Index(i, 6);
            int idx7 = Index(i, 7);
            By[idx0] = By[idx1] + b0 * By[idx2] + b1 * By[idx3] + b2 * By[idx4] + b3 * By[idx5] + b4 * By[idx6] + b5 * By[idx7];
            Bz[idx0] = Bz[idx1] + b0 * Bz[idx2] + b1 * Bz[idx3] + b2 * Bz[idx4] + b3 * Bz[idx5] + b4 * Bz[idx6] + b5 * Bz[idx7];
            Ey[idx0] = Ey[idx1] + b0 * Ey[idx2] + b1 * Ey[idx3] + b2 * Ey[idx4] + b3 * Ey[idx5] + b4 * Ey[idx6] + b5 * Ey[idx7];
            Ez[idx0] = Ez[idx1] + b0 * Ez[idx2] + b1 * Ez[idx3] + b2 * Ez[idx4] + b3 * Ez[idx5] + b4 * Ez[idx6] + b5 * Ez[idx7];
            Ay[idx0] = Ay[idx1] + b0 * Ay[idx2] + b1 * Ay[idx3] + b2 * Ay[idx4] + b3 * Ay[idx5] + b4 * Ay[idx6] + b5 * Ay[idx7];
            Az[idx0] = Az[idx1] + b0 * Az[idx2] + b1 * Az[idx3] + b2 * Az[idx4] + b3 * Az[idx5] + b4 * Az[idx6] + b5 * Az[idx7];

            By[idx1] = By[idx0];
            Bz[idx1] = Bz[idx0];
            Ey[idx1] = Ey[idx0];
            Ez[idx1] = Ez[idx0];
            Ay[idx1] = Ay[idx0];
            Az[idx1] = Az[idx0];
        }
    }
}

void EMFieldSolver::DumpCharge() {
    if (chargeStream == NULL) {
        std::string fileName = "output/charge.txt";
        chargeStream = new std::ofstream(fileName);
        (*chargeStream) << std::scientific << std::setprecision(settings.output.precision);
        DumpCharge();
    } else {

        for (unsigned int k=0; k<settings.q.size(); k++) {

            for (unsigned int i = 0; i < charge.size(); i++) {
                (*chargeStream) << charges[k][i] << " ";
            }

            (*chargeStream) << std::endl;

        }

    }

}

void EMFieldSolver::DumpEnergy() {
    if (energyStreams.size() == 0) {
        for (unsigned int i=0; i<settings.q.size(); i++) {
            std::stringstream fileName;
            fileName << "output/dNdP_" << i << ".txt";
            energyStreams.push_back(std::make_unique<std::ofstream>(fileName.str()));
            (*energyStreams[i]) << std::scientific << std::setprecision(settings.output.precision);
        }
        DumpEnergy();
    } else {
        for (unsigned int k=0; k<settings.q.size(); k++) {
            for (unsigned int i = 0; i < energies[k].size(); i++) {
                (*energyStreams[k]) << energies[k][i] << " ";
            }
            (*energyStreams[k]) << std::endl;
        }
    }
}

void EMFieldSolver::DumpEFieldLongitudinal() {
    if (ELongStream == NULL) {
        std::string fileName = "output/EFieldLong.txt";
        ELongStream = new std::ofstream(fileName);
        (*ELongStream) << std::scientific << std::setprecision(settings.output.precision);
        DumpEFieldLongitudinal();
    } else {
        for (unsigned int i = 0; i < charge.size(); i++) {
            (*ELongStream) << GetEfield(i) << " ";
        }

        (*ELongStream) << std::endl;
    }
}

void EMFieldSolver::DumpEFieldTransverse() {
    if (ETransStream == NULL) {
        std::string fileName = "output/EFieldTrans.txt";
        ETransStream = new std::ofstream(fileName);
        (*ETransStream) << std::scientific << std::setprecision(settings.output.precision);
        DumpEFieldTransverse();
    } else {
        for (unsigned int i = 0; i < n_prepad + x_size + n_postpad; i++) {
            (*ETransStream) << Ey[i] << " ";
        }

        (*ETransStream) << std::endl;

        for (unsigned int i = 0; i < n_prepad + x_size + n_postpad; i++) {
            (*ETransStream) << Ez[i] << " ";
        }

        (*ETransStream) << std::endl;
    }
}

void EMFieldSolver::DumpPotential() {
    if (potentialStream == NULL) {
        std::string fileName = "output/potential.txt";
        potentialStream = new std::ofstream(fileName);
        (*potentialStream) << std::scientific << std::setprecision(settings.output.precision);
        DumpPotential();
    } else {
        for (unsigned int i = 0; i < x_size; i++) {
            (*potentialStream) << PHI[i] << " ";
        }

        (*potentialStream) << std::endl;
    }
}

void EMFieldSolver::DumpBFieldTransverse() {
    if (BStream == NULL) {
        std::string fileName = "output/BFieldTrans.txt";
        BStream = new std::ofstream(fileName);
        (*BStream) << std::scientific << std::setprecision(settings.output.precision);
        DumpBFieldTransverse();
    } else {
        for (unsigned int i = 0; i < n_prepad + x_size + n_postpad; i++) {
            (*BStream) << By[i] << " ";
        }

        (*BStream) << std::endl;

        for (unsigned int i = 0; i < n_prepad + x_size + n_postpad; i++) {
            (*BStream) << Bz[i] << " ";
        }

        (*BStream) << std::endl;
    }
}

void EMFieldSolver::DumpAsqField() {
    if (AsqStream == NULL) {
        std::string fileName = "output/ASquared.txt";
        AsqStream = new std::ofstream(fileName);
        (*AsqStream) << std::scientific << std::setprecision(settings.output.precision);
        DumpAsqField();
    } else {
        for (unsigned int i = 0; i < x_size; i++) {
            (*AsqStream) << GetASquared(i) << " ";
        }

        (*AsqStream) << std::endl;
    }
}

void EMFieldSolver::DumpTime(double time) {
    if (timeStream == NULL) {
        std::string fileName = "output/time.txt";
        timeStream = new std::ofstream(fileName);
        (*timeStream) << std::scientific << std::setprecision(settings.output.precision);
        DumpTime(time);
    } else {
        (*timeStream) << time << std::endl; //Have to flush this one..
    }
}

void EMFieldSolver::RGKCalculateRHS(int step) {
    unsigned int iloop = x_size + n_prepad + n_postpad;
    #pragma ivdep
    for (unsigned int i = 0; i < iloop; i++) {
        int idxs2 = Index(i, step + 2);
        By[idxs2] = 0;
        Bz[idxs2] = 0;
        Ey[idxs2] = 0;
        Ez[idxs2] = 0;
        Ay[idxs2] = 0;
        Az[idxs2] = 0;
    }

    double dx_inv = 1 / settings.GetDx(0);

    {
        int i = 0;
        int idxs2 = Index(i, step + 2);
        int idx1 = Index(i, 1);
        int idxp1 = Index(i + 1, 1);
        By[idxs2] = dx_inv * (Ez[idxp1] - Ez[idx1]);
        Bz[idxs2] = -dx_inv * (Ey[idxp1] - Ey[idx1]);
        Ey[idxs2] = -dx_inv * (Bz[idx1] - settings.GetBZ(0, settings.time)) * eps0_inv * mu_inv;
        Ez[idxs2] = dx_inv * (By[idx1] - settings.GetBY(0, settings.time)) * eps0_inv * mu_inv;
        Ay[idxs2] = -Ey[idx1];
        Az[idxs2] = -Ez[idx1];
    }
    int imax = n_prepad;
    #pragma ivdep
    for (int i = 1; i < imax; i++) {
        int idxs2 = Index(i, step + 2);
        int idx1 = Index(i, 1);
        int idxp1 = Index(i + 1, 1);
        int idxn1 = Index(i - 1, 1);
        By[idxs2] = dx_inv * (Ez[idxp1] - Ez[idx1]);
        Bz[idxs2] = -dx_inv * (Ey[idxp1] - Ey[idx1]);
        Ey[idxs2] = -dx_inv * (Bz[idx1] - Bz[idxn1]) * eps0_inv * mu_inv;
        Ez[idxs2] = dx_inv * (By[idx1] - By[idxn1]) * eps0_inv * mu_inv;
        Ay[idxs2] = -Ey[idx1];
        Az[idxs2] = -Ez[idx1];
    }

    double c1 = -1.0 / 24;
    double c2 = 9.0 / 8.0;
    imax = n_prepad + (int)x_size;
    #pragma ivdep
    for (int i = n_prepad; i < imax; i++) {
        int idxs2 = Index(i, step + 2);
        int idx1 = Index(i, 1);
        int idxp1 = Index(i + 1, 1);
        int idxn1 = Index(i - 1, 1);
        int idxp2 = Index(i + 2, 1);
        int idxn2 = Index(i - 2, 1);
        By[idxs2] = dx_inv * (c1 * (Ez[idxp2] - Ez[idxn1]) + c2 * (Ez[idxp1] - Ez[idx1]));
        Bz[idxs2] = -dx_inv * (c1 * (Ey[idxp2] - Ey[idxn1]) + c2 * (Ey[idxp1] - Ey[idx1]));
        Ey[idxs2] = -dx_inv * (c1 * (Bz[idxp1] - Bz[idxn2]) + c2 * (Bz[idx1] - Bz[idxn1])) * eps0_inv * mu_inv - eps0_inv * J[i - n_prepad] * Ay[idx1];
        Ez[idxs2] = dx_inv * (c1 * (By[idxp1] - By[idxn2]) + c2 * (By[idx1] - By[idxn1])) * eps0_inv * mu_inv - eps0_inv * J[i - n_prepad] * Az[idx1];
        Ay[idxs2] = -Ey[idx1];
        Az[idxs2] = -Ez[idx1];
    }
    imax = n_prepad + (int)x_size + n_postpad - 2;
    #pragma ivdep
    for (int i = n_prepad + (int)x_size; i < imax; i++) {
        int idxs2 = Index(i, step + 2);
        int idx1 = Index(i, 1);
        int idxp1 = Index(i + 1, 1);
        int idxn1 = Index(i - 1, 1);
        By[idxs2] = dx_inv * (Ez[idxp1] - Ez[idx1]);
        Bz[idxs2] = -dx_inv * (Ey[idxp1] - Ey[idx1]);
        Ey[idxs2] = -dx_inv * (Bz[idx1] - Bz[idxn1]) * eps0_inv * mu_inv;
        Ez[idxs2] = dx_inv * (By[idx1] - By[idxn1]) * eps0_inv * mu_inv;
        Ay[idxs2] = -Ey[idx1];
        Az[idxs2] = -Ez[idx1];
    }
}

void EMFieldSolver::InterpolateToFaces() {
    #pragma simd
    for (unsigned int i = 0; i < x_size; i++) {
        double a_y = 0.0;
        {
            double f1 = Ay[Index(n_prepad + i - 2, 1)];
            double f2 = Ay[Index(n_prepad + i - 1, 1)];
            double f3 = Ay[Index(n_prepad + i, 1)];
            double f4 = Ay[Index(n_prepad + i + 1, 1)];

            double fL = (1.0 / 6) * (-f1 + 5 * f2 + 2 * f3);
            double fR = (1.0 / 6) * (2 * f2 + 5 * f3 - f4);

            double AL = f1 - 2 * f2 + f3;
            double BL = f3 - f1;

            double AR = f2 - 2 * f3 + f4;
            double BR = f4 - f2;

            double bL = 4.0 / 3 * (AL * AL) + 0.5 * AL * BL + 0.25 * (BL * BL);
            double bR = 4.0 / 3 * (AR * AR) - 0.5 * AR * BR + 0.25 * (BR * BR);

            double mm(1.0e-10);

            double oL = 0.5 / ((mm + bL) * (mm + bL));
            double oR = 0.5 / ((mm + bR) * (mm + bR));

            double wL = oL / (oL + oR);
            double wR = oR / (oL + oR);

            a_y = wL * fL + wR * fR;
        }
        double a_z = 0.0;
        {
            double f1 = Az[Index(n_prepad + i - 2, 1)];
            double f2 = Az[Index(n_prepad + i - 1, 1)];
            double f3 = Az[Index(n_prepad + i, 1)];
            double f4 = Az[Index(n_prepad + i + 1, 1)];

            double fL = (1.0 / 6) * (-f1 + 5 * f2 + 2 * f3);
            double fR = (1.0 / 6) * (2 * f2 + 5 * f3 - f4);

            double AL = f1 - 2 * f2 + f3;
            double BL = f3 - f1;

            double AR = f2 - 2 * f3 + f4;
            double BR = f4 - f2;

            double bL = 4.0 / 3 * (AL * AL) + 0.5 * AL * BL + 0.25 * (BL * BL);
            double bR = 4.0 / 3 * (AR * AR) - 0.5 * AR * BR + 0.25 * (BR * BR);

            double mm(1.0e-10);

            double oL = 0.5 / ((mm + bL) * (mm + bL));
            double oR = 0.5 / ((mm + bR) * (mm + bR));

            double wL = oL / (oL + oR);
            double wR = oR / (oL + oR);

            a_z = wL * fL + wR * fR;
        }
        a_squared[i] = (a_y * a_y) + (a_z * a_z);
    }

}

void EMFieldSolver::EnforceChargeNeutralization() {
    AssembleRhoAndJ();

    std::fill(J.begin(), J.end(), 0.0);

    for (unsigned int i = 0; i < charge.size(); i++) {
        neutralizationCharge[i] = -charge[i];
    }
}

double EMFieldSolver::EstimateCFLBound() {
    int n = meshes.size();

    std::vector<double> &m = settings.m;
    std::vector<double> &q = settings.q;
    std::vector<double> dps;

    dps.assign(n, 0.0);

    double dpsMax(0.0);

    for (int i = 0; i < n; i++) {
        dps[i] = 1 / settings.GetDp(0, i);
        dpsMax = std::max(dpsMax, std::fabs(q[i]) * dps[i]);
    }

    double pc(0.0);

    for (unsigned int i = 0; i < x_size; i++) {
        double Azv = Az[i];
        double Ayv = Ay[i];
        double As = Ayv * Ayv + Azv * Azv;

        double temp(0.0);
        #pragma novector
        for (int j = 0; j < n; j++) {
            temp = std::max(temp, std::fabs(q[j]) / m[j] / std::sqrt(1 + As / ((m[j] * cs) * (m[j] * cs))) * dps[j]);
        }

        pc = std::max(pc, temp * std::fabs(Ayv * Bz[i] - Azv * By[i]) + dpsMax * std::fabs(GetEfield(i)));
    }

    return 1.0 / std::max((cs / settings.GetDx(0) + pc), 1e-40);
}

double EMFieldSolver::GetMagneticForce(int i) {

    i+=n_prepad;
    i=i>-1?i:0;
    i=i<(x_size+n_prepad+n_postpad)?i:(x_size+n_prepad+n_postpad-1);

    return Az[Index(i,1)]*By[Index(i,1)]-Ay[Index(i,1)]*Bz[Index(i,1)];
}
