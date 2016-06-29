#include "veritas.hpp"
#include "Settings.hpp"
#include "Rectangle.hpp"
#include "BoundaryCondition.hpp"
#include "EMSolver.hpp"

Rectangle::Rectangle(int n_x, int n_p, int x_pos, int p_pos, int depth, Settings &settings, const std::shared_ptr<Rectangle> &bc, bool up, bool down, bool left, bool right, int particleType) : refinementRatio(settings.refinementRatio), depth(depth), particleType(particleType), dp(settings.GetDp(depth, particleType)), dx(settings.GetDx(depth)), m_inv(1 / settings.GetMass(particleType)), settings(settings), EMSolver(NULL), up(up), down(down), left(left), right(right), n_x(n_x), n_p(n_p), x_pos(x_pos), p_pos(p_pos), relativeToBottom(std::pow(refinementRatio, depth)) {

    // Solution
    f.assign(3 * (n_x + 4) * (n_p + 4), 0.0);

    // Edge values
    fp.assign((n_x + 4) * (n_p + 4), 0.0);
    fx.assign((n_x + 4) * (n_p + 4), 0.0);

    // Edge fields.
    ep.assign((n_x + 4) * (n_p + 4), 0.0);
    ex.assign((n_x + 4) * (n_p + 4), 0.0);

    // Fluxes at different times
    FxL.assign(6 * (n_x + 4) * (n_p + 4), 0.0);
    FxH.assign(6 * (n_x + 4) * (n_p + 4), 0.0);
    FpL.assign(6 * (n_x + 4) * (n_p + 4), 0.0);
    FpH.assign(6 * (n_x + 4) * (n_p + 4), 0.0);

    // Fluxes at a given step.
    FxLS.assign((n_x + 4) * (n_p + 4), 0.0);
    FpLS.assign((n_x + 4) * (n_p + 4), 0.0);

    // Fluxes at a given step.
    FxDS.assign((n_x + 4) * (n_p + 4), 0.0);
    FpDS.assign((n_x + 4) * (n_p + 4), 0.0);

    // Limiter
    Rp.assign((n_x + 4) * (n_p + 4), 0.0);
    Rm.assign((n_x + 4) * (n_p + 4), 0.0);
    Cp.assign((n_x + 4) * (n_p + 4), 0.0);
    Cx.assign((n_x + 4) * (n_p + 4), 0.0);

    exterior_boundary_xm.assign(n_p / refinementRatio + 2, bc);
    exterior_boundary_xp.assign(n_p / refinementRatio + 2, bc);
    exterior_boundary_pm.assign(n_x / refinementRatio, bc);
    exterior_boundary_pp.assign(n_x / refinementRatio, bc);

    is_exterrior_boundary_same_level_xm.assign(n_p / refinementRatio + 2, true);
    is_exterrior_boundary_same_level_xp.assign(n_p / refinementRatio + 2, true);
    is_exterrior_boundary_same_level_pm.assign(n_x / refinementRatio, true);
    is_exterrior_boundary_same_level_pp.assign(n_x / refinementRatio, true);

    finer_level.assign((n_x + 4) * (n_p + 4), nullptr);
    finer_level_x.assign((n_x + 4) * (n_p + 4), nullptr);
    finer_level_p.assign((n_x + 4) * (n_p + 4), nullptr);

    is_nested.assign((n_x + 4) * (n_p + 4), false);
    is_interrior_level_boundary_x.assign((n_x + 4) * (n_p + 4), false);
    is_interrior_level_boundary_p.assign((n_x + 4) * (n_p + 4), false);
    is_interpolated.assign((n_x + 4) * (n_p + 4), false);

    chargeRNested.assign(n_x + 4, 0.0);
    currentRNested.assign(n_x + 4, 0.0);
    chargeRNoNested.assign(n_x + 4, 0.0);
    currentRNoNested.assign(n_x + 4, 0.0);

    chargeR.assign(n_x * relativeToBottom, 0.0);
    energyR.assign(n_p * relativeToBottom, 0.0);
    currentR.assign(n_x * relativeToBottom, 0.0);

    errorWeights.assign(5, 0.0);

    errorWeights[0] = 0.5 / settings.GetfMax(particleType) * std::pow(refinementRatio, -0.5 * depth);
    errorWeights[1] = 0.5 / settings.GetfMax(particleType) * std::pow(refinementRatio, -0.5 * depth);
    errorWeights[2] = 0.5 / settings.GetfMax(particleType) * std::pow(refinementRatio, -2.0 * depth);
    errorWeights[3] = 0.5 / settings.GetfMax(particleType) * std::pow(refinementRatio, -2.0 * depth);
    errorWeights[4] = settings.sizeWeight / settings.GetfMax(particleType) * std::pow(refinementRatio, -settings.maxDepth + depth);

    interpolationMatrix.resize(12, 0.0);

    interpolationMatrix[0] = 0.104166666666667;
    interpolationMatrix[1] = -0.708333333333334;
    interpolationMatrix[2] = 0.708333333333334;
    interpolationMatrix[3] = -0.104166666666667;
    interpolationMatrix[4] = 0.117647058823529;
    interpolationMatrix[5] = 0.029411764705882;
    interpolationMatrix[6] = 0.029411764705882;
    interpolationMatrix[7] = 0.117647058823529;
    interpolationMatrix[8] = -0.083333333333333;
    interpolationMatrix[9] = 0.166666666666667;
    interpolationMatrix[10] = -0.166666666666667;
    interpolationMatrix[11] = 0.083333333333334;

    interpolatCoefsREL.resize(3 * relativeToBottom, 0.0);
    interpolatCoefsREF.resize(3 * refinementRatio, 0.0);

    for (int i = 0; i < refinementRatio; i++) {
        double tl = -0.5 + i / (double)refinementRatio;
        double tr = -0.5 + (i + 1.0) / (double)refinementRatio;

        interpolatCoefsREF[3 * i] = (tl + tr) * 0.5;
        interpolatCoefsREF[3 * i + 1] = (tl * tl + tl * tr + tr * tr) / 3.0 - 1.0 / 12;
        interpolatCoefsREF[3 * i + 2] = (tl * tl * tl + tl * tl * tr + tl * tr * tr + tr * tr * tr) * 0.25;
    }

    for (int i = 0; i < relativeToBottom; i++) {
        double tl = -0.5 + i / (double)relativeToBottom;
        double tr = -0.5 + (i + 1.0) / (double)relativeToBottom;

        interpolatCoefsREL[3 * i] = (tl + tr) * 0.5;
        interpolatCoefsREL[3 * i + 1] = (tl * tl + tl * tr + tr * tr) / 3.0 - 1.0 / 12;
        interpolatCoefsREL[3 * i + 2] = (tl * tl * tl + tl * tl * tr + tl * tr * tr + tr * tr * tr) * 0.25;
    }

    if (LOUD) std::cout << n_x << ", " << n_p << ", " << x_pos << ", " << p_pos << ", " << depth << ", " << particleType << ", " << dp << ", " << dx << ", " << relativeToBottom << ", " << up << ", " << down << ", " << left << ", " << right << std::endl;
}

Rectangle::Rectangle() : settings(settings) {
}

Rectangle::~Rectangle() {
}

std::vector<double> Rectangle::GetInterpolantsREF(double f1, double f2, double f3, double f4, double f5) {
    std::vector<double> temp;

    temp.resize(refinementRatio, 0.0);

    f5 -= f3; f4 -= f3; f2 -= f3; f1 -= f3;

    double a1 = interpolationMatrix[0] * f1 + interpolationMatrix[1] * f2 + interpolationMatrix[2] * f4 + interpolationMatrix[3] * f5;
    double a2 = interpolationMatrix[4] * f1 + interpolationMatrix[5] * f2 + interpolationMatrix[6] * f4 + interpolationMatrix[7] * f5;
    double a3 = interpolationMatrix[8] * f1 + interpolationMatrix[9] * f2 + interpolationMatrix[10] * f4 + interpolationMatrix[11] * f5;
    #pragma ivdep
    for (int i = 0; i < refinementRatio; i++) {
        temp[i] = interpolatCoefsREF[3 * i] * a1 + interpolatCoefsREF[3 * i + 1] * a2 + interpolatCoefsREF[3 * i + 2] * a3 + f3;
    }

    return temp;
}

std::vector<double> Rectangle::GetInterpolantsREL(double f1, double f2, double f3, double f4, double f5) {
    std::vector<double> temp;
    int rtb = relativeToBottom;
    temp.resize(rtb, 0.0);

    f5 -= f3; f4 -= f3; f2 -= f3; f1 -= f3;

    double a1 = interpolationMatrix[0] * f1 + interpolationMatrix[1] * f2 + interpolationMatrix[2] * f4 + interpolationMatrix[3] * f5;
    double a2 = interpolationMatrix[4] * f1 + interpolationMatrix[5] * f2 + interpolationMatrix[6] * f4 + interpolationMatrix[7] * f5;
    double a3 = interpolationMatrix[8] * f1 + interpolationMatrix[9] * f2 + interpolationMatrix[10] * f4 + interpolationMatrix[11] * f5;
    #pragma ivdep
    for (int i = 0; i < rtb; i++) {
        temp[i] = interpolatCoefsREL[(int)3 * i] * a1 + interpolatCoefsREL[(int)3 * i + (int)1] * a2 + interpolatCoefsREL[(int)3 * i + (int)2] * a3 + f3;
    }

    return temp;
}

void Rectangle::CalculateRhoAndJ() {
    std::fill(chargeR.begin(), chargeR.end(), 0.0);
    std::fill(currentR.begin(), currentR.end(), 0.0);

    double q = settings.GetCharge(particleType);
    double c1 = m_inv * c_inv;
    double c2 = 1 / c1;
    double c3 = 1 / 48.0;
    int rtb = relativeToBottom;

    #if defined(USINGMKL)
    //unsigned int oldmode = vmlSetMode( VML_EP );
    std::vector<double> loginput;
    std::vector<double> logoutput;
    std::vector<double> a_squared;
    loginput.assign(3*rtb,0.0);
    logoutput.assign(3*rtb,0.0);
    a_squared.assign(rtb,0.0);
    #endif

    for (int i = 0; i < n_x; i++) {
        int j = 0;
        std::vector<double> temp = (GetInterpolantsREL(f[Index3(i - 2, j, 1)], f[Index3(i - 1, j, 1)], f[Index3(i, j, 1)], f[Index3(i + 1, j, 1)], f[Index3(i + 2, j, 1)]));
        std::vector<double> tempm1 = (GetInterpolantsREL(f[Index3(i - 2, j - 1, 1)], f[Index3(i - 1, j - 1, 1)], f[Index3(i, j - 1, 1)], f[Index3(i + 1, j - 1, 1)], f[Index3(i + 2, j - 1, 1)]));
        std::vector<double> tempp1 = (GetInterpolantsREL(f[Index3(i - 2, j + 1, 1)], f[Index3(i - 1, j + 1, 1)], f[Index3(i, j + 1, 1)], f[Index3(i + 1, j + 1, 1)], f[Index3(i + 2, j + 1, 1)]));
        double te = 0.0, tp1 = 0.0, tm1 = 0.0;

        for (int k = 0; k < rtb; k++) {
            te += temp[k];
            tp1 += tempp1[k];
            tm1 += tempm1[k];
        }

        double cor = f[Index3(i, j, 1)] - (1.0 / relativeToBottom) * te;
        double corm1 = f[Index3(i, j - 1, 1)] - (1.0 / relativeToBottom) * tm1;
        double corp1 = f[Index3(i, j + 1, 1)] - (1.0 / relativeToBottom) * tp1;

        for (int k = 0; k < rtb; k++) {
            temp[k] += cor;
            tempm1[k] += corm1;
            tempp1[k] += corp1;
        }

        double gm = 0.0, gm1 = 0.0, gp1 = 0.0;
        for (int j = 0; j < n_p; j++) {
            if (!is_nested[IndexNS(i, j)]) {
                #if defined(USINGMKL)
                #pragma ivdep
                for (int k = 0; k < rtb; k++) {
                    a_squared[k] = ((EMFieldSolver*)EMSolver)->EMFieldSolver::GetCellAverageASquared((i + x_pos) * rtb + k);

                    chargeR[i * rtb + k] += temp[k];

                    double gammam1 = Gamma(Momentum(j - 1), q * q * a_squared[k]);
                    double gamma0 = Gamma(Momentum(j), q * q * a_squared[k]);
                    double gamma1 = Gamma(Momentum(j + 1), q * q * a_squared[k]);
                    double gamma2 = Gamma(Momentum(j + 2), q * q * a_squared[k]);

                    //NOTE: While it vectorises much more efficiently here if we dont have these memory strides (ie three separate variables),
                    //It takes much longer to set up the process. So we save time this way.
                    loginput[k] = (gamma1 + c1 * Momentum(j + 1)) / (gamma0 + c1 * Momentum(j));
                    loginput[k+rtb] = (gamma0 + c1 * Momentum(j)) / (gammam1 + c1 * Momentum(j - 1));
                    loginput[k+2*rtb] = (gamma2 + c1 * Momentum(j + 2)) / (gamma1 + c1 * Momentum(j + 1));
                }
                vdLn(rtb,&loginput[0],&logoutput[0]);
                vdLn(rtb,&loginput[rtb],&logoutput[rtb]);
                vdLn(rtb,&loginput[2*rtb],&logoutput[2*rtb]);
                for (int k = 0; k < rtb; k++) {
                    //No vector performance gains here due to the memory stride.
                    double gm = c2 * logoutput[k];
                    double gm1 = c2 * logoutput[k+rtb];
                    double gp1 = c2 * logoutput[k+2*rtb];

                    currentR[i * rtb + k] += temp[k] * gm + c3 * (gp1 - gm1) * (tempp1[k] - tempm1[k]);
                }
                #else
                for (int k = 0; k < rtb; k++) {
                    double a_squared = ((EMFieldSolver*)EMSolver)->EMFieldSolver::GetCellAverageASquared((i + x_pos) * rtb + k);

                    double gammam1 = Gamma(Momentum(j - 1), q * q * a_squared);
                    double gamma0 = Gamma(Momentum(j), q * q * a_squared);
                    double gamma1 = Gamma(Momentum(j + 1), q * q * a_squared);
                    double gamma2 = Gamma(Momentum(j + 2), q * q * a_squared);

                    gm = c2 * std::log((gamma1 + c1 * Momentum(j + 1)) / (gamma0 + c1 * Momentum(j)));
                    gm1 = c2 * std::log((gamma0 + c1 * Momentum(j)) / (gammam1 + c1 * Momentum(j - 1)));
                    gp1 = c2 * std::log((gamma2 + c1 * Momentum(j + 2)) / (gamma1 + c1 * Momentum(j + 1)));
                }
                for (int k = 0; k < rtb; k++) {
                    chargeR[i * rtb + k] += temp[k];
                    currentR[i * rtb + k] += temp[k] * gm + c3 * (gp1 - gm1) * (tempp1[k] - tempm1[k]);
                }
                #endif
            }

            tempm1.swap(temp);
            temp.swap(tempp1);
            tempp1 = (GetInterpolantsREL(f[Index3(i - 2, j + 2, 1)], f[Index3(i - 1, j + 2, 1)], f[Index3(i, j + 2, 1)], f[Index3(i + 1, j + 2, 1)], f[Index3(i + 2, j + 2, 1)]));

            tp1 = 0.0;

            for (int k = 0; k < rtb; k++) {
                tp1 += tempp1[k];
            }

            corp1 = f[Index3(i, j + 2, 1)] - (1.0 / relativeToBottom) * tp1;

            #pragma ivdep
            for (int k = 0; k < rtb; k++) {
                tempp1[k] += corp1;
            }
        }

        for (int k = 0; k < rtb; k++) {
            chargeR[i * rtb + k] *= dp * q;
            currentR[i * rtb + k] *= -q * q / settings.GetMass(particleType);
        }
    }
}

void Rectangle::CalculateEnergy() {
    std::fill(energyR.begin(), energyR.end(), 0.0);
    int rtb = relativeToBottom;
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_p; j++) {
            if (!is_nested[IndexNS(i, j)]) {
                std::vector<double> tempp = (GetInterpolantsREL(f[Index3(i, j - 2, 1)], f[Index3(i, j - 1, 1)], f[Index3(i, j, 1)], f[Index3(i, j + 1, 1)], f[Index3(i, j + 2, 1)]));
                for (int k = 0; k < rtb; k++) {
                    energyR.at(j * rtb + k) += tempp[k];
                }
            }
        }
    }
    for (int j = 0; j < n_p; j++) {
        for (int k = 0; k < rtb; k++) {
            energyR.at(j * rtb + k) *= dx;
        }
    }
}

double Rectangle::GetValueFromSameLevel(int i, int j,int val) {
    int i_neighbour = i - x_pos;
    int j_neighbour = j - p_pos;

    return f[Index3(i_neighbour, j_neighbour, val)];
}

double Rectangle::GetValueFromFinerLevel(int i, int j,int val) {
    int i_finer = i * refinementRatio - x_pos;
    int j_finer = j * refinementRatio - p_pos;
    double temp(0.0);

    for (int k = 0; k < refinementRatio; k++) {
        for (int l = 0; l < refinementRatio; l++) {
            temp += f[Index3(i_finer + k, j_finer + l, val)];
        }
    }

    temp /= std::pow((double)refinementRatio, 2.0);
    return temp;
}

void Rectangle::UpdateInterriorPoints(int val) {
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_p; j++) {
            if (is_nested[IndexNS(i, j)]) {
                f[Index3(i, j, val)] = finer_level[IndexNS(i, j)]->GetValueFromFinerLevel(x_pos + i, p_pos + j,val);
            }
        }
    }
}

void Rectangle::SetFieldSolver(EMFieldSolver *solver) {
    EMSolver = solver;
}

std::vector<double> Rectangle::GetWenoValueFromCoarseLevel(int i, int j, int d,int val) {
    int i_coarse = i / refinementRatio - x_pos;
    int j_coarse = j / refinementRatio - p_pos;

    std::vector< std::vector<double> > temps;

    for (int k = -2; k < 3; k++) {
        temps.push_back(GetInterpolantsREF(f[Index3(i_coarse - 2, j_coarse + k, val)], f[Index3(i_coarse - 1, j_coarse + k, val)], f[Index3(i_coarse, j_coarse + k, val)], f[Index3(i_coarse + 1, j_coarse + k, val)], f[Index3(i_coarse + 2, j_coarse + k, val)]));
    }

    std::vector<double> interpolants;
    interpolants.resize(refinementRatio * refinementRatio);
    double sum(0.0);

    for (int k = 0; k < refinementRatio; k++) {
        std::vector<double> interpolants_part = GetInterpolantsREF(temps[0][k], temps[1][k], temps[2][k], temps[3][k], temps[4][k]);

        for (int l = 0; l < refinementRatio; l++) {
            interpolants[k * refinementRatio + l] = interpolants_part[l];
            sum += interpolants_part[l];
        }
    }

    double correction = f[Index3(i_coarse, j_coarse, val)] - 1.0 / std::pow(refinementRatio, 2) * sum;

    for (int k = 0; k < refinementRatio; k++) {
        for (int l = 0; l < refinementRatio; l++) {
            interpolants[k * refinementRatio + l] += correction;
        }
    }

    if (d == 0) {
        std::vector<double> interpolants_reduced;
        interpolants_reduced.resize(refinementRatio * 2);

        for (int k = 0; k < refinementRatio; k++) {
            interpolants_reduced[2 * k] = interpolants[refinementRatio * k + refinementRatio - 1];
            interpolants_reduced[2 * k + 1] = interpolants[refinementRatio * k + refinementRatio - 2];
        }

        return interpolants_reduced;
    } else if (d == 1) {
        std::vector<double> interpolants_reduced;
        interpolants_reduced.resize(refinementRatio * 2);

        for (int k = 0; k < refinementRatio; k++) {
            interpolants_reduced[2 * k] = interpolants[k];
            interpolants_reduced[2 * k + 1] = interpolants[refinementRatio + k];
        }

        return interpolants_reduced;
    } else if (d == 2) {
        std::vector<double> interpolants_reduced;
        interpolants_reduced.resize(refinementRatio * 2);

        for (int k = 0; k < refinementRatio; k++) {
            interpolants_reduced[2 * k] = interpolants[refinementRatio * k];
            interpolants_reduced[2 * k + 1] = interpolants[refinementRatio * k + 1];
        }

        return interpolants_reduced;
    } else if (d == 3) {
        std::vector<double> interpolants_reduced;
        interpolants_reduced.resize(refinementRatio * 2);

        for (int k = 0; k < refinementRatio; k++) {
            interpolants_reduced[2 * k] = interpolants[refinementRatio * (refinementRatio - 1) + k];
            interpolants_reduced[2 * k + 1] = interpolants[refinementRatio * (refinementRatio - 2) + k];
        }

        return interpolants_reduced;
    } else {
        return interpolants;
    }
}

void Rectangle::UpdateDifferentLevelBoundaries(int val) {
    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xm.size() - 2; i++) {
        if (!is_exterrior_boundary_same_level_xm[i + 1]) {
            if (auto spt = exterior_boundary_xm[i + 1].lock()) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos - 1, p_pos + i * refinementRatio, 3,val);

                for (int j = 0; j < refinementRatio; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(-2, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xp.size() - 2; i++) {
        if (!is_exterrior_boundary_same_level_xp[i + 1]) {
            if (auto spt = exterior_boundary_xp[i + 1].lock()) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos + n_x, p_pos + i * refinementRatio, 1,val);

                for (int j = 0; j < refinementRatio; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pm.size(); i++) {
        if (!is_exterrior_boundary_same_level_pm[i]) {
            if (auto spt = exterior_boundary_pm[i].lock()) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos + i * refinementRatio, p_pos - 1, 0,val);

                for (int j = 0; j < refinementRatio; j++) {
                    f[Index3(i * refinementRatio + j, -1, val)] = temp[2 * j];
                    f[Index3(i * refinementRatio + j, -2, val)] = temp[2 * j + 1];
                }
            } else {
                std::cerr << "WARNING: exterior_boundary_pm[" << i << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pp.size(); i++) {
        if (!is_exterrior_boundary_same_level_pp[i]) {
            if (auto spt = exterior_boundary_pp[i].lock()) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos + i * refinementRatio, p_pos + n_p, 2,val);
                for (int j = 0; j < refinementRatio; j++) {
                    f[Index3(i * refinementRatio + j, n_p, val)] = temp[2 * j];
                    f[Index3(i * refinementRatio + j, n_p + 1, val)] = temp[2 * j + 1];
                }
            } else {
                std::cerr << "WARNING: exterior_boundary_pp[" << i << "] pointer is expired" << std::endl;
            }
        }
    }

    {
        int i = -1;

        if (auto spt = exterior_boundary_xm[i + 1].lock()) {
            if (!is_exterrior_boundary_same_level_xm[i + 1]) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos - 1, p_pos + i * refinementRatio, 3,val);
                for (int j = refinementRatio - 2; j < refinementRatio; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(-2, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                for (int j = refinementRatio - 2; j < refinementRatio; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 1, p_pos + i * refinementRatio + j,val);
                    f[Index3(-2, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 2, p_pos + i * refinementRatio + j,val);
                }
            }
        } else {
            std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
        }
    }

    {
        int i = is_exterrior_boundary_same_level_xm.size() - 2;

        if (auto spt = exterior_boundary_xm[i + 1].lock()) {
            if (!is_exterrior_boundary_same_level_xm[i + 1]) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos - 1, p_pos + i * refinementRatio, 3,val);
                for (int j = 0; j < 2; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(-2, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                for (int j = 0; j < 2; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 1, p_pos + i * refinementRatio + j,val);
                    f[Index3(-2, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 2, p_pos + i * refinementRatio + j,val);
                }
            }
        } else {
            std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
        }
    }

    {
        int i = -1;

        if (auto spt = exterior_boundary_xp[i + 1].lock()) {
            if (!is_exterrior_boundary_same_level_xp[i + 1]) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos + n_x, p_pos + i * refinementRatio, 1,val);
                for (int j = refinementRatio - 2; j < refinementRatio; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                for (int j = refinementRatio - 2; j < refinementRatio; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x, p_pos + i * refinementRatio + j,val);
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x + 1, p_pos + i * refinementRatio + j,val);
                }
            }
        } else {
            std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
        }
    }

    {
        int i = is_exterrior_boundary_same_level_xp.size() - 2;

        if (auto spt = exterior_boundary_xp[i + 1].lock()) {
            if (!is_exterrior_boundary_same_level_xp[i + 1]) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos + n_x, p_pos + i * refinementRatio, 1,val);
                for (int j = 0; j < 2; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                for (int j = 0; j < 2; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x, p_pos + i * refinementRatio + j,val);
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x + 1, p_pos + i * refinementRatio + j,val);
                }
            }
        } else {
            std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
        }
    }
}

void Rectangle::UpdateSameLevelBoundaries(int val) {
    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xm.size() - 2; i++) {
        if (is_exterrior_boundary_same_level_xm[i + 1]) {
            if (auto spt = exterior_boundary_xm[i + 1].lock()) {
                for (int j = 0; j < refinementRatio; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 1, p_pos + i * refinementRatio + j,val);
                    f[Index3(-2, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 2, p_pos + i * refinementRatio + j,val);
                }
            } else {
                std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xp.size() - 2; i++) {
        if (is_exterrior_boundary_same_level_xp[i + 1]) {
            if (auto spt = exterior_boundary_xp[i + 1].lock()) {
                for (int j = 0; j < refinementRatio; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x, p_pos + i * refinementRatio + j,val);
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x + 1, p_pos + i * refinementRatio + j,val);
                }
            } else {
                std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pm.size(); i++) {
        if (is_exterrior_boundary_same_level_pm[i]) {
            if (auto spt = exterior_boundary_pm[i].lock()) {
                for (int j = 0; j < refinementRatio; j++) {
                    f[Index3(i * refinementRatio + j, -1, val)] = spt->GetValueFromSameLevel(x_pos + i * refinementRatio + j, p_pos - 1,val);
                    f[Index3(i * refinementRatio + j, -2, val)] = spt->GetValueFromSameLevel(x_pos + i * refinementRatio + j, p_pos - 2,val);
                }
            } else {
                std::cerr << "WARNING: exterior_boundary_pm[" << i << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pp.size(); i++) {
        if (is_exterrior_boundary_same_level_pp[i]) {
            if (auto spt = exterior_boundary_pp[i].lock()) {
                for (int j = 0; j < refinementRatio; j++) {
                    f[Index3(i * refinementRatio + j, n_p, val)] = spt->GetValueFromSameLevel(x_pos + i * refinementRatio + j, p_pos + n_p,val);
                    f[Index3(i * refinementRatio + j, n_p + 1, val)] = spt->GetValueFromSameLevel(x_pos + i * refinementRatio + j, p_pos + n_p + 1,val);
                }
            } else {
                std::cerr << "WARNING: exterior_boundary_pp[" << i << "] pointer is expired" << std::endl;
            }
        }
    }
}

void Rectangle::InitializeDistribution() {
    double temp, xp, pp;

    int nvals = std::pow(refinementRatio, depth + settings.quadratureDepth);

    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_p; j++) {
            temp = 0.0;

            for (int k = 0; k < nvals; k++) {
                for (int l = 0; l < nvals; l++) {
                    xp = ((0.5 + k) / nvals + i + x_pos) * dx;
                    pp = Momentum((0.5 + l) / nvals + j);
                    temp += settings.InitialDistribution(xp, pp, particleType);
                }
            }

            f[Index3(i, j, 0)] = temp * (1 / std::pow(nvals, 2.0));
            f[Index3(i, j, 1)] = temp * (1 / std::pow(nvals, 2.0));
        }
    }
}

bool Rectangle::OnLine(int x, int y1, int y2) {
    return ((x - y1) > -1) && ((y2 - x) > -1);
}

void Rectangle::CalculateConnectivitySame(const std::shared_ptr<Rectangle> &rectangle) {
    int n_x_2 = rectangle->n_x;
    int n_p_2 = rectangle->n_p;
    int x_pos_2 = rectangle->x_pos;
    int p_pos_2 = rectangle->p_pos;
    int x_pos_r = x_pos_2 - x_pos;
    int p_pos_r = p_pos_2 - p_pos;

    if (x_pos_r == n_x) {
        int lower = std::max(p_pos_2, p_pos);
        int upper = std::min(p_pos + n_p, p_pos_2 + n_p_2);

        if (OnLine(-1, p_pos_r, p_pos_r + n_p_2 - 1)) {
            exterior_boundary_xp.at(0) = rectangle;
            is_exterrior_boundary_same_level_xp.at(0) = true;
        }

        if (OnLine(n_p, p_pos_r, p_pos_r + n_p_2 - 1)) {
            exterior_boundary_xp[exterior_boundary_xp.size() - 1] = rectangle;
            is_exterrior_boundary_same_level_xp[exterior_boundary_xp.size() - 1] = true;
        }

        while (lower < upper) {
            int i1 = (lower - p_pos) / refinementRatio + 1;
            exterior_boundary_xp.at(i1) = rectangle;
            is_exterrior_boundary_same_level_xp.at(i1) = true;
            lower += refinementRatio;
        }
    }

    if (x_pos_r == -n_x_2) {
        int lower = std::max(p_pos_2, p_pos);
        int upper = std::min(p_pos + n_p, p_pos_2 + n_p_2);

        if (OnLine(-1, p_pos_r, p_pos_r + n_p_2 - 1)) {
            exterior_boundary_xm[0] = rectangle;
            is_exterrior_boundary_same_level_xm[0] = true;
        }

        if (OnLine(n_p, p_pos_r, p_pos_r + n_p_2 - 1)) {
            exterior_boundary_xm[exterior_boundary_xm.size() - 1] = rectangle;
            is_exterrior_boundary_same_level_xm[exterior_boundary_xm.size() - 1] = true;
        }

        while (lower < upper) {
            int i1 = (lower - p_pos) / refinementRatio + 1;
            exterior_boundary_xm.at(i1) = rectangle;
            is_exterrior_boundary_same_level_xm.at(i1) = true;
            lower += refinementRatio;
        }
    }

    if (p_pos_r == n_p) {
        int lower = std::max(x_pos_2, x_pos);
        int upper = std::min(x_pos + n_x, x_pos_2 + n_x_2);

        if (OnLine(-1, x_pos_r, x_pos_r + n_x_2 - 1)) {
            exterior_boundary_xm[exterior_boundary_xm.size() - 1] = rectangle;
            is_exterrior_boundary_same_level_xm[exterior_boundary_xm.size() - 1] = true;
        }

        if (OnLine(n_x, x_pos_r, x_pos_r + n_x_2 - 1)) {
            exterior_boundary_xp[exterior_boundary_xp.size() - 1] = rectangle;
            is_exterrior_boundary_same_level_xp[exterior_boundary_xp.size() - 1] = true;
        }

        while (lower < upper) {
            int i1 = (lower - x_pos) / refinementRatio;
            exterior_boundary_pp.at(i1) = rectangle;
            is_exterrior_boundary_same_level_pp.at(i1) = true;
            lower += refinementRatio;
        }
    }

    if (p_pos_r == -n_p_2) {
        int lower = std::max(x_pos_2, x_pos);
        int upper = std::min(x_pos + n_x, x_pos_2 + n_x_2);

        if (OnLine(-1, x_pos_r, x_pos_r + n_x_2 - 1)) {
            exterior_boundary_xm[0] = rectangle;
            is_exterrior_boundary_same_level_xm[0] = true;
        }

        if (OnLine(n_x, x_pos_r, x_pos_r + n_x_2 - 1)) {
            exterior_boundary_xp[0] = rectangle;
            is_exterrior_boundary_same_level_xp[0] = true;
        }

        while (lower < upper) {
            int i1 = (lower - x_pos) / refinementRatio;
            exterior_boundary_pm.at(i1) = rectangle;
            is_exterrior_boundary_same_level_pm.at(i1) = true;
            lower += refinementRatio;
        }
    }
}

void Rectangle::CalculateConnectivityFromFiner(const std::shared_ptr<Rectangle> &current, const std::shared_ptr<Rectangle> &rectangle) {
    int x_pos_2 = rectangle->x_pos / refinementRatio - x_pos;
    int p_pos_2 = rectangle->p_pos / refinementRatio - p_pos;
    int n_x_2 = rectangle->n_x / refinementRatio;
    int n_p_2 = rectangle->n_p / refinementRatio;
    int lower_x = std::max(0, x_pos_2);
    int upper_x = std::min(n_x, x_pos_2 + n_x_2);
    int lower_p = std::max(0, p_pos_2);
    int upper_p = std::min(n_p, p_pos_2 + n_p_2);

    for (int i = lower_x; i < upper_x; i++) {
        for (int j = lower_p; j < upper_p; j++) {
            is_nested[IndexNS(i, j)] = true;
            finer_level[IndexNS(i, j)] = rectangle;
        }
    }

    if (((x_pos_2 + n_x_2) < (n_x + 1)) && ((x_pos_2 + n_x_2) > -1)) {
        for (int i = lower_p; i < upper_p; i++) {
            finer_level_x[IndexNS(x_pos_2 + n_x_2, i)] = rectangle;
            is_interrior_level_boundary_x[IndexNS(x_pos_2 + n_x_2, i)] = true;
        }
    }

    if (((x_pos_2 + n_x_2) < n_x) && ((x_pos_2 + n_x_2) > -1)) {
        for (int i = lower_p; i < upper_p; i++) {
            int j = i - p_pos_2 + 1;
            rectangle->exterior_boundary_xp[j] = current;
            rectangle->is_exterrior_boundary_same_level_xp[j] = false;
        }

        if (OnLine(p_pos_2 - 1, 0, n_p - 1)) {
            rectangle->exterior_boundary_xp[0] = current;
            rectangle->is_exterrior_boundary_same_level_xp[0] = false;
        }

        if (OnLine(p_pos_2 + n_p_2, 0, n_p - 1)) {
            rectangle->exterior_boundary_xp[rectangle->exterior_boundary_xp.size() - 1] = current;
            rectangle->is_exterrior_boundary_same_level_xp[rectangle->exterior_boundary_xp.size() - 1] = false;
        }
    }

    if ((x_pos_2 < (n_x + 1)) && (x_pos_2 > -1)) {
        for (int i = lower_p; i < upper_p; i++) {
            finer_level_x[IndexNS(x_pos_2, i)] = rectangle;
            is_interrior_level_boundary_x[IndexNS(x_pos_2, i)] = true;
        }
    }

    if ((x_pos_2 < (n_x + 1)) && (x_pos_2 > 0)) {
        for (int i = lower_p; i < upper_p; i++) {
            int j = i - p_pos_2 + 1;
            rectangle->exterior_boundary_xm[j] = current;
            rectangle->is_exterrior_boundary_same_level_xm[j] = false;
        }

        if (OnLine(p_pos_2 - 1, 0, n_p - 1)) {
            rectangle->exterior_boundary_xm[0] = current;
            rectangle->is_exterrior_boundary_same_level_xm[0] = false;
        }

        if (OnLine(p_pos_2 + n_p_2, 0, n_p - 1)) {
            rectangle->exterior_boundary_xm[rectangle->exterior_boundary_xm.size() - 1] = current;
            rectangle->is_exterrior_boundary_same_level_xm[rectangle->exterior_boundary_xm.size() - 1] = false;
        }
    }

    if (((p_pos_2 + n_p_2) < (n_p + 1)) && ((p_pos_2 + n_p_2) > -1)) {
        for (int i = lower_x; i < upper_x; i++) {
            finer_level_p[IndexNS(i, p_pos_2 + n_p_2)] = rectangle;
            is_interrior_level_boundary_p[IndexNS(i, p_pos_2 + n_p_2)] = true;
        }
    }

    if (((p_pos_2 + n_p_2) < n_p) && ((p_pos_2 + n_p_2) > -1)) {
        for (int i = lower_x; i < upper_x; i++) {
            int j = i - x_pos_2;
            rectangle->exterior_boundary_pp[j] = current;
            rectangle->is_exterrior_boundary_same_level_pp[j] = false;
        }
    }

    if ((p_pos_2 < (n_p + 1)) && (p_pos_2 > -1)) {
        for (int i = lower_x; i < upper_x; i++) {
            finer_level_p[IndexNS(i, p_pos_2)] = rectangle;
            is_interrior_level_boundary_p[IndexNS(i, p_pos_2)] = true;
        }
    }

    if ((p_pos_2 < (n_p + 1)) && (p_pos_2 > 0)) {
        for (int i = lower_x; i < upper_x; i++) {
            int j = i - x_pos_2;
            rectangle->exterior_boundary_pm[j] = current;
            rectangle->is_exterrior_boundary_same_level_pm[j] = false;
        }
    }
}

void Rectangle::getError(std::vector<coords> &flaggedCells, int particleType) {
    //Identify error on current rectangle data

    int offset = 0;

    for (int i = 0; i <= settings.maxDepth - depth; i++) {
        offset += std::pow(settings.refinementRatio, i);
    }

    offset *= 3;

    int lower_x = std::max(x_pos, offset);
    int upper_x = std::min(x_pos + n_x, settings.GetXSize(depth) - offset);

    int lower_p = std::max(p_pos, offset);
    int upper_p = std::min(p_pos + n_p, settings.GetPSize(depth, particleType) - offset);

    for (int i = lower_x; i < upper_x; i++) {
        for (int j = lower_p; j < upper_p; j++) {
            if ((ErrorEstimate(i - x_pos, j - p_pos) > settings.refinementCriteria) || (settings.RefinementOverride((i + 0.5) * dx, Momentum(j - p_pos), depth, particleType))) {
                flaggedCells.push_back(std::make_pair(i, j));
            }
        }
    }
}

void Rectangle::GetDataFromCoarseLevelRectangle(const std::shared_ptr<Rectangle> &rectangle) {
    int x_pos_2 = rectangle->x_pos / refinementRatio;
    int p_pos_2 = rectangle->p_pos / refinementRatio;
    int n_x_2 = rectangle->n_x / refinementRatio;
    int n_p_2 = rectangle->n_p / refinementRatio;

    int lower_x = std::max(x_pos, x_pos_2 - 1);
    int upper_x = std::min(x_pos + n_x, x_pos_2 + n_x_2 + 1);

    int lower_p = std::max(p_pos, p_pos_2 - 1);
    int upper_p = std::min(p_pos + n_p, p_pos_2 + n_p_2 + 1);

    for (int i = lower_x; i < upper_x; i++) {
        for (int j = lower_p; j < upper_p; j++) {
            std::vector<double> interpolants = GetWenoValueFromCoarseLevel(refinementRatio * i + 1, refinementRatio * j + 1, -1);

            for (int k = 0; k < refinementRatio; k++) {
                for (int l = 0; l < refinementRatio; l++) {
                    rectangle->f[rectangle->Index3((i - x_pos_2) * refinementRatio + k, (j - p_pos_2) * refinementRatio + l, 0)] = interpolants[l + refinementRatio * k];
                    rectangle->f[rectangle->Index3((i - x_pos_2) * refinementRatio + k, (j - p_pos_2) * refinementRatio + l, 1)] = rectangle->f[rectangle->Index3((i - x_pos_2) * refinementRatio + k, (j - p_pos_2) * refinementRatio + l, 0)];

                    rectangle->is_interpolated[rectangle->IndexNS((i - x_pos_2) * refinementRatio + k, (j - p_pos_2) * refinementRatio + l)] = true;
                }
            }
        }
    }
}

void Rectangle::GetDataFromSameLevelRectangle(const std::shared_ptr<Rectangle> &rectangle) {
    int x_pos_2 = rectangle->x_pos;
    int p_pos_2 = rectangle->p_pos;

    int n_x_2 = rectangle->n_x;
    int n_p_2 = rectangle->n_p;

    int lower_x = std::max(x_pos, x_pos_2 - 2);
    int upper_x = std::min(x_pos + n_x, x_pos_2 + n_x_2 + 2);

    int lower_p = std::max(p_pos, p_pos_2 - 2);
    int upper_p = std::min(p_pos + n_p, p_pos_2 + n_p_2 + 2);

    for (int i = lower_x; i < upper_x; i++) {
        for (int j = lower_p; j < upper_p; j++) {
            rectangle->f[rectangle->Index3(i - x_pos_2, j - p_pos_2, 0)] = GetValueFromSameLevel(i, j);
            rectangle->f[rectangle->Index3(i - x_pos_2, j - p_pos_2, 1)] = rectangle->f[rectangle->Index3(i - x_pos_2, j - p_pos_2, 0)];

            rectangle->is_interpolated[rectangle->IndexNS(i - x_pos_2, j - p_pos_2)] = true;
        }
    }
}

double Rectangle::GetWenoEdgeValue(double f1, double f2, double f3, double f4, bool right) {

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

    double wL0 = wL * (0.75 + wL * (wL - 0.5));
    double wR0 = wR * (0.75 + wR * (wR - 0.5));

    if (right) {
        wL = std::max(wL0, wR0) / (wL0 + wR0);
        wR = 1 - wL;
    } else {
        wL = std::min(wL0, wR0) / (wL0 + wR0);
        wR = 1 - wL;
    }

    return (wL * (fL) + wR * (fR));
}

double Rectangle::RGKGetFluxX(int i, int j) {
    if (!is_interrior_level_boundary_x[IndexNS(i, j)]) {
        double w3 = 1.0 / 48.0;

        double dx_inv = 1 / dx;
        double dp_inv = 1 / dp;

        double q = settings.GetCharge(particleType);
        double cc = cs * cs * settings.GetMass(particleType);
        double as = q * q * EMSolver->GetASquared(GetFinestIndex(i));
        double am =  dp_inv * cc * (Gamma(Momentum(j + 1), as) - Gamma(Momentum(j), as));
        double ap1 = dp_inv * cc * (Gamma(Momentum(j + 2.0), as) - Gamma(Momentum(j + 1.0), as));
        double am1 = dp_inv * cc * (Gamma(Momentum(j), as) - Gamma(Momentum(j - 1.0), as));
        double fm = GetWenoEdgeValue(f[Index3(i - 2, j, 1)], f[Index3(i - 1, j, 1)], f[Index3(i, j, 1)], f[Index3(i + 1, j, 1)], am > 0.0);
        double fp1 = GetWenoEdgeValue(f[Index3(i - 2, j + 1, 1)], f[Index3(i - 1, j + 1, 1)], f[Index3(i, j + 1, 1)], f[Index3(i + 1, j + 1, 1)], ap1 > 0.0);
        double fm1 = GetWenoEdgeValue(f[Index3(i - 2, j - 1, 1)], f[Index3(i - 1, j - 1, 1)], f[Index3(i, j - 1, 1)], f[Index3(i + 1, j - 1, 1)], am1 > 0.0);

        return dx_inv * (fm * am + w3 * (fp1 - fm1) * (ap1 - am1));
    } else {
        return finer_level_x[IndexNS(i, j)]->CalculateFluxToCoarseX(x_pos + i, p_pos + j);
    }
}

double Rectangle::GetEfield(int i) {
    int j = GetFinestIndex(i);

    double temp(0.0);

    for (int k = 0; k < relativeToBottom; k++) {
        temp += EMSolver->GetEfield(j + k);
    }

    temp *= (1.0 / relativeToBottom);

    return temp;
}

double Rectangle::RGKGetFluxP(int i, int j) {
    if (!is_interrior_level_boundary_p[IndexNS(i, j)]) {
        double w3 = 1 / 48.0;

        double dp_inv = 1 / dp;
        double dx_inv = 1 / dx;

        double cc = cs * cs * settings.GetMass(particleType);
        double q = settings.GetCharge(particleType);
        double as_0 = q * q * EMSolver->GetASquared(GetFinestIndex(i - 1));
        double as_1 = q * q * EMSolver->GetASquared(GetFinestIndex(i));
        double as_2 = q * q * EMSolver->GetASquared(GetFinestIndex(i + 1));
        double as_3 = q * q * EMSolver->GetASquared(GetFinestIndex(i + 2));
        double Em = q * GetEfield(i);
        double Ep1 = q * GetEfield(i + 1);
        double Em1 = q * GetEfield(i - 1);
        double momentum = Momentum(j);

        double am =Em - cc * dx_inv * (Gamma(momentum, as_2) - Gamma(momentum, as_1));
        double ap1 = Ep1 - cc * dx_inv * (Gamma(momentum, as_3) - Gamma(momentum, as_2));
        double am1 = Em1 - cc * dx_inv * (Gamma(momentum, as_1) - Gamma(momentum, as_0));
        double fm = GetWenoEdgeValue(f[Index3(i, j - 2, 1)], f[Index3(i, j - 1, 1)], f[Index3(i, j, 1)], f[Index3(i, j + 1, 1)], am > 0.0);
        double fp1 = GetWenoEdgeValue(f[Index3(i + 1, j - 2, 1)], f[Index3(i + 1, j - 1, 1)], f[Index3(i + 1, j, 1)], f[Index3(i + 1, j + 1, 1)], ap1 > 0.0);
        double fm1 = GetWenoEdgeValue(f[Index3(i - 1, j - 2, 1)], f[Index3(i - 1, j - 1, 1)], f[Index3(i - 1, j, 1)], f[Index3(i - 1, j + 1, 1)], am1 > 0.0);

        return dp_inv * (fm * am + w3 * (fp1 - fm1) * (ap1 - am1));
    } else {
        return finer_level_p[IndexNS(i, j)]->CalculateFluxToCoarseP(x_pos + i, p_pos + j);
    }
}

void Rectangle::GetDataFromCoarseNewLevelRectangle(const std::shared_ptr<Rectangle> &rectangle) {
    // Here, we make a restriction to refinementratio 2.

    int x_pos_2 = rectangle->x_pos / refinementRatio;
    int p_pos_2 = rectangle->p_pos / refinementRatio;
    int n_x_2 = rectangle->n_x / refinementRatio;
    int n_p_2 = rectangle->n_p / refinementRatio;

    int lower_x = std::max(x_pos, x_pos_2 - 1);
    int upper_x = std::min(x_pos + n_x, x_pos_2 + n_x_2 + 1);

    int lower_p = std::max(p_pos, p_pos_2 - 1);
    int upper_p = std::min(p_pos + n_p, p_pos_2 + n_p_2 + 1);

    for (int i = lower_x; i < upper_x; i++) {
        for (int j = lower_p; j < upper_p; j++) {
            if (!rectangle->is_interpolated[rectangle->IndexNS(refinementRatio * (i - x_pos_2) + 1, refinementRatio * (j - p_pos_2) + 1)]) {
                std::vector<double> interpolants = GetWenoValueFromCoarseLevel(refinementRatio * i + 1, refinementRatio * j + 1, -1);

                for (int k = 0; k < refinementRatio; k++) {
                    for (int l = 0; l < refinementRatio; l++) {
                        rectangle->f[rectangle->Index3((i - x_pos_2) * refinementRatio + k, (j - p_pos_2) * refinementRatio + l, 0)] = interpolants[l + refinementRatio * k];
                        rectangle->f[rectangle->Index3((i - x_pos_2) * refinementRatio + k, (j - p_pos_2) * refinementRatio + l, 1)] = rectangle->f[rectangle->Index3((i - x_pos_2) * refinementRatio + k, (j - p_pos_2) * refinementRatio + l, 0)];
                    }
                }
            }
        }
    }
}

void Rectangle::UpdateCornerPoints(int val) {
    {
        int i = -1;

        if (auto spt = exterior_boundary_xm[i + 1].lock()) {
            if (!is_exterrior_boundary_same_level_xm[i + 1]) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos - 1, p_pos + i * refinementRatio, 3,val);
                for (int j = refinementRatio - 2; j < refinementRatio; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(-2, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                for (int j = refinementRatio - 2; j < refinementRatio; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 1, p_pos + i * refinementRatio + j,val);
                    f[Index3(-2, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 2, p_pos + i * refinementRatio + j,val);
                }
            }
        } else {
            std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
        }
    }

    {
        int i = is_exterrior_boundary_same_level_xm.size() - 2;

        if (auto spt = exterior_boundary_xm[i + 1].lock()) {
            if (!is_exterrior_boundary_same_level_xm[i + 1]) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos - 1, p_pos + i * refinementRatio, 3,val);
                for (int j = 0; j < 2; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(-2, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                for (int j = 0; j < 2; j++) {
                    f[Index3(-1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 1, p_pos + i * refinementRatio + j,val);
                    f[Index3(-2, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos - 2, p_pos + i * refinementRatio + j,val);
                }
            }
        } else {
            std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
        }
    }

    {
        int i = -1;

        if (auto spt = exterior_boundary_xp[i + 1].lock()) {
            if (!is_exterrior_boundary_same_level_xp[i + 1]) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos + n_x, p_pos + i * refinementRatio, 1,val);
                for (int j = refinementRatio - 2; j < refinementRatio; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                for (int j = refinementRatio - 2; j < refinementRatio; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x, p_pos + i * refinementRatio + j,val);
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x + 1, p_pos + i * refinementRatio + j,val);
                }
            }
        } else {
            std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
        }
    }

    {
        int i = is_exterrior_boundary_same_level_xp.size() - 2;

        if (auto spt = exterior_boundary_xp[i + 1].lock()) {
            if (!is_exterrior_boundary_same_level_xp[i + 1]) {
                std::vector<double> temp = spt->GetWenoValueFromCoarseLevel(x_pos + n_x, p_pos + i * refinementRatio, 1,val);
                for (int j = 0; j < 2; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = temp[2 * j];
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = temp[2 * j + 1];
                }
            } else {
                for (int j = 0; j < 2; j++) {
                    f[Index3(n_x, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x, p_pos + i * refinementRatio + j,val);
                    f[Index3(n_x + 1, i * refinementRatio + j, val)] = spt->GetValueFromSameLevel(x_pos + n_x + 1, p_pos + i * refinementRatio + j,val);
                }
            }
        } else {
            std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
        }
    }
}

double Rectangle::RGKGetFluxPL(int i, int j) {
    if (!is_interrior_level_boundary_p[IndexNS(i, j)]) {
        double dp_inv = 1 / dp;
        double dx_inv = 1 / dx;

        double cc = cs * cs * settings.GetMass(particleType);
        double q = settings.GetCharge(particleType);
        double as_1 = q * q * EMSolver->GetASquared(GetFinestIndex(i));
        double as_2 = q * q * EMSolver->GetASquared(GetFinestIndex(i + 1));
        double Em = q * GetEfield(i);
        double momentum = Momentum(j);

        double bm=q*q*m_inv*GetMagneticForce(i);
        double Am=q*q*GetCellAverageASquared(i);

        double am = Em - cc * dx_inv * (Gamma(momentum, as_2) - Gamma(momentum, as_1));

        return dp_inv * (am > 0.0 ? f[Index3(i,j-1,1)] : f[Index3(i,j,1)] ) * am;
    } else {
        return finer_level_p[IndexNS(i, j)]->CalculateFluxToCoarsePL(x_pos + i, p_pos + j);
    }
}

double Rectangle::RGKGetFluxXL(int i, int j) {
    if (!is_interrior_level_boundary_x[IndexNS(i, j)]) {
        double dx_inv = 1 / dx;
        double dp_inv = 1 / dp;

        double q = settings.GetCharge(particleType);
        double cc = cs * cs * settings.GetMass(particleType);
        double as = q * q * EMSolver->GetASquared(GetFinestIndex(i));
        double am = dp_inv * cc * (Gamma(Momentum(j + 1), as) - Gamma(Momentum(j), as));

        return dx_inv * (am > 0.0 ? f[Index3(i-1,j,1)] : f[Index3(i,j,1)] )* am ;
    } else {
        return finer_level_x[IndexNS(i, j)]->CalculateFluxToCoarseXL(x_pos + i, p_pos + j);
    }
}

void Rectangle::FCTTimeStep(double timestep,int step,int subStep) {

    int xm = left ? 1 : 0;
    int xp = right ? n_x : (n_x + 1);

    int pp = up ? n_p : (n_p + 1);
    int pm = down ? 1 : 0;
    int nx1 = n_x + 1;
    int nx2 = n_x + 2;
    int np1 = n_p + 1;
    int np2 = n_p + 2;

    if (subStep==0) {

        double q = settings.GetCharge(particleType);
        double w3 = 1 / 48.0;
        double dx_inv = 1 / dx;
        double dp_inv = 1 / dp;
        double cc = cs * cs * settings.GetMass(particleType);

        for (int i = -1; i < nx2; i++) {
            double as = q * q * EMSolver->GetASquared(GetFinestIndex(i));
            #pragma ivdep
            for (int j = -1; j < np1; j++) {
                double am =dp_inv * cc * (Gamma(Momentum(j + 1), as) - Gamma(Momentum(j), as));
                ex[IndexNS(i, j)] = am;
            }
        }

        for (int i = 0; i < nx1; i++) {
            for (int j = -1; j < np1; j++) {
                fx[IndexNS(i, j)] = GetWenoEdgeValue(f[Index3(i - 2, j, 1)], f[Index3(i - 1, j, 1)], f[Index3(i, j, 1)], f[Index3(i + 1, j, 1)], ex[IndexNS(i, j)] > 0.0);
            }
        }

        for (int i = -1; i < nx1; i++) {
            double as_1 = q * q * EMSolver->GetASquared(GetFinestIndex(i));
            double as_2 = q * q * EMSolver->GetASquared(GetFinestIndex(i + 1));
            double Em = q * GetEfield(i);
            #pragma ivdep
            for (int j = -1; j < np2; j++) {
                double momentum = Momentum(j);
                double am =Em - cc * dx_inv * (Gamma(momentum, as_2) - Gamma(momentum, as_1));
                ep[IndexNS(i, j)] = am;
            }
        }

        for (int i = -1; i < nx1; i++) {
            for (int j = 0; j < np1; j++) {
                fp[IndexNS(i, j)] = GetWenoEdgeValue(f[Index3(i, j - 2, 1)], f[Index3(i, j - 1, 1)], f[Index3(i, j, 1)], f[Index3(i, j + 1, 1)], ep[IndexNS(i, j)] > 0.0);
            }
        }

        for (int i = 0; i < nx1; i++) {
            for (int j = -1; j < np1; j++) {
                double fluxH(0.0);

                if (!is_interrior_level_boundary_x[IndexNS(i, j)]) {
                    double am = ex[IndexNS(i, j)];
                    double ap1 = ex[IndexNS(i, j + 1)];
                    double am1 = ex[IndexNS(i, j - 1)];
                    double fm = fx[IndexNS(i, j)];
                    double fp1 = fx[IndexNS(i, j + 1)];
                    double fm1 = fx[IndexNS(i, j - 1)];

                    fluxH = dx_inv * (fm * am + w3 * (fp1 - fm1) * (ap1 - am1));

                } else {
                    fluxH = finer_level_x[IndexNS(i, j)]->CalculateFluxToCoarseX(x_pos + i, p_pos + j);
                }

                FxH[Index6(i, j,step)]=fluxH;
            }
        }
        for (int i = 0; i < nx1; i++) {
            for (int j = -1; j < np1; j++) {
                double fluxL(0.0);

                if (!is_interrior_level_boundary_x[IndexNS(i, j)]) {
                    double am = ex[IndexNS(i, j)];

                    fluxL = dx_inv * ((am>0.0 ? f[Index3(i-1,j,1)] : f[Index3(i,j,1)]) * am );

                } else {

                    fluxL = finer_level_x[IndexNS(i, j)]->CalculateFluxToCoarseXL(x_pos + i, p_pos + j);
                }

                FxL[Index6(i, j,step)]=fluxL;
            }
        }

        for (int i = -1; i < nx1; i++) {
            for (int j = 0; j < np1; j++) {
                double fluxH(0.0);

                if (!is_interrior_level_boundary_p[IndexNS(i, j)]) {
                    double am = ep[IndexNS(i, j)];
                    double ap1 = ep[IndexNS(i + 1, j)];
                    double am1 = ep[IndexNS(i - 1, j)];
                    double fm = fp[IndexNS(i, j)];
                    double fp1 = fp[IndexNS(i + 1, j)];
                    double fm1 = fp[IndexNS(i - 1, j)];

                    fluxH = dp_inv * (fm * am + w3 * (fp1 - fm1) * (ap1 - am1));

                } else {
                    fluxH = finer_level_p[IndexNS(i, j)]->CalculateFluxToCoarseP(x_pos + i, p_pos + j);
                }

                FpH[Index6(i, j,step)]=fluxH;

            }
        }

        for (int i =-1; i < nx1; i++) {
            for (int j = 0; j < np1; j++) {
                double fluxL(0.0);

                if (!is_interrior_level_boundary_p[IndexNS(i, j)]) {
                    double am = ep[IndexNS(i, j)];

                    fluxL = dp_inv * ((am>0.0 ? f[Index3(i,j-1,1)] : f[Index3(i,j,1)]) * am );

                } else {

                    fluxL = finer_level_p[IndexNS(i, j)]->CalculateFluxToCoarsePL(x_pos + i, p_pos + j);
                }

                FpL[Index6(i, j,step)]=fluxL;

            }
        }

        if (step == 0) {
            double a0 = 0.5 * timestep;

            for (int i = -2; i < nx2; i++) {
                for (int j = -2; j < np2; j++) {
                    int idx = IndexNS(i, j);
                    int idx0 = Index6(i, j, 0);
                    FxLS[idx] = a0 * FxL[idx0];
                    FpLS[idx] = a0 * FpL[idx0];
                    FxDS[idx] = a0 * FxH[idx0]-FxLS[idx];
                    FpDS[idx] = a0 * FpH[idx0]-FpLS[idx];
                }
            }
        } else if (step == 1) {
            double a0 = 0.221776 * timestep;
            double a1 = 0.110224 * timestep;

            for (int i = -2; i < nx2; i++) {
                for (int j = -2; j < np2; j++) {
                    int idx = IndexNS(i, j);
                    int idx0 = Index6(i, j, 0);
                    int idx1 = Index6(i, j, 1);
                    FxLS[idx] = a0 * FxL[idx0] + a1*FxL[idx1];
                    FpLS[idx] = a0 * FpL[idx0] + a1 * FpL[idx1];
                    FxDS[idx] = a0 * FxH[idx0]+ a1*FxH[idx1]-FxLS[idx];
                    FpDS[idx] = a0 * FpH[idx0]+a1 * FpH[idx1]-FpLS[idx];
                }
            }
        } else if (step == 2) {
            double a0 = -0.04884659515311857 * timestep;
            double a1 = -0.17772065232640102 * timestep;
            double a2 = 0.8465672474795197 * timestep;

            for (int i = -2; i < nx2; i++) {
                for (int j = -2; j < np2; j++) {
                    int idx = IndexNS(i, j);
                    int idx0 = Index6(i, j, 0);
                    int idx1 = Index6(i, j, 1);
                    int idx2 = Index6(i, j, 2);
                    FxLS[idx] = a0 * FxL[idx0] + a1*FxL[idx1]+a2*FxL[idx2];
                    FpLS[idx] = a0 * FpL[idx0] + a1 * FpL[idx1]+ a2 * FpL[idx2];
                    FxDS[idx] = a0 * FxH[idx0]+ a1*FxH[idx1]+a2*FxH[idx2]-FxLS[idx];
                    FpDS[idx] = a0 * FpH[idx0]+a1 * FpH[idx1]+a2*FpH[idx2]-FpLS[idx];
                }
            }
        } else if (step == 3) {
            double a0 = -0.15541685842491548 * timestep;
            double a1 = -0.3567050098221991 * timestep;
            double a2 = 1.0587258798684427 * timestep;
            double a3 = 0.30339598837867193 * timestep;

            for (int i = -2; i < nx2; i++) {
                for (int j = -2; j < np2; j++) {
                    int idx = IndexNS(i, j);
                    int idx0 = Index6(i, j, 0);
                    int idx1 = Index6(i, j, 1);
                    int idx2 = Index6(i, j, 2);
                    int idx3 = Index6(i, j, 3);
                    FxLS[idx] = a0 * FxL[idx0] + a1*FxL[idx1]+a2*FxL[idx2]+a3*FxL[idx3];
                    FpLS[idx] = a0 * FpL[idx0] + a1 * FpL[idx1]+ a2 * FpL[idx2]+a3 * FpL[idx3];
                    FxDS[idx] = a0 * FxH[idx0]+ a1*FxH[idx1]+a2*FxH[idx2]+a3*FxH[idx3]-FxLS[idx];
                    FpDS[idx] = a0 * FpH[idx0]+a1 * FpH[idx1]+a2*FpH[idx2]+a3 * FpH[idx3]-FpLS[idx];
                }
            }
        } else if (step == 4) {
            double a0 = 0.2014243506726763 * timestep;
            double a1 = 0.008742057842904185 * timestep;
            double a2 = 0.15993995707168115 * timestep;
            double a3 = 0.4038290605220775 * timestep;
            double a4 = 0.22606457389066084 * timestep;

            for (int i = -2; i < nx2; i++) {
                for (int j = -2; j < np2; j++) {
                    int idx = IndexNS(i, j);
                    int idx0 = Index6(i, j, 0);
                    int idx1 = Index6(i, j, 1);
                    int idx2 = Index6(i, j, 2);
                    int idx3 = Index6(i, j, 3);
                    int idx4 = Index6(i, j, 4);
                    FxLS[idx] = a0 * FxL[idx0] + a1*FxL[idx1]+a2*FxL[idx2]+a3*FxL[idx3]+a4*FxL[idx4];
                    FpLS[idx] = a0 * FpL[idx0] + a1 * FpL[idx1]+ a2 * FpL[idx2]+a3 * FpL[idx3]+a4 * FpL[idx4];
                    FxDS[idx] = a0 * FxH[idx0]+ a1*FxH[idx1]+a2*FxH[idx2]+a3*FxH[idx3]+a4*FxH[idx4]-FxLS[idx];
                    FpDS[idx] = a0 * FpH[idx0]+a1 * FpH[idx1]+a2*FpH[idx2]+a3 * FpH[idx3]+a4 * FpH[idx4]-FpLS[idx];
                }
            }
        } else if (step==5) {

            double b0 = 0.15791629516167136 * timestep ;
            double b1 = 0.0 * timestep;
            double b2 = 0.18675894052400077 * timestep;
            double b3 = 0.6805652953093346 * timestep;
            double b4 = -0.27524053099500667 * timestep;
            double b5 = 0.25 * timestep ;

            for (int i = -2; i < nx2; i++) {
                for (int j = -2; j < np2; j++) {
                    int idx = IndexNS(i, j);
                    int idx0 = Index6(i, j, 0);
                    int idx1 = Index6(i, j, 1);
                    int idx2 = Index6(i, j, 2);
                    int idx3 = Index6(i, j, 3);
                    int idx4 = Index6(i, j, 4);
                    int idx5 = Index6(i, j, 5);
                    FxLS[idx] = b0 * FxL[idx0] + b1*FxL[idx1]+b2*FxL[idx2]+b3*FxL[idx3]+b4*FxL[idx4]+b5*FxL[idx5];
                    FpLS[idx] = b0 * FpL[idx0] + b1 * FpL[idx1]+ b2 * FpL[idx2]+b3 * FpL[idx3]+b4 * FpL[idx4]+b5 * FpL[idx5];
                    FxDS[idx] = b0 * FxH[idx0]+ b1*FxH[idx1]+b2*FxH[idx2]+b3*FxH[idx3]+b4*FxH[idx4]+b5*FxH[idx5]-FxLS[idx];
                    FpDS[idx] = b0 * FpH[idx0]+b1 * FpH[idx1]+b2*FpH[idx2]+b3 * FpH[idx3]+b4 * FpH[idx4]+b5 * FpH[idx5]-FpLS[idx];
                }
            }

        }

        for (int i=0; i<n_x; i++) {
            for (int j=0; j<n_p; j++) {
                f[Index3(i,j,2)]=f[Index3(i,j,0)];
            }
        }

        for (int i=xm; i<xp; i++) {
            for (int j=pm; j<pp; j++) {
                f[Index3(i,j,2)]+=FxLS[IndexNS(i,j)];
                f[Index3(i-1,j,2)]-=FxLS[IndexNS(i,j)];
                f[Index3(i,j,2)]+=FpLS[IndexNS(i,j)];
                f[Index3(i,j-1,2)]-=FpLS[IndexNS(i,j)];
            }
        }

    } else if (subStep==1) {
        for (int i=-1; i<nx1; i++) {
            for (int j=-1; j<np1; j++) {
                int idx = IndexNS(i,j);
                int idxij2 = Index3(i,j,2);
                double Pp=std::max(0.0,FxDS[idx])-std::min(0.0,FxDS[IndexNS(i+1,j)])+std::max(0.0,FpDS[idx])-std::min(0.0,FpDS[IndexNS(i,j+1)]);
                double Pm=std::max(0.0,FxDS[IndexNS(i+1,j)])-std::min(0.0,FxDS[idx])+std::max(0.0,FpDS[IndexNS(i,j+1)])-std::min(0.0,FpDS[idx]);

                double w1ma=std::max(f[Index3(i,j,0)],f[idxij2]);
                double w2ma=std::max(f[Index3(i+1,j,0)],f[Index3(i+1,j,2)]);
                double w3ma=std::max(f[Index3(i-1,j,0)],f[Index3(i-1,j,2)]);
                double w4ma=std::max(f[Index3(i,j+1,0)],f[Index3(i,j+1,2)]);
                double w5ma=std::max(f[Index3(i,j-1,0)],f[Index3(i,j-1,2)]);

                double wMax=std::max(w1ma,std::max(w2ma,std::max(w3ma,std::max(w4ma,w5ma))));

                double w1mi=std::min(f[Index3(i,j,0)],f[idxij2]);
                double w2mi=std::min(f[Index3(i+1,j,0)],f[Index3(i+1,j,2)]);
                double w3mi=std::min(f[Index3(i-1,j,0)],f[Index3(i-1,j,2)]);
                double w4mi=std::min(f[Index3(i,j+1,0)],f[Index3(i,j+1,2)]);
                double w5mi=std::min(f[Index3(i,j-1,0)],f[Index3(i,j-1,2)]);

                double wMin=std::min(w1mi,std::min(w2mi,std::min(w3mi,std::min(w4mi,w5mi))));

                double Qm=-wMin+f[idxij2];
                double Qp=wMax-f[idxij2];

                Rp[idx]=Pp>0.0?std::min(1.0,Qp/Pp):0.0;
                Rm[idx]=Pm>0.0?std::min(1.0,Qm/Pm):0.0;
            }
        }

        std::fill(Cp.begin(), Cp.end(), 1.0);
        std::fill(Cx.begin(), Cx.end(), 1.0);

        for (int i=1; i<n_x; i++) {
            for (int j=0; j<n_p; j++) {
                Cx[IndexNS(i,j)]=FxDS[IndexNS(i,j)]>0.0?std::min(Rp[IndexNS(i,j)],Rm[IndexNS(i-1,j)]):std::min(Rp[IndexNS(i-1,j)],Rm[IndexNS(i,j)]);
            }
        }

        for (int i=0; i<n_x; i++) {
            for (int j=1; j<n_p; j++) {
                Cp[IndexNS(i,j)]=FpDS[IndexNS(i,j)]>0.0?std::min(Rp[IndexNS(i,j)],Rm[IndexNS(i,j-1)]):std::min(Rp[IndexNS(i,j-1)],Rm[IndexNS(i,j)]);
            }
        }
    } else if (subStep==2) {
        for (int i=0; i<n_x; i++) {
            for (int j=0; j<n_p; j++) {
                f[Index3(i,j,1)]=f[Index3(i,j,2)];
            }
        }

        for (int i=xm; i<xp; i++) {
            for (int j=pm; j<pp; j++) {

                    f[Index3(i,j,1)]+=Cx[IndexNS(i,j)]*FxDS[IndexNS(i,j)];
                    f[Index3(i-1,j,1)]-=Cx[IndexNS(i,j)]*FxDS[IndexNS(i,j)];
                    f[Index3(i,j,1)]+=Cp[IndexNS(i,j)]*FpDS[IndexNS(i,j)];
                    f[Index3(i,j-1,1)]-=Cp[IndexNS(i,j)]*FpDS[IndexNS(i,j)];

            }
        }

    } else if (subStep==3) {
        for (int i=-2; i<nx2; i++) {
            for (int j=-2; j<np2; j++) {
                f[Index3(i,j,0)]=f[Index3(i,j,1)];
            }
        }
    }
}

void Rectangle::SetCFromDifferentLevel(int i, int j,Rectangle *rectangle,int t) {

    /* o is rectangle

     0 x | o
     1 o | x
     2 o / x
     3 x / o

     */

    int i_caller=i-rectangle->x_pos;
    int j_caller=j-rectangle->p_pos;

    int i_coarse=i/refinementRatio-x_pos;
    int j_coarse=j/refinementRatio-p_pos;

    if (t==0) {
        double cx=Cx[IndexNS(i_coarse,j_coarse)];

        for (int k=0; k<refinementRatio; k++) {
            cx=std::min(cx,rectangle->FxDS[rectangle->IndexNS(i_caller,j_caller+k)]>0.0?rectangle->Rp[rectangle->IndexNS(i_caller,j_caller+k)]:rectangle->Rm[rectangle->IndexNS(i_caller,j_caller+k)]);
        }

        cx=std::min(cx,FxDS[rectangle->IndexNS(i_coarse,j_coarse)]>0.0?Rm[IndexNS(i_coarse-1,j_coarse)]:Rp[IndexNS(i_coarse-1,j_coarse)]);

        for (int k=0; k<refinementRatio; k++) {
            rectangle->Cx[rectangle->IndexNS(i_caller,j_caller+k)]=cx;
        }

        Cx[IndexNS(i_coarse,j_coarse)]=cx;

    } else if(t==1){
        double cx=Cx[IndexNS(i_coarse,j_coarse)];

        for (int k=0; k<refinementRatio; k++) {
            cx=std::min(cx,rectangle->FxDS[rectangle->IndexNS(i_caller,j_caller+k)]>0.0?rectangle->Rm[rectangle->IndexNS(i_caller-1,j_caller+k)]:rectangle->Rp[rectangle->IndexNS(i_caller-1,j_caller+k)]);
        }

        cx=std::min(cx,FxDS[rectangle->IndexNS(i_coarse,j_coarse)]>0.0?Rp[IndexNS(i_coarse,j_coarse)]:Rm[IndexNS(i_coarse,j_coarse)]);

        for (int k=0; k<refinementRatio; k++) {
            rectangle->Cx[rectangle->IndexNS(i_caller,j_caller+k)]=cx;
        }

        Cx[IndexNS(i_coarse,j_coarse)]=cx;

    } else if(t==2){
        double cp=Cp[IndexNS(i_coarse,j_coarse)];

        for (int k=0; k<refinementRatio; k++) {
            cp=std::min(cp,rectangle->FpDS[rectangle->IndexNS(i_caller+k,j_caller)]>0.0?rectangle->Rp[rectangle->IndexNS(i_caller+k,j_caller)]:rectangle->Rm[rectangle->IndexNS(i_caller+k,j_caller)]);
        }

        cp=std::min(cp,FpDS[IndexNS(i_coarse,j_coarse)]>0.0?Rm[IndexNS(i_coarse,j_coarse-1)]:Rp[IndexNS(i_coarse,j_coarse-1)]);

        for (int k=0; k<refinementRatio; k++) {
            rectangle->Cp[rectangle->IndexNS(i_caller+k,j_caller)]=cp;
        }

        Cp[IndexNS(i_coarse,j_coarse)]=cp;

    } else if(t==3){
        double cp=Cp[IndexNS(i_coarse,j_coarse)];

        for (int k=0; k<refinementRatio; k++) {
            cp=std::min(cp,rectangle->FpDS[rectangle->IndexNS(i_caller+k,j_caller)]>0.0?rectangle->Rm[rectangle->IndexNS(i_caller+k,j_caller-1)]:rectangle->Rp[rectangle->IndexNS(i_caller+k,j_caller-1)]);
        }

        cp=std::min(cp,FpDS[IndexNS(i_coarse,j_coarse)]>0.0?Rp[IndexNS(i_coarse,j_coarse)]:Rm[IndexNS(i_coarse,j_coarse)]);

        for (int k=0; k<refinementRatio; k++) {
            rectangle->Cp[rectangle->IndexNS(i_caller+k,j_caller)]=cp;
        }

        Cp[IndexNS(i_coarse,j_coarse)]=cp;

    }

}

void Rectangle::CalculateDifferentBoundaryC() {

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xm.size() - 2; i++) {

        if (!is_exterrior_boundary_same_level_xm[i + 1]) {
            if (auto spt = exterior_boundary_xm[i + 1].lock()) {

                spt->SetCFromDifferentLevel(x_pos,p_pos+refinementRatio*i,this,0);

            } else {
                std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xp.size() - 2; i++) {


        if (!is_exterrior_boundary_same_level_xp[i + 1]) {
            if (auto spt = exterior_boundary_xp[i + 1].lock()) {

                spt->SetCFromDifferentLevel(x_pos+n_x,p_pos+refinementRatio*i,this,1);

            } else {
                std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pm.size(); i++) {

        if (!is_exterrior_boundary_same_level_pm[i]) {
            if (auto spt = exterior_boundary_pm[i].lock()) {

                spt->SetCFromDifferentLevel(x_pos+refinementRatio*i,p_pos,this,2);

            } else {

                std::cerr << "WARNING: exterior_boundary_pm[" << i << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pp.size(); i++) {

        if (!is_exterrior_boundary_same_level_pp[i]) {

            if (auto spt = exterior_boundary_pp[i].lock()) {

                spt->SetCFromDifferentLevel(x_pos+refinementRatio*i,p_pos+n_p,this,3);

            } else {

                std::cerr << "WARNING: exterior_boundary_pp[" << i << "] pointer is expired" << std::endl;
            }
        }
    }
}

double Rectangle::GetMagneticForce(int i) {

    int j=GetFinestIndex(i);

    double temp(0.0);

    for (int k=0;k<relativeToBottom; k++) {
        temp+=EMSolver->GetMagneticForce(j+k);
    }

    temp*=(1.0/relativeToBottom);

    return temp;
}

double Rectangle::GetCellAverageASquared(int i) {

    int j=GetFinestIndex(i);

    double temp(0.0);

    for (int k=0;k<relativeToBottom; k++) {
        temp+=EMSolver->GetCellAverageASquared(j+k);
    }

    temp*=(1.0/relativeToBottom);

    return temp;
}

void Rectangle::SetCFromSameLevel(int i, int j,Rectangle *rectangle,int t) {
    int i_neighbour = i - x_pos;
    int j_neighbour = j - p_pos;

    int i_caller=i-rectangle->x_pos;
    int j_caller=j-rectangle->p_pos;


    if (t==0) {
        for (int k=0; k<refinementRatio; k++) {
            double cx=FxDS[IndexNS(i_neighbour,j_neighbour+k)]>0.0?std::min(rectangle->Rp[rectangle->IndexNS(i_caller,j_caller+k)],Rm[IndexNS(i_neighbour-1,j_neighbour+k)]):std::min(Rp[IndexNS(i_neighbour-1,j_neighbour+k)],rectangle->Rm[rectangle->IndexNS(i_caller,j_caller+k)]);

            Cx[IndexNS(i_neighbour,j_neighbour+k)]=cx;
            rectangle->Cx[rectangle->IndexNS(i_caller,j_caller+k)]=cx;
        }
    } else if (t==1) {
        for (int k=0; k<refinementRatio; k++) {
            double cx=FxDS[IndexNS(i_neighbour,j_neighbour+k)]>0.0?std::min(rectangle->Rm[rectangle->IndexNS(i_caller-1,j_caller+k)],Rp[IndexNS(i_neighbour,j_neighbour+k)]):std::min(Rm[IndexNS(i_neighbour,j_neighbour+k)],rectangle->Rp[rectangle->IndexNS(i_caller-1,j_caller+k)]);

            Cx[IndexNS(i_neighbour,j_neighbour+k)]=cx;
            rectangle->Cx[rectangle->IndexNS(i_caller,j_caller+k)]=cx;
        }
    } else if (t==2) {
        for (int k=0; k<refinementRatio; k++) {
            double cp=FpDS[IndexNS(i_neighbour+k,j_neighbour)]>0.0?std::min(rectangle->Rp[rectangle->IndexNS(i_caller+k,j_caller)],Rm[IndexNS(i_neighbour+k,j_neighbour-1)]):std::min(Rp[IndexNS(i_neighbour+k,j_neighbour-1)],rectangle->Rm[rectangle->IndexNS(i_caller+k,j_caller)]);

            Cp[IndexNS(i_neighbour+k,j_neighbour)]=cp;
            rectangle->Cp[rectangle->IndexNS(i_caller+k,j_caller)]=cp;
        }
    } else if (t==3) {
        for (int k=0; k<refinementRatio; k++) {
            double cp=FpDS[IndexNS(i_neighbour+k,j_neighbour)]>0.0?std::min(rectangle->Rm[rectangle->IndexNS(i_caller+k,j_caller-1)],Rp[IndexNS(i_neighbour+k,j_neighbour)]):std::min(Rm[IndexNS(i_neighbour+k,j_neighbour)],rectangle->Rp[rectangle->IndexNS(i_caller+k,j_caller-1)]);

            Cp[IndexNS(i_neighbour+k,j_neighbour)]=cp;
            rectangle->Cp[rectangle->IndexNS(i_caller+k,j_caller)]=cp;
        }
    }

}

void Rectangle::CalculateSameBoundaryC() {

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xm.size() - 2; i++) {
        if (is_exterrior_boundary_same_level_xm[i + 1]) {
            if (auto spt = exterior_boundary_xm[i + 1].lock()) {
                spt->SetCFromSameLevel(x_pos,p_pos+i*refinementRatio,this,0);
            } else {
                std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xp.size() - 2; i++) {
        if (is_exterrior_boundary_same_level_xp[i + 1]) {
            if (auto spt = exterior_boundary_xp[i + 1].lock()) {
                spt->SetCFromSameLevel(x_pos+n_x,p_pos+i*refinementRatio,this,1);
            } else {
                std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pm.size(); i++) {
        if (is_exterrior_boundary_same_level_pm[i]) {
            if (auto spt = exterior_boundary_pm[i].lock()) {
                spt->SetCFromSameLevel(x_pos+i*refinementRatio,p_pos,this,2);
            } else {
                std::cerr << "WARNING: exterior_boundary_pm[" << i << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pp.size(); i++) {
        if (is_exterrior_boundary_same_level_pp[i]) {
            if (auto spt = exterior_boundary_pp[i].lock()) {
                spt->SetCFromSameLevel(x_pos+i*refinementRatio,p_pos+n_p,this,3);
            } else {
                std::cerr << "WARNING: exterior_boundary_pp[" << i << "] pointer is expired" << std::endl;
            }
        }
    }

}

void Rectangle::UpdateCFromSameLevel(int i, int j,Rectangle *rectangle,int t) {
    int i_neighbour = i - x_pos;
    int j_neighbour = j - p_pos;

    int i_caller=i-rectangle->x_pos;
    int j_caller=j-rectangle->p_pos;


    if (t==0) {
        for (int k=0; k<refinementRatio; k++) {
            double cx=std::min(Cx[IndexNS(i_neighbour,j_neighbour+k)],rectangle->Cx[rectangle->IndexNS(i_caller,j_caller+k)]);

            Cx[IndexNS(i_neighbour,j_neighbour+k)]=cx;
            rectangle->Cx[rectangle->IndexNS(i_caller,j_caller+k)]=cx;
        }
    } else if (t==1) {
        for (int k=0; k<refinementRatio; k++) {
            double cx=std::min(Cx[IndexNS(i_neighbour,j_neighbour+k)],rectangle->Cx[rectangle->IndexNS(i_caller,j_caller+k)]);

            Cx[IndexNS(i_neighbour,j_neighbour+k)]=cx;
            rectangle->Cx[rectangle->IndexNS(i_caller,j_caller+k)]=cx;
        }
    } else if (t==2) {
        for (int k=0; k<refinementRatio; k++) {
            double cp=std::min(Cp[IndexNS(i_neighbour+k,j_neighbour)],rectangle->Cp[rectangle->IndexNS(i_caller+k,j_caller)]);

            Cp[IndexNS(i_neighbour+k,j_neighbour)]=cp;
            rectangle->Cp[rectangle->IndexNS(i_caller+k,j_caller)]=cp;
        }
    } else if (t==3) {
        for (int k=0; k<refinementRatio; k++) {
            double cp=std::min(Cp[IndexNS(i_neighbour+k,j_neighbour)],rectangle->Cp[rectangle->IndexNS(i_caller+k,j_caller)]);

            Cp[IndexNS(i_neighbour+k,j_neighbour)]=cp;
            rectangle->Cp[rectangle->IndexNS(i_caller+k,j_caller)]=cp;
        }
    }

}

void Rectangle::UpdateSameBoundaryC() {

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xm.size() - 2; i++) {
        if (is_exterrior_boundary_same_level_xm[i + 1]) {
            if (auto spt = exterior_boundary_xm[i + 1].lock()) {
                spt->UpdateCFromSameLevel(x_pos,p_pos+i*refinementRatio,this,0);
            } else {
                std::cerr << "WARNING: exterior_boundary_xm[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_xp.size() - 2; i++) {
        if (is_exterrior_boundary_same_level_xp[i + 1]) {
            if (auto spt = exterior_boundary_xp[i + 1].lock()) {
                spt->UpdateCFromSameLevel(x_pos+n_x,p_pos+i*refinementRatio,this,1);
            } else {
                std::cerr << "WARNING: exterior_boundary_xp[" << i + 1 << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pm.size(); i++) {
        if (is_exterrior_boundary_same_level_pm[i]) {
            if (auto spt = exterior_boundary_pm[i].lock()) {
                spt->UpdateCFromSameLevel(x_pos+i*refinementRatio,p_pos,this,2);
            } else {
                std::cerr << "WARNING: exterior_boundary_pm[" << i << "] pointer is expired" << std::endl;
            }
        }
    }

    for (unsigned int i = 0; i < is_exterrior_boundary_same_level_pp.size(); i++) {
        if (is_exterrior_boundary_same_level_pp[i]) {
            if (auto spt = exterior_boundary_pp[i].lock()) {
                spt->UpdateCFromSameLevel(x_pos+i*refinementRatio,p_pos+n_p,this,3);
            } else {
                std::cerr << "WARNING: exterior_boundary_pp[" << i << "] pointer is expired" << std::endl;
            }
        }
    }
}
