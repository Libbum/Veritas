#include "veritas.hpp"
#include "gitsha1.hpp"
#include "Settings.hpp"

Settings::Settings(const Input &grid, const Particles &particles, const Output &out) {
    output = out;
    int maxleft = 29, maxright = 25;

    //Set input values to settings
    dx = grid.dx;
    dp = particles.dp;
    maxDepth = grid.Lfinest - 1;
    x_size = grid.nx;
    p_size = particles.np;
    refinementRatio = grid.r;
    refinementCriteria = grid.refinementCriteria;
    minEfficiency = grid.minEfficiency;
    time = 0.0;
    preLength = grid.preLength;
    postLength = grid.postLength;
    cfl = grid.cfl;
    sizeWeight = grid.sizeWeight;
    loadBalanceLength = grid.loadBalanceLength;
    m = particles.mass;
    q = particles.charge;
    tempEM = grid.tempEM;
    temp = particles.misc;
    pmin = particles.pmin;

    //Set threading counters
    unsigned int maxThreads = omp_get_max_threads();
    numMeshes = std::min((unsigned)q.size(), maxThreads);
    numThreads = std::max(maxThreads / numMeshes, (unsigned)1);

    //User may override some of these
    settingsOverride();

    //Hard failures
    if (refinementRatio % 2) {
        std::cerr << "ERROR: Mesh refinement ratio 'refinementRatio' must be an even value." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (loadBalanceLength % 2) {
        std::cerr << "ERROR: 'loadBalanceLength' must be an even value." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (x_size % refinementRatio != 0) {
        std::cerr << "ERROR: 'x_size' must be completely divisible by the Mesh refinement ratio 'refinementRatio'." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (unsigned int p = 0; p < p_size.size(); p++) {
        if (p_size[p] % refinementRatio != 0) {
            std::cerr << "ERROR: 'p_size' at index (" << p << ") must be completely divisible by the Mesh refinement ratio 'refinementRatio'." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    //Extras or fixes after overrides.
    x_size_finest = x_size * std::pow(refinementRatio, maxDepth);
    for (unsigned int i=0; i<p_size.size();i++) {
        p_size_finest.push_back(p_size[i] * std::pow(refinementRatio, maxDepth));
    }
    fMax.resize(q.size(), 0.0);
    DetermineMaximum();

    //Output run information
    title("VERITAS. Current build SHA1: " + std::string(GIT_SHA1), ' ');
    std::cout << std::endl;

    title(" Configuration Parameters ", '=');

    if (std::any_of(particles.pmin.begin(), particles.pmin.end(), [](int i){return i>=0;})) maxleft--;

    std::cout << std::setfill(' ') << std::setprecision(5) << std::setw(maxleft) << std::left << "Minimum Efficiency: " << std::setw(OUTW / 2 - maxleft) << minEfficiency;
    std::cout << std::setw(maxright) << "Mesh refinement ratio: " << std::setw(OUTW / 2 - maxright) << refinementRatio << std::endl;

    std::cout << std::setw(maxleft) << "Maximum number of levels: " << std::setw(OUTW / 2 - maxleft) << maxDepth + 1;
    std::cout << std::setw(maxright) << "Error tolerance: " << std::setw(OUTW / 2 - maxright) << refinementCriteria << std::endl;

    std::cout << std::setw(maxleft) << "nx: " << std::setw(OUTW / 2 - maxleft) << x_size;
    std::cout << std::setw(maxright) << "PreLength: " << std::setw(OUTW / 2 - maxright) << preLength << std::endl;


    std::cout << std::setw(maxleft) << "dx: " << std::setw(OUTW / 2 - maxleft) << dx;
    std::cout << std::setw(maxright) << "PostLength: " << std::setw(OUTW / 2 - maxright) << postLength << std::endl;

    std::cout << std::setw(maxleft) << "Mesh Threads: " << std::setw(OUTW / 2 - maxleft) << numMeshes;
    std::cout << std::setw(maxright) << "Compute Threads: " << std::setw(OUTW / 2 - maxright) << numThreads << std::endl;

    title(" Particle Properties ", '=');
    unsigned int numSpecies = m.size();

    for (unsigned int i = 0; i < numSpecies; ++i) {
        std::cout << std::setfill(' ') << std::setw(maxleft) << "np: " << std::setw(OUTW / 2 - maxleft) << p_size[i];
        std::cout << std::setw(maxright) << "Mass: " << std::setw(OUTW / 2 - maxright) << m[i] << std::endl;

        std::cout << std::setw(maxleft) << "dp: " << std::setw(OUTW / 2 - maxleft) << dp[i];
        std::cout << std::setw(maxright) << "Charge: " << std::setw(OUTW / 2 - maxright) << q[i] << std::endl;

        std::cout << std::setw((particles.pmin[i] >= 0) ? maxleft : maxleft - 1) << "Momentum Cutoff: " << std::setw(OUTW / 2 - maxleft) << pmin[i] << std::endl;

        unsigned int mcount = temp[i].size();
        std::cout << std::setw(maxleft) << "Misc values: ";

        for (unsigned int j = 0; j < mcount; ++j) {
            if (j < mcount - 1) {
                std::cout << temp[i][j] << "; ";
            } else {
                std::cout << temp[i][j] << std::endl;
            }
        }

        if ((numSpecies > 1) && (i < numSpecies-1)) {
            std::cout << std::setfill('-') << std::setw(80) << "-" << std::endl;
        }
    }

    title(" Output Settings ", '=');
    std::cout << std::setfill(' ') << std::setw(maxleft) << "Rectangle Data: " << std::setw(OUTW / 2 - maxleft) << std::boolalpha << output.rectangleData;
    std::cout << std::setw(maxright) << "Charge: " << std::setw(OUTW / 2 - maxright) << std::boolalpha << output.charge << std::endl;

    std::cout << std::setw(maxleft) << "Longitudinal ElectricField: " << std::setw(OUTW / 2 - maxleft) << std::boolalpha << output.EFieldLongitudinal;
    std::cout << std::setw(maxright) << "Energy: " << std::setw(OUTW / 2 - maxright) << std::boolalpha << output.energy << std::endl;

    std::cout << std::setw(maxleft) << "Transverse ElectricField: " << std::setw(OUTW / 2 - maxleft) << std::boolalpha << output.EFieldTransverse;
    std::cout << std::setw(maxright) << "Electrostatic Potential: " << std::setw(OUTW / 2 - maxright) << std::boolalpha << output.potential << std::endl;

    std::cout << std::setw(maxleft) << "Transverse MagniticField: " << std::setw(OUTW / 2 - maxleft) << std::boolalpha << output.BFieldTransverse;
    std::cout << std::setw(maxright) << "Squared Vector Potental: " << std::setw(OUTW / 2 - maxright) << std::boolalpha << output.AFieldSquared << std::endl;

    std::cout << std::setw(maxleft) << "Output Precision: " << std::setw(OUTW / 2 - maxleft) << output.precision << std::endl;

    std::cout << std::setfill('=') << std::setw(80) << "=" << std::endl;
    std::cout << std::setfill(' ') << std::endl << std::setprecision(6);
}

void Settings::title(std::string output, char spacer) {
    int space = OUTW - output.length();

    if (space % 2) {
        std::cout << std::setfill(spacer) << std::setw(floor(space / 2)) << spacer << output << std::setw(floor(space / 2) + 1) << spacer << std::endl;
    } else {
        std::cout << std::setfill(spacer) << std::setw(space / 2) << spacer << output << std::setw(space / 2) << spacer << std::endl;
    }
}

double Settings::GetDp(int level, int particleType) {
    return std::pow(refinementRatio, level - maxDepth) * dp.at(particleType);
}

double Settings::GetDx(int level) {
    return std::pow(refinementRatio, level - maxDepth) * dx;
}

int Settings::GetXSize(int level) {
    return x_size * std::pow(refinementRatio, -level + maxDepth);
}

int Settings::GetPSize(int level, int particleType) {
    return p_size.at(particleType) * std::pow(refinementRatio, -level + maxDepth);
}

double Settings::GetMass(int i) {
    return m.at(i);
}

double Settings::GetCharge(int i) {
    return q.at(i);
}

double Settings::GetfMax(int i) {
    return fMax.at(i);
}

void Settings::UpdateTime(int step, double dt) {
    if (step == 0) {
    } else if (step == 1) {
        time += (0.5 * dt);
    } else if (step == 2) {
        time += (0.332 - 0.5) * dt;
    } else if (step == 3) {
        time += (0.62 - 0.332) * dt;
    } else if (step == 4) {
        time += (0.85 - 0.62) * dt;
    } else if (step == 5) {
        time += (1.0 - 0.85) * dt;
    }
}

void Settings::DetermineMaximum() {
    double dx = GetDx(0);

    double mValue(0.0);

    for (unsigned int j = 0; j < q.size(); j++) {
        for (int i = 0; i < GetXSize(0); i++) {
            double temp(InitialDistribution((i + 0.5) * dx, 0, j));

            mValue = temp > mValue ? temp : mValue;
        }

        fMax[j] = mValue;
    }
}
