#include "veritas.hpp"
#include "Mesh.hpp"
#include "Settings.hpp"
#include "Level.hpp"
#include "Rectangle.hpp"
#include "BoundaryCondition.hpp"

//Override std::ostream for coordinates
std::ostream & operator<<(std::ostream& o, const coords& p) {
    return o << "(" << p.first << ", " << p.second << ")";
}

//Override std::ostream for rectangles
std::ostream& operator<<(std::ostream& o, const rect& r) {
    return o << r.first.first << " " << r.first.second << " " << r.second.first << " " << r.second.second;
}

//Override std::ostream for rectangle data
std::ostream& operator<<(std::ostream& o, const rectData& data) {
    int one = std::get<0>(data), two = std::get<1>(data);
    double three = std::get<2>(data);
    std::stringstream output;

    output << one << " " << two << " " << std::scientific << three;
    return o << output.str();
}

Mesh::Mesh(int particleType, Settings &settings) : settings(settings), EMSolver(NULL), particleType(particleType) {
    bc = std::make_shared<BoundaryCondition>();
    std::vector<coords> flaggedCells;
    rect extrema;
    level identified;    //identified rectangles and list of flagged coords for each rectangle

    for (int l = 0; l <= settings.maxDepth; l++) {
        //For each level, generate flagged cells from error count.
        if (l > 0) {
            //We're looking at a non ground level in the hierarchy.
            if (flaggedCells.size() > 0) {
                //Need an l+1 level.
                getExtrema(extrema, flaggedCells);
                identified = splitRectangle(extrema, flaggedCells, settings.minEfficiency);
                interpRectanglesUp(identified, l - 1); //Where l-1 is the level id that identfied the rectangles
                hierarchy.push_back(identified);
            }             //else, nothing needs to happen here as we are creating. On update this needs to kill the level.

            //get error for next level (if needed)
            if (l <= settings.maxDepth) {
                flaggedCells.clear();
                getError(l, true, flaggedCells);
            }
        } else {
            //We're on the ground. One rectangle, which is ultimately the extrema and all the contained points.
            int xmax = settings.GetXSize(settings.maxDepth), pmax = settings.GetPSize(settings.maxDepth, particleType);

            extrema = std::make_pair(std::make_pair(0, 0), std::make_pair(xmax, pmax));

            identified.push_back(extrema);
            loadBalanceRectangles(identified);
            hierarchy.push_back(identified);

            getError(l, true, flaggedCells); //Flag for new level.
        }
    }

    promoteHierarchyToMesh(true);

    for (int i = 0; i < settings.maxDepth; i++) {
        PushData();
        updateHierarchy(true);
    }
}

void Mesh::InterpolateRhoAndJToFinestMesh(std::vector<double> &charge, std::vector<double> &J) {
    for (auto &lvl: levels) {
        lvl->InterpolateRhoAndJToFinestMesh(charge, J);
    }
}

void Mesh::InterpolateEnergyToFinestMesh(std::vector<double> &energy) {
    for (auto &lvl: levels) {
        lvl->InterpolateEnergyToFinestMesh(energy);
    }
}

void Mesh::Advance(double timeStep, int step) {

    for (int i = levels.size() - 1; i > -1; i--) {
        levels[i]->FCTTimeStep(timeStep, step,0);
    }

    PushData(2);

    for (int i = levels.size() - 1; i > -1; i--) {
        levels[i]->FCTTimeStep(timeStep, step,1);
    }

    PushBoundaryC();

    for (int i = levels.size() - 1; i > -1; i--) {
        levels[i]->FCTTimeStep(timeStep, step,2);
    }

    PushData(1);

    if (step==5) {
        for (int i = levels.size() - 1; i > -1; i--) {
            levels[i]->FCTTimeStep(timeStep, step,3);
        }
    }

}

void Mesh::PushData(int val) {

    levels.at(0)->PushData(1,val);

    for (unsigned int i = 1; i < levels.size(); i++) {
        levels.at(i)->PushData(0,val);
        levels.at(i)->PushData(1,val);
    }

    levels.at(levels.size() - 1)->PushData(3,val);

    for (unsigned int i = levels.size() - 1; i > 0; i--) {
        levels.at(i - 1)->PushData(2,val);
    }

}

void Mesh::SetFieldSolver(EMFieldSolver *solver) {
    EMSolver = solver;

    for (auto &lvl : levels) {
        lvl->SetFieldSolver(EMSolver);
    }
}

void Mesh::InterMeshDataTransfer(const std::vector<std::unique_ptr<Level> > &levels_n) {
    for (unsigned int i = 0; i < levels_n.size() - 1; i++) {
        levels[i]->GetDataFromCoarserLevel(levels_n[i + 1]);
        levels[i]->GetDataFromSameLevel(levels_n[i]);
    }

    if (levels_n.size() > 1) {
        int i = levels_n.size() - 1;
        levels[i]->GetDataFromSameLevel(levels_n[i]);
    }

    for (int i = levels.size() - 1; i > 0; i--) {
        levels[i - 1]->GetDataFromCoarseNewLevel(levels[i]);
    }
}

void Mesh::updateHierarchy(bool init) {
    //Finest level first. Estimate error -> flag points. if flagged, finest+1 will be needed ¤
    //Estimate error on finest-1 -> flag points. if flagged, finest will need to be regridded, if any finest+1, they will need to fit inside new finest
    //etc until base.

    std::vector<coords> flaggedCells, flaggedl2, flaggedl1, flaggedtmp;
    rect extrema;
    level identified;

    if (hierarchy.empty()) {
        std::cerr << "Exception Occured: hierarchy is empty, cannot update it." << std::endl;
        exit(EXIT_FAILURE);
    }

    int lloop = (int)hierarchy.size() - 1;

    for (int l = lloop; l > -1; l--) {    //From finest to coarsest.
        //No need to update level 0: it's always the same, but we must check it for level 1 flaggedCells..
        if (l < (int)hierarchy.size() - 2) {
            //We need to make sure that all l+2 rectangles reside inside the current flag set
            //Flag all l+2 data points on the current grid.
            flaggedl1.clear();
            flaggedl2.clear();

            if (hierarchy.size() > 2) {
                for (rect &r: hierarchy.at(l + 2)) {
                    flaggedtmp.clear();
                    mergeDownFlaggedData(l, r, flaggedtmp);

                    if (!flaggedtmp.empty()) {
                        flaggedl2.insert(flaggedl2.end(), flaggedtmp.begin(), flaggedtmp.end());
                    }
                }
            }

            getError(l, false, flaggedl1);
            sort(flaggedl1.begin(), flaggedl1.end());
            flaggedCells.reserve(flaggedl1.size() + flaggedl2.size());
            set_union(flaggedl1.begin(), flaggedl1.end(), flaggedl2.begin(), flaggedl2.end(), back_inserter(flaggedCells));

            if (NOISY) {
                std::cout << "Flagged Cells (before union): " << flaggedl1.size() << std::endl;
                std::cout << "Flagged Cells (after union): " << flaggedCells.size() << std::endl;
            }
        } else {
            getError(l, false, flaggedCells);
        }

        if (flaggedCells.size() > 0) {
            getExtrema(extrema, flaggedCells);
            identified = splitRectangle(extrema, flaggedCells, settings.minEfficiency);

            if (l == lloop) {
                if (l < settings.maxDepth) {
                    //Need a new l+1 level.
                    interpRectanglesUp(identified, l);
                    loadBalanceRectangles(identified);
                    hierarchy.push_back(identified);

                    if (LOUD) std::cout << "Constructed new level (" << settings.maxDepth - (l + 1) << ") in hierarchy. Rectangles: " << identified.size() << ", Flagged Cells: " << flaggedCells.size() << std::endl;
                }
            } else {
                interpRectanglesUp(identified, l);
                loadBalanceRectangles(identified);
                hierarchy.at(l + 1) = identified;

                if (LOUD) std::cout << "Regridded level " << settings.maxDepth - (l + 1) << ". Rectangles: " << identified.size() << ", Flagged Cells: " << flaggedCells.size() << std::endl;
            }
        } else {
            if (l < lloop) {
                hierarchy.erase(hierarchy.begin() + l + 1);

                if (LOUD) std::cout << "Deleted level " << settings.maxDepth - (l + 1) << "." << std::endl;
            }
        }

        flaggedCells.clear();
    }

    promoteHierarchyToMesh(init);
}

void Mesh::loadBalanceRectangles(level &identified) {
    level rlevel;
    coords left, right;
    double length = settings.loadBalanceLength;
    double minSize = std::pow(length, 2.0);

    for (rect &r: identified) {
        double size = (r.second.first - r.first.first) * (r.second.second - r.first.second);
        int n_x = std::ceil((r.second.first - r.first.first) / length);
        int n_p = std::ceil((r.second.second - r.first.second) / length);

        if (size > minSize) {
            for (int i = 0; i < n_x; i++) {
                for (int j = 0; j < n_p; j++) {
                    left = std::make_pair(r.first.first + i * (int)length, r.first.second + j * (int)length);
                    right = std::make_pair(std::min(r.first.first + (i + 1) * (int)length, r.second.first.m), std::min(r.first.second + (j + 1) * (int)length, r.second.second.m));
                    rlevel.push_back(std::make_pair(left, right));
                }
            }
        } else {
            rlevel.push_back(r);
        }
    }

    identified = rlevel;
}

void Mesh::getError(const int &lvl, bool init, std::vector<coords> &flaggedCells) {
    //lvl is current level. Will generate flaggedCells for lvl+1
    if (init) {
        //Generate initial data on all levels
        double offset = 0.0;
        std::vector<double> errorWeights;

        errorWeights.assign(2, 0.5 / settings.GetfMax(particleType));

        for (int i = 0; i <= lvl; i++) {
            offset += std::pow(settings.refinementRatio, i);
        }

        offset *= 3;

        int maxx = settings.GetXSize(settings.maxDepth - lvl) - offset;
        int maxp = settings.GetPSize(settings.maxDepth - lvl,  particleType) - offset;
        double dp = settings.GetDp(settings.maxDepth - lvl, particleType), dx = settings.GetDx(settings.maxDepth - lvl);
        double temp, xp, pp, error;

        for (int i = offset; i < maxx; i++) { //Don't allow new levels to have the same bounds (ensure nesting)
            for (int j = offset; j < maxp; j++) {
                xp = (0.5 + i) * dx;
                pp = settings.pmin[particleType] + dp * (0.5 + j);
                temp = settings.InitialDistribution(xp, pp, particleType);
                error = errorWeights[0] * std::fabs(settings.InitialDistribution(xp + dx, pp, particleType) - settings.InitialDistribution(xp - dx, pp, particleType)) + errorWeights[1] * std::fabs(settings.InitialDistribution(xp, pp + dp, particleType) - settings.InitialDistribution(xp, pp - dp, particleType));

                if (error > settings.refinementCriteria) {
                    flaggedCells.push_back(std::make_pair(i, j));
                }
            }
        }
    } else {
        //We're advancing time, use current data on rectangles
        for (auto &r: levels.at(settings.maxDepth - lvl)->rectangles) {
            r->getError(flaggedCells, particleType);
        }
    }
}

void Mesh::interpRectanglesUp(level &identified, const int &lvl) {
    //Rectangle coords are on lvl but are actually needed on lvl+1. We need to re-index
    int nmax = settings.GetXSize(settings.maxDepth - lvl);
    int pmax = settings.GetPSize(settings.maxDepth - lvl, particleType);
    int nmax1 = settings.GetXSize(settings.maxDepth - (lvl + 1));
    int pmax1 = settings.GetPSize(settings.maxDepth - (lvl + 1), particleType);

    for (rect &r: identified) {
        r.second.first.m += 1;
        r.second.second.m += 1;
        r.first.first.m *= (nmax1 / double(nmax));
        r.first.second.m *= (pmax1 / double(pmax));
        r.second.first.m *= (nmax1 / double(nmax));
        r.second.second.m *= (pmax1 / double(pmax));
    }
}

void Mesh::mergeDownFlaggedData(const int &lvl, const rect &r, std::vector<coords> &foundCells) {
    //We want to flag all points on l which reside inside the provided l+2 rectangle
    std::vector<coords>::iterator it;
    int ic, jc;

    int nmax = settings.GetXSize(settings.maxDepth - lvl);
    int pmax = settings.GetPSize(settings.maxDepth - lvl, particleType);
    int nmaxl2 = settings.GetXSize(settings.maxDepth - (lvl + 2));
    int pmaxl2 = settings.GetPSize(settings.maxDepth - (lvl + 2), particleType);

    for (int i = r.first.first - 2 * (int)settings.refinementRatio; i <= r.second.first + 2 * (int)settings.refinementRatio + 1; i++) {
        for (int j = r.first.second - 2 * (int)settings.refinementRatio; j <= r.second.second + 2 * (int)settings.refinementRatio; j++) {
            ic = i * (nmax / double(nmaxl2));
            jc = j * (pmax / double(pmaxl2));

            if ((ic > -1) && (ic < nmax) && (jc > -1) && (jc < pmax)) {            //Ensure proper nesting.
                foundCells.push_back(std::make_pair(ic, jc));
            }
        }
    }

    sort(foundCells.begin(), foundCells.end());
    it = unique(foundCells.begin(), foundCells.end());
    foundCells.resize(distance(foundCells.begin(), it));
}

bool Mesh::choose_second(const coords &lhs, const coords &rhs) {
    return lhs.second < rhs.second;
}

bool Mesh::choose_first(const coords &lhs, const coords &rhs) {
    return lhs.first < rhs.first;
}

int Mesh::sgn(int &x) {
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

void Mesh::printRectangle(rect &r) {
    //Just prints a rectangle to screen. Usefull only for debugging. No newline so you can put multiple rectangles on one line.
    //Using this rather than the std::ostream override so we can print directly to file with that.
    std::cout << "(" << r.first.first.m << ", " << r.first.second.m << ")->(" << r.second.first.m << ", " << r.second.second.m << ")";
}

void Mesh::printRectangle(rect &r, int lvl) {
    //Just prints a rectangle to screen. Usefull only for debugging. No newline so you can put multiple rectangles on one line.
    //Using this rather than the std::ostream override so we can print directly to file with that.
    std::cout << "(" << ceil(r.first.first / std::pow(2, lvl)) << ", " << ceil(r.first.second / std::pow(2, lvl)) << ")->(" << ceil(r.second.first / std::pow(2, lvl)) << ", " << ceil(r.second.second / std::pow(2, lvl)) << ")";
}

int Mesh::countCells(const rect &span) {
    // (xmax-xmin)*(pmax-pmin)
    return (span.second.first - span.first.first + 1) * (span.second.second - span.first.second + 1);
}

void Mesh::isFlaggedInside(std::vector<coords> &flagged, std::vector<bool> &isInside, rect &r) {
    //returns true if flagged value is inside a rectangle.
    for (unsigned int i = 0; i < flagged.size(); i++) {
        isInside.at(i) = ((flagged.at(i).first >= r.first.first) && (flagged.at(i).first <= r.second.first) && (flagged.at(i).second >= r.first.second) && (flagged.at(i).second <= r.second.second)) == 1;
    }
}

void Mesh::getSigs(rect &rectangle, std::vector<coords> &flagged, std::vector<int> &sigX, std::vector<int> &sigP) {
    //Needed to separate this because it's wanted by computeSignatures and shrinkBorders

    for (coords &xp: flagged) {
        sigX.at(xp.first - rectangle.first.first)++;
        sigP.at(xp.second - rectangle.first.second)++;
    }
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int> > Mesh::computeSignatures(rect &rectangle, std::vector<coords> &flagged) {
    double lengthx, lengthp;

    lengthx = rectangle.second.first - rectangle.first.first;
    lengthp = rectangle.second.second - rectangle.first.second;

    std::vector<int> sigX(lengthx + 1), sigP(lengthp + 1);
    std::vector<int> delSigX, delSigP;     //Don't initialise as length can be less than 2.

    //Signatures are of size of the limits we have in the current rectangle, but the rectangle probably doesn't have a lower
    //bound of (0,0). So the for example, a rectangle of (20,5), (50, 17) has sigX(31), sigY(13), but are indexed from 0 still, so
    //in the loop below the ll coords are taken out. we will therefore need to add the ll coords to correctly index the signatures.

    getSigs(rectangle, flagged, sigX, sigP);

    if (NOISY) {
        std::cout << "sigX: ";
        for (int i: sigX) {
            std::cout << i << " ";
        }
        std::cout << std::endl << "sigP: ";
        for (int i: sigP) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }

    //Lapliacians are stored in a similar manner and need to be indexed by ll+1 for obvious reasons.
    // We only need a 1D, finite element laplacian so:

    if (lengthx >= 2) {
        delSigX.resize(lengthx - 1);
        for (unsigned int h = 1; h < sigX.size() - 1; h++) {
            delSigX.at(h - 1) = sigX.at(h - 1) + sigX.at(h + 1) - 2 * sigX.at(h);
        }

        if (NOISY) {
            std::cout << "delSigX: ";
            for (int i: delSigX) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }
    }

    if (lengthp >= 2) {
        delSigP.resize(lengthp - 1);

        for (unsigned int h = 1; h < sigP.size() - 1; h++) {
            delSigP.at(h - 1) = sigP.at(h - 1) + sigP.at(h + 1) - 2 * sigP.at(h);
        }
        if (NOISY) {
            std::cout << "delSigY: ";
            for (int i: delSigP) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }
    }

    return std::make_tuple(sigX, delSigX, sigP, delSigP);
}

coords Mesh::identifyInflection(rect &rectangle, std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int> > &signatures) {
    int inflectX = 0, inflectP = 0, inflect = 0, indexX = 0, indexP = 0, sizeSigX, sizeSigP, max, countmax = 0;
    std::vector<int> delSigX, delSigP, inflectIdx, inflectSize;
    std::vector<int>::iterator maxit;
    std::pair<bool, bool> bisect = std::make_pair(false, false);

    tie(std::ignore, delSigX, std::ignore, delSigP) = signatures;     //Using std::ignore we can pass big things by reference and pull out only the pieces we need. Neat.

    sizeSigX = delSigX.size();
    sizeSigP = delSigP.size();

    //Find largest inflection point.
    if (sizeSigX > 2) {
        for (unsigned int i = 0; i < delSigX.size() - 1; i++) {
            if ((delSigX.at(i) != 0) && (delSigX.at(i + 1) != 0) && (sgn(delSigX.at(i)) != sgn(delSigX.at(i + 1)))) {
                inflectIdx.push_back(i);
                inflectSize.push_back(abs(delSigX.at(i + 1) - delSigX.at(i)));
            }
        }

        inflect = inflectIdx.size();

        if (!inflectSize.empty()) {
            maxit = max_element(inflectSize.begin(), inflectSize.end());
            max = *maxit;
            countmax = count(inflectSize.begin(), inflectSize.end(), max);
        }

        if (inflect == 1) {
            //Only one inflection, cut here
            indexX = inflectIdx.at(0);
            inflectX = inflectSize.at(0);

            if (LOUD) std::cout << "X: One inflection " << indexX << ":" << inflectX << std::endl;
        } else if (inflect > 1) {
            //More than one inflection.
            if (countmax > 1) {
                //we have two maximum count values, take the one closest to the center.

                int center = sizeSigX / 2;         //we're fine to round here, we need an index.
                int range = std::numeric_limits<int>::max(); //NOTE: We will have to go to long ints (in general I mean) if we need grids above this, but I highly doubt that will happen.
                int chosen = 0;

                for (unsigned int i = 0; i < inflectIdx.size(); i++) {
                    if (range > abs(inflectIdx.at(i) - center)) {
                        range = abs(inflectIdx.at(i) - center);
                        chosen = i;
                    }
                }

                indexX = inflectIdx.at(chosen);
                inflectX = inflectSize.at(chosen);

                if (LOUD) std::cout << "X: Take Max in center " << indexX << ":" << inflectX << std::endl;
            } else {
                //Take the max value.
                indexX = inflectIdx.at(maxit - inflectSize.begin());
                inflectX = inflectSize.at(maxit - inflectSize.begin());

                if (LOUD) std::cout << "X: Take largest " << indexX << ":" << inflectX << std::endl;
            }
        } else {
            //We have either something highly symmetric or something pretty dicy.
            //Bisect in the long direction as a last resort.
            if (LOUD) std::cout << "X: Bisect on long edge" << std::endl;

            bisect.first = true;
        }
    }

    inflectIdx.clear();
    inflectSize.clear();

    if (sizeSigP > 2) {
        for (unsigned int i = 0; i < delSigP.size() - 1; i++) {
            if ((delSigP.at(i) != 0) && (delSigP.at(i + 1) != 0) && (sgn(delSigP.at(i)) != sgn(delSigP.at(i + 1)))) {
                inflectIdx.push_back(i);
                inflectSize.push_back(abs(delSigP.at(i + 1) - delSigP.at(i)));
            }
        }

        inflect = inflectIdx.size();

        if (!inflectSize.empty()) {
            maxit = max_element(inflectSize.begin(), inflectSize.end());
            max = *maxit;
            countmax = count(inflectSize.begin(), inflectSize.end(), max);
        }

        if (inflect == 1) {
            //Only one inflection, cut here
            indexP = inflectIdx.at(0);
            inflectP = inflectSize.at(0);

            if (LOUD) std::cout << "P: One inflection " << indexP << ":" << inflectP << std::endl;
        } else if (inflect > 1) {
            //More than one inflection.
            if (countmax > 1) {
                //we have two maximum count values, take the one closest to the center.

                int center = sizeSigP / 2;         //we're fine to round here, we need an index.
                int range = 32767;                 //MAX_INT
                int chosen = 0;

                for (unsigned int i = 0; i < inflectIdx.size(); i++) {
                    if (range > abs(inflectIdx.at(i) - center)) {
                        range = abs(inflectIdx.at(i) - center);
                        chosen = i;
                    }
                }

                indexP = inflectIdx.at(chosen);
                inflectP = inflectSize.at(chosen);

                if (LOUD) std::cout << "P: Take Max in center " << indexP << ":" << inflectP << std::endl;
            } else {
                //Take the max value.
                indexP = inflectIdx.at(maxit - inflectSize.begin());
                inflectP = inflectSize.at(maxit - inflectSize.begin());

                if (LOUD) std::cout << "P: Take largest " << indexP << ":" << inflectP << std::endl;
            }
        } else {
            //We have either something highly symmetric or something pretty dicy.
            //Bisect in the long direction as a last resort.
            if (LOUD) std::cout << "P: Bisect on long edge" << std::endl;

            bisect.second = true;
        }
    }

    if (LOUD) std::cout << "Inflect vals: " << inflectX << " " << inflectP << " (" << bisect.first << "," << bisect.second << ")" << std::endl;

    if (bisect.first && bisect.second) {
        //We have to bisect.
        return std::make_pair(2, 0);
    } else if (((sizeSigX <= 2) && (sizeSigP <= 2)) || ((inflectX == 0) && (inflectP == 0))) {
        //Rectangle is too small, don't cut
        return std::make_pair(3, 0);
    } else if (inflectX > inflectP) {
        indexX += rectangle.first.first + 1;       //Offset index with ll+1 of rectangle
        return std::make_pair(0, indexX);         //0 denotes a cut in X
    } else {
        indexP += rectangle.first.second + 1;
        return std::make_pair(1, indexP);         //1 denotes a cut in P
    }
}

std::tuple<bool, int, int> Mesh::hasHole(std::vector<int> &sig) {
    //TODO: It may be possible to implement this a bit more cleanly now I know we can send variables to lambda functions.
    //Search a signature for a hole
    int start = 0, finish = 0;
    bool found;
    std::vector<int>::iterator it;

    it = search_n(sig.begin(), sig.end(), 1, 0);     //Only one now.
    found = it != sig.end();

    if (found) {
        //We've found at least one hole.
        std::vector<int>::iterator it2;

        //we may be cutting on the boundary, if so don't -1.
        it == sig.begin() ? start = it - sig.begin() : start = it - sig.begin() - 1;

        //The job here is not to be optimal and identify ALL holes, just the one. We can cut this again later by calling a cut on the suboptimal second piece.
        //However, we can check to see if we have a longer hole length.
        it2 = find_if_not(it, sig.end(), [](int i) {
            return i == 0;
        });
        //we may be cutting on the boundary, if so take the end of it..
        it2 == sig.end() ? finish = it2 - sig.begin() - 1 : finish = it2 - sig.begin();
    }

    return std::make_tuple(found, start, finish);
}

void Mesh::getExtrema(rect &extrema, const std::vector<coords> &flaggedCells) {
    coords ll, ur;     //these aren't really required, but the code gets a little terse if we don't use them.

    auto limx = minmax_element(flaggedCells.begin(), flaggedCells.end(), choose_first);
    auto limp = minmax_element(flaggedCells.begin(), flaggedCells.end(), choose_second);

    ll = std::make_pair(std::get<0>(*limx.first), std::get<1>(*limp.first));     //Lower left coordinates of the inital bounds
    ur = std::make_pair(std::get<0>(*limx.second), std::get<1>(*limp.second));     //Upper right coordinates of the inital bounds

    extrema = std::make_pair(ll, ur);

    if (LOUD) {
        std::cout << "Extrema: ";
        printRectangle(extrema);
        std::cout << std::endl;
    }
}

void Mesh::saveRectangleData(std::ofstream &os, rect &rectangle, std::vector<rectData> &output) {
    //rectangle bounds at top, followed by coords and data
    os << "r " << rectangle << std::endl;
    std::ostream_iterator<rectData> out(os, "\n");
    copy(output.begin(), output.end(), out);
}

level Mesh::splitRectangle(rect &rectangle, std::vector<coords> &flagged, const double &minEfficiency) {
    double efficiencyRatio;
    std::vector<rect> resultRect;

    efficiencyRatio = flagged.size() / double(countCells(rectangle));

    if (LOUD) {
        std::cout << "Eff for ";
        printRectangle(rectangle);
        std::cout << ". Flagged: " << flagged.size() << ", Cells: " << countCells(rectangle) << " (" << efficiencyRatio << ")" << std::endl;
    }

    if (efficiencyRatio > minEfficiency) {
        //We're done. Head back up the chain.
        if (LOUD) std::cout << "RESULT: (" << rectangle.first.first.m << ", " << rectangle.first.second.m << "), (" << rectangle.second.first.m << ", " << rectangle.second.second.m << "): " << flagged.size() << std::endl;

        resultRect.push_back(rectangle);
        return resultRect;
    } else {
        std::tuple<bool, int, int> holePosX, holePosP;
        rect pieceOne, pieceTwo;
        std::vector<coords> flaggedOne, flaggedTwo;
        std::vector<bool> isInside;
        std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int> > signatures;
        std::vector<rect> resultRectTwo, resultMerged;
        //Get signatures
        signatures = computeSignatures(rectangle, flagged);

        //Find Islands
        holePosX = hasHole(std::get<0>(signatures));         //sigX
        holePosP = hasHole(std::get<2>(signatures));         //sigP
        //Offset index with ll of rectangle
        std::get<1>(holePosX) += rectangle.first.first;
        std::get<2>(holePosX) += rectangle.first.first;
        std::get<1>(holePosP) += rectangle.first.second;
        std::get<2>(holePosP) += rectangle.first.second;

        if ((std::get<0>(holePosX)) && (std::get<0>(holePosP))) {
            //Holes in both directions.
            if (LOUD) {
                std::cout << "Holes in both directions." << std::endl;
                std::cout << "X. start: " << std::get<1>(holePosX) << ", finish: " << std::get<2>(holePosX) << std::endl;
                std::cout << "P. start: " << std::get<1>(holePosP) << ", finish: " << std::get<2>(holePosP) << std::endl;
            }

            //Unfortunately we can't cut in both dirs unless we branch out by 4 here.
            //It's too much work, so just cut the biggest distance first.
            if (abs(std::get<2>(holePosX) - std::get<1>(holePosX)) > abs(std::get<2>(holePosP) - std::get<1>(holePosP))) {
                //Cut X
                if (LOUD) std::cout << "Cutting X" << std::endl;

                pieceOne = std::make_pair(rectangle.first, coords(std::get<1>(holePosX), rectangle.second.second)); //left
                pieceTwo = std::make_pair(coords(std::get<2>(holePosX), rectangle.first.second), rectangle.second); //right
            } else {
                if (LOUD) std::cout << "Cutting P" << std::endl;

                pieceOne = std::make_pair(rectangle.first, coords(rectangle.second.first, std::get<1>(holePosP)));  //bottom
                pieceTwo = std::make_pair(coords(rectangle.first.first, std::get<2>(holePosP)), rectangle.second);  //top
            }
        } else if (std::get<0>(holePosX)) {
            //Hole in X
            if (LOUD) std::cout << "Cut X. start: " << std::get<1>(holePosX) << ", finish: " << std::get<2>(holePosX) << std::endl;

            pieceOne = std::make_pair(rectangle.first, coords(std::get<1>(holePosX), rectangle.second.second));     //left
            pieceTwo = std::make_pair(coords(std::get<2>(holePosX), rectangle.first.second), rectangle.second);     //right
        } else if (std::get<0>(holePosP)) {
            //Hole in P
            if (LOUD) std::cout << "Cut P. start: " << std::get<1>(holePosP) << ", finish: " << std::get<2>(holePosP) << std::endl;

            pieceOne = std::make_pair(rectangle.first, coords(rectangle.second.first, std::get<1>(holePosP)));      //bottom
            pieceTwo = std::make_pair(coords(rectangle.first.first, std::get<2>(holePosP)), rectangle.second);      //top
        } else {
            //No Holes, choose inflection point.
            coords splitLocation;                                                                                   //.first = 0, then cut in X, = 1, cut in P at index .second, = 2, bisect long, = 3, don't cut at all.

            if (LOUD) std::cout << "Inflection cut" << std::endl;

            splitLocation = identifyInflection(rectangle, signatures);

            if (splitLocation.first == 0) {
                if (LOUD) std::cout << "Cut X: " << splitLocation.second.m << std::endl;

                pieceOne = std::make_pair(rectangle.first, coords(splitLocation.second, rectangle.second.second));                //left
                pieceTwo = std::make_pair(coords(splitLocation.second + 1, rectangle.first.second), rectangle.second);            //right
            } else if (splitLocation.first == 1) {
                if (LOUD) std::cout << "Cut P: " << splitLocation.second.m << std::endl;

                pieceOne = std::make_pair(rectangle.first, coords(rectangle.second.first, splitLocation.second));                 //bottom
                pieceTwo = std::make_pair(coords(rectangle.first.first, splitLocation.second + 1), rectangle.second);             //top
            } else if (splitLocation.first == 2) {
                // Have to bisect.
                double cutPosition, lengthx, lengthp;

                lengthx = rectangle.second.first - rectangle.first.first;
                lengthp = rectangle.second.second - rectangle.first.second;

                if (lengthx > lengthp) {
                    cutPosition = floor((lengthx / 2) + rectangle.first.first);
                    pieceOne = std::make_pair(rectangle.first, coords(cutPosition, rectangle.second.second));                    //left
                    pieceTwo = std::make_pair(coords(cutPosition + 1, rectangle.first.second), rectangle.second);                //right

                    if (LOUD) std::cout << "Bisect in x. Cut at " << cutPosition << ", L:" << lengthx << ", O:" << rectangle.first.first.m << ", T:" << rectangle.second.first.m << std::endl;
                } else {
                    cutPosition = floor((lengthp / 2) + rectangle.first.second);
                    pieceOne = std::make_pair(rectangle.first, coords(rectangle.second.first, cutPosition));                     //bottom
                    pieceTwo = std::make_pair(coords(rectangle.first.first, cutPosition + 1), rectangle.second);                 //top

                    if (LOUD) std::cout << "Bisect in p. Cut at " << cutPosition << ", L:" << lengthp << ", O:" << rectangle.first.second.m << ", T:" << rectangle.second.second.m << std::endl;
                }
            } else {
                //Can't cut, rectangle is too small. Head back up.
                if (LOUD) std::cout << "RESULT: (" << rectangle.first.first.m << ", " << rectangle.first.second.m << "), (" << rectangle.second.first.m << ", " << rectangle.second.second.m << "): " << flagged.size() << std::endl;

                resultRect.push_back(rectangle);
                return resultRect;
            }
        }

        //Pull out flagged coords in new rectangles.
        isInside.resize(flagged.size());
        isFlaggedInside(flagged, isInside, pieceOne);

        for (unsigned int i = 0; i < flagged.size(); i++) {
            if (isInside.at(i)) {
                flaggedOne.push_back(flagged.at(i));
            } else {
                flaggedTwo.push_back(flagged.at(i));
            }
        }

        if (LOUD || NOISY) std::cout << "Found in each piece: " << flaggedOne.size() << ", " << flaggedTwo.size() << std::endl;

        if (NOISY) {
            for (const coords &f: flaggedOne) {
                std::cout << f << ", ";
            }

            std::cout << "\b\b  " << std::endl;

            for (const coords &f: flaggedTwo) {
                std::cout << f << ", ";
            }

            std::cout << "\b\b  " << std::endl;
        }

        //So long as we still have flagged cells, we must go deeper.
        if (flaggedOne.size() > 0) {
            resultRect = splitRectangle(pieceOne, flaggedOne, minEfficiency);
        }

        if (flaggedTwo.size() > 0) {
            resultRectTwo = splitRectangle(pieceTwo, flaggedTwo, minEfficiency);
        }

        // preallocate memory
        resultMerged.reserve(resultRect.size() + resultRectTwo.size());
        //Merge the two resultant sets of data from the recursion.
        resultMerged.insert(resultMerged.end(), resultRect.begin(), resultRect.end());
        resultMerged.insert(resultMerged.end(), resultRectTwo.begin(), resultRectTwo.end());
        return resultMerged;
    }
}

void Mesh::promoteHierarchyToMesh(bool init) {
    //Once a hierarchy is created we need to generate appropreate Level and Rectangle instances for it

    int xmax, pmax, N = hierarchy.size() - 1;
    int empty = settings.maxDepth - N;
    bool up, down, left, right;
    std::vector<std::unique_ptr<Level> > levels_n;

    if (!init) { //Copy old level sructure.
        for (auto &lvl: levels) {
            levels_n.push_back(std::move(lvl));
        }
    }

    levels.clear();

    for (int i = 0; i < empty; i++) { //Generates empty levels
        levels.push_back(std::make_unique<Level>(particleType, i, settings));
    }

    for (int lvl = N; lvl >= 0; lvl--) { //Levels AND Rectangles are indexed finest = 0, coarsest = N.
        xmax = settings.GetXSize(N - lvl);
        pmax = settings.GetPSize(N - lvl, particleType);
        levels.push_back(std::make_unique<Level>(particleType, N - lvl + empty, settings));

        //Rectangles.
        for (rect &r: hierarchy.at(lvl)) {
            //Get boundaries
            up = r.second.second == pmax;
            down = r.first.second == 0;
            left = r.first.first == 0;
            right = r.second.first == xmax;

            levels.at(N - lvl + empty)->rectangles.push_back(std::make_shared<Rectangle>(r.second.first - r.first.first, r.second.second - r.first.second, r.first.first, r.first.second, N - lvl + empty, settings, bc, up, down, left, right, particleType));
        }
    }

    if (init) {
       for (auto &lvl: levels) {
           #pragma omp parallel for schedule(dynamic,1)
           for (unsigned int i = 0; i < lvl->rectangles.size(); i++) {
               lvl->rectangles[i]->InitializeDistribution();
           }
       }
    }

    for (unsigned int l = 0; l < levels.size(); l++) {
        auto &lvl = levels.at(l);
        for (auto &rl: lvl->rectangles) {
            if (l > 0) {
                for (auto &rf: levels.at(l - 1)->rectangles) {
                    rl->CalculateConnectivityFromFiner(rl, rf);
                }
            }
        }
    }

    for (unsigned int l = 0; l < levels.size(); l++) {
        auto &lvl = levels.at(l);
        for (auto &rl: lvl->rectangles) {
            for (auto &rr: lvl->rectangles) {
                if (rl != rr) {
                    rl->CalculateConnectivitySame(rr);
                }
            }
        }
    }

    for (unsigned int l = 0; l < levels.size(); l++) {
        auto &lvl = levels.at(l);
        #pragma omp parallel for schedule(dynamic,1)
        for (int i=0;i<lvl->rectangles.size();i++) {
            if (init) lvl->rectangles[i]->InitializeDistribution();
        }
    }


    SetFieldSolver(EMSolver);

    if (!init) InterMeshDataTransfer(levels_n);

    PushData();

    for (unsigned int l = 0; l < levels.size(); l++) {
        auto &lvl = levels.at(l);
        for (auto &rl: lvl->rectangles) {
            rl->FCTTimeStep(0,-1,3);
        }
    }


}

void Mesh::outputRectangleData(double tidx) {
    rect trect;
    std::ofstream os;
    std::vector<rectData> output;

    for (int l = 0; l <= settings.maxDepth; l++) {
        if (!levels.at(l)->rectangles.empty()) {
            std::stringstream fileName;
            fileName << "output/rectangleData/rectangleData_p" << particleType << "_l" << l << "_t" << std::scientific << tidx << ".txt";
            os.open(fileName.str());

            for (unsigned int r = 0; r < levels.at(l)->rectangles.size(); r++) {
                auto &rectangle = levels.at(l)->rectangles.at(r);

                for (int i = 0; i < rectangle->n_x; i++) {
                    for (int j = 0; j < rectangle->n_p; j++) {
                        output.push_back(std::make_tuple(rectangle->x_pos + i, rectangle->p_pos + j, rectangle->f[rectangle->Index3(i, j, 0)]));
                    }
                }

                trect = std::make_pair(std::make_pair(rectangle->x_pos, rectangle->p_pos), std::make_pair(rectangle->n_x, rectangle->n_p));
                saveRectangleData(os, trect, output);
                output.clear();
            }

            os.close();
        }
    }
}

void Mesh::PushBoundaryC() {
    for (unsigned int i = 0; i < levels.size(); i++) {
        levels.at(i)->PushData(4);
    }

    for (unsigned int i = 0; i < levels.size(); i++) {
        levels.at(i)->PushData(5);
    }

    for (unsigned int i = 0; i < levels.size(); i++) {
        levels.at(i)->PushData(6);
    }

}
