#ifndef __veritas_hpp__
#define __veritas_hpp__

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <tuple>
#include <memory>
#include <sstream>
#include <vector>
#include <omp.h>
#if defined(USINGMKL)
#include <mkl.h>
#endif

#define eps0_inv 1.1294e+11
#define eps0 8.854187817e-12
#define mu 1.2566370614e-6
#define mu_inv 795774.715482
#define cs 299792458.0
#define c_inv 3.33564095e-9
#define eps 0
#define DPI 6.28318530718
#define OUTW 80

class Rectangle;
class Level;
class Mesh;
class Parameters;
class EMFieldSolver;
class Settings;

extern bool LOUD;
extern bool NOISY;

//Using this overloaded int allows us to use the >> and << overloads outside of std for coords and rect types (which would otherwise be illegal).
struct oint {
    int m;
    oint() : m() {}
    oint(int p) : m(p) {}
    operator int() const { return m; }
};

typedef std::pair<oint, oint>            coords;
typedef std::pair<coords, coords>        rect;
typedef std::vector<rect>                level;
typedef std::tuple<oint, oint, double>   rectData;

std::ostream& operator<<(std::ostream& o, const coords& p);
std::ostream& operator<<(std::ostream& o, const rect& r);
std::ostream& operator<<(std::ostream& o, rectData& data);

#endif /* __veritas_hpp__ */
