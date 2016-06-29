#include "veritas.hpp"
#include "BoundaryCondition.hpp"
BoundaryCondition::~BoundaryCondition() {
}

double BoundaryCondition::GetValueFromSameLevel(int i, int j,int val) {
    return 0.0;
}
