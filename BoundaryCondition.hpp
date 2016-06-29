#ifndef __BoundaryCondition_hpp__
#define __BoundaryCondition_hpp__

#include "Rectangle.hpp"
#include "Settings.hpp"

class BoundaryCondition : public Rectangle {
public:
~BoundaryCondition();
virtual double GetValueFromSameLevel(int i, int j,int val=1.0);
};

#endif /* __BoundaryCondition_hpp__ */
