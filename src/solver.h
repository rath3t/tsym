#ifndef TSYM_SOLVER_H
#define TSYM_SOLVER_H

#include <cassert>
#include <chrono>
#include "var.h"
#include "vector.h"
#include "matrix.h"
#include "logging.h"
#include "printer.h"

namespace tsym {
class Solver
{
    public:

    Solver(){};
    virtual ~Solver() {};
    virtual void solve(Matrix& a, Array1D<double>& u, Array1D<double>& f) = 0;
};

}
#include "solver_LU.h"

#endif 



