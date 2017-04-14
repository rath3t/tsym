#ifndef TSYM_SOLVER_LU_H
#define TSYM_SOLVER_LU_H

namespace tsym {
class Solver_LU : public Solver
{
public:
  void solve(tsym::Matrix& a, tsym::Vector& u, tsym::Vector& f);
  
};
}

#endif 