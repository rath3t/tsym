#include "Solver.h"

tsym::Vector tsym::Solver_LU::solve(const tsym::Matrix& A, const tsym::Vector& rhs) const
{
    if (!isSquare())
        TSYM_ERROR("Matrix (%zu, %zu) isn't square!", nRow , nCol);
    else if (rhs.size() != nRow)
        TSYM_ERROR("Matrix dimension %zu doesn't match vector size %zu", nRow, rhs.dim);
    else if (nRow == 0)
        TSYM_ERROR("Matrix and vector with zero dimension can't be solved!");
    else
        return solveChecked(A, rhs);

    TSYM_ERROR("Return vector with zero dimension.");

    return Vector();
}

tsym::Vector tsym::Solver_LU::solveChecked(const tsym::Matrix& A, const tsym::Vector& rhs)
{
    auto ts = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds ms;
    unsigned nPivotSwaps;
    decltype(ts) te;
    Matrix PLU(A);
    Vector b(rhs);
    Vector x(nRow);

    nPivotSwaps = PLU.compPartialPivots(&b);
    factorizeLU(&PLU);

    if (detFromLU(&PLU, nPivotSwaps).isZero()) {
        TSYM_WARNING("Can't solve system of equations with singular coefficient matrix!");
        x = Vector();
    } else {
        compXFromLU(&PLU, x, b);

        te = std::chrono::high_resolution_clock::now();
        ms = std::chrono::duration_cast<std::chrono::microseconds>(te - ts);

        TSYM_INFO("Solved %zu-dim. system of equations in %.2f ms.", nRow,
                static_cast<float>(ms.count())/1000.0);
    }

    return x;

}

void tsym::Matrix::factorizeLU(tsym::Matrix& A)
{
    for (size_t j = 0; j + 1 < nCol; ++j) {
        for (size_t i = j + 1; i < nRow; ++i) {
            A.data[i][j] /= A.data[j][j];
                for (size_t k = j + 1; k < nCol; ++k)
                    A.data[i][k] -= A.data[i][j]*A.data[j][k];
        }
    }
}

void tsym::Matrix::compXFromLU(tsym::Matrix& A, tsym::Vector& x, tsym::Vector& b) const
{
    for (size_t i = 0; i < nRow; ++i)
        for (size_t j = 0; j < i; ++j)
            b.data[i] -= A.data[i][j]*b.data[j];

    assert(areAllItemsZero(x));

    for (size_t i = nRow - 1; i + 1 > 0; --i) {
        for (size_t j = i + 1; j < nCol; ++j) {
            x.data[i] -= A.data[i][j]*x.data[j];
        }

        x.data[i] = ((b.data[i] + x.data[i])/A.data[i][i]).normal();
    }
}

tsym::Var tsym::Matrix::detFromLU(tsym::Matrix& A, unsigned nPivotSwaps) const
{
    Var det(nPivotSwaps % 2 == 0 ? 1 : -1);

    for (size_t i = 0; i < nRow; ++i)
        det *= A.data[i][i];

    det = det.normal();

    return det;
}


