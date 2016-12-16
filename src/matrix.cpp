
#include <cassert>
#include <ctime>
#include "matrix.h"
#include "numeric.h"
#include "logging.h"
#include "printer.h"

namespace tsym {
    namespace {
        const Var& constZero()
        {
            static const Var zero(Numeric::zero());

            return zero;
        }

        Var& zero()
        {
            static Var nonConstZero(Numeric::zero());

            return nonConstZero;
        }

#ifndef NDEBUG
        bool areAllItemsZero(const Vector& vec)
        {
            for (size_t i = 0; i < vec.size(); ++i)
                if (!vec(i).isZero())
                    return false;

            return true;
        }
#endif
    }
}

tsym::Matrix::Matrix() :
    data(NULL),
    nRow(0),
    nCol(0)
{}

tsym::Matrix::Matrix(size_t nRow, size_t nCol) :
    data(NULL),
    nRow(nRow),
    nCol(nCol)
{
    allocateMem();
}

tsym::Matrix::Matrix(const Matrix& other) :
    nRow(other.nRow),
    nCol(other.nCol)
{
    allocateMem();
    copyValuesFromMatrix(other);
}

tsym::Matrix& tsym::Matrix::operator = (const Matrix& rhs)
{
    if (this == &rhs)
        return *this;

    if (nRow != rhs.nRow || nCol != rhs.nCol) {
        deleteMem();
        nRow = rhs.nRow;
        nCol = rhs.nCol;
        allocateMem();
    }

    copyValuesFromMatrix(rhs);

    return *this;
}

tsym::Matrix::~Matrix()
{
    deleteMem();
}

void tsym::Matrix::allocateMem()
{
    data = new Var*[nRow];

    for (size_t i = 0; i < nRow; ++i)
        data[i] = new Var[nCol];
}

void tsym::Matrix::copyValuesFromMatrix(const Matrix& other)
{
    /* nCol and nRow of both matrices have been checked, they are equal. */
    for (size_t i = 0; i < nRow; ++i)
        for (size_t j = 0; j < nCol; ++j)
            data[i][j] = other.data[i][j];
}

void tsym::Matrix::deleteMem()
{
    if (nRow == 0 && nCol == 0)
        return;

    for (size_t i = 0; i < nRow; i++)
        delete[] data[i];

    delete[] data;
}

tsym::Var& tsym::Matrix::operator () (size_t i, size_t j)
{
    if (nRow == 0 && nCol == 0)
        logging::error() << "Matrix has zero column/row size! Return zero.";
    else if (i > nRow || j > nCol)
        logging::error() << "Matric indices (" << i << ", " << j << ") out of bounds! Return zero.";
    else
        return data[i][j];

    /* A non-const reference to a static object is returned, it could possibly be modified from
     * the outside. So initialize it to zero before returning is necessary. */
    zero() = constZero();

    return zero();
}

const tsym::Var& tsym::Matrix::operator () (size_t i, size_t j) const
{
    if (nRow == 0 && nCol == 0)
        logging::error() << "Matrix has zero column/row size! Return zero.";
    else if (i > nRow || j > nCol)
        logging::error() << "Matric indices (" << i << ", " << j << ") out of bounds! Return zero.";
    else
        return data[i][j];

    return constZero();
}

tsym::Matrix& tsym::Matrix::operator += (const Matrix& rhs)
{
    if (nRow == rhs.nRow && nCol == rhs.nCol) {
        for (size_t i = 0; i < nRow; ++i)
            for (size_t j = 0; j < nCol; ++j)
                data[i][j] += rhs.data[i][j];
    } else
        logging::error() << "Matrix dimensions " << nRow << " and " << rhs.nRow <<
            " don't match! Return unmodified left hand side.";

    return *this;
}

tsym::Matrix& tsym::Matrix::operator -= (const Matrix& rhs)
{
    if (nRow == rhs.nRow && nCol == rhs.nCol) {
        for (size_t i = 0; i < nRow; ++i)
            for (size_t j = 0; j < nCol; ++j)
                data[i][j] -= rhs.data[i][j];
    } else
        logging::error() << "Matrix dimensions " << nRow << " and " << rhs.nRow <<
            " don't match! Return unmodified left hand side.";

    return *this;
}

tsym::Matrix& tsym::Matrix::operator *= (const Matrix& rhs)
{
    if (nCol == rhs.nRow)
        multiplyChecked(rhs);
    else
        logging::error() << "Matrix dimensions " << nRow << " and " << rhs.nRow <<
            " don't match! Return matrix with zero entries.";

    return *this;
}

void tsym::Matrix::multiplyChecked(const Matrix& rhs)
{
    const Matrix copy(*this);

    deleteMem();
    nCol = rhs.nCol;
    allocateMem();

    for (size_t i = 0; i < copy.nRow; ++i)
        for (size_t j = 0; j < rhs.nCol; ++j)
            for (size_t k = 0; k < copy.nCol; ++k)
                data[i][j] += copy.data[i][k]*rhs.data[k][j];
}

tsym::Matrix& tsym::Matrix::operator *= (const Var& rhs)
{
    for (size_t i = 0; i < nRow; ++i)
        for (size_t j = 0; j < nCol; ++j)
            data[i][j] *= rhs;

    return *this;
}

tsym::Vector tsym::Matrix::operator * (const Vector& rhs) const
{
    Vector result(nRow);

    if (nCol == rhs.dim)
        for (size_t i = 0; i < nRow; ++i)
            for (size_t j = 0; j < nCol; ++j)
                result.data[i] += data[i][j]*rhs.data[j];
    else
        logging::error() << nCol << " matrix columns don't match vector size (" << rhs.dim <<
            ")! Return vector with zero entries.";

    return result;
}

const tsym::Matrix& tsym::Matrix::operator + () const
{
    return *this;
}

tsym::Matrix tsym::Matrix::operator - () const
{
    Matrix res(nRow, nCol);

    for (size_t i = 0; i < nRow; ++i)
        for (size_t j = 0; j < nCol; ++j)
            res.data[i][j] = -data[i][j];

    return res;
}

tsym::Matrix tsym::Matrix::transpose() const
{
    Matrix res(nCol, nRow);

    for (size_t i = 0; i < nCol; ++i)
        for (size_t j = 0; j < nRow; ++j)
            res.data[i][j] = data[j][i];

    return res;
}

tsym::Vector tsym::Matrix::solve(const Vector& rhs) const
{
    if (!isSquare())
        logging::error() << "Matrix (" << nRow << ", " << nCol << ") isn't square!";
    else if (rhs.size() != nRow)
        logging::error() << "Matrix dim. " << nRow << " doesn't match vector size " << rhs.dim;
    else if (nRow == 0)
        logging::error() << "Matrix and vector with zero dimension can't be solved!";
    else
        return solveChecked(rhs);

    logging::error() << "Return vector with zero dimension.";

    return Vector();
}

tsym::Vector tsym::Matrix::solveChecked(const Vector& rhs) const
{
    const time_t start = time(NULL);
    Matrix PLU(*this);
    Vector b(rhs);
    Vector x(nRow);

    PLU.compPartialPivots(&b);
    PLU.factorizeLU();

    if (PLU.detFromLU().isZero()) {
        logging::warning() << "Can't solve system of equations with singular coefficient matrix!";
        x = Vector();
    } else {
        PLU.compXFromLU(x, b);
        logging::info() << "Solved " << nRow << "-dim. system of equations in " <<
            difftime(time(NULL), start) << " s.";
    }

    return x;
}

void tsym::Matrix::compPartialPivots(Vector *b)
{
    for (size_t j = 0; j + 1 < nCol; ++j) {
        if (!data[j][j].isZero())
            continue;

        for (size_t i = j + 1; i < nRow; ++i)
            if (!data[i][j].isZero()) {
                swapRows(i, j);
                if (b != NULL)
                    std::swap(b->data[j], b->data[i]);
                break;
            }
    }
}

void tsym::Matrix::swapRows(size_t index1, size_t index2)
{
    for (size_t j = 0; j < nCol; ++j)
        std::swap(data[index1][j], data[index2][j]);
}

void tsym::Matrix::factorizeLU()
{
    for (size_t j = 0; j + 1 < nCol; ++j) {
        for (size_t i = j + 1; i < nRow; ++i) {
            data[i][j] = data[i][j]/data[j][j];
                for (size_t k = j + 1; k < nCol; ++k)
                    data[i][k] = data[i][k] - data[i][j]*data[j][k];
        }
    }
}

void tsym::Matrix::compXFromLU(Vector& x, Vector& b) const
{
    for (size_t i = 0; i < nRow; ++i)
        for (size_t j = 0; j < i; ++j)
            b.data[i] -= data[i][j]*b.data[j];

    assert(areAllItemsZero(x));

    for (size_t i = nRow - 1; i + 1 > 0; --i) {
        for (size_t j = i + 1; j < nCol; ++j) {
            x.data[i] -= data[i][j]*x.data[j];
        }

        x.data[i] = ((b.data[i] + x.data[i])/data[i][i]).normal();
    }
}

tsym::Matrix tsym::Matrix::inverse() const
{
    if (!isSquare())
        logging::error() << "Inversion for " << nRow << "x" << nCol << " maxtrix impossible!";
    else if (det() == 0)
        logging::error() << "Matrix is singular, no inversion possible!";
    else
        return checkedInverse();

    logging::error() << "Return zero dimension matrix.";

    return Matrix();
}

tsym::Matrix tsym::Matrix::checkedInverse() const
{
    Matrix inverse(nRow, nRow);
    Vector inverseCol;
    Vector b(nRow);

    for (size_t i = 0; i < nRow; ++i, b.data[i - 1] = 0) {
        b.data[i] = 1;

        inverseCol = solve(b);

        for (size_t j = 0; j < nRow; ++j)
            inverse.data[j][i] = inverseCol.data[j];
    }

    return inverse;
}

tsym::Var tsym::Matrix::det() const
{
    if (nRow == nCol && nRow != 0)
        return checkedDet();

    logging::error() << "Illegal determinant request for " << nRow << "x" << nCol << " matrix!";
    logging::error() << "Return zero determinant.";

    return constZero();
}

tsym::Var tsym::Matrix::checkedDet() const
{
    Matrix PLU(*this);

    PLU.compPartialPivots(NULL);
    PLU.factorizeLU();

    return PLU.detFromLU();
}

tsym::Var tsym::Matrix::detFromLU() const
{
    Var det(1);

    for (size_t i = 0; i < nRow; ++i)
        det *= data[i][i];

    det = det.normal();

    return det;
}

size_t tsym::Matrix::rowSize() const
{
    return nRow;
}

size_t tsym::Matrix::colSize() const
{
    return nCol;
}

bool tsym::Matrix::isSymmetric() const
{
    if (nRow != nCol)
        return false;

    for (size_t i = 1; i < nRow; ++i)
        for (size_t j = 0; j < i; ++j)
            if (data[i][j] != data[j][i])
                return false;

    return true;
}

bool tsym::Matrix::isSquare() const
{
    return nRow == nCol;
}

bool tsym::Matrix::equal(const Matrix& other) const
{
    if (nRow != other.nRow || nCol != other.nCol)
        return false;

    for (size_t i = 0; i < nRow; ++i)
        for (size_t j = 0; j < nCol; ++j)
            if (data[i][j] != other.data[i][j])
                return false;

    return true;
}

bool tsym::operator == (const Matrix& lhs, const Matrix& rhs)
{
    return lhs.equal(rhs);
}

bool tsym::operator != (const Matrix& lhs, const Matrix& rhs)
{
    return !lhs.equal(rhs);
}

tsym::Matrix tsym::operator + (Matrix lhs, const Matrix& rhs)
{
    lhs += rhs;

    return lhs;
}

tsym::Matrix tsym::operator - (Matrix lhs, const Matrix& rhs)
{
    lhs -= rhs;

    return lhs;
}

tsym::Matrix tsym::operator * (Matrix lhs, const Matrix& rhs)
{
    lhs *= rhs;

    return lhs;
}

tsym::Matrix tsym::operator * (Matrix lhs, const Var& rhs)
{
    lhs *= rhs;

    return lhs;
}

tsym::Matrix tsym::operator * (const Var& lhs, Matrix rhs)
{
    rhs *= lhs;

    return rhs;
}

std::ostream& tsym::operator << (std::ostream& stream, const Matrix& m)
{
    Printer printer(m);

    printer.print(stream);

    return stream;
}
