#include <cmath>

#include <omp.h>

#include "rmatrixpreconditioner.h"

void RMatrixPreconditioner::_init(const RMatrixPreconditioner *pMatrixPreconditioner)
{
    if (pMatrixPreconditioner)
    {
        this->matrixPreconditionerType = pMatrixPreconditioner->matrixPreconditionerType;
        this->data = pMatrixPreconditioner->data;
    }
}

RMatrixPreconditioner::RMatrixPreconditioner(const RSparseMatrix &matrix, RMatrixPreconditionerType matrixPreconditionerType, unsigned int blockSize)
    : matrixPreconditionerType(matrixPreconditionerType)
{
    this->_init();

    switch (matrixPreconditionerType)
    {
        case R_MATRIX_PRECONDITIONER_NONE:
            break;
        case R_MATRIX_PRECONDITIONER_JACOBI:
            this->constructJacobi(matrix);
            break;
        case R_MATRIX_PRECONDITIONER_BLOCK_JACOBI:
            this->constructBlockJacobi(matrix,blockSize);
            break;
        default:
            throw RError(RError::Type::Application,R_ERROR_REF,"Invalid matrix preconditioner type \'%d\'",matrixPreconditionerType);
    }
}

RMatrixPreconditioner::RMatrixPreconditioner(const RMatrixPreconditioner &matrixPreconditioner)
{
    this->_init(&matrixPreconditioner);
}

RMatrixPreconditioner::~RMatrixPreconditioner()
{
}

RMatrixPreconditioner &RMatrixPreconditioner::operator =(const RMatrixPreconditioner &matrixPreconditioner)
{
    this->_init(&matrixPreconditioner);
    return (*this);
}

void RMatrixPreconditioner::compute(const RRVector &x, RRVector &y) const
{
    switch (matrixPreconditionerType)
    {
        case R_MATRIX_PRECONDITIONER_NONE:
        {
#pragma omp single
            y.resize(x.getNRows());
#pragma omp barrier
#pragma omp for
            for (int64_t i=0;i<int64_t(x.getNRows());i++)
            {
                y[uint(i)] = x[uint(i)];
            }
            break;
        }
        case R_MATRIX_PRECONDITIONER_JACOBI:
            this->computeJacobi(x,y);
            break;
        case R_MATRIX_PRECONDITIONER_BLOCK_JACOBI:
            this->computeBlockJacobi(x,y);
            break;
        default:
            throw RError(RError::Type::Application,R_ERROR_REF,"Invalid matrix preconditioner type \'%d\'",matrixPreconditionerType);
    }
}

void RMatrixPreconditioner::constructJacobi(const RSparseMatrix &matrix)
{
    unsigned int nRows = matrix.getNRows();

    this->data.resize(nRows,1);
    this->data.fill(0.0);

    for (unsigned int i=0;i<nRows;i++)
    {
        unsigned int columnPosition = 0;
        if (matrix.findColumnPosition(i,i,columnPosition))
        {
            double value = matrix.getValue(i,columnPosition);
            this->data[i][0] = (value == 0.0) ? 0.0 : 1.0 / value;
        }
    }
}

void RMatrixPreconditioner::constructBlockJacobi(const RSparseMatrix &matrix, unsigned int blockSize)
{
    unsigned int width = blockSize;

    // Matrix mus be dividable by block size.
    unsigned int nRows = matrix.getNRows();
    while (nRows % width != 0)
    {
        width --;
    }

    this->data.resize(nRows,width,0.0);

    unsigned int nBlocks = nRows/width;

#pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(nBlocks);i++)
    {
        for (unsigned int k=0;k<width;k++)
        {
            unsigned int m = uint(i)*width + k;
            for (unsigned int l=0;l<width;l++)
            {
                unsigned int n = uint(i)*width + l;
                unsigned int columnPosition = 0;
                if (matrix.findColumnPosition(m,n,columnPosition))
                {
                    this->data[m][l] = matrix.getValue(m,columnPosition);
                }
            }
        }

        RRMatrix inverse(width,width,0.0);
        for (unsigned int k=0;k<width;k++)
        {
            unsigned int m = uint(i)*width + k;
            for (unsigned int l=0;l<width;l++)
            {
                inverse[k][l] = this->data[m][l];
            }
        }
        inverse.invert();

        for (unsigned int k=0;k<width;k++)
        {
            unsigned int m = uint(i)*width + k;
            for (unsigned int l=0;l<width;l++)
            {
                this->data[m][l] = inverse[k][l];
            }
        }
    }
}

void RMatrixPreconditioner::computeJacobi(const RRVector &x, RRVector &y) const
{
    unsigned int nRows = this->data.getNRows();

#pragma omp single
    y.resize(nRows);
#pragma omp barrier

#pragma omp for
    for (int64_t i=0;i<int64_t(nRows);i++)
    {
        y[uint(i)] = this->data[uint(i)][0] * x[uint(i)];
    }
}

void RMatrixPreconditioner::computeBlockJacobi(const RRVector &x, RRVector &y) const
{
    unsigned int nRows = this->data.getNRows();
    unsigned int width = this->data.getNColumns();

#pragma omp single
    y.resize(nRows);
#pragma omp barrier

    unsigned int nBlocks = nRows/width;

#pragma omp for
    for (int64_t i=0;i<int64_t(nBlocks);i++)
    {
        for (unsigned int k=0;k<width;k++)
        {
            unsigned int m = uint(i)*width + k;
            double value = 0.0;
            for (unsigned int l=0;l<width;l++)
            {
                value += this->data[m][l] * x[uint(i)*width + l];
            }
            y[m] = value;
        }
    }
}
