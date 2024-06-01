#ifndef JACOBI_H
#define JACOBI_H

#include "LinearSystemSolver/types.h"

namespace Jacobi
{
    struct Result
    {
        Eigen::VectorXf x;
        Eigen::VectorXf residual;
    };

    /*
     * INPUT
     * matrix A
     * vector b (termine noto)
     * vector xk (our temp solution)
     * tolerance
     */
    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance);

    Eigen::VectorXf getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk);

}


#endif //JACOBI_H
