#ifndef JACOBI_H
#define JACOBI_H

#include "types.h"

namespace Jacobi
{
    struct Result
    {
        Eigen::VectorXf xk;
        Eigen::VectorXf residual;
        int numIter;
        float residualNorm;
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
