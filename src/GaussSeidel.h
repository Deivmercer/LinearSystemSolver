//
// Created by Davide Costantini on 17/05/24.
//

#ifndef GAUSS_SEIDEL_GAUSSSEIDEL_H
#define GAUSS_SEIDEL_GAUSSSEIDEL_H

#include "types.h"
#include "libs/unsupported/Eigen/SparseExtra"

namespace GaussSeidel
{
    struct Result
    {
        Eigen::VectorXf x;
        Eigen::VectorXf residual;
    };

    Matrix getP(const Matrix& matrix);

    Matrix getN(const Matrix& matrix);

    Result getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk);

    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance);

    bool thresholdReached(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance);

    const int MAX_ITERATIONS = 1000;
}

#endif //GAUSS_SEIDEL_GAUSSSEIDEL_H
