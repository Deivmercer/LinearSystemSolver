//
// Created by Davide Costantini on 19/05/24.
//

#ifndef GRADIENTMETHOD_H
#define GRADIENTMETHOD_H

#include "types.h"

namespace GradientMethod
{
    struct Result
    {
        Eigen::VectorXf x;
        Eigen::VectorXf residual;
    };

    Result getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk);

    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance);

    const int MAX_ITERATIONS = 1000;
};

#endif //GRADIENTMETHOD_H
