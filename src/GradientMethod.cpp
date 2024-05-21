//
// Created by Davide Costantini on 19/05/24.
//

#include "GradientMethod.h"
#include "Utils.h"

#ifndef NDEBUG
#include <iostream>
#endif

namespace GradientMethod
{
    Result getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk)
    {
        Eigen::VectorXf r = b - (A * xk);
        Eigen::VectorXf y = A * r;
        float a = r.transpose() * r;
        float bk = r.transpose() * y;
        float alpha = a / bk;
        Eigen::VectorXf nextXk = xk + (alpha * r);

        return {nextXk, r};

    }

    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance)
    {
        int i = 0;
        Result r = {xk, xk};

        while (i < MAX_ITERATIONS && !Utils::thresholdReached(A, b, r.x, tolerance))
        {
            r = getNextXk(A, b, r.x);

#ifndef NDEBUG
            std::cout << "[GradientMethod::solve] Iteration: " << i + 1 << std::endl;
            std::cout << "[GradientMethod::solve] x: " << r.x.transpose() << std::endl;
            std::cout << "[GradientMethod::solve] residual: " << r.residual.transpose() << std::endl << std::endl;
#endif

            ++i;
        }

        return r;
    }
}