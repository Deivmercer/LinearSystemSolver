#include "Jacobi.h"
#include "LinearSystemSolver/Utils.h"

#ifndef NDEBUG
#include <iostream>
#endif

namespace Jacobi
{
    Result getNextXk(const Matrix &A, const Eigen::VectorXf &b, const Eigen::VectorXf &xk, Eigen::VectorXf invP)
    {

        //x^(k+1) = x^k + P^(-1)(Ax^(k) - b)
        //rk = b − Ax^(k)
        Eigen::VectorXf residual = b - (A * xk);
        Eigen::VectorXf rk = invP.asDiagonal() * residual;
        Eigen::VectorXf nextXk = xk + rk;

        return {nextXk, residual};
    }

    Result solve(const Matrix &A, const Eigen::VectorXf &b, const Eigen::VectorXf &xk, const float tolerance)
    {
        Result result = {xk, xk};
        Eigen::VectorXf invP = Utils::invertDiagonal(A);

        int i = 0;

        while (i < 10000 && !(Utils::thresholdReached(A, b, result.x, tolerance)))
        {
            result = getNextXk(A, b, result.x, invP);
            ++i;
        }

#ifndef NDEBUG
        std::cout << "[Jacobi::solve] Iteration: " << i + 1 << std::endl;
        std::cout << "[Jacobi::solve] x: " << result.x.transpose() << std::endl;
        std::cout << "[Jacobi::solve] residual: " << result.residual.transpose() << std::endl << std::endl;
#endif

        return result;
    }
}
