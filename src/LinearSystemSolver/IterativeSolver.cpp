//
// Created by David on 31/05/2024.
//

#include "IterativeSolver.h"
#include "Utils.h"

#ifndef NDEBUG
#include <iostream>
#endif

IterativeSolver::IterativeSolver() : n_iteration(0)
{

}

void IterativeSolver::solve(const Matrix &A, const Eigen::VectorXf &b, const Eigen::VectorXf &xk, float tolerance)
{
    x = xk;

    while (n_iteration < MAX_ITERATIONS && !Utils::thresholdReached(A, b, x, tolerance))
    {
        getNextXk(A, b, x);

#ifndef NDEBUG
        std::cout << "[solve] Iteration: " << n_iteration + 1 << std::endl;
        std::cout << "[solve] x: " << x.transpose() << std::endl;
        std::cout << "[solve] residual: " << residual.transpose() << std::endl << std::endl;
#endif

        ++n_iteration;
    }
}

Eigen::VectorXf IterativeSolver::getX() const
{
    return x;
}

Eigen::VectorXf IterativeSolver::getResidual() const
{
    return residual;
}

int IterativeSolver::getIterations() const
{
    return n_iteration;
}