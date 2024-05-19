//
// Created by Davide Costantini on 17/05/24.
//

#include "GaussSeidel.h"
#include "BackwardSubstitution.h"

#ifndef NDEBUG
#include <iostream>
#endif

namespace GaussSeidel
{
    Matrix getLowerMatrix(const Matrix& matrix)
    {
        return matrix.triangularView<Eigen::Lower>();
    }

    Matrix getUpperMatrix(const Matrix& matrix)
    {
        return matrix.triangularView<Eigen::Upper | Eigen::ZeroDiag>() * -1;
    }

    Result getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk)
    {
        // TODO: controllare teorema 4.3
        // Extracting the upper and lower triangular components

        Matrix P = GaussSeidel::getLowerMatrix(A);
        Matrix N = GaussSeidel::getUpperMatrix(A);

        Eigen::VectorXf r = b - (A * xk);
        Eigen::VectorXf y = BackwardSubstitution::solve(P, r);
        Eigen::VectorXf nextXk = xk + y;

        return {nextXk, r};
    }

    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const float tolerance)
    {
        int i = 0;
        Result r = {xk, xk};

        while (i < MAX_ITERATIONS && !thresholdReached(A, b, r.x, tolerance))
        {
            r = getNextXk(A, b, r.x);

#ifndef NDEBUG
            std::cout << "[GaussSeidel::solve] Iteration: " << i + 1 << std::endl;
            std::cout << "[GaussSeidel::solve] x: " << r.x.transpose() << std::endl;
            std::cout << "[GaussSeidel::solve] residual: " << r.residual.transpose() << std::endl << std::endl;
#endif

            ++i;
        }

        return r;
    }

    bool thresholdReached(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const float tolerance)
    {
        // ||b - Axk|| / ||b|| < tol
        Eigen::VectorXf axk = A * xk;
        Eigen::VectorXf normDistance = b - axk;

        return (normDistance.squaredNorm() / b.squaredNorm()) < tolerance;
    }
}
