//
// Created by Davide Costantini on 17/05/24.
//

#include "GaussSeidel.h"

namespace GaussSeidel
{
    Matrix getP(const Matrix& matrix)
    {
        return matrix.triangularView<Eigen::Lower>();
    }

    Matrix getN(const Matrix& matrix)
    {
        Matrix N(matrix.triangularView<Eigen::Upper>());
        N *= -1;
        // TODO: rimuovere da N i valori sulla diagonale

        return N;
    }

    Result getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk)
    {
        // TODO: controllare teorema 4.3

        Matrix P = GaussSeidel::getP(A);
        Matrix N = GaussSeidel::getN(A);

        Eigen::VectorXf r = b - (A * xk);
        // TODO: y = P.inverse() * rk
        Eigen::VectorXf nextXk = xk; // + y

        return {nextXk, r};
    }

    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const float tolerance)
    {
        int i = 0;
        Result r = {xk, xk};

        while (i < MAX_ITERATIONS && !thresholdReached(A, b, r.x, tolerance))
        {
            r = getNextXk(A, b, r.x);
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
