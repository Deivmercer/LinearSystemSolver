#ifndef JACOBI_H
#define JACOBI_H

#include "IterativeSolver.h"
#include "Utils.h"

class Jacobi : public IterativeSolver
{
    Eigen::VectorXf invP;

    void getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk) override;

public:
    void solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance)
    {
        invP = Utils::invertDiagonal(A);
        // Chiamare il metodo solve della superclasse
        IterativeSolver::solve(A, b, xk, tolerance);
    }
};

#endif //JACOBI_H