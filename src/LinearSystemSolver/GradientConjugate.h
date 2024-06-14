#ifndef GRADIENTCONJUGATE_H
#define GRADIENTCONJUGATE_H

#include "IterativeSolver.h"
#include "Utils.h"

class GradientConjugate : public IterativeSolver
{
    Eigen::VectorXf direction;

    void getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk) override;

public:
    void solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance)
    {
        //INIT residual and direction
        //r^(0) = b − Ax^(0) e d^(0) = r^(0)
        residual = b - (A * xk);
        direction = residual;

        IterativeSolver::solve(A, b, xk, tolerance);
    }

    Eigen::VectorXf getDirection ()
    {
        return direction;
    }

    bool checkConvergence(const Matrix& A) override;
};

#endif //GRADIENTCONJUGATE_H