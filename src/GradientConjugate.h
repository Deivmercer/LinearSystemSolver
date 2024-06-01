#ifndef GRADIENTCONJUGATE_H
#define GRADIENTCONJUGATE_H

#include "LinearSystemSolver/types.h"

namespace GradientConjugate
{
    struct Result
    {
        Eigen::VectorXf x;
        Eigen::VectorXf residual;
        Eigen::VectorXf direction;
    };

    /*
     * INPUT
     * matrix A
     * vector b (termine noto)
     * vector xk (our temp solution)
     * residual rk
     * direction dk
     * tolerance
     */

    // chiamo questo per iniziare ed inizializzare rk e dk
    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance);

    // effettivo solver
    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const Eigen::VectorXf& rk, const Eigen::VectorXf& dk, float tolerance);

    Eigen::VectorXf getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const Eigen::VectorXf& rk, const Eigen::VectorXf& dk);

    Eigen::VectorXf getNextDk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const Eigen::VectorXf& rk, const Eigen::VectorXf& dk);

};

#endif //GRADIENTCONJUGATE_H
