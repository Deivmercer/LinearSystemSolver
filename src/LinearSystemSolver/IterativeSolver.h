
#ifndef ITERATIVESOLVER_H
#define ITERATIVESOLVER_H

#include "types.h"

class IterativeSolver
{
    protected:
        Eigen::VectorXf x;
        Eigen::VectorXf residual;
        int n_iteration;

    public:
        IterativeSolver();
        const int MAX_ITERATIONS = 20000;
        void solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance);
        virtual void getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk) = 0;
        Eigen::VectorXf getX() const;
        Eigen::VectorXf getResidual() const;
        int getIterations() const;
        virtual bool checkConvergence(const Matrix& A) = 0;
};

#endif //ITERATIVESOLVER_H
