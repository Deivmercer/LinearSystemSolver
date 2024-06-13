#ifndef GRADIENTCONJUGATE_H
#define GRADIENTCONJUGATE_H

#include "IterativeSolver.h"
#include "Utils.h"

class GradientConjugate : public IterativeSolver
{
    Eigen::VectorXf direction;

    void getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk) override;

public:

    Eigen::VectorXf getDirection ()
    {
        return direction;
    }

    bool checkConvergence(const Matrix& A) override;
};

#endif //GRADIENTCONJUGATE_H