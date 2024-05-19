#include "Jacobi.h"
#include "Utils.h"

namespace Jacobi
{
    Eigen::VectorXf getNextXk(Matrix &A, const Eigen::VectorXf &b, const Eigen::VectorXf &xk)
    {
        return {0, 0};
    }

    Result solve(Matrix &A, const Eigen::VectorXf &b, const Eigen::VectorXf &xk, const float tolerance)
    {
        Eigen::VectorXf invP = Utils::invertDiagonal(A);

        Result r = {xk, xk, 0, 0};
        return r;
    }
}
