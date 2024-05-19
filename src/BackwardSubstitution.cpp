//
// Created by Davide Costantini on 18/05/24.
//

#include "BackwardSubstitution.h"

#ifndef NDEBUG
#include <iostream>
#endif

namespace BackwardSubstitution
{
    Eigen::VectorXf solve(const Matrix& U, const Eigen::VectorXf& b)
    {
        Eigen::Index n = U.rows() - 1;
        Eigen::VectorXf x(b.size());
        x.setZero();

        assert(U.coeff(n, n) != 0);
        x.coeffRef(n) = b.coeff(n) / U.coeff(n, n);

        for (Eigen::Index i = n - 1; i > -1; --i)
        {
            assert(U.coeff(i, i) != 0);
            x.coeffRef(i) = (b.coeff(i) - (U.row(i) * x)) / U.coeff(i, i);
        }

#ifndef NDEBUG
        std::cout << "[BackwardSubstitution::solve] x: " << x.transpose() << std::endl << std::endl;
#endif

        return x;
    }
}