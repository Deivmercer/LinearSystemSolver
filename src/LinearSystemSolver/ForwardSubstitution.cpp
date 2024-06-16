//
// Created by Davide Costantini on 18/05/24.
//

#include "ForwardSubstitution.h"

#ifndef NDEBUG
#include <iostream>
#endif

Eigen::VectorXf ForwardSubstitution::solve(const Matrix& L, const Eigen::VectorXf& b)
{
    Eigen::VectorXf x(b.size());
    x.setZero();

    assert(L.coeff(0, 0) != 0);
    x.coeffRef(0) = b.coeff(0) / L.coeff(0, 0);

    for (Eigen::Index i = 1; i < L.rows(); ++i)
    {
        assert(L.coeff(i, i) != 0);
        x.coeffRef(i) = (b.coeff(i) - (L.row(i) * x)) / L.coeff(i, i);
    }

#ifndef NDEBUG
    std::cout << "[ForwardSubstitution::solve] x: " << x.transpose() << std::endl << std::endl;
#endif

    return x;
}