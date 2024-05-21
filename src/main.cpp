#include <iostream>
#include "types.h"
#include "GaussSeidel.h"
#include "Jacobi.h"
#include "GradientMethod.h"
#include "libs/unsupported/Eigen/SparseExtra"
#include "Utils.h"

int main()
{
    // Parsing the MatrixMarket file
    Matrix matrix;
    Eigen::loadMarket(matrix, "../sample-mtx/spa1.mtx");
    if(!Utils::checkSize(matrix))
    {
        throw std::invalid_argument("Matrix not valid");
    }
    // Initializing the initial guess
    Eigen::VectorXf x0;
    x0.resize(matrix.cols());
    x0.setZero();

    // Initializing the b vector
    Eigen::VectorXf b;
    b.resize(matrix.cols());
    b.setOnes();

    GaussSeidel::Result rgs = GaussSeidel::solve(matrix, b, x0, 0.001f);
    std::cout << "[GaussSeidel::solve] x: " << rgs.x.transpose() << std::endl;
    std::cout << "[GaussSeidel::solve] residual: " << rgs.residual.transpose() << std::endl << std::endl;

    GradientMethod::Result rgm = GradientMethod::solve(matrix, b, x0, 0.001f);
    std::cout << "[GradientMethod::solve] x: " << rgm.x.transpose() << std::endl;
    std::cout << "[GradientMethod::solve] residual: " << rgm.residual.transpose() << std::endl << std::endl;

    Jacobi::Result r = Jacobi::solve(matrix, b, x0, 0.001f);

    return 0;
}
