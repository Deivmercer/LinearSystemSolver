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
    Eigen::VectorXf x;
    x.resize(matrix.cols());
    x.setOnes();
    Eigen::VectorXf b;
    b = matrix * x;

    GaussSeidel::Result rgs = GaussSeidel::solve(matrix, b, x0, 0.001f);
    std::cout << "[GaussSeidel::solve] x: " << rgs.x.transpose() << std::endl;
    std::cout << "[GaussSeidel::solve] residual: " << rgs.residual.transpose() << std::endl;
    float gserror = Utils::euclideanNorm(rgs.x, x) / Utils::euclideanNorm(x);
    std::cout << "[GaussSeidel::solve] Relative error: " << gserror << std::endl << std::endl;

    GradientMethod::Result rgm = GradientMethod::solve(matrix, b, x0, 0.001f);
    std::cout << "[GradientMethod::solve] x: " << rgm.x.transpose() << std::endl;
    std::cout << "[GradientMethod::solve] residual: " << rgm.residual.transpose() << std::endl;
    float gmerror = Utils::euclideanNorm(rgm.x, x) / Utils::euclideanNorm(x);
    std::cout << "[GaussSeidel::solve] Relative error: " << gmerror << std::endl << std::endl;

    Jacobi::Result rjm = Jacobi::solve(matrix, b, x0, 0.001f);
    std::cout << "[JacobiMethod::solve] x: " << rjm.x.transpose() << std::endl;
    std::cout << "[JacobiMethod::solve] residual: " << rjm.residual.transpose() << std::endl;
    float jmerror = Utils::euclideanNorm(rjm.x, x) / Utils::euclideanNorm(x);
    std::cout << "[JacobiMethod::solve] Relative error: " << jmerror << std::endl << std::endl;
    return 0;
}
