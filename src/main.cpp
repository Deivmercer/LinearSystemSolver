#include <iostream>
#include "main.h"
#include "LinearSystemSolver/GaussSeidel.h"
#include "LinearSystemSolver/GradientConjugate.h"
#include "LinearSystemSolver/GradientMethod.h"
#include "LinearSystemSolver/Jacobi.h"
#include "LinearSystemSolver/libs/unsupported/Eigen/SparseExtra"
#include "LinearSystemSolver/Utils.h"

int main()
{
    DemoMatrix demo = DEMO_MATRICES[0];

    // Parsing the MatrixMarket file
    Matrix matrix;
    Eigen::loadMarket(matrix, demo.path());
    if(!Utils::checkSize(matrix))
    {
        throw std::invalid_argument("Matrix not valid");
    }

    std::cout << "Matrix: " << demo.name << std::endl;
    std::cout << "Tolerance: " << demo.tolerance << std::endl << std::endl;

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




    GradientConjugate gc;
    gc.solve(matrix, b, x0, demo.tolerance);
    std::cout << "[GradientConjugate::solve] Iterations: " << gc.getIterations() << std::endl;
    std::cout << "[GradientConjugate::solve] x: " << gc.getX().transpose() << std::endl;
    std::cout << "[GradientConjugate::solve] direction: " << gc.getDirection().transpose() << std::endl;
    std::cout << "[GradientConjugate::solve] residual: " << gc.getResidual().transpose() << std::endl << std::endl;
    float gcerror = Utils::euclideanNorm(gc.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GradientConjugate::solve] Relative error: " << gcerror << std::endl << std::endl;

    return 0;
}
