#include <iostream>
#include "main.h"
#include "LinearSystemSolver/GaussSeidel.h"
#include "GradientConjugate.h"
#include "Jacobi.h"
#include "LinearSystemSolver/GradientMethod.h"
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

    GaussSeidel gs;
    gs.solve(matrix, b, x0, demo.tolerance);
    std::cout << "[GaussSeidel::solve] Iterations: " << gs.getIterations() << std::endl;
    std::cout << "[GaussSeidel::solve] x: " << gs.getX().transpose() << std::endl;
    std::cout << "[GaussSeidel::solve] residual: " << gs.getResidual().transpose() << std::endl << std::endl;
    float gserror = Utils::euclideanNorm(gs.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GaussSeidel::solve] Relative error: " << gserror << std::endl << std::endl;

    GradientMethod gm;
    gm.solve(matrix, b, x0, demo.tolerance);
    std::cout << "[GaussSeidel::solve] Iterations: " << gm.getIterations() << std::endl;
    std::cout << "[GradientMethod::solve] x: " << gm.getX().transpose() << std::endl;
    std::cout << "[GradientMethod::solve] residual: " << gm.getResidual().transpose() << std::endl << std::endl;
    float gmerror = Utils::euclideanNorm(gm.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GradientMethod::solve] Relative error: " << gmerror << std::endl << std::endl;

    Jacobi::Result rjm = Jacobi::solve(matrix, b, x0, demo.tolerance);
    std::cout << "[JacobiMethod::solve] x: " << rjm.x.transpose() << std::endl;
    std::cout << "[JacobitMethod::solve] residual: " << rjm.residual.transpose() << std::endl << std::endl;
    float jmerror = Utils::euclideanNorm(rjm.x, x) / Utils::euclideanNorm(x);
    std::cout << "[JacobiMethod::solve] Relative error: " << jmerror << std::endl << std::endl;


    GradientConjugate::Result rgc = GradientConjugate::solve(matrix, b, x0, demo.tolerance);
    std::cout << "[GradientConjugate::solve] x: " << rgc.x.transpose() << std::endl;
    std::cout << "[GradientConjugate::solve] direction: " << rgc.direction.transpose() << std::endl;
    std::cout << "[GradientConjugate::solve] residual: " << rgc.residual.transpose() << std::endl << std::endl;
    float gcerror = Utils::euclideanNorm(rjm.x, x) / Utils::euclideanNorm(x);
    std::cout << "[GradientConjugate::solve] Relative error: " << gcerror << std::endl << std::endl;

    return 0;
}
