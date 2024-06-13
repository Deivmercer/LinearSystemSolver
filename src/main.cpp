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
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;

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

    Jacobi jm;
    start = std::chrono::steady_clock::now();
    jm.solve(matrix, b, x0, demo.tolerance);
    end = std::chrono::steady_clock::now();
    std::cout << "[Jacobi] Iterations: " << jm.getIterations() << std::endl;
    std::cout << "[Jacobi] x: " << jm.getX().transpose() << std::endl;
    std::cout << "[Jacobi] residual: " << jm.getResidual().transpose() << std::endl;
    float jmerror = Utils::euclideanNorm(jm.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[Jacobi] Relative error: " << jmerror << std::endl;
    std::cout << "[Jacobi] Execution time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl << std::endl;

    GaussSeidel gs;
    start = std::chrono::steady_clock::now();
    gs.solve(matrix, b, x0, demo.tolerance);
    end = std::chrono::steady_clock::now();
    std::cout << "[GaussSeidel] Iterations: " << gs.getIterations() << std::endl;
    std::cout << "[GaussSeidel] x: " << gs.getX().transpose() << std::endl;
    std::cout << "[GaussSeidel] residual: " << gs.getResidual().transpose() << std::endl;
    float gserror = Utils::euclideanNorm(gs.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GaussSeidel] Relative error: " << gserror << std::endl;
    std::cout << "[GaussSeidel] Execution time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl << std::endl;

    GradientMethod gm;
    start = std::chrono::steady_clock::now();
    gm.solve(matrix, b, x0, demo.tolerance);
    end = std::chrono::steady_clock::now();
    std::cout << "[GradientMethod] Iterations: " << gm.getIterations() << std::endl;
    std::cout << "[GradientMethod] x: " << gm.getX().transpose() << std::endl;
    std::cout << "[GradientMethod] residual: " << gm.getResidual().transpose() << std::endl;
    float gmerror = Utils::euclideanNorm(gm.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GradientMethod] Relative error: " << gmerror << std::endl;
    std::cout << "[GradientMethod] Execution time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl << std::endl;

    GradientConjugate gc;
    start = std::chrono::steady_clock::now();
    gc.solve(matrix, b, x0, demo.tolerance);
    end = std::chrono::steady_clock::now();
    std::cout << "[GradientConjugate] Iterations: " << gc.getIterations() << std::endl;
    std::cout << "[GradientConjugate] x: " << gc.getX().transpose() << std::endl;
    std::cout << "[GradientConjugate] direction: " << gc.getDirection().transpose() << std::endl;
    std::cout << "[GradientConjugate] residual: " << gc.getResidual().transpose() << std::endl;
    float gcerror = Utils::euclideanNorm(gc.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GradientConjugate] Relative error: " << gcerror << std::endl;
    std::cout << "[GradientConjugate] Execution time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl << std::endl;

    return 0;
}
