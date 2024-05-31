#include "GradientConjugate.h"
#include "Utils.h"

#ifndef NDEBUG
#include <iostream>
#endif

namespace GradientConjugate
{

    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance)
    {
        //r^(0) = b − Ax^(0) e d^(0) = r^(0)
        Eigen::VectorXf r0 = b - (A * xk);
        Eigen::VectorXf d0 = r0;

        return solve(A,b, xk, r0, d0,  tolerance);
    }

    Result solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const Eigen::VectorXf& rk, const Eigen::VectorXf& dk, float tolerance)
    {
        Result result = {xk, rk, dk};
        int i = 0;

        while (i < 10000 && !(Utils::thresholdReached(A, b, result.x, tolerance)))
        {
            result.x = getNextXk(A, b, result.x, result.residual, result.direction);
            result.residual = b - (A * result.x);
            result.direction = getNextDk(A, b, result.x, result.residual, result.direction);;
            ++i;
        }

#ifndef NDEBUG
        std::cout << "[GradientConjugate::solve] Iteration: " << i + 1 << std::endl;
        std::cout << "[GradientConjugate::solve] x: " << result.x.transpose() << std::endl;
        std::cout << "[GradientConjugate::solve] direction: " << result.direction.transpose() << std::endl << std::endl;
        std::cout << "[GradientConjugate::solve] residual: " << result.residual.transpose() << std::endl << std::endl;
#endif

        return result;
    }

    /*
     *1: r^(k) = b − Ax^(k); rk -> residuo
     *2: y^(k) = A*d^(k); yk
     *3: z^(k) = A*r^(k); zk
     *4: αk = (d^(k)· r^(k)) / (d^(k)· y^(k)); ak -> coefficiente, esso richiede un prodotto scalar
     *5: x^(k+1) = x^(k) + αkd^(k);
     */
    Eigen::VectorXf getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const Eigen::VectorXf& rk, const Eigen::VectorXf& dk)
    {
        Eigen::VectorXf rk_new = b - (A * xk);
        Eigen::VectorXf yk = A * dk;
        Eigen::VectorXf zk = A * rk_new;
        float ak = (dk.dot(rk)) / (dk.dot(yk));
        Eigen::VectorXf xk_next = xk + (ak * dk);



        //Return xk_next e dk_next (iterata successiva e direzione)
        return xk_next;
    }

    /*
    * 6: r^(k+1) = b − Ax^(k+1);
    * 7: w^(k) = Ar^(k+1);
    * 8: βk = (d^(k)· w^(k))/(d^(k)· y^(k));
    * 9: d^(k+1) = r^(k+1) − βkd^(k);
    */
    Eigen::VectorXf getNextDk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const Eigen::VectorXf& rk, const Eigen::VectorXf& dk)
    {
        Eigen::VectorXf yk = A * dk; //Mi sere per calcolare yk. C'è un modo per avere metodi divisi ma mantenere valore o ricalcolo sempre?

        Eigen::VectorXf rk_next = b - (A * xk);
        Eigen::VectorXf wk = A * rk_next;
        float bk = (dk.dot(wk)) / (dk.dot(yk));
        Eigen::VectorXf dk_next = rk_next - (bk * dk);

        return dk_next;
    }

}