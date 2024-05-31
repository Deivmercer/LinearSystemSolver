
#ifndef ITERATIVESOLVER_H
#define ITERATIVESOLVER_H



class IterativeSolver {

    private:
     Eigen::VectorXf x;
     Eigen::VectorXf residual;
     int n_iteration;

    public:
     void solve(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, float tolerance);

    virtual void getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk);
};



#endif //ITERATIVESOLVER_H
