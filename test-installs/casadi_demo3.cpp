/**
   Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
   @author: Joel Andersson, KU Leuven 2013
*/

#include <pinocchio/autodiff/casadi.hpp>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

// A function working with Eigen::Matrix'es parameterized by the Scalar type
template <typename Scalar, typename T1, typename T2, typename T3, typename T4>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
eigenFun(Eigen::MatrixBase<T1> const &A,
         Eigen::MatrixBase<T2> const &a,
         Eigen::MatrixBase<T3> const &B,
         Eigen::MatrixBase<T4> const &b)
{
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> c(4);
  c.segment(1, 3) = A * a.segment(1, 3) - B.transpose() * b;
  c[0] = 0.;

  return c;
}

int main(int argc, char *argv[])
{
  // Declare casadi symbolic matrix arguments
  casadi::SX cs_a = casadi::SX::sym("a", 4);
  casadi::SX cs_b = casadi::SX::sym("b", 3);

  // Declare Eigen matrices
  Eigen::Matrix<casadi::SX, 3, 3> A, B;
  Eigen::Matrix<casadi::SX, -1, 1> a(4), c(4);
  Eigen::Matrix<casadi::SX, 3, 1> b;

  std::cout << "A = \n" << A << std::endl;
  std::cout << "a = \n" << a << std::endl;
  std::cout << "b = \n" << b << std::endl;

  // Let A, B be some numeric matrices
  for (Eigen::Index i = 0; i < A.rows(); ++i)
  {
    for (Eigen::Index j = 0; j < A.cols(); ++j)
    {
      A(i, j) = 10. * static_cast<double>(i) + static_cast<double>(j);
      B(i, j) = -10. * static_cast<double>(i) - static_cast<double>(j);
    }
  }

  std::cout << "A = \n" << A << std::endl;

  // Let a, b be symbolic arguments of a function
  pinocchio::casadi::copy(cs_b, b);
  pinocchio::casadi::copy(cs_a, a);

  // Call the function taking Eigen matrices
  c = eigenFun<casadi::SX>(A, a, B, b);

  // Casadi dense matrix to Eigen sparse matrix
  casadi::DM Matrx = casadi::DM::eye(7);
  casadi::Sparsity SpA = Matrx.get_sparsity();

  std::vector<casadi_int> output_row, output_col;
  SpA.get_triplet(output_row, output_col);
  std::vector<double> values = Matrx.get_nonzeros();

  using T = Eigen::Triplet<double>;
  std::vector<T> TripletList;
  TripletList.resize(values.size());
  for(int k = 0; k < values.size(); ++k)
      TripletList[k] = T(output_row[k], output_col[k], values[k]);

  Eigen::SparseMatrix<double> SpMatrx(Matrx.size1(), Matrx.size2());
  SpMatrx.setFromTriplets(TripletList.begin(), TripletList.end());
  std::cout << "SpMatrx = \n" << SpMatrx << std::endl;

  // Eigen sparse matrix to casadi dense matrix
  Eigen::MatrixXd dMat = Eigen::MatrixXd(SpMatrx);
  size_t rows = dMat.rows();
  size_t cols = dMat.cols();

  Matrx.resize(rows,cols);
  Matrx = casadi::DM::zeros(rows,cols);

  std::memcpy(Matrx.ptr(), dMat.data(), sizeof(double)*rows*cols); 
  std::cout << "Matrx = \n" << Matrx << std::endl;

  // Copy the result from Eigen matrices to casadi matrix
  casadi::SX cs_c = casadi::SX(casadi::Sparsity::dense(c.rows(), 1));
  pinocchio::casadi::copy(c, cs_c);

  // Display the resulting casadi matrix
  std::cout << "c = \n" << cs_c << std::endl;

  // Do some AD
  casadi::SX dc_da = jacobian(cs_c, cs_a);

  // Display the resulting jacobian
  std::cout << "dc/da = " << dc_da << std::endl;

  // Create a function which takes a, b and returns c and dc_da
  casadi::Function fun("fun", casadi::SXVector{cs_a, cs_b}, casadi::SXVector{cs_c, dc_da});
  std::cout << "fun = \n" << fun << std::endl;

  // Evaluate the function
  casadi::DMVector res = fun(casadi::DMVector{std::vector<double>{1., 2., 3., 4.}, std::vector<double>{-1., -2., -3.}});
  std::cout << "fun(a, b)= \n" << res << std::endl;

  return 0;
}
