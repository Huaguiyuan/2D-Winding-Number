#include "twodwind.h"
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <complex>
#include <eigen3/Eigen/Eigenvalues>

sys::sys(double pX_, double pY_, double mu_, double delta_)
{
    std::complex<double> eps( mu_ - (cos( pX_ ) + cos( pY_)), 0.0 );
    std::complex<double> gap( delta_ * sin( pX_ ), delta_ * sin( pY_ ) );

    hamiltonian(0,0) = 0.5 * eps;
    hamiltonian(0,1) = 0.5 * gap;
    hamiltonian(1,0) = 0.5 * std::conj(gap);
    hamiltonian(1,1) = -0.5 * eps;

    Eigen::ComplexEigenSolver<Eigen::Matrix2cd> ces(hamiltonian);

    Eigen::Vector2cd eigenvaluetemp;
    eigenvaluetemp = ces.eigenvalues();

    Eigen::Vector2cd temp;

    if (std::real(eigenvaluetemp[0]) < std::real(eigenvaluetemp[1]))
    {
        eigenvalue[0] = eigenvaluetemp[0];
        eigenvalue[1] = eigenvaluetemp[1];
        temp = ces.eigenvectors().col(0);
    }
    else
    {
        eigenvalue[0] = eigenvaluetemp[1];
        eigenvalue[1] = eigenvaluetemp[0];
        temp = ces.eigenvectors().col(1);
    }

    eigenvector = temp;

    sVec(0) = std::real(temp.dot(pauli(1).output() * temp));
    sVec(1) = std::real(temp.dot(pauli(2).output() * temp));
    sVec(2) = std::real(temp.dot(pauli(3).output() * temp));

    gapdiff = std::real(eigenvalue[1] - eigenvalue[0]);
   // std::cout << gapdiff << " ";
}

