#ifndef TWODWIND_H
#define TWODWIND_H

#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>

class pauli //defines a pauli matrix based on pauli_type_ = {1,2,3}
{
    private:
    Eigen::Matrix2cd emptyp;

    public:
    pauli( int pauli_type_ );
    Eigen::Matrix2cd output() { return emptyp;}
};

class sys //defines the hamiltonian and computes the eigenvectors, eigenvalues and winding vector
{
    private:
    Eigen::Matrix2cd hamiltonian;
    Eigen::Vector3d sVec;
    Eigen::Vector2cd eigenvalue;
    Eigen::Vector2cd eigenvector;
    double gapdiff;


    public:
    sys(double pX_, double pY_, double mu, double delta);
    Eigen::Vector3d output() { return sVec; }
    double getgap() { return gapdiff; }
    Eigen::Vector2cd vecoutput() { return eigenvector;}
};

class wind //calculates the winding number and the gap of the system
{
    private:
    std::vector<std::vector<sys>> momSpace;
    std::vector<std::vector<sys>> momSpaceShiftX;
    std::vector<std::vector<sys>> momSpaceShiftY;
    std::vector<double> pXVec;
    std::vector<double> pXVecShift;
    std::vector<double> pYVec;
    std::vector<double> pYVecShift;
    std::vector<double> gaps;
    double windN;
    double mingap;


    public:
    wind( double mu_, double delta_ , int pint_);
    void linspace(double pmin, double pmax, int pint);
    void printWind() {std::cout << windN;}
    double getWind() {return windN;}
    double getgapmin() { return mingap; }
};

class chern
{
    private:
    std::vector<std::vector<sys>> momSpace;
    std::vector<double> pVec;
    double ChernN;

    public:
    chern( double mu_, double delta_, int pint_ );
    void linspace(double pmin, double pmax, int pint);
    double getChern() { return ChernN; }
};

class phasespace
{
    private:
    std::vector<std::vector<wind>> space;
    std::vector<double> muVec;
    std::vector<double> deltaVec;

    public:
    phasespace(double mumin_, double mumax_, int muint_, double deltamin_, double deltamax_, int deltaint_, int pint_);
    void mulinspace(double mumin_, double mumax_, int muint_);
    void deltalinspace(double deltamin_, double deltamax_, int deltaint_);
    double getmu( int munum_ ) {return muVec[munum_];}
    double getdelta( int deltanum_ ) {return deltaVec[deltanum_];}
    double getWinding(int munum_, int deltanum_);
};

#endif
