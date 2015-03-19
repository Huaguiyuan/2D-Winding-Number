#include "twodwind.h"
#include <iostream>

int main()
{
    double mumin_   = -3.0;     //minimum chemical potential
    double mumax_   = 3.0;      //maximum chemical potential
    int muint_      = 21;       //descritisation of chemical potential
    double deltamin_= -3.0;     //minimum pairing coupling
    double deltamax_= 3.0;      //maximum pairing coupling
    int deltaint_   = 21;       //descretisation of pairing coupling
    int pint_       = 31;      //descretisation of momentum space


    //wind test(1.0,1.0,150);
    //std::cout << test.getWind();

    //chern test(-1.0,3.0,200);
    //std::cout << test.getChern();

    phasespace PS(mumin_, mumax_, muint_, deltamin_, deltamax_, deltaint_, pint_);
    return 0;
}
