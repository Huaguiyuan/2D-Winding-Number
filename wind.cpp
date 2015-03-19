#include "twodwind.h"
#include <iostream>

#define PI 3.141592

void wind::linspace(double pmin, double pmax, int pint)
{
    double step = (pmax-pmin) / (pint-1);
    double shift = 1e-2;

    while (!(pmin > pmax))
    {
        pXVec.push_back(pmin);
        pXVecShift.push_back(pmin+shift);
        pYVec.push_back(pmin);
        pYVecShift.push_back(pmin+shift);
        pmin += step;
        //std::cout << pmin << " ";
    }
}


wind::wind( double mu_, double delta_, int pint_ )
{
    int i,j;
    linspace(-PI,PI,pint_);

    momSpace.resize( pint_ );
    momSpaceShiftX.resize( pint_ );
    momSpaceShiftY.resize( pint_ );

    for ( auto& row : momSpace ) {
        row.reserve( pint_ );
    }

    for ( auto& row : momSpaceShiftX ) {
        row.reserve( pint_ );
    }

    for ( auto& row : momSpaceShiftY ) {
        row.reserve( pint_ );
    }

    for(i = 0; i < pint_; i++)
    {
        for( j = 0; j < pint_; j++ )
        {
            momSpace[i][j] = sys(pXVec[i],pYVec[j], mu_, delta_);
            momSpaceShiftX[i][j] = sys(pXVecShift[i],pYVec[j], mu_, delta_);
            momSpaceShiftY[i][j] = sys(pXVec[i],pYVecShift[j], mu_, delta_);
        }

    }


    double shift = 1e-2;
    windN = 0;
    Eigen::Vector3d temp, tempShiftX, tempShiftY, dsXVec, dsYVec, crossed;
    double windPart;
    double windTemp = 0;
    double areaelement = (2*PI)/pint_;
    double gaptemp = 100.0;

    for(i = 0; i < pint_-1; i++)
    {

        for(j = 0; j < pint_-1; j++)
        {
            temp = momSpace[i][j].output();
            tempShiftX = momSpaceShiftX[i][j].output();
            tempShiftY = momSpaceShiftY[i][j].output();
            dsXVec = (tempShiftX - temp)/shift;
            dsYVec = (tempShiftY - temp)/shift;
            crossed = dsYVec.cross(dsXVec);
            windPart = temp.dot(crossed);
            //std::cout << std::endl;
            //std::cout << dsXVec << std::endl;
            //std::cout << temp[0] << "\t" << temp[1] << "\t" << temp [2] << std::endl;
            //std::cout << dsXVec[0] << "\t" << dsXVec[1] << "\t" << dsXVec [2] << std::endl;
            //std::cout << crossed[0] << "\t" << crossed[1] << "\t" << crossed[2] << std::endl;
            windTemp += windPart;

            //std::cout << windTemp << "\t" << windPart << std::endl;
        }


        //std::cout << std::endl;
        //std::cout << tempShift << '\t' << temp;
    }

    windN = (1/(4*PI)) * windTemp * areaelement * areaelement;

    for (i = 0; i < pint_-1; i++)
    {
        for( j = 0; j<pint_ -1; j++ )
        {

            if ( momSpace[i][j].getgap() < gaptemp )
            {
                gaptemp = momSpace[i][j].getgap();
            }
        }
    }

    mingap = gaptemp;
}

