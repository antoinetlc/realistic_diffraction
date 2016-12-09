/**
  * Practical acquisition and rendering of diffraction effects in surface reflectance.
  * Antoine Toisoul, Abhijeet Ghosh.
  * Imperial College London, December 2016.
  */

#include "mathfunctions.h"

using namespace std;
using namespace cv;

void numericalIntegration_onLambda(int colorChannel, vector<float *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult)
{
    for(int u = 0 ; u<height ; u++)
    {

        #pragma omp parallel for num_threads(omp_get_max_threads())
        for(int v = 0 ; v<width ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                float currentValue = integrand_forAllUVLambda[w][u*width+v];
                float nextValue = integrand_forAllUVLambda[w+1][u*width+v];

                integrationResult.at<Vec3f>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}
