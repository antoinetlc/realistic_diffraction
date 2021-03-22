/**
  * Practical acquisition and rendering of diffraction effects in surface reflectance.
  * Antoine Toisoul, Abhijeet Ghosh.
  * Imperial College London, December 2016.
  */

#include <iostream>
#include <math.h>
#include <string.h>

#include "PFMReadWrite.h"
#include "huygens.h"
#include "convolution.h"

using namespace std;
using namespace cv;

int main(int argc, char *argv[])
{
    //Exemple of usage
    //Load a HDR diffraction pattern captured with a spectral filter
    string path = string("path/to/pattern.pfm");
    Mat diffractionPattern = loadPFM(path);

    //Pattern measured with green spectral filter.
    float lambdaMeasurement = 0.53;
    float colorChannel = 1;

    //In centimeters
    float distanceLightSourceCm = 300.0;

    //Width and height in order to have a scale between centimeters and pixels
    //Necessary for the sine calculations
    float widthObjectCm= 7.46;
    float heightObjectCm = 14.63;

    float widthObjectPx = 1326.0;
    float heightObjectPx = 2393.0;

    //Center of the specular lobe : center of the pixels with the highest intensity
    Point2f center(1105, 943);

    //Load your light spectrum (full white here)
    vector<float> spectralPowerDistribution(81, 1.0);

    //Fresnel at normal incidence
    float F0 = 0.04;

    //Compute and display the table
    computeDiffractionTableFromSpectralMeasurement(diffractionPattern, center, F0, widthObjectCm, heightObjectCm, widthObjectPx, heightObjectPx,
                                         lambdaMeasurement, colorChannel, distanceLightSourceCm, spectralPowerDistribution, 1024, 1024, 5.0, 81);

    //Preconvolution of Environment map
    //The diffraction lookup table must be in RGB and sampled linearly
    Mat diffractionPattern = loadPFM("/path/to/diffractionPattern.pfm");

    Mat environmentMap = loadPFM("/path/to/environemntMap.pfm");
    Mat filteredEM = convolveEnvironmentMap(em, diffractionPattern, 1.0, 1.0);

    return 0;
}

