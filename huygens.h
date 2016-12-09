/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#ifndef DIFFRACTIONSHADERS_H
#define DIFFRACTIONSHADERS_H

#define M_PI 3.14159265358979323846
#define EPSILON 0.000001

#include "mathfunctions.h"
#include "spectrum.h"
#include "PFMReadWrite.h"

/*---- Standard library ----*/
#include <iostream>
#include <cmath>
#include <string>

/*---- OpenCV ----*/
#include <opencv2/core/core.hpp>
#include <opencv/highgui.h>

#include <omp.h>

/**
 * Function to compute the diffraction look up table given the diffraction pattern measured with a spectral filter.
 * @brief computeDiffractionTable_Measurement2D
 * @param diffractionPattern : image of the diffraction pattern with the spectral filter.
 * @param center : location of the specular lobe (order 0 of diffraction).
 * @param F0 : Fresnel value at normal incidence of the sample.
 * @param widthObjectCm : width in centimeters of the object.
 * @param heightObjectCm : height in centimeters of the object.
 * @param widthObjectPx : width in pixels of the object.
 * @param heightObjectPx : height in pixels of the object.
 * @param lambdaMeasurement : wavelength at which the diffraction pattern was measured.
 * @param colorChannel : color channel on which the intensity of the Mat diffractionPattern will be read (0 is blue, 1 is green, 2 is red).
 * @param distanceLightSource : distance between the sample and the light source during the measurement.
 * @param spectrumName : name of the spectrum that will be loaded (full white by default).
 * @param widthTable : width of the final diffraction look up table.
 * @param heightTable : height of the final diffraction look up table.
 * @param power : use and odd power (e.g, 3.0, 5.0) for non linear sampling.
 * @param numberOfWavelengths : Number of wavelengths used in the sampling. Note the the CIE color matching functions are given every 5 nanometers which corresponds to 81 wavelengths in the range 380 780 nanometers.
 */
void computeDiffractionTable_Measurement2D(cv::Mat const &diffractionPattern, cv::Point2f const &center, float F0,
                                           float widthObjectCm, float heightObjectCm, float widthObjectPx, float heightObjectPx,
                                           float lambdaMeasurement, int colorChannel,
                                           float distanceLightSource,
                                           std::string spectrumName, int widthTable, int heightTable, float power = 1.0, int numberOfWavelengths = 81);


cv::Mat XYZToRGB_sRGB(const cv::Mat &imageXYZ);

//Values for the CIE colour matching functions sampled ever 5 nanometers
static float cie_colour_matching_function[81][3] = {
       {(float) 0.0014,(float) 0.0000,(float) 0.0065}, {(float) 0.0022,(float) 0.0001,(float) 0.0105}, {(float) 0.0042,(float) 0.0001,(float) 0.0201},
       {(float) 0.0076,(float) 0.0002,(float) 0.0362}, {(float) 0.0143,(float) 0.0004,(float) 0.0679}, {(float) 0.0232,(float) 0.0006,(float) 0.1102},
       {(float) 0.0435,(float) 0.0012,(float) 0.2074}, {(float) 0.0776,(float) 0.0022,(float) 0.3713}, {(float) 0.1344,(float) 0.0040,(float) 0.6456},
       {(float) 0.2148,(float) 0.0073,(float) 1.0391}, {(float) 0.2839,(float) 0.0116,(float) 1.3856}, {(float) 0.3285,(float) 0.0168,(float) 1.6230},
       {(float) 0.3483,(float) 0.0230,(float) 1.7471}, {(float) 0.3481,(float) 0.0298,(float) 1.7826}, {(float) 0.3362,(float) 0.0380,(float) 1.7721},
       {(float) 0.3187,(float) 0.0480,(float) 1.7441}, {(float) 0.2908,(float) 0.0600,(float) 1.6692}, {(float) 0.2511,(float) 0.0739,(float) 1.5281},
       {(float) 0.1954,(float) 0.0910,(float) 1.2876}, {(float) 0.1421,(float) 0.1126,(float) 1.0419}, {(float) 0.0956,(float) 0.1390,(float) 0.8130},
       {(float) 0.0580,(float) 0.1693,(float) 0.6162}, {(float) 0.0320,(float) 0.2080,(float) 0.4652}, {(float) 0.0147,(float) 0.2586,(float) 0.3533},
       {(float) 0.0049,(float) 0.3230,(float) 0.2720}, {(float) 0.0024,(float) 0.4073,(float) 0.2123}, {(float) 0.0093,(float) 0.5030,(float) 0.1582},
       {(float) 0.0291,(float) 0.6082,(float) 0.1117}, {(float) 0.0633,(float) 0.7100,(float) 0.0782}, {(float) 0.1096,(float) 0.7932,(float) 0.0573},
       {(float) 0.1655,(float) 0.8620,(float) 0.0422}, {(float) 0.2257,(float) 0.9149,(float) 0.0298}, {(float) 0.2904,(float) 0.9540,(float) 0.0203},
       {(float) 0.3597,(float) 0.9803,(float) 0.0134}, {(float) 0.4334,(float) 0.9950,(float) 0.0087}, {(float) 0.5121,(float) 1.0000,(float) 0.0057},
       {(float) 0.5945,(float) 0.9950,(float) 0.0039}, {(float) 0.6784,(float) 0.9786,(float) 0.0027}, {(float) 0.7621,(float) 0.9520,(float) 0.0021},
       {(float) 0.8425,(float) 0.9154,(float) 0.0018}, {(float) 0.9163,(float) 0.8700,(float) 0.0017}, {(float) 0.9786,(float) 0.8163,(float) 0.0014},
       {(float) 1.0263,(float) 0.7570,(float) 0.0011}, {(float) 1.0567,(float) 0.6949,(float) 0.0010}, {(float) 1.0622,(float) 0.6310,(float) 0.0008},
       {(float) 1.0456,(float) 0.5668,(float) 0.0006}, {(float) 1.0026,(float) 0.5030,(float) 0.0003}, {(float) 0.9384,(float) 0.4412,(float) 0.0002},
       {(float) 0.8544,(float) 0.3810,(float) 0.0002}, {(float) 0.7514,(float) 0.3210,(float) 0.0001}, {(float) 0.6424,(float) 0.2650,(float) 0.0000},
       {(float) 0.5419,(float) 0.2170,(float) 0.0000}, {(float) 0.4479,(float) 0.1750,(float) 0.0000}, {(float) 0.3608,(float) 0.1382,(float) 0.0000},
       {(float) 0.2835,(float) 0.1070,(float) 0.0000}, {(float) 0.2187,(float) 0.0816,(float) 0.0000}, {(float) 0.1649,(float) 0.0610,(float) 0.0000},
       {(float) 0.1212,(float) 0.0446,(float) 0.0000}, {(float) 0.0874,(float) 0.0320,(float) 0.0000}, {(float) 0.0636,(float) 0.0232,(float) 0.0000},
       {(float) 0.0468,(float) 0.0170,(float) 0.0000}, {(float) 0.0329,(float) 0.0119,(float) 0.0000}, {(float) 0.0227,(float) 0.0082,(float) 0.0000},
       {(float) 0.0158,(float) 0.0057,(float) 0.0000}, {(float) 0.0114,(float) 0.0041,(float) 0.0000}, {(float) 0.0081,(float) 0.0029,(float) 0.0000},
       {(float) 0.0058,(float) 0.0021,(float) 0.0000}, {(float) 0.0041,(float) 0.0015,(float) 0.0000}, {(float) 0.0029,(float) 0.0010,(float) 0.0000},
       {(float) 0.0020,(float) 0.0007,(float) 0.0000}, {(float) 0.0014,(float) 0.0005,(float) 0.0000}, {(float) 0.0010,(float) 0.0004,(float) 0.0000},
       {(float) 0.0007,(float) 0.0002,(float) 0.0000}, {(float) 0.0005,(float) 0.0002,(float) 0.0000}, {(float) 0.0003,(float) 0.0001,(float) 0.0000},
       {(float) 0.0002,(float) 0.0001,(float) 0.0000}, {(float) 0.0002,(float) 0.0001,(float) 0.0000}, {(float) 0.0001,(float) 0.0000,(float) 0.0000},
       {(float) 0.0001,(float) 0.0000,(float) 0.0000}, {(float) 0.0001,(float) 0.0000,(float) 0.0000}, {(float) 0.0000,(float) 0.0000,(float) 0.0000}
};

#endif
