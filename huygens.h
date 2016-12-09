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
#include <time.h>
#include <string>
#include <fstream>
#include <sstream>
#include <complex>

/*---- Qt ----*/
#include <QObject>

/*---- OpenCV ----*/
#include <opencv2/core/core.hpp>
#include <opencv/highgui.h>

#include <omp.h>

using namespace std;
using namespace cv;

/**
 * Function to compute the diffraction look up table using the spectral picture of the pattern.
 * @brief computeDiffractionTable_Measurement2D
 * @param object
 * @param spectrumName
 * @param width
 * @param height
 * @param numberOfWavelengths
 * @param power
 */
void computeDiffractionTable_Measurement2D(Mat const &diffractionPattern, Point2f const &center, float F0,
                                           float widthObjectCm, float heightObjectCm, float widthObjectPx, float heightObjectPx,
                                           float lambdaMeasurement, int colorChannel,
                                           float distanceLightSource,
                                           string spectrumName, int widthTable, int heightTable, int numberOfWavelengths, float power);
#endif
