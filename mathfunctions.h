/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#ifndef MATHFUNCTIONS
#define MATHFUNCTIONS

#define M_PI 3.14159265358979323846

#include <iostream>
#include <cmath>
#include <vector>

/*---- OpenCV ----*/
#include <opencv2/core/core.hpp>

#include <omp.h>

void numericalIntegration_onLambda(int colorChannel, std::vector<float *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, float samplingDistanceWavelength, cv::Mat& integrationResult);

#endif // MATHFUNCTIONS

