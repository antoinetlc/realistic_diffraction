/**
  * Practical acquisition and rendering of diffraction effects in surface reflectance.
  * Antoine Toisoul, Abhijeet Ghosh.
  * Imperial College London, December 2016.
  */

#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#define M_PI 3.14159265358979323846

#include <cmath>
#include <vector>

/*---- OpenCV ----*/
#include <opencv2/core/core.hpp>

#include <omp.h>

/**
 * Numerical integration on the wavelength using the trapezoidal rule.
 * @brief numericalIntegration_onLambda
 * @param colorChannel
 * @param integrand_forAllUVLambda
 * @param width
 * @param height
 * @param numberOfWavelengths
 * @param samplingDistanceWavelength
 * @param integrationResult
 */
void numericalIntegration_onLambda(int colorChannel, std::vector<float *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, float samplingDistanceWavelength, cv::Mat& integrationResult);

#endif // MATHFUNCTIONS_H

