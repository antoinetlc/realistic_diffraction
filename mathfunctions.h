/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#ifndef MATHFUNCTIONS
#define MATHFUNCTIONS

#define M_PI 3.14159265358979323846

#include "coloursystem.h"

#include <iostream>
#include <cmath>
#include <vector>

#include <QMatrix3x3>
#include <QVector3D>

/*---- OpenCV ----*/
#include <opencv2/core/core.hpp>

/*----- Dlib -----*/
#include <dlib/optimization.h>

#include <omp.h>

using namespace std;
using namespace cv;

void invertMatrix3x3(float matrix[][3], float invertedMatrix[][3]);
void printMatrix3x3(float matrix[][3]);

/**
 * Function that converts cartesian coordinates to spherical coordinates
 * theta is in the range [0:Pi]
 * phi is in the range [0:2Pi]
 * @param INPUT : x is a const float& corresponding to the x coordinate in the cartesian coordinate system.
 * @param INPUT : y is a const float& corresponding to the y coordinate in the cartesian coordinate system
 * @param INPUT : z is a const float& corresponding to the z coordinate in the cartesian coordinate system
 * @param OUTPUT : r is a float& corresponding to the r coordinate in the spherical coordinate system
 * @param OUTPUT : theta is a float& corresponding to the polar angle theta in the spherical coordinate system
 * @param OUTPUT : phi is a float& corresponding to the azimuthal angle phi in the spherical coordinate system
 * @param INPUT : offset is a const float that adds an offset to the angle phi.
 */
void cartesianToSpherical(const float &x, const float &y, const float &z, float &r, float &theta, float &phi);

/**
* Function that converts cartesian coordinates to spherical coordinates for an entire vector
* theta is in the range [0:Pi]
* phi is in the range [0:2Pi]
* @param INPUT : Each element of cartesianVector is a vector that contains 3 elements : x,y,z
* @param OUTPUT : Each element of sphericalVector is a vector that contains 3 elements : r, phi, theta
* @param INPUT : offset is a const float that adds an offset to the angle phi.
*/
void cartesianToSphericalVector(const vector<vector<float> > &cartesianVector, vector<vector<float> > &sphericalVector);

/**
* Calculate number mod(modulo) where modulo and number are real numbers
* @param INPUT : number is a float
* @param INPUT : modulo is a float
* @return outputs the result of number mod(modulo)
*/
float moduloRealNumber(float number, float modulo);

/**
 * Clamp a number so it is in the range [inf ; sup].
 * If number<inf, then the function returns inf.
 * If number>sup, then the function returns sup.
 * If number is in the range [inf ; sup] then the function returns the number.
 * @param INPUT : value is a float which will be clamped.
 * @param INPUT : inf is the lower bound of the range.
 * @param INPUT : sup is the upper bound of the range.
 * @return The clamped value in the range [inf ; sup]. If sup < inf then the function returns the number.
 */
float clamp(float value, float inf, float sup);

double sinc(double value);

double sincModel(double x, const dlib::matrix<double,1,1>& parameter);
double residual (const std::pair<double, double>& data, const dlib::matrix<double,1,1>& parameter);
double residual_derivative(const pair<double, double>& data, const dlib::matrix<double,1,1>& parameter);


double sincModel2(double x, const dlib::matrix<double,2,1> &parameter);
double residual2(const pair<double, double>& data, const dlib::matrix<double,2,1> &parameter);
dlib::matrix<double,2,1> residual_derivative2(const pair<double, double>& data, const dlib::matrix<double,2,1> &parameter);

/**
 * Compute a least square optimisation to fit a function on the data : (x, data). The result of the fitted function is in the array result.
 * offsetOrder0 is the distance around the central peak that will not be taken into account for the optimisation (avoid an overfitting on the central peak).
 * parameter is an ouput that contains the fitted parameter of the model function.
 * @brief leastSquaresOptimisation
 * @param x
 * @param data
 * @param size
 * @param result
 * @param offsetOrder0
 * @param parameter
 */
void leastSquaresOptimisation(double *x, double *data, int size, int offsetOrder0, double* result, dlib::matrix<double,1,1>& parameter);

/**
 * Function to calculate the integral of a sampled function. size is the number of elements in the array. The sampling is done every samplingDistance.
 * @brief numbericalIntegration
 * @param functionSamples
 * @param samplingDistance
 * @param inf
 * @param sup
 * @return
 */
float numericalIntegration(float functionSamples[], int size, float samplingDistance);

/**
 * Numerical integration over the visible spectrum using the trapezoidal rule.
 * @brief numericalIntegration_onLambda
 * @param p
 * @param current_n
 * @param current_m
 * @param colorChannel
 * @param factN
 * @param factM
 * @param integrand_forAllUVLambda
 * @param numberOfUV
 * @param numberOfWavelengths
 * @param samplingDistance
 * @param Ip
 */
void numericalIntegration_onLambda(int colorChannel, vector<float *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult);

void numericalIntegration_onLambda(int colorChannel, vector<float *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult);

void numericalIntegration_onLambda(int colorChannel, vector<double *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, double samplingDistanceWavelength, Mat& integrationResult);

void numericalIntegration_onLambda(int colorChannel, vector<double *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, double samplingDistanceWavelength, Mat& integrationResult);

/**
 * Numerical integration over the visible spectrum using the trapezoidal rule. No optimisations
 * @brief numericalIntegration_onLambda_non_optimised
 * @param p
 * @param colorChannel
 * @param factN
 * @param factM
 * @param integrand_forAllUVLambda
 * @param numberOfUV
 * @param numberOfWavelengths
 * @param samplingDistance
 * @param Ip
 */
void numericalIntegration_onLambda_non_optimised(int p, int colorChannel, double factN, double factM, vector<float *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, float samplingDistance, Mat Ip[]);


/**
 * Function to compute the factorial of n. No overflow for n<100 (at least).
 * @brief factorial
 * @param n
 * @return
 */
double factorial(double n);

/**
 * Given a 2D discrete function described as a matrix, returns function2D^power
 * @brief pow2D
 * @param function2D
 * @param power
 * @return
 */
Mat pow2D(Mat& function2D, unsigned int power);

/**
 * Returns the sign of the value x.
 * @brief sign
 * @param x
 * @return
 */
float sign(float x);

void preComputeGaussian(int numberOfUV, int size, float lambdaMin,
              int numberOfWavelength, float samplingDistanceWavelength, vector<float *> &result);


QVector3D multiplyQMatrixVector3x3(QMatrix3x3 matrix, QVector3D vector);

/**
 * Sets the norm of the vector to 1
 * @brief normalizeVector
 * @param vector
 */
void normalizeVector(Mat &vector);

/**
 * Sets the norm of the vector to 1.
 * @brief normalizeVector
 * @param x
 * @param y
 * @param z
 */
void normalizeVector(float &x, float &y, float &z);

Mat makeRotationMatrix(Point3f axis, float sin, float cos);

/**
 * Return the value of the Gauss error function (erf) with a standard deviation sigma.
 * @brief erf
 * @param x
 * @param sigma
 * @return
 */
float erfSigma(float x, float sigma);

/**
 * Return the value of the Gauss error function (erf) with a standard deviation sigma.
 * @brief erfSigma
 * @param x
 * @param sigma
 * @return
 */
double erfSigma(double x, double sigma);
#endif // MATHFUNCTIONS

