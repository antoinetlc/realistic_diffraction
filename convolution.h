/**
  * Practical acquisition and rendering of diffraction effects in surface reflectance.
  * Antoine Toisoul, Abhijeet Ghosh.
  * Imperial College London, December 2016.
  */
#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "PFMReadWrite.h"
#include "coordinates.h"

#include <QCoreApplication>
#include <QVector3D>

#include <iostream>
#include <string>

#include <omp.h>

#include <opencv2/core/core.hpp>

/**
 * Convolves the environment map with the diffraction pattern assuming normal = viewing vector = reflection vector.
 * The diffraction pattern lookup table must be a RGB lookup table linearly sampled.
 * The x and y ranges that the lookup table spans must be given as inputs.
 * @brief convolveEnvironmentMap
 * @param environmentMap
 * @param diffractionPattern
 * @param xRange
 * @param yRange
 */
cv::Mat convolveEnvironmentMap(cv::Mat &environmentMap, cv::Mat &diffractionPattern, float xRange, float yRange);

#endif // CONVOLUTION_H
