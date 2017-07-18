/**
  * Practical acquisition and rendering of diffraction effects in surface reflectance.
  * Antoine Toisoul, Abhijeet Ghosh.
  * Imperial College London, December 2016.
  */

#ifndef COORDINATES_H
#define COORDINATES_H

#define M_PI 3.14159265358979323846

#include <iostream>
#include <cmath>
#include <vector>

/*---- OpenCV ----*/
#include <opencv2/core/core.hpp>

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
 * Conversion from spherical to cartesian coordinates using the y axis as the vertical axis.
 * @brief sphericalToCartesian
 * @param r
 * @param theta
 * @param phi
 * @param x
 * @param y
 * @param z
 */
void sphericalToCartesian(float r, float theta, float phi, float &x, float &y, float &z);

#endif // COORDINATES_H

