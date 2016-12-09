/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#ifndef COLOURSYSTEM_H
#define COLOURSYSTEM_H

#include "mathfunctions.h"

#include <iostream>
#include <string>
#include <cmath>

#include <QVector>
#include <QVector3D>

#include <opencv2/core/core.hpp>
#include <opencv/highgui.h>
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
using namespace cv;
/**
 * @brief The colourSystem class
 * Values for the colour system and the CIE matching functions taken from http://www.fourmilab.ch/documents/specrend/
 * sRGB : http://en.wikipedia.org/wiki/SRGB
 * XYZ to RGB matrix computation : http://www.brucelindbloom.com/index.html?WorkingSpaceInfo.html
 */
class ColourSystem
{
    public:
        ColourSystem(string colourSystemName);
        ~ColourSystem();

        /**
         * Compute the RGB color of each wavelength of the spectrum
         * Assume that the spectrum of a wavelength is a Dirac function (1.0 if the wavelength 0.0 otherwise)
         * @brief computeWavelengthRGB
         * @param numberOfWavelengths
         * @param lambdaMin
         * @param rgbColors
         */
        void computeWavelengthRGB(int numberOfWavelengths, float lambdaMin, QVector<QVector3D>& rgbColors);

        /**
         * Function used for debugging. Display an image 401*512 of the 401 wavelengths between 380 and 780 nm.
         * @brief displaySpectrum
         * @param rgbColors
         */
        void displaySpectrum(QVector<QVector3D>& rgbColors);

        /**
         * Method to convert a XYZ matrix to RGB in the sRGB system.
         * In matrix the encoding is XYZ = RGB.
         * @brief ColourSystem::XYZToRGB
         * @param matrix
         * @return
         */
        Mat XYZToRGB(const Mat &matrix);

        /**
         * Method to convert a RGB matrix to XYZ in the sRGB system.
         * Assume linear RGB (no gamma correction).
         * In matrix the encoding is XYZ = RGB.
         * @brief RGBToXYZ
         * @param matrix
         * @return
         */
        Mat RGBToXYZ(const Mat &matrix);

        /**
         * Method to convert a XYZ matrix to Lab in the sRGB system.
         * In matrix the encoding is Lab = RGB.
         * @brief ColourSystem::XYZToLAB
         * @param matrix
         * @return
         */
        Mat XYZToLAB(const Mat &matrix);

    private:
        string m_colourSystemName;
        float m_xRed;
        float m_yRed;
        float m_xBlue;
        float m_yBlue;
        float m_xGreen;
        float m_yGreen;
        float m_xWhitePoint;
        float m_yWhitePoint;
        float m_zWhitePoint;
        float m_gamma;
        float m_RGBToXYZMatrix[3][3];
};




#endif // COLOURSYSTEM_H
