/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#include "huygens.h"

void computeDiffractionTable_Measurement2D(Mat const &diffractionPattern, Point2f const &center, float F0,
                                           float widthObjectCm, float heightObjectCm, float widthObjectPx, float heightObjectPx,
                                           float lambdaMeasurement, int colorChannel,
                                           float distanceLightSource,
                                           string spectrumName, int widthTable, int heightTable, int numberOfWavelengths, float power)
{
    /*--- Initialisation ---*/
    //integrand_forAllUVLambdaX/Y/Z contains the value of the integrand for all (u,v,lambda).
    //At the end it will be integrated over lambda to get the final lookup table
    vector<float*> integrand_forAllUVLambdaX(numberOfWavelengths);
    vector<float*> integrand_forAllUVLambdaY(numberOfWavelengths);
    vector<float*> integrand_forAllUVLambdaZ(numberOfWavelengths);

    //For each lambda, we have a widthTable*heightTable table
    for(int k = 0 ; k<numberOfWavelengths ; k++)
    {
        integrand_forAllUVLambdaX[k] = new float[widthTable*heightTable];
        integrand_forAllUVLambdaY[k] = new float[widthTable*heightTable];
        integrand_forAllUVLambdaZ[k] = new float[widthTable*heightTable];
    }

    //Load the spectrum
    Spectrum spectrum(spectrumName);

    //Compute the extreme values of sines in the X and Y directions (Sines at each corner of the image)
    //This will be used later to avoid computing the table outside the range of the measured diffraction pattern
    int widthPattern = diffractionPattern.cols;
    int heightPattern = diffractionPattern.rows;

    float x = 0.0;
    float y = 0.0;

    //Top Left
    x = (-center.x)*widthObjectCm/widthObjectPx;
    y = (-center.y)*heightObjectCm/heightObjectPx;
    float sinXTopLeft = x/(sqrt(distanceLightSource*distanceLightSource+x*x));
    float sinYTopLeft = y/(sqrt(distanceLightSource*distanceLightSource+y*y));

    //Bottom Left
    x = (-center.x)*widthObjectCm/widthObjectPx;
    y = (heightPattern-1-center.y)*heightObjectCm/heightObjectPx;
    float sinXBottomLeft = x/(sqrt(distanceLightSource*distanceLightSource+x*x));
    float sinYBottomLeft = y/(sqrt(distanceLightSource*distanceLightSource+y*y));

    //Bottom Right
    x = (widthPattern-1-center.x)*widthObjectCm/widthObjectPx;
    y = (heightPattern-1-center.y)*heightObjectCm/heightObjectPx;
    float sinXBottomRight = x/(sqrt(distanceLightSource*distanceLightSource+x*x));
    float sinYBottomRight = y/(sqrt(distanceLightSource*distanceLightSource+y*y));

    //Top Right
    x = (widthPattern-1-center.x)*widthObjectCm/widthObjectPx;
    y = (-center.y)*heightObjectCm/heightObjectPx;
    float sinXTopRight = x/(sqrt(distanceLightSource*distanceLightSource+x*x));
    float sinYTopRight = y/(sqrt(distanceLightSource*distanceLightSource+y*y));

    /*--- Compute the look up table ---*/
    Mat result = Mat::zeros(widthTable,heightTable,CV_32FC3);
    float halfSizeX = widthTable/2.0;
    float halfSizeY = heightTable/2.0;

    //For each values of v and u
    for(int v = 0 ; v<widthTable ; v++)
    {
        #pragma omp parallel for num_threads(omp_get_max_threads())
        for(int u = 0 ; u<heightTable ; u++)
        {
            //For each wavelength
            for(int w = 0 ; w<numberOfWavelengths ; w++)
            {
                //Current wavelength
                float currentLambda = 0.38+0.005*w;

                //Current sines
                float  sinX = 2.0*pow((v-halfSizeX)/(widthTable), power);
                float  sinY = 2.0*pow((u-halfSizeY)/(heightTable), power);

                //The sines are scaled depending on the wavelength (Equation 15 in the paper)
                //For currentLambda = lambdaMeasurement the value of the current intensity is directly given by the measured diffraction pattern
                //For currentLambda != lambdaMeasurement the value of the current intensity is slightly offset by lambdaMeasurement/currentLambda
                sinX *= lambdaMeasurement/currentLambda;
                sinY *= lambdaMeasurement/currentLambda;

                float intensity = 0.0;

                //Do not look for a value of sin_x, sin_y if the current value is not in the range of the diffraction pattern picture
                bool isOutsideRange = (sinXBottomLeft>sinX && sinYBottomLeft>sinY)
                                        || (sinXBottomRight<sinX && sinYBottomRight>sinY)
                                        || (sinXTopLeft>sinX && sinYTopLeft<sinY)
                                        || (sinXTopRight<sinX && sinYTopRight<sinY);

                //Do not look for a value of sin_x, sin_y if the current value is not in the range of the diffraction pattern picture
                if(!isOutsideRange)
                {

                    //The following code finds the correct index ()
                    int kCurrentSinY = 0;
                    int lCurrentSinX = 0;

                    //currentX and current Y must be in the same quadrant as sinX and sinY.
                    float currentX = 0.0;
                    float currentY = 0.0;

                    //Find the X and Y coordinates in centimeters
                    if(fabs(sinX)>EPSILON && fabs(sinY)>EPSILON)
                    {
                        currentX = distanceLightSource*sinX/(sqrt(1-sinX*sinX));
                        currentY = distanceLightSource*sinY/(sqrt(1-sinY*sinY));
                    }
                    else if(fabs(sinX)>EPSILON && fabs(sinY)<EPSILON)
                    {
                        currentX = distanceLightSource*sinX/(sqrt(1-sinX*sinX));
                        currentY = 0.0;
                    }
                    else if(fabs(sinX)<EPSILON && fabs(sinY)>EPSILON)
                    {
                        currentX =0.0;
                        currentY = distanceLightSource*sinY/(sqrt(1-sinY*sinY));
                    }
                    else if(fabs(sinX)<EPSILON && fabs(sinY)<EPSILON)
                    {
                        currentX = 0.0;
                        currentY = 0.0;
                    }

                    //Convert back to pixels with indices in the correct range of the image
                    currentX *= widthObjectPx/widthObjectCm;
                    currentY *= heightObjectPx/heightObjectCm;

                    currentX += center.x;
                    currentY += center.y;

                    lCurrentSinX = floor(currentX);
                    kCurrentSinY = floor(currentY);

                    //Make sure we are inside the image
                    if(kCurrentSinY>=0 && lCurrentSinX>=0 && kCurrentSinY<diffractionPattern.rows && lCurrentSinX<diffractionPattern.cols)
                        intensity = diffractionPattern.at<Vec3f>(kCurrentSinY,lCurrentSinX).val[colorChannel];
                    else
                        intensity = 0.0;
                }
                else //else do not look for values
                {
                     intensity = 0.0;
                }

                float spectralPower =  spectrum.getSpectralPowerDistribution(w);

                //Normalisation factor due to the measurement at normal incidence (equation 15 paper)
                float normalizationFactor = pow(lambdaMeasurement/currentLambda,2.0)/(4.0*M_PI*M_PI*F0*F0);

                //Store the 3 values of the integrand
                integrand_forAllUVLambdaX[w][u*widthTable+v] = normalizationFactor*intensity*spectralPower*cie_colour_matching_function[w][0];
                integrand_forAllUVLambdaY[w][u*widthTable+v] = normalizationFactor*intensity*spectralPower*cie_colour_matching_function[w][1];
                integrand_forAllUVLambdaZ[w][u*widthTable+v] = normalizationFactor*intensity*spectralPower*cie_colour_matching_function[w][2];
             }//End for wavelengths
            }//End for v
        }//End for u

    //Integrate using the Trapezoidal rule
    numericalIntegration_onLambda(2, integrand_forAllUVLambdaX, widthTable, numberOfWavelengths, 0.005, result);
    numericalIntegration_onLambda(1, integrand_forAllUVLambdaY, widthTable, numberOfWavelengths, 0.005, result);
    numericalIntegration_onLambda(0, integrand_forAllUVLambdaZ, widthTable, numberOfWavelengths, 0.005, result);

    //Convert XYZ to RGB
    ColourSystem cs("sRGB");
    Mat RGB = cs.XYZToRGB(result);
    savePFM(RGB, qApp->applicationDirPath().toStdString()+ "diffractionTable.pfm");

    //Free the memory
    for(int k = 0 ; k<numberOfWavelengths ; k++)
    {
        delete[] integrand_forAllUVLambdaX[k];
        delete[] integrand_forAllUVLambdaY[k];
        delete[] integrand_forAllUVLambdaZ[k];
    }
}
