/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#include "coloursystem.h"

ColourSystem::ColourSystem(string colourSystemName)
{
}

ColourSystem::~ColourSystem()
{

}

Mat ColourSystem::XYZToRGB(const Mat &matrix)
{  //Computation for the sRGB system
    //float XYZToRGBmatrix[3][3] = {{3.1338561, -1.6168667, -0.4906146}, {-0.9787684, 1.9161415, 0.0334540}, {0.0719453, -0.2289914, 1.4052427}};
    float XYZToRGBmatrix[3][3] = {{3.2404542, -1.5371385, -0.4985314}, {-0.9692660, 1.8760108, 0.0415560}, {0.0556434 , -0.2040259, 1.0572252}};

    Mat matrixFloat = matrix.clone();
    matrixFloat.convertTo(matrixFloat, CV_32FC3);

    float X = 0.0, Y = 0.0, Z = 0.0;
    float R = 0.0, G = 0.0, B = 0.0;

    float weight = 0.0;

    for(int i = 0 ; i< matrix.rows ; i++)
    {
        for(int j = 0 ; j< matrix.cols ; j++)
        {
            //XYZ = RGB (BGR in openCV)
            X = matrixFloat.at<Vec3f>(i,j).val[2];
            Y = matrixFloat.at<Vec3f>(i,j).val[1];
            Z = matrixFloat.at<Vec3f>(i,j).val[0];

            //Calculate the RGB values for the wavelength from the XYZ
            R = XYZToRGBmatrix[0][0]*X+XYZToRGBmatrix[0][1]*Y+XYZToRGBmatrix[0][2]*Z;
            G = XYZToRGBmatrix[1][0]*X+XYZToRGBmatrix[1][1]*Y+XYZToRGBmatrix[1][2]*Z;
            B = XYZToRGBmatrix[2][0]*X+XYZToRGBmatrix[2][1]*Y+XYZToRGBmatrix[2][2]*Z;


            matrixFloat.at<Vec3f>(i,j).val[2] = R;
            matrixFloat.at<Vec3f>(i,j).val[1] = G;
            matrixFloat.at<Vec3f>(i,j).val[0] = B;


             if(R<0 || G<0 || B<0)
            {
                weight = -min(R,min(G,B)); //We set the most negative one to zero and add its opposite value to the others
            }

            R += weight;
            G += weight;
            B += weight;

            matrixFloat.at<Vec3f>(i,j).val[2] = R;
            matrixFloat.at<Vec3f>(i,j).val[1] = G;
            matrixFloat.at<Vec3f>(i,j).val[0] = B;
/*
            currentValue = max(R, max(G, B));

            if(currentValue >maximumValue)
            {
                maximumValue = currentValue;
            }*/
        }
    }

   // matrixFloat /= maximumValue;

    //Gamma correction
   /* for(int i = 0 ; i< matrix.rows ; i++)
    {
        for(int j = 0 ; j< matrix.cols ; j++)
        {
            //XYZ = RGB (BGR in openCV)
            R = matrixFloat.at<Vec3f>(i,j).val[2];
            G = matrixFloat.at<Vec3f>(i,j).val[1];
            B = matrixFloat.at<Vec3f>(i,j).val[0];

            matrixFloat.at<Vec3f>(i,j).val[2] = pow(R, 1.0/2.2);
            matrixFloat.at<Vec3f>(i,j).val[1] = pow(G, 1.0/2.2);
            matrixFloat.at<Vec3f>(i,j).val[0] = pow(B, 1.0/2.2);
        }
    }
*/

    return matrixFloat;
}
