/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#include "coloursystem.h"

ColourSystem::ColourSystem(string colourSystemName)
{

    //Load the chromaticities depending on the system
    if(colourSystemName == "NTSC")
    {
        m_colourSystemName = string("NTSC");
        m_xRed = (float) 0.67;
        m_yRed = (float) 0.33;
        m_xGreen = (float) 0.21;
        m_yGreen = (float) 0.71;
        m_xBlue = (float) 0.14;
        m_yBlue = (float) 0.08;
        m_xWhitePoint = (float) 0.3101;
        m_yWhitePoint = (float) 0.3162;
        m_gamma = 0;
    }
    else if(colourSystemName == "EBU (PAL/SECAM)")
    {
        m_colourSystemName = string("EBU (PAL/SECAM)");
        m_xRed = (float) 0.64;
        m_yRed = (float) 0.33;
        m_xGreen = (float) 0.29;
        m_yGreen = (float) 0.60;
        m_xBlue = (float) 0.15;
        m_yBlue = (float) 0.06;
        m_xWhitePoint = (float) 0.3127;
        m_yWhitePoint = (float) 0.3291;
        m_gamma = 0;
    }
    else if(colourSystemName == "SMPTE")
    {
        m_colourSystemName = string("SMPTE");
        m_xRed = (float) 0.630;
        m_yRed = (float) 0.340;
        m_xGreen = (float) 0.310;
        m_yGreen = (float) 0.595;
        m_xBlue = (float) 0.155;
        m_yBlue = (float) 0.070;
        m_xWhitePoint = (float) 0.3127;
        m_yWhitePoint = (float) 0.3291;
        m_gamma = 0;
    }
    else if(colourSystemName == "HDTV")
    {
      m_colourSystemName = string("HDTV");
      m_xRed = (float) 0.670;
      m_yRed = (float) 0.330;
      m_xGreen = (float) 0.210;
      m_yGreen = (float) 0.710;
      m_xBlue = (float) 0.150;
      m_yBlue = (float) 0.060;
      m_xWhitePoint = (float) 0.3127;
      m_yWhitePoint = (float) 0.3291;
      m_gamma = 0;
    }
    else if(colourSystemName == "CIE")
    {
        m_colourSystemName = string("CIE");
        m_xRed = (float) 0.7355;
        m_yRed = (float) 0.2645;
        m_xGreen = (float) 0.2658;
        m_yGreen = (float) 0.7243;
        m_xBlue = (float) 0.1669;
        m_yBlue = (float) 0.0085;
        m_xWhitePoint = (float) 0.33333333;
        m_yWhitePoint = (float) 0.33333333;
        m_gamma = 0;
    }
    else if(colourSystemName == "CIE REC 709")
    {
        m_colourSystemName = string("CIE REC 709");
        m_xRed = (float) 0.64;
        m_yRed = (float) 0.33;
        m_xGreen = (float) 0.30;
        m_yGreen = (float) 0.60;
        m_xBlue = (float) 0.15;
        m_yBlue = (float) 0.06;
        m_xWhitePoint = (float) 0.33333333;
        m_yWhitePoint = (float) 0.33333333;
        m_gamma = 0;
    }
    else if(colourSystemName == "sRGB")
    {
        m_colourSystemName = string("sRGB");
        m_xRed = (float) 0.64;
        m_yRed = (float) 0.33;
        m_xGreen = (float) 0.30;
        m_yGreen = (float) 0.60;
        m_xBlue = (float) 0.15;
        m_yBlue = (float) 0.06;
        m_xWhitePoint = (float) 0.3127;
        m_yWhitePoint = (float) 0.3290;
        m_zWhitePoint = (float) 0.3583;
        m_gamma = 0;
    }

    /*--------- Computation of the RGB to XYZ matrix --------*/
    //YWhite assumed to be 1.0
    float Sr = 0.0, Sg = 0.0, Sb =0.0;

    float YWhite = 1.0;

    float whitePointMatrix[3][3] ={{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    float whitePointMatrixInverted[3][3] ={{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};

    //xyY to XYZ for the white point
    float XWhite = YWhite*m_xWhitePoint/m_yWhitePoint;
    float ZWhite = YWhite*(1.0-m_xWhitePoint-m_yWhitePoint)/m_yWhitePoint;

    whitePointMatrix[0][0] = m_xRed/m_yRed;            whitePointMatrix[0][1] = m_xGreen/m_yGreen;               whitePointMatrix[0][2] =  m_xBlue/m_yBlue;
    whitePointMatrix[1][0] = 1.0;                      whitePointMatrix[1][1] = 1.0;                             whitePointMatrix[1][2] = 1.0;
    whitePointMatrix[2][0] = (1-m_xRed-m_yRed)/m_yRed; whitePointMatrix[2][1] = (1-m_xGreen-m_yGreen)/m_yGreen;  whitePointMatrix[2][2] = (1-m_xBlue-m_yBlue)/m_yBlue;

    invertMatrix3x3(whitePointMatrix, whitePointMatrixInverted);

    //Calculate Sr, Sg, Sb necessary for the calculation of RGB to XYZ matrix
    Sr = whitePointMatrixInverted[0][0]*XWhite+whitePointMatrixInverted[0][1]*YWhite+whitePointMatrixInverted[0][2]*ZWhite;
    Sg = whitePointMatrixInverted[1][0]*XWhite+whitePointMatrixInverted[1][1]*YWhite+whitePointMatrixInverted[1][2]*ZWhite;
    Sb = whitePointMatrixInverted[2][0]*XWhite+whitePointMatrixInverted[2][1]*YWhite+whitePointMatrixInverted[2][2]*ZWhite;

    m_RGBToXYZMatrix[0][0] = Sr*whitePointMatrix[0][0];     m_RGBToXYZMatrix[0][1] = Sg*whitePointMatrix[0][1];     m_RGBToXYZMatrix[0][2] = Sb*whitePointMatrix[0][2];
    m_RGBToXYZMatrix[1][0] = Sr*whitePointMatrix[1][0];     m_RGBToXYZMatrix[1][1] = Sg*whitePointMatrix[1][1];     m_RGBToXYZMatrix[1][2] = Sb*whitePointMatrix[1][2];
    m_RGBToXYZMatrix[2][0] = Sr*whitePointMatrix[2][0];     m_RGBToXYZMatrix[2][1] = Sg*whitePointMatrix[2][1];     m_RGBToXYZMatrix[2][2] = Sb*whitePointMatrix[2][2];
}

ColourSystem::~ColourSystem()
{

}

void ColourSystem::computeWavelengthRGB(int numberOfWavelengths, float lambdaMin, QVector<QVector3D>& rgbColors)
{
    //For each avelength the spectrum is assumed to be a Dirac function (energy 1.0 at the wavelength and 0.0 anywhere else)
    float currentWavelength = 0.0;
    int colorMatchingFunctionIndex = 0;

    //xyz_bar : valuesof the colour matching function for a specific wavelength
    float x_bar = 0.0, y_bar = 0.0, z_bar = 0.0;
    float X = 0.0, Y = 0.0, Z = 0.0;    //CIE XYZ value for the current wavelength

    //Compute the XYZ to RGB matrix
    float XYZToRGBmatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    invertMatrix3x3(m_RGBToXYZMatrix, XYZToRGBmatrix);

    //For each wavelength compute the value of X,Y and Z
    QVector<QVector3D> XYZColours;
    for(int i = 0 ; i< numberOfWavelengths ; i++)
    {
        //Wavelength for which X,Y and Z are computed
        currentWavelength = lambdaMin + (float) i;

        //Find the values of the color matching function for that wavelength
        //The values of the colour matching function are given every 5 nanometers
        //In between a linear interpolation is used y = (1-alpha)y1 + alpha*y2

        float alpha = (currentWavelength-lambdaMin)/5.0-floor((currentWavelength-lambdaMin)/5.0);
        colorMatchingFunctionIndex = floor((currentWavelength-lambdaMin)/5.0);

        if(colorMatchingFunctionIndex<80)
        {
            x_bar = (1.0-alpha)*cie_colour_matching_function[colorMatchingFunctionIndex][0]+alpha*cie_colour_matching_function[colorMatchingFunctionIndex+1][0];
            y_bar =(1.0-alpha)*cie_colour_matching_function[colorMatchingFunctionIndex][1]+alpha*cie_colour_matching_function[colorMatchingFunctionIndex+1][1];
            z_bar = (1.0-alpha)*cie_colour_matching_function[colorMatchingFunctionIndex][2]+alpha*cie_colour_matching_function[colorMatchingFunctionIndex+1][2];
        }
        else
        {
            x_bar = cie_colour_matching_function[colorMatchingFunctionIndex][0];
            y_bar = cie_colour_matching_function[colorMatchingFunctionIndex][1];
            z_bar = cie_colour_matching_function[colorMatchingFunctionIndex][2];
        }

        //Compute the CIE XYZ (= matching funtions here due to specific spectrum)
        X = x_bar;
        Y = y_bar;
        Z = z_bar;

        QVector3D currentColor = QVector3D(X,Y,Z);
        XYZColours.push_back(currentColor);
    }


    //The colours have been computed in XYZ space
    //Now set the luminance of each component to the range [0 ; 1] for a correct display
    float maximumValueY = 0.0;
    for(int i = 0 ; i< numberOfWavelengths ; i++)
    {
        if(XYZColours[i].y()>maximumValueY)
        {
            maximumValueY = XYZColours[i].y();
        }
    }

    for(int i = 0 ; i< numberOfWavelengths ; i++)
    {
       XYZColours[i].setX(XYZColours[i].x()/maximumValueY);
       XYZColours[i].setY(XYZColours[i].y()/maximumValueY);
       XYZColours[i].setZ(XYZColours[i].z()/maximumValueY);
    }

    //Now calculate the RGB values
    //Calculate the RGB values for the wavelength from the XYZ
    float R = 0.0, G = 0.0, B = 0.0; //RGB for the current wavelength

    for(int i = 0 ; i< numberOfWavelengths ; i++)
    {
        X = XYZColours[i].x();
        Y = XYZColours[i].y();
        Z = XYZColours[i].z();

        R = XYZToRGBmatrix[0][0]*X+XYZToRGBmatrix[0][1]*Y+XYZToRGBmatrix[0][2]*Z;
        G = XYZToRGBmatrix[1][0]*Y+XYZToRGBmatrix[1][1]*Y+XYZToRGBmatrix[1][2]*Z;
        B = XYZToRGBmatrix[2][0]*Z+XYZToRGBmatrix[2][1]*Y+XYZToRGBmatrix[2][2]*Z;

        //R,G,B can be negative (outside the RGB gamut)
        //Add a positive weight such as the three components are positive
        float weight = 0.0;

        if(R<0 || G<0 || B<0)
        {
            weight = -min(R,min(G,B)); //We set the most negative one to zero and add its opposite value to the others
        }

        R += weight;
        G += weight;
        B += weight;

        QVector3D currentColor = QVector3D(R,G,B);
        rgbColors.push_back(currentColor);
    }


    /*-------Apply a gamma correction-------*/
    //2.2 gamma correction
    for(int k = 0 ; k<rgbColors.size() ; k++)
    {
       R = rgbColors[k].x();
       G = rgbColors[k].y();
       B = rgbColors[k].z();

       R = pow(R, 1.0/2.2);
       G = pow(G, 1.0/2.2);
       B = pow(B, 1.0/2.2);

       rgbColors[k].setX(R);
       rgbColors[k].setY(G);
       rgbColors[k].setZ(B);
    }

}

void ColourSystem::displaySpectrum(QVector<QVector3D>& rgbColors)
{
    Mat colors = Mat::zeros(Size(401,512), CV_32FC3);
    Mat channel[3];

    split(colors, channel);

    //OpenCv uses BGR
    for(int i = 0 ; i<512 ; i++)
    {
        for(int j = 0 ; j<401 ; j++)
        {
            channel[0].at<float>(i,j) = 255.0*rgbColors[j].z();
            channel[1].at<float>(i,j) = 255.0*rgbColors[j].y();
            channel[2].at<float>(i,j) = 255.0*rgbColors[j].x();
        }
    }

    merge(channel, 3 , colors);
    colors.convertTo(colors, CV_8UC3);

    imshow("Colors", colors);
    waitKey(0);

}

Mat ColourSystem::XYZToRGB(const Mat &matrix)
{  //Computation for the sRGB system
    //float XYZToRGBmatrix[3][3] = {{3.1338561, -1.6168667, -0.4906146}, {-0.9787684, 1.9161415, 0.0334540}, {0.0719453, -0.2289914, 1.4052427}};
    float XYZToRGBmatrix[3][3] = {{3.2404542, -1.5371385, -0.4985314}, {-0.9692660, 1.8760108, 0.0415560}, {0.0556434 , -0.2040259, 1.0572252}};

    Mat result = Mat(matrix.rows, matrix.cols, CV_32FC3);
    Mat matrixFloat = matrix.clone();
    matrixFloat.convertTo(matrixFloat, CV_32FC3);

    float X = 0.0, Y = 0.0, Z = 0.0;
    float R = 0.0, G = 0.0, B = 0.0;

    float weight = 0.0;
    float currentValue = 0.0, maximumValue = 0.0;

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


      /*       if(R<0 || G<0 || B<0)
            {
                weight = -min(R,min(G,B)); //We set the most negative one to zero and add its opposite value to the others
            }

            R += weight;
            G += weight;
            B += weight;

            matrixFloat.at<Vec3f>(i,j).val[2] = R;
            matrixFloat.at<Vec3f>(i,j).val[1] = G;
            matrixFloat.at<Vec3f>(i,j).val[0] = B;

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

Mat ColourSystem::RGBToXYZ(const Mat &matrix)
{
    //Computation for the sRGB system
    //Assume that the gamma has already been removed : linear RGB
    float RGBToXYZMatrix[3][3] = {{0.4360747, 0.3850649, 0.1430804}, {0.2225045, 0.7168786, 0.0606169}, {0.0139322 , 0.0971045, 0.7141733}};

    Mat result = Mat(matrix.rows, matrix.cols, CV_32FC3);

    float X = 0.0, Y = 0.0, Z = 0.0;
    float R = 0.0, G = 0.0, B = 0.0;

    float weight = 0.0;
    float currentValue = 0.0, maximumValue = 0.0;

    for(int i = 0 ; i< matrix.rows ; i++)
    {
        for(int j = 0 ; j< matrix.cols ; j++)
        {
            //XYZ = RGB (BGR in openCV)
            R = matrix.at<Vec3f>(i,j).val[2];
            G = matrix.at<Vec3f>(i,j).val[1];
            B = matrix.at<Vec3f>(i,j).val[0];

            //Calculate the RGB values for the wavelength from the XYZ
            X = RGBToXYZMatrix[0][0]*R+RGBToXYZMatrix[0][1]*G+RGBToXYZMatrix[0][2]*B;
            Y = RGBToXYZMatrix[1][0]*R+RGBToXYZMatrix[1][1]*G+RGBToXYZMatrix[1][2]*B;
            Z = RGBToXYZMatrix[2][0]*R+RGBToXYZMatrix[2][1]*G+RGBToXYZMatrix[2][2]*B;

            //Encoding : RGB = XYZ
            result.at<Vec3f>(i,j).val[2] = X;
            result.at<Vec3f>(i,j).val[1] = Y;
            result.at<Vec3f>(i,j).val[0] = Z;
        }
    }

    return result;
}

Mat ColourSystem::XYZToLAB(const Mat &matrix)
{
    //Computation for the sRGB system
    //CIE standards
    float epsilon = 0.008856;
    float kappa = 903.3;

    float XWhite = m_xWhitePoint/m_yWhitePoint;
    float YWhite = 1.0;
    float ZWhite = (1.0-m_xWhitePoint-m_yWhitePoint)/m_yWhitePoint;

    Mat result = Mat(matrix.rows, matrix.cols, CV_32FC3);

    for(int i = 0 ; i< matrix.rows ; i++)
    {
        for(int j = 0 ; j< matrix.cols ; j++)
        {

            float X = matrix.at<Vec3f>(i,j).val[2];
            float Y = matrix.at<Vec3f>(i,j).val[1];
            float Z = matrix.at<Vec3f>(i,j).val[0];

            float xr = X/XWhite;
            float yr = Y/YWhite;
            float zr = Z/ZWhite;

            float fx = 0.0, fy = 0.0, fz = 0.0;
            float L = 0.0, a = 0.0, b = 0.0;

            if(xr>epsilon)
            {
                fx = pow(xr, 1.0/3.0);
            }
            else
            {
                fx = (kappa*xr+16.0)/116.0;
            }

            if(yr>epsilon)
            {
                fy = pow(yr, 1.0/3.0);
            }
            else
            {
                fy = (kappa*yr+16.0)/116.0;
            }

            if(zr>epsilon)
            {
                fz = pow(zr, 1.0/3.0);
            }
            else
            {
                fz = (kappa*zr+16.0)/116.0;
            }

            L = 116.0*fy-16.0;
            a = 500.0*(fx-fy);
            b = 200.0*(fy-fz);

            //Lab = RGB (BGR in openCV)
            result.at<Vec3f>(i,j).val[2] = L;
            result.at<Vec3f>(i,j).val[1] = a;
            result.at<Vec3f>(i,j).val[0] = b;
        }
    }

    return result;
}
