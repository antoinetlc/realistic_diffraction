/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/
#include <QApplication>

#include <iostream>
#include <math.h>
#include <string.h>

#include "huygens.h"

using namespace std;
using namespace cv;


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

   float widthObjectCm= 7.46;
   float heightObjectCm = 14.63;

   //Green
   float lambdaMeasurement = 0.53;
   float colorChannel = 1;

   float L = 300.0;
   float widthObjectPx = 1326.0;
   float heightObjectPx = 2393.0;
   Point2f center(1105, 943);

   string path = string("/PhoneLG/lg_pattern");

   vector<float> spectralPowerDistribution(81, 1.0);

   Mat diffractionPattern = loadPFM(qApp->applicationDirPath().toStdString() + path + ".pfm");
   computeDiffractionTable_Measurement2D(diffractionPattern, center, 0.04, widthObjectCm, heightObjectCm, widthObjectPx, heightObjectPx,
                                         lambdaMeasurement, colorChannel, L, spectralPowerDistribution, 1024, 1024, 5.0, 81);

    return 0;
}

