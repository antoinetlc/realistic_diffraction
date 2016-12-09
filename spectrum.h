/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <QVector>
#include <QVector3D>
#include <QColor>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

/*---- Qt ----*/
#include <QApplication>

using namespace std;

class Spectrum
{
    public:
        Spectrum();
        Spectrum(string spectrumName);
        ~Spectrum();

        void loadSPectrum();
        void normaliseSPD();

        float getLambdaMin() const;
        float getLambdaMax() const;
        float getSamplingDistance() const;
        int getNumberOfWavelength() const;
        QColor getRGBColor() const;
        float getSpectralPowerDistribution(int index);

    private:

        string m_spectrumName;
        float m_lambdaMin;
        float m_lambdaMax;
        float m_samplingDistance;
        int m_numberOfWavelengths;
        QColor m_RGBColor;
        float *m_spectralPowerDistribution;

};

#endif // SPECTRUM_H
