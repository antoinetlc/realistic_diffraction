/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#include "spectrum.h"

Spectrum::Spectrum(): m_spectrumName(""), m_lambdaMin(0.38), m_lambdaMax(0.78), m_samplingDistance(0.005), m_numberOfWavelengths(81),
    m_RGBColor(QColor())
{

}

Spectrum::Spectrum(string spectrumName): m_spectrumName(spectrumName), m_lambdaMin(0.38), m_lambdaMax(0.78), m_samplingDistance(0.005), m_numberOfWavelengths(81),
    m_RGBColor(QColor())
{

    m_spectralPowerDistribution = new float[m_numberOfWavelengths];


    this->loadSPectrum();

}

Spectrum::~Spectrum()
{
    delete[] m_spectralPowerDistribution;
}

void Spectrum::loadSPectrum()
{

    ifstream SPDFile;

    if(m_spectrumName == "fullWhite")
    {
        for(int i = 0 ; i<m_numberOfWavelengths ; i++)
        {
           m_spectralPowerDistribution[i] = 1.0;
        }
    }
    else if(m_spectrumName == "incandescent")
    {
        SPDFile = ifstream(qApp->applicationDirPath().toStdString() + "/Data/illuminantA.txt", ios::in);

        /*--Read SPD values---*/
        string wavelength, SPDValue;

        for(int k = 0 ; k<m_numberOfWavelengths ; k++)
        {
            //First line is parallel and second is cross polarised
            SPDFile >> wavelength >> SPDValue;
            m_spectralPowerDistribution[k] = atof(SPDValue.c_str());
        }
    }
    else if(m_spectrumName == "incandescentBlue")
    {
        SPDFile = ifstream(qApp->applicationDirPath().toStdString() + "/Data/illuminantA.txt", ios::in);

        /*--Read SPD values---*/
        string wavelength, SPDValue;

        for(int k = 0 ; k<m_numberOfWavelengths ; k++)
        {
            //First line is parallel and second is cross polarised
            SPDFile >> wavelength >> SPDValue;
            m_spectralPowerDistribution[m_numberOfWavelengths-1-k] = atof(SPDValue.c_str());
        }
    }
    else if(m_spectrumName == "F1")
    {
        SPDFile = ifstream(qApp->applicationDirPath().toStdString() + "/Data/illuminantF1.txt", ios::in);

        /*--Read SPD values---*/
        string wavelength, SPDValue;

        for(int k = 0 ; k<m_numberOfWavelengths ; k++)
        {
            //First line is parallel and second is cross polarised
            SPDFile >> wavelength >> SPDValue;
            m_spectralPowerDistribution[k] = atof(SPDValue.c_str());
            cout << m_spectralPowerDistribution[k] << endl;
        }
    }
    else if(m_spectrumName == "custom")
    {
        SPDFile = ifstream(qApp->applicationDirPath().toStdString() + "/Data/custom.txt", ios::in);

        /*--Read SPD values---*/
        string wavelength, SPDValue;

        for(int k = 0 ; k<m_numberOfWavelengths ; k++)
        {
            //First line is parallel and second is cross polarised
            SPDFile >> wavelength >> SPDValue;
            m_spectralPowerDistribution[k] = atof(SPDValue.c_str());
            cout << m_spectralPowerDistribution[k] << endl;
        }
    }
    else
    {
            cout << m_spectrumName << " does not exist ! " << endl;
    }

    this->normaliseSPD();

}


void Spectrum::normaliseSPD()
{
    //Find Maximum of SPD
    float maximum = 0.0;

    for(int k = 0 ; k<m_numberOfWavelengths ; k++)
    {
       if(m_spectralPowerDistribution[k]>maximum)
       {
           maximum = m_spectralPowerDistribution[k];
       }
    }

    //Normalise in [0;1] range
    for(int k = 0 ; k<m_numberOfWavelengths ; k++)
    {
       m_spectralPowerDistribution[k] /= maximum;
    }
}

float Spectrum::getLambdaMin() const
{
    return m_lambdaMin;
}

float Spectrum::getLambdaMax() const
{
    return m_lambdaMax;
}

float Spectrum::getSamplingDistance() const
{
    return m_samplingDistance;
}

int Spectrum::getNumberOfWavelength() const
{
    return m_numberOfWavelengths;
}

QColor Spectrum::getRGBColor() const
{
    return m_RGBColor;
}


float Spectrum::getSpectralPowerDistribution(int index)
{
    return m_spectralPowerDistribution[index];
}
