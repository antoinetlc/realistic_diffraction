/*******************************************************
 * Copyright (C) 2015 Antoine Toisoul <antoine.toisoul@telecom-paristech.org>
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Antoine Toisoul <antoine.toisoul@telecom-paristech.org>, November 2015
 *******************************************************/

#include "mathfunctions.h"

void invertMatrix3x3(float matrix[][3], float invertedMatrix[][3])
{
    float determinant = matrix[0][0]*matrix[1][1]*matrix[2][2]+
            matrix[1][0]*matrix[2][1]*matrix[0][2]+
            matrix[0][1]*matrix[1][2]*matrix[2][0]-
            matrix[0][2]*matrix[1][1]*matrix[2][0]-
            matrix[0][1]*matrix[1][0]*matrix[2][2]-
            matrix[1][2]*matrix[2][1]*matrix[0][0];

    if(determinant != 0.0)
    {
        invertedMatrix[0][0] = (matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2])/determinant;
        invertedMatrix[0][1] = (matrix[0][2]*matrix[2][1]-matrix[2][2]*matrix[0][1])/determinant;
        invertedMatrix[0][2] = (matrix[0][1]*matrix[1][2]-matrix[1][1]*matrix[0][2])/determinant;

        invertedMatrix[1][0] = (matrix[1][2]*matrix[2][0]-matrix[2][2]*matrix[1][0])/determinant;
        invertedMatrix[1][1] = (matrix[0][0]*matrix[2][2]-matrix[2][0]*matrix[0][2])/determinant;
        invertedMatrix[1][2] = (matrix[0][2]*matrix[1][0]-matrix[1][2]*matrix[0][0])/determinant;

        invertedMatrix[2][0] = (matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1])/determinant;
        invertedMatrix[2][1] = (matrix[0][1]*matrix[2][0]-matrix[2][1]*matrix[0][0])/determinant;
        invertedMatrix[2][2] = (matrix[0][0]*matrix[1][1]-matrix[1][0]*matrix[0][1])/determinant;
    }

}

void printMatrix3x3(float matrix[][3])
{
    for(int i = 0 ; i<3 ; i++)
    {
        for(int j = 0 ; j<3 ; j++)
        {
            cout << matrix[i][j] <<  "   ";
        }
        cout << endl;
    }
}

void cartesianToSpherical(const float &x, const float &y, const float &z, float &r, float &theta, float &phi)
{
    r = sqrt(x*x + y*y + z*z);
    //theta = M_PI-acos(y/r);
    theta = acos(y/r);
    //atan2 compute the the correct value for phi depending on the quadrant (x,z) are in
    //The value return is in the range -Pi ; Pi. Hence Pi is added for the result to be in the range 0 2Pi.
    phi = moduloRealNumber(atan2(x, z), (float) 2.0*M_PI); //Phi is between 0 and 2*PI
}


void cartesianToSphericalVector(const vector<vector<float> > &cartesianVector, vector<vector<float> > &sphericalVector)
{

    float r = 0.0;
    float theta = 0.0;
    float phi = 0.0;

    sphericalVector.clear();
    for(unsigned int i = 0 ; i<cartesianVector.size() ; i++)
    {
        vector<float> currentSpherical;
        cartesianToSpherical(cartesianVector[i][0],cartesianVector[i][1],cartesianVector[i][2],r,theta,phi);
        currentSpherical.push_back(r);
        currentSpherical.push_back(phi);
        currentSpherical.push_back(theta);
        sphericalVector.push_back(currentSpherical);
    }
}

float moduloRealNumber(float number, float modulo)
{
    int quotient = 0;
    float result = 0.0;

    quotient = (int) floor(number/modulo);
    result = number-quotient*modulo;

    return result;
}


float clamp(float value, float inf, float sup)
{
    float result = 0.0f;

    if(inf<=sup)
    {
        if(value < inf)
            result = inf;
        else if(value > sup)
            result = sup;
        else
            result = value;
    }
    else
    {
        cerr << "sup must be greater than inf !" << endl;
        return value;
    }

    return result;
}

double sinc(double value)
{
    double epsilon = 0.0000001;

    if(abs(value)< epsilon)
    {
        //Taylor series sin(x)/x = 1-x^2/6
        return 1-value*value/6.0;
    }
    else
    {
        return sin(value)/value;
    }
}

double sincModel(double x, const dlib::matrix<double,1,1> &parameter)
{   
    const double value = pow(sinc(parameter*x),2.0);

    return value;
}

double residual(const pair<double, double>& data, const dlib::matrix<double,1,1> &parameter)
{
   //  cout << data.first(0) << endl;
    return sincModel(data.first, parameter) - data.second;
}

double residual_derivative(const pair<double, double>& data, const dlib::matrix<double,1,1> &parameter)
{
    double der;


    const double x = data.first;

    der = 2*sinc(parameter*x)*(cos(parameter*x)-sinc(parameter*x))/parameter;

    return der;
}

double sincModel2(double x, const dlib::matrix<double,2,1> &parameter)
{
    const double value = parameter(0)*pow(sinc(parameter(1)*x),2.0);

    return value;
}

double residual2(const pair<double, double>& data, const dlib::matrix<double,2,1> &parameter)
{
    return sincModel2(data.first, parameter) - data.second;
}

dlib::matrix<double,2,1> residual_derivative2(const pair<double, double>& data, const dlib::matrix<double,2,1> &parameter)
{
    dlib::matrix<double,2,1> der;

    const double x = data.first;

    der(0) = pow(sinc(parameter(1)*x) ,2.0);
    der(1) = 2*parameter(0)*sinc(parameter(1)*x)*(cos(parameter(1)*x)-sinc(parameter(1)*x))/parameter(1);

    return der;
}

void leastSquaresOptimisation(double *x, double *data, int size, int offsetOrder0, double* result, dlib::matrix<double,1,1>& parameter)
{

    //The data used for the optimisation method is different
    //The central peak is removed so that the function is not overfitted to the central peak
    std::vector<std::pair<double, double> > data_samples;

    int center = size/2;

    for (int i=0; i<size; i++)
    {
        //Delete the central peak
        if(i>center+offsetOrder0 || i<center-offsetOrder0)
        {
            data_samples.push_back(std::make_pair(x[i], data[i]));
        }

    }


    //Least squares minimisation with a parameter of 1.0 to start
    parameter = 1.0;

    //The next line checks if the analytic expression of the derivative of the residual is correct
    //The value has to be very small
    cout << "Derivative error: " << length(residual_derivative(data_samples[0], parameter) - derivative(residual)(data_samples[0], parameter)) << endl;

    //Optimisation
    dlib::solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7).be_verbose(), residual, residual_derivative, data_samples, parameter);

    cout << "The infered parameter is : "<< trans(parameter) << endl;

    //Generate the fitted function
    for (int i=0; i<size; i++)
    {
        result[i] = pow(sinc(parameter(0)*x[i]),2.0);
    }

}


float numericalIntegration(float functionSamples[], int size, float samplingDistance)
{
    //Assumption (sup-inf)/samplingDistance) is an integer
    float integral = 0.0;

    //Size-1 trapezes
    for(int i = 0 ; i<(size-1) ; i++)
    {
        //Add the area of a trapeze
        integral += samplingDistance*(functionSamples[i+1]+functionSamples[i])/2.0;
    }

    return integral;
}

void numericalIntegration_onLambda(int colorChannel, vector<float *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult)
{
   // float currentValue = 0.0;
    //float nextValue = 0.0;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int u = 0 ; u<numberOfUV ; u++)
    {

        for(int v = 0 ; v<numberOfUV ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                float currentValue = integrand_forAllUVLambda[w][u*numberOfUV+v];
                float nextValue = integrand_forAllUVLambda[w+1][u*numberOfUV+v];

                integrationResult.at<Vec3f>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}


void numericalIntegration_onLambda(int colorChannel, vector<float *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult)
{
   // float currentValue = 0.0;
    //float nextValue = 0.0;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int u = 0 ; u<height ; u++)
    {

        for(int v = 0 ; v<width ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                float currentValue = integrand_forAllUVLambda[w][u*width+v];
                float nextValue = integrand_forAllUVLambda[w+1][u*width+v];

                integrationResult.at<Vec3f>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}

void numericalIntegration_onLambda(int colorChannel, vector<double *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, double samplingDistanceWavelength, Mat& integrationResult)
{
   // double currentValue = 0.0;
    //double nextValue = 0.0;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int u = 0 ; u<height ; u++)
    {

        for(int v = 0 ; v<width ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                double currentValue = integrand_forAllUVLambda[w][u*width+v];
                double nextValue = integrand_forAllUVLambda[w+1][u*width+v];

                integrationResult.at<Vec3f>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}


void numericalIntegration_onLambda(int colorChannel, vector<double *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, double samplingDistanceWavelength, Mat& integrationResult)
{
   // float currentValue = 0.0;
    //float nextValue = 0.0;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int u = 0 ; u<numberOfUV ; u++)
    {
        for(int v = 0 ; v<numberOfUV ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                double currentValue = integrand_forAllUVLambda[w][u*numberOfUV+v];
                double nextValue = integrand_forAllUVLambda[w+1][u*numberOfUV+v];

                integrationResult.at<Vec3d>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}

void numericalIntegration_onLambda_non_optimised(int p, int colorChannel, double factN, double factM, vector<float *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, float samplingDistance, Mat Ip[])
{
    for(int u = 0 ; u<numberOfUV ; u++)
    {
        for(int v = 0 ; v<numberOfUV ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                float currentValue = integrand_forAllUVLambda[w][u*numberOfUV+v];
                float nextValue = integrand_forAllUVLambda[w+1][u*numberOfUV+v];

                Ip[p].at<Vec3f>(u,v).val[colorChannel] += static_cast<float>(1.0/(factN*factM))*samplingDistance*(currentValue+nextValue)/2.0;
            }
        }
    }
}



double factorial(double n)
{
    if(n==0)
        return 1;
    else
        return n*factorial(n-1);
}


Mat pow2D(Mat& function2D, unsigned int power)
{
    int width = function2D.cols;
    int height = function2D.rows;
    Mat result = Mat::ones(width, height, CV_32FC1);

    function2D.convertTo(function2D, CV_32FC1);

    if(power >0)
    {
        result = function2D.clone();
        for(unsigned int k = 2 ; k<=power ; k++)
        {
           result = result.mul(function2D);
        }
    }

    return result;
}

float sign(float x)
{
    if(x<0)
        return -1.0;
    else
        return 1.0;
}

void preComputeGaussian(int numberOfUV, int size, float lambdaMin,
              int numberOfWavelength, float samplingDistanceWavelength, vector<float*>& result)
{
    int halfSizeUV = (numberOfUV-1)/2;
    float x = 0.0;
    float u = 0.0;
    float ul = 0.0;
    float a = 3.0;
    float sigmaS = 65.0/4.0;//1118.3;//65.0/4.0;//10.0;//SigmaS was 10.0 //The height map represents a size of 65 micrometer^2 in real life
    float sigmaF = 1.0/(2.0*M_PI*sigmaS);

    //Samplings
    //Xenopeltis_256 afm crop : 0.27968
    //Xenopeltis_zoom_256 afm crop : 0.0271875
    //CD_256 afm crop : 0.030438
    //Blazed grating 256 afm crop : 0.2859
    //Elaphe 128 0.0698
    float samplingDistanceInGrating = 0.083;//81.25;//0.08125;//0.8266;//;gratings 0.08125;//0.5;//maxres 0.850;//1024 0.919; //micrometers

    float currentLambda = 0.0;

    for(int k = -halfSizeUV ; k<=halfSizeUV ; k++)
    {
        ul = static_cast<float>(k)/static_cast<float>(halfSizeUV);
        u = 2.0*pow(ul, a);

        for(int t = 0 ; t<size ; t++)
        {
            for(int l = 0 ; l<numberOfWavelength ; l++)
            {
                currentLambda = lambdaMin + (float)l*samplingDistanceWavelength;
                x = u/currentLambda-(t-(float)size/2.0)/(float)(size*samplingDistanceInGrating);

                result[l][(k+halfSizeUV)*size+t] = static_cast<float>(exp(-x*x/(2.0*sigmaF*sigmaF)));
            }
        }
    }
}


QVector3D multiplyQMatrixVector3x3(QMatrix3x3 matrix, QVector3D vector)
{
    const float *matrixData = matrix.constData();

    //First row
    float x = matrixData[0]*vector.x()+matrixData[1]*vector.y()+matrixData[2]*vector.z();
    //Second row
    float y = matrixData[3]*vector.x()+matrixData[3+1]*vector.y()+matrixData[3+2]*vector.z();
    //Third row
    float z = matrixData[6]*vector.x()+matrixData[6+1]*vector.y()+matrixData[6+2]*vector.z();
cout << x << " " << y << " " << z << endl;
    return QVector3D(x, y ,z);
}

void normalizeVector(Mat &vector)
{
    float x = vector.at<float>(0,0);
    float y = vector.at<float>(1,0);
    float z = vector.at<float>(2,0);
    float norm = sqrt(x*x+y*y+z*z);

    vector.at<float>(0,0) = x/norm;
    vector.at<float>(1,0) = y/norm;
    vector.at<float>(2,0) = z/norm;

}



void normalizeVector(float &x, float &y, float &z)
{
    float norm = sqrt(x*x+y*y+z*z);

    if(norm>0)
    {
        x /= norm;
        y /= norm;
        z /= norm;
    }
}

Mat makeRotationMatrix(Point3f axis, float sin, float cos)
{
    Mat rotationMatrix = Mat::zeros(3,3, CV_32FC1);

    Mat Q = Mat::zeros(3,3, CV_32FC1);

    Q.at<float>(0,1) = -axis.z;
    Q.at<float>(0,2) = axis.y;

    Q.at<float>(1,0) = axis.z;
    Q.at<float>(1,2) = -axis.x;

    Q.at<float>(2,0) = axis.y;
    Q.at<float>(2,1) = -axis.x;

    rotationMatrix += Mat::eye(3,3, CV_32FC1);
    rotationMatrix += sin*Q;
    rotationMatrix += (1.0-cos)*Q*Q;

    return rotationMatrix;
}

float erfSigma(float x, float sigma)
{
    if(x>0.0 || x<0.0)
        return static_cast<float>(erf((double)x/(sqrt(2)*sigma)));
    else
        return 0.0;
}

double erfSigma(double x, double sigma)
{
    if(x>0.0 || x<0.0)
        return erf(x/(sqrt(2)*sigma));
    else
        return 0.0;
}
