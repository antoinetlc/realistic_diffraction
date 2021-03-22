/**
  * Practical acquisition and rendering of diffraction effects in surface reflectance.
  * Antoine Toisoul, Abhijeet Ghosh.
  * Imperial College London, December 2016.
  */

#include "convolution.h"

using namespace std;
using namespace cv;

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
Mat convolveEnvironmentMap(Mat &environmentMap, Mat &diffractionPattern, float xRange, float yRange)
{
    int heightEM = environmentMap.rows;
    int widthEM = environmentMap.cols;

    int widthPattern = diffractionPattern.cols;
    int heightPattern = diffractionPattern.cols;

    //Mat to store result
    Mat envMapFiltered = Mat::zeros(heightEM, widthEM, CV_32FC3);


    //Convolves each pixel of the environment map with the diffraction pattern
    for(int i = 0 ; i<heightEM ; ++i)
    {
        #pragma omp parallel for num_threads(omp_get_max_threads())
        for(int j = 0 ; j<widthEM ; ++j)
        {
                //Calculate reflection vector at the current pixel
                float r = 1.0;
                float phi(2.0*M_PI*(float)(j-(float)widthEM/2.0)/(float)widthEM);//Phi is 0 at the center of the environment map
                float theta(M_PI*(float)i/(float)heightEM);

                //Convert it to cartesian
                float x(0.0), y(0.0), z(0.0);
                sphericalToCartesian(r, theta, phi, x, y ,z);

                //Assumption that the reflection vector is equal to the viewing direction and the normal
                //See epic games :http://blog.selfshadow.com/publications/s2013-shading-course/karis/s2013_pbs_epic_notes_v2.pdf
                QVector3D reflectionVector = QVector3D(x,y,z);
                reflectionVector.normalize();

                QVector3D viewingVector = QVector3D(x,y,z);
                viewingVector.normalize();

                QVector3D normal = QVector3D(x,y,z);
                normal.normalize();

                float totalWeight[3] = {0.0, 0.0, 0.0};

                //Compute the convolution for each half vector of the diffraction lookup table
                for(int m = 0 ; m<heightPattern ; ++m)
                {
                    for(int n = 0 ; n<widthPattern ; ++n)
                    {
                        //For each color channel
                        for(int k = 0 ; k<3 ; ++k)
                        {

                            float intensityLobe = diffractionPattern.at<Vec3f>(m,n).val[k];

                            if(intensityLobe>0.0)
                            {
                                //For a given sample on the diffraction distribution
                                //Calculate the corresponding half vector
                                float x(0.0), y(0.0), z(0.0);

                                //(i,j) are indices of the image sampled in a linear range.
                                x = xRange*2.0*((float)n-(float)widthPattern/2.0)/(float)widthPattern;
                                y = yRange*2.0*((float)heightPattern/2.0-(float)m)/(float)heightPattern;

                                z = 1.0+sqrt(1.0-x*x-y*y);

                                QVector3D halfVector(x, y ,z);
                                halfVector.normalize();

                                //Find the light direction
                                QVector3D lightDirection = 2.0*QVector3D::dotProduct(viewingVector, halfVector)*halfVector-viewingVector;
                                lightDirection.normalize();

                                //Convert the light direction to spherical coordinates
                                cartesianToSpherical(lightDirection.x(), lightDirection.y(), lightDirection.z(),
                                                         r, theta, phi);


                                //Find the pixel the light direction corresponds to
                                int iLightIndex = floor((M_PI-theta)/M_PI*(height-1));
                                int jLightIndex = floor(width-1-(phi+M_PI)/(2.0*M_PI)*(width-1));

                                envMapFiltered.at<Vec3f>(i,j).val[k] += intensityLobe*environmentMap.at<Vec3f>(iLightIndex,jLightIndex).val[k]*sin(theta);
                                totalWeight[k] += intensityLobe;
                            }
                        }
                    }
                }

                float maxWeight = max(totalWeight[2], max(totalWeight[1], totalWeight[0]));
                envMapFiltered.at<Vec3f>(i,j).val[2] /= maxWeight;
                envMapFiltered.at<Vec3f>(i,j).val[1] /= maxWeight;
                envMapFiltered.at<Vec3f>(i,j).val[0] /= maxWeight;

        }
    }

    return envMapFiltered;
}
