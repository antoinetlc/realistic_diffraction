# Realistic Diffraction

Practical Acquisition and Rendering of Diffraction Effects in Surface Reflectance. Antoine Toisoul and Abhijeet Ghosh. ACM Transactions on Graphics, 2017 (presentation at SIGGRAPH 2017).


[![[SIGGRAPH 2017] - Practical acquisition and rendering of diffraction effects in surface reflectance](https://img.youtube.com/vi/DNOej4KgtJY/0.jpg)](https://www.youtube.com/watch?v=DNOej4KgtJY)

#### Description

Realistic diffraction is an implementation of the paper *Practical Acquisition and Rendering of Diffraction Effects in Surface Reflectance. Antoine Toisoul and Abhijeet Ghosh. ACM Transactions on Graphics, 2017.* 
It allows : 
* the computation of diffraction lookup tables given the HDR photograph of a diffraction pattern measured with a spectral filter
* the prefiltering of an environment map with a diffraction lookup table for real-time rendering of diffraction effects in arbitrary environments

The software requires OpenCV and has been tested with OpenCV 2.4.11.

A video of the renderings is available [here](https://www.youtube.com/watch?v=DNOej4KgtJY)

#### Usage

An example of how to compute a diffraction lookup table and the convolution are given in `main.cpp`. The function that computes the lookup table is given in `huygens.cpp`. The function that computes the prefiltering is given in `convolution.cpp`.

#### Citation

Please cite the following paper when using Realistic Diffraction :

    Practical Acquisition and Rendering of Diffraction Effects in Surface Reflectance. 
    Antoine Toisoul and Abhijeet Ghosh. ACM Transactions on Graphics,36(5), 2017.

#### Credits and License

Realistic Diffraction. Author : Antoine TOISOUL. Copyright © 2016 Antoine TOISOUL, Imperial College London. All rights reserved.

Realistic Diffraction is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. Realistic Diffraction is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. You should have received a copy of the GNU Lesser General Public License along with this program. If not, see http://www.gnu.org/licenses/.

#### Contact

Antoine Toisoul : antoine.toisoul 'AT' telecom-paristech.org
