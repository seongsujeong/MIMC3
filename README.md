MIMC3
=====

Code by Seongsu Jeong at the Ohio State University (now at University of California, Irvine)

* Algorithm
  * The algorithm implemented in here is described in a paper "Improved Multiple Matching Method for Observing Glacier Motion With Repeat Image Feature Tracking" by [Jeong et al in 2017](https://ieeexplore.ieee.org/document/7827084) at IEEE Transactions on Geoscience and Remote Sensing (DOI: 10.1109/TGRS.2016.2643699)
  * The video clip below is about how DLC works to find the matching point.
  ![DLC](docs/DLC_slow.mov)

* Prerequisites and Dependency:
  * You will need a C compiler that supports OpenMP. GCC should do the compilation done.
  * libtiff is necessary in order to build the MIMC3. The source code along with the instruction can be found at http://www.libtiff.org

* Compilation command example:
  * In shell:

  ```bash
  gcc -fopenmp -O3 -g -o MIMC3 -I$HOME/miniconda3/include/ MIMC_main.c GMA.c georefimg.c MIMC_module.c MIMC_misc.c $HOME/miniconda3/lib/libtiff.dylib
  ```

* Usage:
  * Add the path of libtiff.dylib (in case of mac) or libtiff.so (linux) to LD_LIBRARY_PATH in environment setting.
  * MIMC3 [first image file] [second image file] [xyuvav file] [output directory]
  * Note that the dimension of the first and the second image has to be the same.
  * xyuvav is n-by-6 matrix that contains the grid location (image & map) as well as the *a priori* velocity information, that MIMC3 refers.

* Version history:
  * 1.0 : MIMC algorithm first proposed in a paper "Efficient automated glacier surface velocity measurement from repeat images using multi-image/multichip and null exclusion feature tracking" by Ahn and Howat in IEEE Transactions on Geoscience and Remote Sensing (DOI: 10.1109/TGRS.2011.2114891). Implementation was in MATLAB associated with C-MEX
  * 2.0 : The algorithm was revamped by DLC and QM. Implementation was in MATLAB, with utilization of C-MEX and parallel toolbox. This is the first implementation of the algorithms presented at [Jeong et al. (2017)](https://ieeexplore.ieee.org/document/7827084).
  * **3.0 : Codes in this repo.** Minor revision of the algorithm implemented in MIMC 2.0. Fully implemented in C program with utilization of OpenMP.
