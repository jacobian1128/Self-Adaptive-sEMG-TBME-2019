# Self-Adaptive-sEMG-TBME-2019

These are C++ examples implementing the publication below.

## Publication
Seongsik Park, Wan Kyun Chung and Keehoon Kim, "Training-free Bayesian self-adaptive classification for sEMG pattern recognition including motion transition," *IEEE Transactions on Biomedical Engineering* ([DOI:10.1109/TBME.2019.2947089](http://doi.org/10.1109/TBME.2019.2947089)).

## Environment and Dependency
Library, SDK or API that each example requires

### CUDA and MFC
* [CUDA toolkit 10.2](https://developer.nvidia.com/cuda-toolkit) for GPU parallel process
* [VTK 9.0.0](https://vtk.org/) for visualization built by VS2019

### Thalmic MYO armband
* [Eigen](http://eigen.tuxfamily.org) to handle vectors and matrices
* [Myo Windows SDK](https://support.getmyo.com/hc/en-us/articles/360018409792-Myo-Connect-SDK-and-firmware-downloads) to utilize sEMG inside MYO armband

### MATLAB mex function
* [Eigen](http://eigen.tuxfamily.org) to handle vectors and matrices
* [MATLAB C++ MEX API](https://www.mathworks.com/help/matlab/matlab_external/cpp-mex-api.html) to realize mex function for MATLAB
