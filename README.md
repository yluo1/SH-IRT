# Spherical Harmonic - Impulse Response Tools (SH-IRT)

This open-source library provides tools for modifying, generating, and analyzing spatial room impulse responses (SRIRs) in the spherical harmonics (SH) domain. SH-IRT supports the following features:

* **Operations**:  Rotation, reflection, addition, multiplication, conjugation, magnitude square, convolution, normalization, translation, complex-real conversion
* **Evaluations**: SH basis, expansion, integral (sum, power), covariance
* **Fitting**: Least squares (truncated singular value decomposition, Tiknonov regularization, magnitude squared, sum-of-magnitude squared, mix-of-magnitude squared, mixture power), radial basis functions
* **Encoding**: Projections, radial basis functions, spherical baffled piston
* **Filtering**: Time / frequency domain LTI convolution, time-varying exponentiating and recursive convolution, Tiknonov least squares FIR fitting, iterative re-weighted least squares IIR fitting
* **Density Modeling**: Maximum likelihood fitting, inverse sampling, scattering transport, spherical sliced Wasserstein optimal transport
* **SRIR Generation**: Cross-directivity SH image-source model (ISM) expansion, spherical Poisson-process
* **T60 Modeling**: Spherical frequency decay time Gaussian process (GP) regression and sampling
* **Misc**: Projection plotting, SH-AmbiX conversion, spherical coordinate grid generation (uniform spherical Fibonacci, Platonic-solid vertices, polar sampling)

## Table of Contents

* [Installation](#installation)
* [Guides and Tutorials](#guides-and-tutorials)
* [References](#references)

## Installation
Add the following sub-folders to your search path by calling the Matlab command-line function:
```
startup
```

| Folder | Description |
| :--- | :--- | 
| sh | Spherical Harmonics |
| ft | Filter Toolbox|
| gp |  Gaussian Processes |
|bf| Beamformer Designs |
|tst| Unit Testing |
| plt | Plotting Examples |
| err | Error Metrics|
| figs | Figures | 

### Software Requirements
* Implemented in Matlab (2024b) with the following toolboxes:
  ```
  --- Required MathWorks Products ---
  MATLAB (Version: 24.2)
  Optimization Toolbox (Version: 24.2)
  Signal Processing Toolbox (Version: 24.2)
  Symbolic Math Toolbox (Version: 24.2)
  Statistics and Machine Learning Toolbox (Version: 24.2)
  Global Optimization Toolbox (Version: 24.2)
  ```
*  Semi-definite programs require the [convex optimization library](https://cvxr.com/cvx/). If you’re installing on Apple silicon, you can install the pre-built binaries following the [instructions](https://ask.cvxr.com/t/announcement-cvx-for-apple-silicon/12280/).

## Guides and Tutorials

* [Functions Guide](guide_func.md)
* [Spherical Harmonic Operators Tutorial](tutorial_sh_ops.md)
* [Spherical Harmonic Basis Fitting Tutorial](tutorial_sh_fit.md)
* [Spatial Room Impulse Response (SRIR) T60 Augmentation Tutorial](tutorial_t60_aug.md)

## References

This library was developed from methods in the following works:

> Luo, Y. (2026). [Fast Time-Varying Exponentiated Convolution Methods for Generative Direction Dependent Reverberation](https://arxiv.org/abs/2609.24809), In AES Show 2026 Nashville. Audio Engineering Society; 2026.pp. 59:1-10.
>
> Luo, Y. (2026). [Spherical Harmonic Sliced Wasserstein Displacement Interpolation for Acoustic Source and Reflection Density Modeling](https://arxiv.org/abs/2609.22028). 	arXiv:2609.22028.
> 
> Luo, Y. (2024, August). [Constant Directivity Loudspeaker Beamforming](https://ieeexplore.ieee.org/abstract/document/10715159). In 2024 32nd European Signal Processing Conference (EUSIPCO) (pp. 246-250). IEEE.
> 
> Luo, Y. (2021). [Spherical Harmonic Covariance and Magnitude Function Encodings for Beamformer Design](https://link.springer.com/article/10.1186/s13636-021-00230-7). EURASIP Journal on Audio, Speech, and Music Processing, 2021(1), 41.
> 
> Luo, Y., Kim, W. (2020). [Fast Source-Room-Receiver Acoustics Modeling](https://ieeexplore.ieee.org/abstract/document/9287377). In 2020 28th European Signal Processing Conference (EUSIPCO) (pp. 51-55). IEEE.

 

