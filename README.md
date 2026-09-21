# Spherical Harmonic - Impulse Response Tools (SH-IRT)

This open-source library provides tools for modifying, generating, and analyzing spatial room impulse responses (SRIRs) in the spherical harmonics (SH) domain. SH-IRT supports the following features:

* **Operators**:  Rotation, reflection, convolution, multiplication, conjugation, square magnitude, integration, complex-real conversion
* **Fitting**: Least squares, magnitude squared least squares, sum-of-magnitude squared least squares
* **Filtering**: Time / frequency domain LTI convolution, time-varying exponentiating and recursive convolution
* **Density Modeling**: Maximum likelihood fitting, inverse sampling, scattering transport, spherical sliced Wasserstein optimal transport
* **Encoding**: Projections, radial basis functions, spherical piston
* **Generating**: Cross-directivity SH image-source model (ISM) expansion, spherical Poisson-process
* **T60 Modeling**: Spherical frequency decay time Gaussian process (GP) regression and sampling
* **Misc**: Projection plotting, SH-AmbiX conversion, uniform spherical Fibonacci and Platonic-solid point sampling

## Content

* [Installation](#installation)
* [Guides and Tutorials](#guides-and-tutorials)
* [References](#references)

## Installation
Add the following sub-folders to your search path by calling the Matlab command-line function:
```
startup
```

| Folder | Description |
| --- | --- | 
| sh | Spherical harmonics functions |
| ft | Time-varying filtering|
| gp | Spherical x frequency Gaussian processes |
| err | Error functions|
| plt | Plotting examples |
|tst| Unit testing |

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
*  Semi-definite programs require the [convex optimization library](https://cvxr.com/cvx/). If you’re installing on Apple silicon, you can use the pre-built binaries following the [instructions](https://ask.cvxr.com/t/announcement-cvx-for-apple-silicon/12280/).

## Guides and Tutorials

[Functions Guide](guide_func.md)

[Spatial Room Impulse Response (SRIR) T60 Augmentation Tutorial](tutorial_t60_aug.md)

## References

This library was developed from methods in the following works:

> Luo, Y. (2026). "Fast Time-Varying Exponentiated Convolution Methods for Generative Direction Dependent Reverberation", Proceedings of the 161th Audio Engineering Society Convention.
>
> Luo, Y. (2026). [Spherical Harmonic Sliced Wasserstein Displacement Interpolation for Acoustic Source and Reflection Density Modeling](https://arxiv.org/abs/2609.22028). 	arXiv:2609.22028
> 
> Luo, Y. (2021). [Spherical harmonic covariance and magnitude function encodings for beamformer design](https://link.springer.com/article/10.1186/s13636-021-00230-7). EURASIP Journal on Audio, Speech, and Music Processing, 2021(1), 41.
> 
> Luo, Y., Kim, W. (2020). [Fast source-room-receiver acoustics modeling](https://ieeexplore.ieee.org/abstract/document/9287377). In 2020 28th European Signal Processing Conference (EUSIPCO) (pp. 51-55). IEEE.

 

