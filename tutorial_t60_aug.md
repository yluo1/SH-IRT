# Tutorial: T60 Reverberation Augmentation



## Sampling T60 Functions from Gaussian Processes
We can sample smooth T60 functions $T_{60}(\omega, \theta, \phi)$ of frequency and spherical coordinates from Gaussian processes distributions.

We can construct positive semi-definite covariance functions from chordal distances of spherical coordinates $(\theta, \phi)$, $(\theta’, \phi’)$ and there corresponding unit-directions  $\bf{v}$, $\bf{v}'$  on the unit sphere given by 

$$d(\theta, \phi, \theta', \phi') = 2 \sin \left ( \frac{ \left | \cos^{-1} (\bf{v}^T \bf{v}') \right |   }{2} \right ).$$

The squared exponential of chordal distance and non-stationary frequency is positive semi-definite and given by
$$k(\bf{x},\bf{x}') = \sigma^2 \sqrt{\frac{2 \lambda \lambda'}{\lambda + (\lambda')^2}} \exp \left (  - \frac{d^2(\theta, \phi, \theta', \phi') }{(\lambda^2 + (\lambda')^2) / 2}  \right ), $$

where $\lambda = \ell / f^{\gamma}$ is the scaled wavelength for velocity $\ell$ (m/s), ordinary frequency $f = \angle \omega / (2 \pi)$ (Hz), and power $\gamma$. Increasing the hyper-parameter $\gamma$ decreases the covariance between low and high frequencies of the T60 functions. Therefore, larger $\gamma$ decreases smoothness in high frequency as shown in the following figure:

```
gp_plt('cov_sqx_chw_ns');
```
<img src="./figs/figs_t60/cov_sqx_chw_ns_1.png" alt="Squared Exponential Chordal Distance with Non-stationary Frequency" width="1200"/>

where the maximum covariances occur at $\lambda = \ell / f_0^{\gamma}$ for varying $\lambda' = \ell / f_1^{\gamma}$.