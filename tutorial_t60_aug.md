# Tutorial: Spatial Room Impulse Response T60 Augmentation

In this tutorial, we cover sound decay time augmentation of room impulse responses (RIRs) from our paper 
>Yuancheng Luo, "Fast Time-Varying Exponentiated Convolution Methods for Generative Direction Dependent Reverberation", Proceedings of the 160th Audio Engineering Society Convention.

The high-level steps are as follows:
* [Sampling T60 from Gaussian Processes (GPs):](#sampling-t60-functions-from-gaussian-processes) Sample or design a sound decay time function $T_{60}(\omega, \theta, \phi)$ of angular frequency $\omega$, and spherical coordinates $(\theta, \phi)$ co-latitude and azimuth respectively
* [Exponentiating FIR Optimization:](#fitting-fir-to-t60-functions) Fit a short finite impulse response (FIR) exponentiating filter $\bf{g}$ to sampled $T_{60}(\omega, \theta, \phi)$ functions
* [Specifying or generating colorless RIR:](#generating-colorless-room-impulse-responses) Specify an input RIR or generate a colorless (constant T60) RIR $\bf{h}$
* [Augmenting RIR:](#applying-time-varying-exponentiated-convolution) Apply time-varying exponentiated convolution $\textbf{f} = f(\bf{g}, \bf{h})$

## Sampling T60 Functions from Gaussian Processes
We can sample smooth T60 functions $T_{60}(\omega, \theta, \phi)$ of frequency and spherical coordinates from Gaussian processes distributions defined by prior mean and covariance functions.

### Mean Function Prior Specifications

A simple T60 prior mean is the power law function given by

$$\overline{T}_{60}(\omega) =  \alpha  \max \left ( 1, \nu(\omega)  \right ) ^{-\beta}, \quad \nu(\omega) = \frac{\omega}{2 \pi}, $$

where $\alpha$ is the T60 (seconds) at 0 Hz, $\beta$ is the decay exponent, and $\nu(\omega)$ is ordinary frequency. Increasing the decay rate $\beta$ inflects the T60 curve towards 0 seconds for larger frequencies in the following figure: 

```
gp_plt('mu_pow');
```
<img src="./figs/figs_t60/mu_pow_1.png" alt="Power Function" width="600"/>

We can introduce a 3rd shape parameter to control the tilt following an idealized low-pass filter response design given by

$$\overline{T}_{60}(\omega) =  \frac{\alpha}{1 + \left ( \frac{\nu(\omega)}{f_c}  \right )^{\beta}}, $$

where $f_c$ is the cross-over frequency relative to the flat T60 curve $(\beta = 0)$ shown in the following figure:

```
gp_plt('mu_lpf');
```
<img src="./figs/figs_t60/mu_lpf_1.png" alt="Power Function" width="900"/>

The GP prior mean function can be evaluated via `gp_mu.m` by specifying its options struct `gp_mu_opts.m` and plotted:
```
omega = 2 * pi * logspace(log10(20), log10(24000), 256)'; % Angular frequencies

mu_power = gp_mu(omega, gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.5, 'enable_disp', true));
mu_lpf   = gp_mu(omega, gp_mu_opts('mu_func', 'LPF',   'mu_alpha', 1, 'mu_beta', 0.5, 'mu_fc', 500, 'enable_disp', true));
```

### Covariance Function Prior Specifications

We can construct positive semi-definite covariance functions from chordal distances of spherical coordinates $(\theta, \phi)$, $(\theta’, \phi’)$ and there corresponding unit-directions  $\bf{v}$, $\bf{v}'$  on the unit sphere in $\mathbb{R}^3$ given by 

$$d(\theta, \phi, \theta', \phi') = 2 \sin \left ( \frac{ \left | \, \cos^{-1} (\bf{v}^T \bf{v}') \right |   }{2} \right ) = \left \| \bf{v} - \bf{v}'  \right \|_2 .$$

The squared exponential of chordal distance with non-stationary kernel wavelengths is therefore a non-stationary Gaussian kernel[^PACIOREK_NS] in 3-dimensions where $\bf{\Sigma} = \lambda^2 \bf{I} \in \mathbb{R}^{3 \times 3}$, and is positive semi-definite following  

$$k(\bf{x},\bf{x}') = \sigma^2 \left ( \frac{2 \lambda \lambda'}{\lambda^2 + (\lambda')^2} \right )^{3/2} \exp \left (  - \frac{d^2(\theta, \phi, \theta', \phi') }{(\lambda^2 + (\lambda')^2) / 2}  \right ),$$

where $\lambda = \ell / f^{\gamma}$ is the scaled wavelength for velocity $\ell$ (m/s), ordinary frequency $f = \angle \omega / (2 \pi)$ (Hz), and power $\gamma$. Increasing the hyper-parameter $\gamma$ decreases the covariance between low and high frequencies of the T60 functions. Therefore, larger $\gamma$ decreases smoothness in high frequency as shown in the following figure:

```
gp_plt('cov_sqx_chw_ns');
```
<img src="./figs/figs_t60/cov_sqx_chw_ns_1.png" alt="Squared Exponential Chordal Distance with Non-stationary Frequency" width="1200"/>

where the maximum covariances occur at $\lambda = \ell / f_0^{\gamma}$ for varying $\lambda' = \ell / f_1^{\gamma}$. The covariance function’s shape follows

* T60 at frequency $f_0$ covaries more with lower frequencies than with higher frequencies given the same angular separation.
* T60s at lower frequencies are smoother than at higher frequencies.


We compare our non-stationary covariance with the stationary product of squared exponential of chordal distance and squared exponential of log-frequency distances given by

$$k(\bf{x},\bf{x}') = \sigma^2 \exp \left (  - \frac{d^2(\theta, \phi, \theta', \phi') }{2 \ell_c^2}  \right )  \exp \left (  - \frac{(log(f) - log(f'))^2 }{2 \ell_f^2}  \right ) ,$$

where $\ell_c$, $\ell_f$ are length-scale hyper-parameters of the chordal and log-frequency distances respectively. Therefore, the covariance function is stationary w.r.t. the chordal and log-frequency distances as shown in the following figure:

```
gp_plt('cov_sqx_chw_sqx');
```
<img src="./figs/figs_t60/cov_sqx_chw_sqx_1.png" alt="Squared Exponential Chordal Distance x  Squared Exponential of Log-Frequency" width="1200"/>

where frequency $f = f_0$ for varying $f' = f_1$. Unlike the non-stationary covariance, the low and high frequencies covary only by their octave separation, and are independent of the absolute frequency.

The GP prior covariance function can be evaluated via `gp_cov.m` by specifying its options struct `gp_cov_opts.m` and plotted:

```
N_B = 16;     % Number of log-uniform angular frequencies
N_E = 20;     % Number of uniform spherical coordinates

omega         = 2 * pi * logspace(log10(20), log10(24000), N_B)';
[theta, phi]  = sh_fib(N_E);
freq          = max(1, omega / (2 * pi)); %[N_B x 1]

freq_grid     = repmat(freq, [1, N_E]);
theta_grid    = repmat(theta', [N_B, 1]);
phi_grid      = repmat(phi', [N_B, 1]);

X = [2 * pi * freq_grid(:), theta_grid(:), phi_grid(:)]; %[N_B * N_E x 1] Input vector

[K, dK_mat_list, dK_name_list] = gp_cov(X, X, gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', 1, 'cov_ell', 343, 'cov_gamma', 1, 'enable_disp', true));
[K, dK_mat_list, dK_name_list] = gp_cov(X, X, gp_cov_opts('cov_func', 'cov_sqx_chw_sqx', 'cov_sigma', 1, 'cov_ell_c', 1, 'cov_ell_f', 2, 'enable_disp', true));
```

### Gaussian Process Specifications

### Sampling from the Prior and Posterior Distributions

## Fitting FIR to T60 Functions

## Generating Colorless Room Impulse Responses

## Applying Time-Varying Exponentiated Convolution

[^PACIOREK_NS]: Paciorek, C., & Schervish, M. (2003). Nonstationary covariance functions for Gaussian process regression. Advances in neural information processing systems, 16.

