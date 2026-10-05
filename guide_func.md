# Functions Guide

Table of contents:
* [Spherical Harmonics](#spherical-harmonic-functions)
  * [Operators](#sh-operators)
  * [Evaluation](#sh-evaluations)
  * [SRIR Generation](#sh-spatial-room-impulse-response-generation)
  * [Basis Fitting](#sh-basis-fitting)
  * [Encodings](#sh-function-encodings)
  * [Density Function](#sh-probability-density-functions)
  * [Grid Generation](#sh-spherical-coordinate-grids)
  * [Miscellaneous](#sh-miscellaneous)
* [Filter Toolbox](#filter-toolbox)
  * [Exponentiating Filtering](#ft-time-varying-exponentiation)
  * [Filter Fitting](#ft-filter-fitting)
* [Gaussian Process](#gaussian-processes)

# Spherical Harmonic Functions

In this library, we default to the complex spherical harmonic (SH) basis functions given by

$$Y_l^m (\theta, \phi) = \sqrt{\frac{(2l + 1)}{4 \pi} \frac{(l-m)!}{(l+m)!} } P_l^m(\cos \theta ) e^{j m \phi},$$

where $P_l^m(\cos \theta)$ are the associated Legendre polynomials of degree $l$ and order $m$ that includes the Condon-Shortley phase $(-1)^m$ term, and spherical coordinates $(\theta, \phi)$ are co-latitude and azimuth respectively.

Functions are represented by their truncated SH expansion coefficients given by

$$f_C(\theta, \phi) = \sum_{l=0}^{L_C} \sum_{m=-l}^l Y_l^m (\theta, \phi) C_l^m, $$

and vectorized by increasing degree and order $\bf{C} = \left [ C_0^0 , C_1^{-1}, C_1^{0}, C_1^{1}, … , C_{L_C}^{L_C}   \right ] \in \mathbb{C}^{(L_C + 1)^2}$.

Some library functions support expansion coefficients belonging to the real spherical harmonics form following:

$$ Y_{lm}(\theta, \phi) = \left \lbrace \begin{array}{cc}
(-1)^m \sqrt{2} \sqrt{\frac{(2l + 1)}{4 \pi} \frac{(l-|m|)!}{(l+|m|)!} } P_l^{|m|}(\cos \theta ) \sin (m |\phi|), & m < 0\\
%%%%%%%%%
\sqrt{\frac{2l + 1}{4 \pi} } P_l^0 (\cos \theta) , &  m = 0 \\
%%%%%%%%
(-1)^m \sqrt{2} \sqrt{\frac{(2l + 1)}{4 \pi} \frac{(l-m)!}{(l+m)!} } P_l^m(\cos \theta ) \cos (m \phi), & m > 0
\end{array} \right . ,$$

whereby a function’s truncated real SH expansion follows

$$f_C(\theta, \phi) = \sum_{l=0}^{L_C} \sum_{m=-l}^l Y_{lm} (\theta, \phi) C_l^m, $$ 

and is vectorized in the same format as the complex case. In such instances, the library function has an `is_real` input argument which can be set to `true` for the real case, and `false` for the complex case.


## SH Operators
| File | Description |
| :--- | :--- | 
| sh_rot.m | Rotation | 
| sh_refl.m | Reflection | 
| sh_mul.m | Multiplication | 
| sh_conv.m | Convolution | 
| sh_nrm.m | Normalizations | 
| sh_int.m | Integration | 
| sh_msq.m | Magnitude squared| 
| sh_sqrt.m | Square root magnitude approximation | 
| sh_conj.m | Conjugation | 
| sh_cpx2re.m | Complex to real form | 
| sh_re2cpx.m | Real to complex form | 
| sh_resize.m | Truncation and zero-padding | 

## SH Evaluations
| File | Description |
| :--- | :--- | 
|sh_val.m| Evaluate SH bases |
|sh_dec.m| Evaluate SH expansions |

## SH Spatial Room Impulse Response Generation
| File | Description |
| :--- | :--- | 
|sh_ism.m| Image-source model |
|sh_rand.m| Random field |
|sh_rand_pp.m| Poisson process |

## SH Basis Fitting
| File | Description |
| :--- | :--- | 
|sh_fit_svd.m | Truncated singular value decomposition least squares: Maximum singular value fraction, total variance fraction, bottom fraction, discrete Picard criterion |
|sh_fit_tr.m | Tikhonov regularized least squares: l2 norm, quadratic degree penalty, discrete Picard criterion |
|sh_fit_rbf.m | Radial basis function fitting: Squared exponential, Matérn, exponential |
|sh_fit_msq.m | Magnitude squared least-squares: Magnitude squared, sum-of-magnitude squared, mix-of-magnitude squared, mixture power|


## SH Filtering
| File | Description |
| :--- | :--- | 
|sh_filter.m | Filter time domain SH expansion |
|sh_filter_freq.m | Filter frequency domain SH expansion |
|sh_filter_delay.m | Add sample delay to SH expansion at angular frequencies |
|sh_filter_dir.m | Directional filtering |
|sh_filter_val.m | Evaluate FIR fitted SH expansions at angular frequencies for resampling SH expansions |
|sh_freqz.m | Frequency response of SH expansion |
|sh_exp_conv.m | Direction independent time-varying exponentiating convolution|
|sh_exp_conv_gp.m | Direction dependent T60 time-varying exponentiating convolution |
|sh_rec_conv.m | Time-varying recursive convolution |

## SH Function Encodings
| File | Description |
| :--- | :--- | 
|sh_enc_proj.m | Dirac-delta projection into spherical harmonics |
|sh_enc_proj_msq.m | Dirac-delta projection magnitude square into spherical harmonics |
|sh_enc_uni.m | Constant function |
|sh_enc_rbf.m | Radial basis functions: Squared exponential, Matérn, exponential, sinc |
|sh_enc_pist_sphere.m | External piston on sphere frequency responses|

## SH Probability Density Functions
| File | Description |
| :--- | :--- | 
|sh_pdf_fit.m | Density function fitting|
|sh_pdf_dist.m | Density function distances: Spherical sliced Wasserstein, Kullback-Lieber divergence|
|sh_pdf_sample.m | Sampling spherical coordinates from density|
|sh_pdf_scatter.m | Transport density function towards uniform density|
|sh_pdf_transport.m | Transport between density functions|
|sh_pdf_preset.m | Preset density functions|
|sh_cdf_inv_theta.m | Inverse sampling marginal cumulative distribution function over co-latitude|
|sh_cdf_inv_phi_cond.m | Inverse sampling cumulative distribution function over azimuth given co-latitude|

## SH Spherical Coordinate Grids
| File | Description |
| :--- | :--- | 
|sh_grd_fib.m| Spherical Fibonacci points |
|sh_grd_plat.m| Platonic solid vertex points |
|sh_grd_caps.m| Uniform points along co-latitude and azimuth (polar concentrated)|
|sh_grd_rand.m | Uniform random spherical coordinates |

## SH Miscellaneous
| File | Description |
| :--- | :--- | 
|sh_plt.m| Plot SH expansion|
|sh2ambx.m| SH to AmbiX format |
|ambx2sh.m| AmbiX to SH format |


# Filter Toolbox

## FT Time-varying Exponentiation
| File | Description |
| :--- | :--- | 
|ft_exp_conv_opt.m| Exponentiated convolution optimized |
|ft_exp_conv_direct.m| Exponentiated convolution direct |
|ft_exp_design.m| Exponentiating FIR filter design |
|ft_two_tap_FIR.m| Two-tap FIR filter design |
|ft_rec_conv_opt.m| Recursive convolution optimized |
|ft_rec_conv_direct.m| Recursive convolution direct |

## FT Filter Fitting
| File | Description |
| :--- | :--- | 
|ft_fit_bnd_minphase.m| Magnitude bounded minimum phase FIR fit to magnitude targets|
|ft_fit_tr.m | Tiknonov regularized least squares FIR fit to complex targets: l2 norm, exponential and circular exponential window, discrete Picard criterion |
|ft_fit_sh.m| Tiknonov regularized least squares FIR fit to SH coefficients: l2 norm, exponential and circular exponential window|
|ft_freq_wt.m | Frequency weighting functions|

# Gaussian Processes

| File | Description |
| :--- | :--- | 
|gp_mu.m| Prior mean function selector|
|gp_cov.m| Covariance function selector |
|gp_t60_optimize.m | Hyper-parameter optimization of T60 model|
|gp_t60_sample.m | T60 function sampling from GP prior or posterior|
|gp_plt.m| Plotting |
|mu_pow.m| Power-law prior mean function |
|mu_lpf.m| Low-pass prior mean function |
|cov_sqx_chw_ns.m| Squared exponential chordal x non-stationary frequency covariance function |
|cov_sqx_chw_sqx.m| Squared exponential chordal x squared exponential frequency covariance function |
