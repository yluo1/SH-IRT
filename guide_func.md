# Spherical Harmonic (SH) Functions

In this library, we default to the complex spherical harmonic basis functions given by

$$Y_l^m (\theta, \phi) = \sqrt{\frac{(2l + 1)}{4 \pi} \frac{(l-m)!}{(l+m)!} } P_l^m(\cos \theta ) e^{j m \phi},$$

where $P_l^m(\cos \theta)$ are the associated Legendre polynomials of degree $l$ and order $m$ that includes the Condon-Shortley phase $(-1)^m$ term, and spherical coordinates $(\theta, \phi)$ are co-latitude and azimuth respectively.

Functions are represented by their truncated SH expansion coefficients given by

$$f_C(\theta, \phi) = \sum_{l=0}^{L_C} \sum_{m=-l}^l Y_l^m (\theta, \phi) C_l^m, $$

and vectorized by increasing degree and order $\bf{C} = \left [ C_0^0 , C_1^{-1}, C_1^{0}, C_1^{1}, … , C_{L_C}^{L_C}   \right ] \in \mathbb{C}^{(L_C + 1)^2}$.

Some library functions support expansion coefficients belonging to the real spherical harmonics form following:

$$ Y_{lm}(\theta, \phi) = \left \lbrace \begin{array}{cc}
(-1)^m \sqrt{2} \sqrt{\frac{(2l + 1)}{4 \pi} \frac{(l-|m|)!}{(l+|m|)!} } P_l^{|m|}(\cos \theta ) \sin (m |\phi|), & m < 0\\[6pt]
%%%%%%%%%
\sqrt{\frac{2l + 1}{4 \pi} } P_l^0 (\cos \theta) , &  m = 0 \\[6pt]
%%%%%%%%
(-1)^m \sqrt{2} \sqrt{\frac{(2l + 1)}{4 \pi} \frac{(l-m)!}{(l+m)!} } P_l^m(\cos \theta ) \cos (m \phi), & m > 0
\end{array} \right . ,$$

whereby a function’s truncated real SH expansion follows

$$f_C(\theta, \phi) = \sum_{l=0}^{L_C} \sum_{m=-l}^l Y_{lm} (\theta, \phi) C_l^m, $$ 

and is vectorized in the same format as the complex case. In such instances, the library function has an `is_real` input argument which can be set to `true` for the real case, and `false` for the complex case.


## Operators
| File | Description |
| --- | --- | 
| sh_rot.m | Rotation | 
| sh_refl.m | Reflection | 
| sh_mul.m | Multiplication | 
| sh_conv.m | Convolution | 
| sh_nrm.m | Normalizations | 
| sh_int.m | Integration | 
| sh_msq.m | Squared magnitude | 
| sh_sqrt.m | Square root magnitude approximation | 
| sh_conj.m | Conjugation | 
| sh_cpx2re.m | Complex to real form | 
| sh_re2cpx.m | Real to complex form | 
| sh_resize.m | Truncation and zero-padding | 

## Evaluations
| File | Description |
| --- | --- | 
|sh_val.m| Evaluate SH bases |
|sh_dec.m| Evaluate SH expansions |

## Spatial Room Impulse Response Generation
| File | Description |
| --- | --- | 
|sh_ism.m| Image-source model |
|sh_rand.m| Random field |
|sh_rand_pp.m| Poisson-process |

## Empirical Fitting
| File | Description |
| --- | --- | 
|sh_fit_svd.m | Truncated singular value decomposition least-squares |
|sh_fit_msq.m | Magnitude squared least-squares |
|sh_fit_rbf.m | Radial basis function fitting |


## Filtering
| File | Description |
| --- | --- | 
|sh_filter.m | Filter time domain SH expansion |
|sh_filter_freq.m | Filter frequency domain SH expansion |
|sh_filter_dir.m | Directional filtering |
|sh_freqz.m | Frequency response of SH expansion |
|sh_exp_conv.m | Direction independent time-varying exponentiating convolution|
|sh_exp_conv_gp.m | Direction dependent T60 time-varying exponentiating convolution |
|sh_rec_conv.m | Time-varying recursive convolution |

## Function Encodings
| File | Description |
| --- | --- | 
|sh_enc_proj.m | Dirac-delta projection into spherical harmonics |
|sh_enc_proj_msq.m | Dirac-delta projection magnitude square into spherical harmonics |
|sh_enc_uni.m | Constant function |
|sh_enc_rbf.m | Radial basis functions |
|sh_enc_pist_sphere.m | External piston on sphere frequency responses|

## Probability Density Functions
| File | Description |
| --- | --- | 
|sh_pdf_fit.m | Density function fitting|
|sh_pdf_sample.m | Sampling spherical coordinates from density|
|sh_pdf_scatter.m | Transport density function towards uniform density|
|sh_pdf_transport.m | Transport between density functions|
|sh_pdf_preset.m | Preset density functions|
|sh_cdf_inv_theta.m | Inverse sampling marginal cumulative distribution function over co-latitude|
|sh_cdf_inv_phi_cond.m | Inverse sampling cumulative distribution function over azimuth given co-latitude|

## Miscellaneous
| File | Description |
| --- | --- | 
|sh_plt.m| Plot SH expansion|
|sh_fib.m| Generate spherical Fibonacci points |
|sh_plat.m| Generate platonic solid vertex points |
|sh_rand_unis.m | Uniform spherical coordinate sampling|
|sh2ambx.m| SH to AmbiX format |
|ambx2sh.m| AmbiX to SH format |

