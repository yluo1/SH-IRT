# Tutorial: Spherical Harmonic (SH) Arrays (Part I) Setup, Specifications, and Capon Beamformer Designs

Beamformer designs in the spherical harmonic (SH) domain modifies array elements' frequency responses in the far-field as to optimize spatial objectives subject to application dependent constraints. Array elements consist of transducers that convert digital and electrical signals between spatial-analog domains such as acoustic (microphone/loudspeaker), electromagnetic (antenna), seismic (seismometer) fields. Therefore, a transducer’s far-field frequency responses can be characterized by its analog/digital transfer functions over spherical coordinates, and continuously represented via SH basis functions.

In part I of this tutorial, we cover formulations of several array designs with SH transfer functions, steering and spatial covariance matrix specifications, and extensions of the classic Capon[^CAPON] beamformer in the spherical harmonics domain. The high-level topics are as follows:

* [Spherical Harmonic Arrays](#spherical-harmonic-arrays)
  * [Free-field Phased Arrays](#free-field-phased-array)
  * [Piston on Spherical Baffle Arrays](#piston-on-spherical-baffle-arrays)
  * [Empirical SH Expansion Arrays](#empirical-spherical-harmonic-arrays)
* [Minimum Variance Distortionless Response (MVDR)](#minimum-variance-distortionless-response)
  * [Covariance Matrix Condition Number Regularization](#covariance-matrix-condition-number-regularization)
  * [Spherical Harmonic MVDR Beamforming](#spherical-harmonic-mvdr-beamforming)
  


## Spherical Harmonic Arrays
 
Arrays of both analytic and fitted SH expansions of various transfer function models are assembled in the function [sh_enc_arr.m](sh/sh_enc_arr.m). Common array element topologies include uniform line, circular, and spherical placements. Custom topologies are supported via optional inputs. We cover several transfer function model cases in the following sub-sections.

### Free-field Phased Array

In a phased array, the transfer function $H_n(\omega, \bf{u} )$ at angular frequency $\omega = 2 \pi f$ between a planewave with incidence unit-direction $\bf{u} \in \mathbb{R}^{3 \times 1}$, and the $n^{th}$ array element at the cartesian coordinates $\bf{x}_n \in \mathbb{R}^{3 \times 1}$ is given by

$$ H_n(\omega, \bf{u} ) = e^{ j \omega \frac{\bf{u}^T \bf{x}_n}{c} }, $$

where $c$ is the velocity of sound (m/s). For convenience, we can omit to omit $c$ by expressing the transfer function $ H_n(\kappa, \bf{u} ) = e^{j \kappa \bf{u}^T \bf{x}_n } $ in terms of angular wavenumber $\kappa = \omega / c$. Array elements orthogonal to the incident wave direction $\bf{u}$ have zero-delay. Elements in the direction of the planewave propagation are delayed by $\lVert \bf{x}_n \rVert$. For example, the $N=3$ element circular planar array topology and each element’s frequency responses over the horizontal plane are shown below:

 | Topology|  Element 1 |  Element 2  |  Element 3 |
 | :---: |  :---: |  :---:  |  :---: |
  | <img src="./figs/figs_bf/arr/uni_circ_N_3_kappa_Inf_1.png" width="400"/> |  <img src="./figs/figs_bf/arr/uni_circ_N_3_kappa_Inf_2.png" width="400"/> |  <img src="./figs/figs_bf/arr/uni_circ_N_3_kappa_Inf_3.png" width="400"/>  |  <img src="./figs/figs_bf/arr/uni_circ_N_3_kappa_Inf_4.png" width="400"/> |

where the +x axis direction corresponds to the spherical coordinates $(\theta = 90^{\circ}, \phi = 0)$ plane wave direction. In the far-field, the magnitude is assumed to be constant for all plane wave incident directions as the distance from a planewave source to any two elements at distance $r$ and $r + \delta$ has unity relative magnitude ratio $\lim_{r \rightarrow \infty} \frac{r}{r+\delta} = 1$ at the limit. Moreover, if the element is centered at the origin, then the phase responses become constant ($0$  radians). We simply use the 0th order SH basis to encode this omni-directional directivity. Therefore, frequency responses of off-centered elements rotate only the phase ponses w.r.t. the planewave direction in the function [sh_translate.m](sh/sh_translate.m), and therefore require higher-order SH expansions.

In practice, we see that SH expansions of functions near their physical center locations require lower-order expansions. Subsequent calculations of an array’s spatial covariance, steering vectors, and downstream beamformer designs can be managed without explicit SH translation of each element’s location. The following example specifies two MVDR beamformers over the phased array:


```
% Sample usage: Design MVDR beamformer over circular omni-directivity transducer array
% Note: Uniform line and circular arrays yield for phased arrays have real responses for all and even N respectively

Fs = 16000;
M = 256;
N = 3;
radius = 0.1;
delay_dist = 1;
max_odr = 36;
mode = 'uni_circ';
freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);

theta = pi/2;
phi = 0;
kappa = inf;
mode_mvdr = 'uniform';

C_arr = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
bf_mvdr_sh(C_arr, false, freq, theta, phi, kappa, mode_mvdr, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);

[C_arr_noncenter, center] = sh_enc_arr(0, mode, N, radius, delay_dist, freq, 'noncentered', true, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
bf_mvdr_sh(C_arr_noncenter, false, freq, theta, phi, kappa, mode_mvdr, 'center', center, 'wavnum', 2 * pi * freq / 343, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);
```
where `C_arr` and  ``C_arr_noncenter`` are $36^{th}$ and $0^{th}$ max-order SH expansions coefficients respectively. Reducing the max-order expansion in the former introduces artifacts in high-frequency due to spatial aliasing.

### Piston on Spherical Baffle Arrays

While free-field phased-array elements are idealized point sources and receivers, real-world world arrays may be baffled and their elements occupy a physical space. Realistic arrays require numerical simulations or measurements to compute their transfer functions. One exception is a piston on a spherical baffle[^ANSOL] which has an analytic solution in the spherical harmonic domain that we implement in the function [sh_enc_pist_sphere.m](sh/sh_enc_pist_sphere.m). An exterior piston spans the positive half-space intersection of a great circle slice of a rigid sphere. The negative half-space acts as the spherical baffle.

Multiple pistons (2.54 cm, 4 mm diameters) can be arrayed over a spherical baffle (10 cm radius) as seen in the following loudspeaker and microphone array examples respectively:

```
Fs = 24000;
M = 512;
N = 3;
radius = 0.1;
delay_dist = 1;
freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
max_odr = 50;

% Circular large piston array on the horizontal plane
spis_arc_len = 0.0254;

[C_arr, center, h_arr_figs] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, ...
    'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);


% Circlar small piston array near the pole
spis_arc_len = 0.004;
custom_theta_phi = [deg2rad(10) * ones(1,N); deg2rad(0:(360/N):(360 - 360/N))];

[C_arr, center, h_arr_figs] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, ...
    'spis_arc_len', spis_arc_len, 'custom_theta_phi', custom_theta_phi, ...
    'enable_disp', true, 'disp_dB_lim', [-36, 12]);
```

 | Topology|  Element 1 |  Element 2  |  Element 3 |
 | :---: |  :---: |  :---:  |  :---: |
  | <img src="./figs/figs_bf/arr/spis_circ_N_3_kappa_Inf_1.png" width="400"/> |  <img src="./figs/figs_bf/arr/spis_circ_N_3_kappa_Inf_2.png" width="400"/> |  <img src="./figs/figs_bf/arr/spis_circ_N_3_kappa_Inf_3.png" width="400"/>  |  <img src="./figs/figs_bf/arr/spis_circ_N_3_kappa_Inf_4.png" width="400"/> |
  | <img src="./figs/figs_bf/arr/spis_custom_N_3_kappa_Inf_1.png" width="400"/> |  <img src="./figs/figs_bf/arr/spis_custom_N_3_kappa_Inf_2.png" width="400"/> |  <img src="./figs/figs_bf/arr/spis_custom_N_3_kappa_Inf_3.png" width="400"/>  |  <img src="./figs/figs_bf/arr/spis_custom_N_3_kappa_Inf_4.png" width="400"/> |

where larger piston sizes are more directive in high frequency. Unlike the free-field phased-arrays, high max-order SH expansions are required regardless of centering the expansion centers at the piston locations.


### Empirical Spherical Harmonic Arrays

The [SH fitting](tutorial_sh_fit.md) tutorial covers several methods for fitting SH basis functions to empirical complex frequency responses. Let us consider a sample SH fitting of a loudspeaker piston mounted on a rotationally invariant baffle that is resting on a stand, with expansion center coincident to the piston surface. We wish to array several such pistons under rotation on the horizontal plane as to steer a beam in frontal half-space. The first step is to translate the expansion center to a rotationally invariant coordinate within the baffle’s interior as in the following code snippet:
```
dat = load('dat/sample_exp_centered.mat');

is_real = false;
N_C = size(dat.C, 1);
P = sqrt(N_C) - 1;
dB_lim = [-40, 20];

% Plot SH expansion centered on piston surface
h_figs_cp = sh_plt(dat.C(:, 3:end), 'horizontal', is_real, 'hz', dat.freq(3:end), 'dB_lim', dB_lim);

% Translate the expansion center to the baffle centroid
omega = 2 * pi * dat.freq;
[C] = sh_translate(dat.C, omega, is_real, dat.loc, 'max_odr_fit', P + 1);
h_figs_cb = sh_plt(C(:, 3:end), 'horizontal', is_real, 'hz', dat.freq(3:end), 'dB_lim', dB_lim);
```
 | Expansion Center at Piston's Location|  Expansion Center at Baffle's Interior |
 | :---: |  :---: |  
  | <img src="./figs/figs_bf/arr/fit_pist_center_pist.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_center_baffle.png" width="400"/> |  

where the SH expansion centered on the piston surface is translated to its physical cartesian coordinates w.r.t. a new center of rotation within the baffle’s interior. As a result, the translated expansion’s phase responses over the spherical coordinates increase in variance as the piston’s radiating components have moved away from the new origin. The magnitude responses remain unchanged as the frequency responses are distance invariant in the far-field.

Next, we setup an array of $N=3$ pistons oriented at `[-60, 0, 60]` degrees azimuth on the horizontal plane via the `fit_rot_custom` mode in the function [sh_enc_arr.m](sh/sh_enc_arr.m):

```
% Duplicate and array three expansions rotated in the frontal half-space
delay_dist = 0;
radius = abs(dat.loc(1));
custom_theta_phi = deg2rad([90, -60; 90, 60; 90 0])';
custom_rot_zyx   = deg2rad([-60 0 0; 60 0 0; 0 0 0])';
C_arr = sh_enc_arr(P+1, 'fit_rot_custom', 3, radius, delay_dist, dat.freq, ...
    'custom_C', C,  'custom_theta_phi', custom_theta_phi,  'custom_rot_zyx', custom_rot_zyx, ...
    'enable_disp', true, 'disp_dB_lim', dB_lim);
```

 | Topology|  Element 1 |  Element 2  |  Element 3 |
 | :---: |  :---: |  :---:  |  :---: |
  | <img src="./figs/figs_bf/arr/fit_pist_arr_N_3_topology.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_arr_N_3_1.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_arr_N_3_2.png" width="400"/>  |  <img src="./figs/figs_bf/arr/fit_pist_arr_N_3_3.png" width="400"/> |

whereby rotating the SH expansions via Tait-Bryan intrinsic rotation angles in `custom_rot_zyx` coincides with rotating the +x axis facing piston on the baffle w.r.t. a common origin in the baffle’s interior. The spherical coordinates of the rotated pistons are specified in the `custom_theta_phi` variable.

Next, we cover sampling and integration methods of the array element's frequency responses over the spherical coordinates for the Capon beamformer.


## Minimum Variance Distortionless Response

Capon’s minimum variance distortionless response (MVDR) beamformer finds the weight vector $\bf{w} \in \mathbb{C}^{N \times 1} $ that minimizes the power response over a spatial domain with covariance matrix $\bf{R} \in \mathbb{C}^{N \times N}$  given by

$$ \bf{w}_{\ast} = \textrm{arg min}_{\bf{w}} \bf{w}^H \bf{R} \bf{w}, \quad s.t. \quad \bf{w}^H \bf{v} = 1, $$

where the linear equality constraint ensures that the response is undistorted w.r.t. the steering vector $\bf{v} \in \mathbb{C}^{N \times 1}$. If the steering vector is a function $\bf{v}(\bf{u})$ of an array's far-field response along a unit-direction $\bf{u}$, then we can specify a far-field spatial covariance matrix given by

$$\bf{R} = \mathbb{E}_{\bf{u} \sim f_R(\bf{u}) } \left [ \bf{v}(\bf{u})  \bf{v}^H(\bf{u})  \right ],$$

which represents the expected covariance between array element’s responses integrated over a user defined density function $`f_R(\bf{u})`$ on the unit-sphere or in spherical coordinates. Common density functions such as the squared exponential and Matérn family are covariance functions that can be expanded and mixed along the SH bases[^LUO_COV], and further manipulated via the [SH operators](tutorial_sh_ops.md). A preset of SH density function expansions are available in the function [sh_pdf_preset.m](sh/sh_pdf_preset.m), and SH densities can be fitted to empirical distributions following [SH basis fit tutorial](tutorial_sh_fit.md). Discrete densities for finite directions are also supported. 

Once the steering direction $\bf{u}$ in spherical coordinates, and density function $`f_R(\bf{u})`$ are specified, the analytic MVDR solution or minimizer yields the conjugated beamformer weights given by 

$$ \bf{w}_{\ast} = \frac{\bf{R}^{-1} \bf{v} }{ \bf{v}^H \bf{R}^{-1} \bf{v}},$$

which typically forms a hyper-directive beam in the direction of $\bf{u}$ whilst simultaneously rejecting the output domain power in the directions sampled from $`f_R(\bf{u})`$. Such beamformer designs are well studied for phased-array topologies which yield solutions with useful properties. For example, MVDR beamformer responses for lined and circular uniform spaced phase arrays are strictly real for arbitrary and even order $N$ respectively. The following MVDR beam patterns for spatial covariances sampled over a uniform density are computed in the function [plot_bf_mvdr_sh_example.m](plt/plot_bf_mvdr_sh_example.m):

```
plot_bf_mvdr_sh_example('uni_line', 3);
plot_bf_mvdr_sh_example('uni_circ', 3);
plot_bf_mvdr_sh_example('uni_circ', 4);
```

| Topology|  Beamformer Response |  Beamformer Statistics  |  
| :---: |  :---: |  :---:  | 
| <img src="./figs/figs_bf/arr/uni_line_N_3_kappa_Inf_1.png" width="320"/> |  <img src="./figs/figs_bf/arr/uni_line_N_3_kappa_Inf_5.png" width="400"/> |  <img src="./figs/figs_bf/arr/uni_line_N_3_kappa_Inf_6.png" width="400"/>
| <img src="./figs/figs_bf/arr/uni_circ_N_3_kappa_Inf_1.png" width="320"/> |  <img src="./figs/figs_bf/arr/uni_circ_N_3_kappa_Inf_5.png" width="400"/> |  <img src="./figs/figs_bf/arr/uni_circ_N_3_kappa_Inf_6.png" width="400"/>
| <img src="./figs/figs_bf/arr/uni_circ_N_4_kappa_Inf_1.png" width="320"/> |  <img src="./figs/figs_bf/arr/uni_circ_N_4_kappa_Inf_6.png" width="400"/> |  <img src="./figs/figs_bf/arr/uni_circ_N_4_kappa_Inf_7.png" width="400"/>  |

where we see that $N=3$ line array and $N=4$ circular array’s main lobe are in-phase, and the side-lobes alternate in phase. For circular uniform arrays of odd number of elements such as, the phase responses are not as regular as seen in the $N=3$ element example.

The beamformer statistics plot the directivity index $\textrm{DI}(\bf{w})$, and electrical gain $\textrm{G}(\bf{w})$  given by
$$\textrm{DI}(\bf{w}) = \frac{\lVert \bf{w}^H \bf{v} \rVert^2}{ \bf{w}^H \bf{R}_{omni} \bf{w}}, \quad
\textrm{G}(\bf{w}) =\bf{w}^H \bf{w},$$

where MVDR solutions maximize DI as the numerator term is unity due to the distortionless constraint, and the denominator term is minimized via the objective function. Examining lower frequencies where a wavelength falls below the array size, the DI plateaus as further side-lobe attenuation is inefficient as all the element responses highly covary, requiring increasing electrical gain to satisfy the on-axis distortionless constraint. The latter is typically remedied by regularizing the spatial covariance matrix via diagonal loading $\bf{R} + \sigma^2 \bf{I}$ to also penalize electrical power. Regularizing the covariance matrix introduces a necessary bias to the MVDR solution w.r.t. the original problem. We consider several regularization methods that bound the covariance matrix’s so-called condition number.

### Covariance Matrix Condition Number Regularization

Covariance matrices $\bf{R}$ for array designs are typically ill-conditioned at lower frequencies. Measurements across small element spacings yield small differences in the phase at longer wavelengths, thereby increasing its sensitivity to noise. This is quantified by a the covariance matrix’s condition number $\kappa$ or ratio of the largest and smallest eigenvalues, which becomes large and can grow to infinity. As a result, small eigenvalues explode the MVDR solution’s numerator $\bf{R}^{-1}\bf{v}$ vector when scaled by the inverse eigenvalues. Application of the Kantorovich inequality gives a conservative upper-bound[^MITTAL] on the electrical gain $`\textrm{G}(\bf{w}_{\ast})`$ of MVDR solution $`\bf{w}_*`$  following

$$\textrm{G}(\bf{w}_{\ast}) = \bf{w}_{\ast}^H \bf{w}_{\ast} = \frac{ \bf{v}^H \bf{R}^{-2} \bf{v}}{\lvert  \bf{v}^H \bf{R}^{-1} \bf{v} \rvert^2}, \qquad
\frac{1}{\bf{v}^H \bf{v}} \leq \textrm{G}(\bf{w}_{\ast}) \leq  
\frac{1}{\bf{v}^H \bf{v}} \frac{(\kappa + 1)^2}{4\kappa},$$

whereby the upper-bound of the overall array gain $`1 \leq \textrm{AG} = (\bf{v}^H \bf{v} ) \textrm{G}(\bf{w}_{\ast})`$ is ensured by solving for the minimum condition number given by

$$\kappa_{min} = 2 \textrm{AG} - 1 + \sqrt{\textrm{AG} (\textrm{AG} - 1) }.$$ 

Several regularization methods that precisely constrain the condition number $\kappa$ of the regularized covariance matrix $\hat{\bf{R}} = \bf{R}(\kappa)$  are considered in the MVDR design function [bf_mvdr.m](bf/bf_mvdr.m):

The first three regularizers (diagonal, threshold, and exponential) apply shape-preserving modifications to the spectra of $\bf{R}$.
* `diag`: 	Minimally shift all eigenvalues upwards by  $`\frac{\lambda_{max} - \lambda_{min} }{\kappa - 1}`$
* `thres`: 	Replace all eigenvalues below threshold  $`\lambda_{\kappa} =  \frac{\lambda_{max} }{ \kappa }`$ with $`\lambda_{\kappa}`$ 
* `exp` 	Exponentiate all eigenvalues by $`\frac{\log \kappa }{ \log \lambda_{max} - \log \lambda_{min} }`$

The remaining three regularizers (sum-of-squared, sum-of-absolute, max absolute spectral errors) finds new eigenvalues that minimize various matrix norms via linear and quadratic programming:
* `sse`: Minimize Frobenius norm given by $`\min_{\hat{\lambda}}  \sum_{n=1}^N \lvert \hat{\lambda}_n - \lambda_n \rvert^2`$
* `ase`: Minimize nuclear norm given by  $`\min_{\hat{\lambda}}  \sum_{n=1}^N \lvert \hat{\lambda}_n - \lambda_n \rvert`$
* `minimax`: Minimize infinity norm given by $`\textrm{min max}_{\hat{\lambda}}  \sum_{n=1}^N \lvert \hat{\lambda}_n - \lambda_n \rvert`$

Larger $\kappa$ reduces the amount of regularization at the cost of larger electrical and array gains. Smaller $\kappa$ introduces bias but reduces the latter. Therefore, it is desirable to minimize the amount of regularization (maximize $\kappa$) subject to a design’s array gain upper-bound $`\textrm{AG} \leq \textrm{AG}_{max}`$. For example, we can initialize $\kappa$ to $`\kappa_{min}`$, compute one of the several regularized MVDR solution weights $`\bf{w}_{\ast}`$, increase $\kappa$ and then re-compute the weights until $`\textrm{AG} \rightarrow \textrm{AG}_{max}`$. This is implemented in the spherical harmonic MVDR function [bf_mvdr_sh.m](bf/bf_mvdr_sh.m) where it’s desirable that the  solutions also smoothly vary over frequency. The following script computes beamformer weights for the previous $N=4$ element uniform circular phased-array with $6$ and $24$ dB max array gain power under `diag` regularization respectively:

```
plot_bf_mvdr_sh_example('uni_circ', 4, inf, 'diag', 'max_arr_gain_pow_dB', 6);
plot_bf_mvdr_sh_example('uni_circ', 4, inf, 'diag', 'max_arr_gain_pow_dB', 24);
```
| Max Array Gain Power (dB)|  Beamformer Response |  Beamformer Statistics  |  
| :---: |  :---: |  :---:  | 
| 6|  <img src="./figs/figs_bf/arr/uni_circ_N_4_kappa_Inf_arrgain_6_6.png" width="400"/> |  <img src="./figs/figs_bf/arr/uni_circ_N_4_kappa_Inf_arrgain_6_7.png" width="400"/>  |
| 24|  <img src="./figs/figs_bf/arr/uni_circ_N_4_kappa_Inf_arrgain_24_6.png" width="400"/> |  <img src="./figs/figs_bf/arr/uni_circ_N_4_kappa_Inf_arrgain_24_7.png" width="400"/>  |

whereby decreasing the max array gain flattens the low frequency electrical gains, and reverts the directivity to being omni-directional compared to constant directivity in the unconstrained array gain solution. At 6 dB array gain upper-bound, the beam pattern widens below 550 Hz; increasing upper-bound to 24 dB lowers the widening frequency to 170 Hz. The other non-diagonal regularization methods also produce solutions of similar DI with varied smoothness over frequency:

```
enable_export = false;
plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'diag', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'thres', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'exp', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'sse', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'ase', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'minimax', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
```

| `diag` |  `thres` |  `exp`  | 
| :---: |  :---: |  :---:  |
|  <img src="./figs/figs_bf/arr/uni_circ_N_5_diag_kappa_Inf_arrgain_6_8.png" width="400"/>  | <img src="./figs/figs_bf/arr/uni_circ_N_5_thres_kappa_Inf_arrgain_6_8.png" width="400"/>  | <img src="./figs/figs_bf/arr/uni_circ_N_5_exp_kappa_Inf_arrgain_6_8.png" width="400"/>  |

| `sse`  |  `ase`  |  `minimax` |
| :---: |  :---: | :---: |
<img src="./figs/figs_bf/arr/uni_circ_N_5_sse_kappa_Inf_arrgain_6_8.png" width="400"/>  | <img src="./figs/figs_bf/arr/uni_circ_N_5_ase_kappa_Inf_arrgain_6_8.png" width="400"/>  | <img src="./figs/figs_bf/arr/uni_circ_N_5_minimax_kappa_Inf_arrgain_6_8.png" width="400"/>  | 

* `diag`: Good compromise between low bias (large DI) and smooth electrical channel powers.
* `exp`: Produces high bias (lowest DI) and smoothest electrical channel power until the lowest frequencies where the smallest eigenvalues are sensitive to noise.
*  `thres`, `sse`, `ase`, `minimax`: Similar DI trend to `diag` but exhibits a visible kink. Electrical gains are piece-wise smooth with discontinuous regimes over frequency. `minimax` exaggerates the kink and is therefore least-smooth.

### Spherical Harmonic MVDR Beamforming

In far-field MVDR beamformer designs, the steering vectors are normally sampled from the array element’s far-field responses in spherical coordinate directions $(\theta, \phi)$, which are already parameterized in the output matrix `C` of the SH expansions of [sh_enc_arr.m](sh/sh_enc_arr.m). The far-field SH MVDR function [bf_mvdr_sh.m](bf/bf_mvdr_sh.m) solves for the beamformer weights that steers a beam in the aforementioned `theta`, `phi` input direction, while minimizing the output power in the domain as defined by the spatial covariance matrix $\bf{R}$. The latter constructs `R` via different sampling schemes for integration:

* `uniform`: Uniform density over spherical coordinates maximizes DI.
* `density`: SH density expansion `options.R_pdf` over spherical coordinates  attenuates broad spatial direction regions (e.g. side wall, ceiling, floor, back-wall reflection windows )
* `direction`: Specifies a discrete list of spherical coordinates in `options.R_theta_phi` for steering sharp nulls. Requires finite max-condition constraint `kappa` or array gain if the number of directions is less than $N$.
* `custom`: Specifies non far-field covariance in `options.R_custom` (e.g. near-field covariance, electrical power)

and where the input maximum condition number `kappa` regularization of `R`, and optional array gain upper-bound `options.max_arr_gain` control for excessive electrical and white-noise gains.

Let us consider different spatial covariances modes for MVDR designs on an $N=3$ element micro [piston array](#piston-on-spherical-baffle-arrays) near the pole of the spherical baffle:

| Sampling Method |  Beamformer Response |  Beamformer Statistics  |  
| :---: |  :---: |  :---:  | 
| `uniform` (on-axis steering direction): Smooth broadband DI ~ 6 dB main-lobe visible below 10 kHz | <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_onaxis_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_onaxis_bs.png" width="400"/>  |
| `uniform` (60 degree azimuth steering direction): Worst-case steering angle between two transducers slightly degrades DI by 1/2 dB. Beam pattern rotates that of on-axis steering angle. | <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_60deg_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_60deg_bs.png" width="400"/>  |
| `density` (reject broad 90 degree azimuth window): Shallow null emerges near 90 degree azimuth| <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_estbrd_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_estbrd_bs.png" width="400"/>  |
| `density` (narrow 90 degree azimuth window): Dual deep nulls emerge near 90 degree azimuth. | <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_estnrw_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_estnrw_bs.png" width="400"/>  |
| `density` (narrow cross-talk window): Singular null steered towards the steering direction. Narrow-band lobes emerge in mid-band.| <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_contra_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_contra_bs.png" width="400"/>  |
| `direction` (90 degree azimuth, heavy regularization `kappa = 1e3`): Steers a deep and singular null in the mid-high band. Low-frequencies revert to omni-directional with 10 dB reduction to electrical gain. | <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_arrg1e2_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_arrg1e2_bs.png" width="400"/>  |
| `direction` (90 degree azimuth, light regularization `kappa = 1e8`): Steers a   broad-band deep and singular null. | <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_arrg1e8_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_arrg1e8_bs.png" width="400"/>  |


MVDR beamformers are suitable for hyper-directive when the steering vectors' frequency responses are relatively flat and omni-directional. Otherwise, satisfying the distortionless response constraint requires excessive electrical gain, which boosts any input white-noise input in the signal. For example, a piston’s output power rolls-off outside its operating frequency range, and increases in directivity for wavelengths less than its size. MVDR beamformers designed for [loudspeaker pistons](#empirical-spherical-harmonic-arrays) require excessive electrical power gain below its operating frequency as seen in the aforementioned $N=3$ array design. We can indirectly constrain the optimization’s electrical gain by imposing a maximum array gain and pre-normalizing the steering vectors’ power to unity:

```
% Design MVDR beam with and without array gain constraint and steering vector normalization
theta_steer = deg2rad(90);
phi_steer = 0;
kappa = inf;
mode_mvdr = 'uniform';

[~, ~, ~, ~, h_figs_MVDR_unconstr] = bf_mvdr_sh(C_arr(:, :, 4:end), false, dat.freq(4:end), theta_steer, phi_steer, kappa, mode_mvdr, ...
    'enable_disp', true, 'disp_dB_lim', dB_lim);

[~, ~, ~, ~, h_figs_MVDR_constr] = bf_mvdr_sh(C_arr(:, :, 4:end), false, dat.freq(4:end), theta_steer, phi_steer, kappa, mode_mvdr, ...
    'max_arr_gain', db2pow(12), 'norm_steer_vec', true, ...
    'enable_disp', true, 'disp_dB_lim', dB_lim);
```
 | MVDR|  Beamformer Response | Beamformer Statistics|
 | :---: |  :---: |   :---: |  
| Unconstrained | <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_unconstr_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_unconstr_bs.png" width="400"/> | 
| Max Array Gain (12 dB), normalized steering vector | <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_constr_bp.png" width="400"/> |  <img src="./figs/figs_bf/arr/fit_pist_MVDR_N_3_constr_bs.png" width="400"/> |

However, the approach doesn’t lower the total electrical gain (including steering vector normalization) required to satisfy the distortionless response constraint; excessive makeup gains below the piston’s operating range (4 kHz) must compensate for its natural frequency roll-off behavior. 

Above the roll-off from 4-13 kHz, we see a narrowing main lobe flanked by increasing number of undesirable side-lobes generated from the mixture of center and phantom-center beams. This combination yields a more directive beam pattern, but at the cost of greater variance in the off-axis responses. At 13 kHz, we observe a strong narrow-band break-up mode where a sharp notch in the axial frequency response causes the DI to become negative (DI of Ch. 3). The beamformer design mitigates this deficiency by forming a narrow-band phantom-center beam with the left and right transducers (electrical gain peak at 13 kHz for Chs. 1-2). Above 13 kHz, the beamformer rolls-off the left and right transducer responses as the phantom-center contributes less to beampattern’s main lobe than the single use of the center transducer.

In the next tutorial, we cover heterogeneous arrays with mixed operating frequency ranges, and generalized Rayleigh quotient beamformer designs.


[^CAPON]: J. Capon, "High-resolution frequency-wavenumber spectrum analysis," in Proceedings of the IEEE, vol. 57, no. 8, pp. 1408-1418, Aug. 1969, doi: 10.1109/PROC.1969.7278.

[^MITTAL]: Mittal, M., Corey, R. M., Buck, J. R., & Singer, A. C. (2026). Adaptive Diagonal Loading for Norm Constrained Beamforming. arXiv preprint arXiv:2605.04342.

[^ANSOL]: Advanced Numerical Solutions. (n.d.). Exterior piston baffled in a sphere radiation problem (Coustyx Validation Manual). Ansol. ansol.us
http://ansol.us/Products/Coustyx/Validation/Indirect/BoundaryExcitation/Sphere/SphericalBaffle/Downloads/dataset_description.pdf

[^LUO_COV]: Luo, Y. (2021). [Spherical Harmonic Covariance and Magnitude Function Encodings for Beamformer Design](https://link.springer.com/article/10.1186/s13636-021-00230-7). EURASIP Journal on Audio, Speech, and Music Processing, 2021(1), 41.