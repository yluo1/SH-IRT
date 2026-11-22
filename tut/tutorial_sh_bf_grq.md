# Tutorial: Spherical Harmonic (SH) Arrays (Part II) Generalized Rayleigh Quotient Beamformer Designs for Heterogeneous Arrays

In the previous tutorial on [Spherical Harmonic Array Designs (Part I)](tutorial_sh_bf_arr.md), we introduced our setup of spherical harmonic (SH) arrays, specification and regularization methods of spatial covariance matrices, and SH integration with classic Capon MVDR beamformer designs. This allowed us to array homogeneous transducers, generate and steer hyper-directive beams with and without null placements, and control the array-gain at a cost to the array’s directivity via conditioning the spatial covariance matrix. As a result, design choices centered around trade-offs between input side constraints and specifications (white-noise-gain, electrical gain, array-gain, spatial covariance’s density function), and the judgement of the output beam pattern and its directivity.

In part II of the array design tutorial, we consider beamformer designs that directly control for an array's directivity using constrained optimization techniques from the following work:

> Luo, Y. (2024, August). [Constant Directivity Loudspeaker Beamforming](https://ieeexplore.ieee.org/abstract/document/10715159). In 2024 32nd European Signal Processing Conference (EUSIPCO) (pp. 246-250). IEEE.

which first extends an array’s directivity index using a generalized Rayleigh quotient (GRQ), and then constrains the GRQ while optimizing for different design objectives. We then introduce a channel-frequency dependent regularization method of the GRQ to support heterogeneous transducer arrays designs for loudspeaker arrays. The latter is a common fixture of multi-way systems where different transducer types have varied operating frequency ranges (e.g. woofer, midrange, tweeter), and are arrayed to form virtual full-range loudspeakers. The high-level topics are as follows:

* [Generalized Rayleigh Quotient based Directivity](#generalized-rayleigh-quotient-based-directivity)
  * [Loudspeaker Directivity and Array Example](#loudspeaker-directivity-and-array-example)
  * [Channel and Frequency Weighted Directivity](#channel-and-frequency-weighted-directivity)
* [Constant Directivity Beamformer Designs](#constant-directivity-beamformer-designs)
  * [Minimum Spectral Norm](#minimum-spectral-norm)
  * [Minimum Norm Projection](#minimum-norm-projection)
  * [Maximum Sensitivity](#maximum-sensitivity)
  * [Maximum Efficiency](#maximum-efficiency)

## Generalized Rayleigh Quotient based Directivity

The Rayleigh quotient (RQ) $`\frac{\bf{x}^H \bf{M} \bf{x}}{ \bf{x}^H  \bf{x}}`$ is a recurring quantity in array designs that measure an array’s efficiency in terms of its output / input power ratio. For an $N$ element (channel) array, we wish to determine its weight vector $\bf{x} \in \mathbb{C}^{N \times 1}$ (beamforming weights), that satisfy some measure of efficiency as defined by the design choice of the measurement's covariance matrix $`\bf{M} \in \mathbb{C}^{N \times N}`$, as well as its bounds (extrema). For example, hyper-directive Capon beamformers for the far-field would specify covariance matrix $\bf{M}$ by integrating the steering-vectors uniformly over spherical coordinates. In loudspeaker array designs, the beamformer's RQ efficiency can be directly interpreted in terms of the acoustic sound power level per unit of input power in watts. Its RQ bounds can be determined from the spectra of the covariance matrix $\bf{M}$ under the following change-of-variables:

$$\textrm{RQ}(\bf{M}, \bf{x}) = \frac{\bf{x}^H \bf{M} \bf{x}}{\bf{x}^H \bf{x}} 
= \frac{ \bf{y}^H \bf{\Lambda} \bf{y}}{\bf{y}^H \bf{y}} 
= \frac{\sum_{n=1}^N \lambda_n \lvert y_n \rvert ^2 }{ \sum_{n=1}^N  \lvert y_n \rvert ^2 }
 , \qquad \bf{x} = \bf{V}^H \bf{y}, $$

where $`\bf{M} = \bf{V} \bf{\Lambda} \bf{V}^H \in \mathbb{C}^{N \times N}`$ is diagonalized by eigenvectors in the column matrix $`\bf{V} \in \mathbb{C}^{N \times N}`$ and eigenvalues in $`\bf{\Lambda} = \textrm{diag}(\bf\lambda{}) \in \mathbb{R}^{N \times N}`$. The change-of-variables to vector $\bf{y}$ allows us to factor the numerator and then cancel the denominator after substituting the smallest and largest eigenvalues $`\lambda_{min}, \lambda_{max}`$ into $`\lambda_n`$. This yields the RQ’s bounds given by

$$\lambda_{min} \leq \textrm{RQ}(\bf{M}, \bf{x}) \leq  \lambda_{max}, $$

where the RQ’s extrema are obtained when $\bf{x}$ are eigenvectors corresponding to the smallest or largest eigenvalues upto arbitrary phase rotation. 

Next, we introduce the so-called generalized Rayleigh quotient (GRQ) $`\frac{\bf{x}^H \bf{A} \bf{x}}{ \bf{x}^H \bf{R} \bf{x}}`$, which evaluates the array’s efficiency between measurement across two domains in the numerator and denominator by their covariance matrices $`\bf{A}, \bf{R} \in \mathbb{C}^{N \times N}`$ . This allows us to generalize far-field directivity by specifying both domains to be the acoustic sound power level averaged over two different spherical coordinate fields (density functions). Thus, GRQ encompasses both measures of conventional directivity and acoustic cancellation efficiency depending on the choice of the density functions. For positive-definite $\bf{R}$, the GRQ’s upper bounds can be determined by transforming GRQ to the standard RQ for via a change-of-variables given by

$$\textrm{GRQ}(\bf{A}, \bf{R}, \bf{w}) 
= \frac{\bf{w}^H \bf{A} \bf{w}}{\bf{w}^H \bf{R} \bf{w}} 
= \frac{\bf{x}^H \bf{M} \bf{x}  }{\bf{x}^H \bf{x}} = \textrm{RQ}(\bf{M}, \bf{x}), \qquad  
\bf{M} = \bf{R}^{-\frac{1}{2} } \bf{A} \bf{R}^{-\frac{H}{2}}, \quad
\bf{w} = \bf{R}^{-\frac{H}{2}} \bf{x},$$

where $\bf{R} = \bf{R}^{\frac{1}{2}} \bf{R}^{\frac{H}{2}}$ is the Cholesky decomposition. Therefore, if $\bf{x}$ is an eigenvector of $\bf{M}$, then we can show $\bf{w}$ to be the generalized eigenvector of $\bf{A}, \bf{R}$ corresponding to the same eigenvalue $\lambda$ as follows:

$$\bf{M} \bf{x} = \bf{R}^{-\frac{1}{2} } \bf{A} \bf{R}^{-\frac{H}{2}} \bf{x} = \lambda \bf{x} \quad \Rightarrow \quad
\bf{A}  \bf{R}^{-\frac{H}{2}} \bf{x} =  \lambda  \bf{R}^{\frac{1}{2}} \bf{x} \quad \Rightarrow \quad
\bf{A} \bf{w} = \lambda \bf{R} \bf{w},
$$

which allows the GRQ extrema and spectra to be efficiently found via power iterations. Note that we also require $\bf{R}$ to also be well-conditioned, which becomes relevant to regularizing heterogenous arrays where channel-frequency dependent output power can vary considerably. In the remaining tutorial and subsequent far-field GRQ array designs, we restrict matrices $\bf{A}$ and $\bf{R}$ to spatial covariances between far-field steering vectors $\bf{a}(\bf{\Omega}) \in \mathbb{C}^{N \times 1}$ in the direction of spherical coordinates $\Omega = (\textrm{co-latitude } \theta, \textrm{azimuth } \phi)$ given by

$$\bf{A} = \mathbb{E}_{\bf{\Omega} \sim f_A(\bf{\Omega}) } \left [ \bf{a}(\bf{\Omega})  \bf{a}^H(\bf{\Omega})  \right ], \quad
\bf{R} = \mathbb{E}_{\bf{\Omega} \sim f_R(\bf{\Omega}) } \left [ \bf{a}(\bf{\Omega})  \bf{a}^H(\bf{\Omega})  \right ],$$

which are integrated over spatial density functions $f_A(\bf{\Omega})$, and $f_R(\bf{\Omega})$ respectively. The former defines the directions of `acceptance`, and the latter gives the directions of `rejection`, which interprets the GRQ to be a ratio of the array’s output power in directions that boost and attenuate respectively. Thus, we obtain a  flexible index for directivity that is defined by $f_A(\bf{\Omega})$, and $f_R(\bf{\Omega})$, which can vary across domain-specific contexts. Next, we give an example of a directivity centric loudspeaker array design for pistons mounted on a shared spherical baffle:

### Loudspeaker Directivity and Array Example

In loudspeaker transducer and array characterizations, the anechoic acoustic power responses are commonly averaged across different spatial windows[^CTA_2034_A] in standard `spinorama` calculations:

<img align="right" src="../figs/figs_grd/grd_cea2034a_ER.png" width="250"/>
<img align="right" src="../figs/figs_grd/grd_cea2034a_lis.png" width="250"/>

* On-axis: Single direction $0^{\circ}$ horizontal and vertical 
* Listening window: Between ± 30º horizontal in 10º increments, and ± 10º vertical
* Early-reflection window:
  * Floor Bounce: 20º, 30º, 40º down
  * Ceiling Bounce: 40º, 50º, 60º up
  * Front Wall Bounce: 0º, ± 10º, ± 20º, ± 30º horizontal
  * Side Wall Bounces: ± 40º, ± 50º, ± 60º, ± 70º, ± 80º horizontal
  * Rear Wall Bounces: 90º to 270º horizontal in 10º increments
* Sound power window: All directions

Two directivity indices are specified:
* Sound Power Directivity Index (SPDI): Listening window / Sound power window
* Early Reflections Directivity Index (ERDI): Listening window / Early-reflection window

which are pertinent to predicting and equalizing a loudspeaker’s acoustic responses intended for a listener in the listening window within a room; wall-reflections produce multiple sound propagation paths that contribute to impulse responses. SPDI approximates the relative contribution of energy between the direct and reverberant paths of acoustic impulse response; low-values predict degradation to clarity and a recessed sound-stage due to acoustic masking by the reverberation tail. ERDI approximates the contribution between direct and early reflection paths; low-values predict stronger comb-filter colorations due to strong early reflections.

Both directivity indexes can be specified by the GRQ of discrete density functions $`f_A(\bf{\Omega})`$, and $`f_R(\bf{\Omega})`$ of directions from the function [sh_grd_cea.m](../sh/sh_grd_cea.m) belonging to each spatial window. If more measurements are available, or the far-field responses are already parameterized via the spherical harmonics, then we can increase the accuracy of the indices when the spatial window and density functions are continuous. This is useful in high-frequency and ultra-sound ranges where the wavelengths are small, and for integrating over the resulting spatially narrow acoustic lobes. The latter also occurs in beam patterns of loudspeaker arrays, which introduces side-lobes below wavelengths that fit between the transducer spacings of an array. 

Consider the following piston array specified in the function [sh_enc_arr.m](../sh/sh_enc_arr.m):

```
Fs = 24000;
M = 512;
N = 4;
radius = 0.0254 * 3; % 3 inch
delay_dist = 1;
freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
is_real = false;

spis_arc_len = 0.0254;
max_odr = 30;

% Generate array
[C_arr, ~, h_arr_figs] = sh_enc_arr(max_odr, 'spis_circ', N, radius, delay_dist, freq, ...
    'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);

% Couple transducers 2, 4
C_arr = cat(2, C_arr(:,1,:), C_arr(:,2,:)+C_arr(:,4,:),  C_arr(:,3,:));
```

Four pistons are equally spaced in a circular array on a shared spherical baffle, and facing outwards. Their frequency responses over the horizontal plane are rotated. Pistons 2 and 4 are coupled into a single element in the subsequent center-channel beamformer design:

 | Topology|  Element 1 |  Element 2  |  Element 3 |  Element 4 |
 | :---: |  :---: |  :---:  |  :---: |  :---: |
 | <img src="../figs/figs_bf/grq_max/arr_1.png" width="400"/>  | <img src="../figs/figs_bf/grq_max/arr_2.png" width="400"/>   | <img src="../figs/figs_bf/grq_max/arr_3.png" width="400"/> | <img src="../figs/figs_bf/grq_max/arr_4.png" width="400"/> | <img src="../figs/figs_bf/grq_max/arr_5.png" width="400"/>|

We can specify continuous density functions over the spherical coordinates from a list of presets in the function [sh_pdf_preset](../sh/sh_pdf_preset.m), which are represented by their SH expansion coefficients `A_pdf`, and `R_pdf` for computing the acceptance and rejection covariance respectively. In this example, the acceptance density function is a forward facing listening window, and the rejection density function broadly covers side-wall reflections shown below:

```
% Generate accept density for spatial covariance A, R
[A_pdf, h_A_pdf_fig] = sh_pdf_preset('FwdCenter', 'SqExp', 8, is_real, 'ell', 0.5, 'enable_disp', true);
[R_pdf, h_R_pdf_fig] = sh_pdf_preset('FwdCenter', 'SqExpLeftRight', 8, is_real, 'ell', 0.75, 'enable_disp', true);
```

| Acceptance Density |  Rejection Density | 
| :---: |  :---: | 
| <img src="../figs/figs_bf/grq_max/A_pdf.png" width="400"/>  | <img src="../figs/figs_bf/grq_max/R_pdf.png" width="400"/>   | 

The beamformer's input density functions, target GRQ, and objective function are then specified in the design function [bf_grq_sh.m](../bf/bf_grq_sh.m). Setting the target GRQ to be `inf` constrains the beamformer’s weights to that of the GRQ maximizer, which maximally boosts the acoustic power in the listening window, and attenuates the acoustic power along side-wall directions.

```
varargin = {'A_pdf', A_pdf, 'R_pdf', R_pdf, 'enable_disp', true, 'disp_dB_lim', [-48, inf], 'disp_tau_dB_lim', [-inf, inf], 'disp_DI_dB_lim', [-inf, inf], 'disp_EP_dB_lim', [-48, inf], 'disp_spin', true};
[~, f_fig_bf] = bf_grq_sh(C_arr, false, freq, db2pow(inf), 'sensitivity', varargin{:});
```

The resulting beam pattern yields a regular main lobe along the listening window from low to mid-frequencies, and attenuated the side-lobes via acoustic null placement towards sides of the array via the non-forward facing pistons. The beam pattern then flares out in the upper-mid frequencies before narrowing in high-frequency due to the increased directivity of forward-facing piston's responses as shown below:

| GRQ Range |  Beamformer Response  | 
| :---: |  :---: | 
| <img src="../figs/figs_bf/grq_max/bf_sensitivity_tau_Inf_tau_tgt.png" width="400"/>  | <img src="../figs/figs_bf/grq_max/bf_sensitivity_tau_Inf_bp.png" width="400"/>   | 

Note that the beamformer solver’s method is set to `sensitivity`, which additionally constrains the on-axis response to be distortionless. Thus, we can evaluate the tradeoff between the directivity index and electrical gain similar to the Capon’s beamformer designs as follows:

| Beamformer Statistics |  Spinorama | 
| :---: |  :---: | 
| <img src="../figs/figs_bf/grq_max/bf_sensitivity_tau_Inf_bs.png" width="400"/>   | <img src="../figs/figs_bf/grq_max/bf_sensitivity_tau_Inf_spin.png" width="400"/>   | 

The beamformer's directivity index is smooth and near constant at low-frequency, briefly drops in the mid-frequencies, and then rapidly increasing in high-frequency. The GRQ yields a flatter index as it averages over broader accept region, and narrower reject region. We can interpret the relation between electrical gain, directivity, and the spinorama’s window power averages as follows:

* The low-frequency regime used all three channels and required more electrical power to both place spatial side-nulls and equalize to a distortionless response at a cost to acoustic / electric domain efficiency. The mid-frequency regime transitions and briefly swaps between the use of the side and rear-firing pistons. The high-frequency regime switches to the dominant use of the forward facing piston by attenuating the non-forward pistons as indicated in the channel dependent electrical gains.

* Each piston’s directivity index w.r.t. the array’s steering direction can be sorted: Piston 1 is front-facing and therefore the most efficient in both GRQ and acoustic/electric domains. Piston 2 is coupled with 4 is less-efficient as it’s side-firing, but can be used to narrow the main lobe in low-mid frequencies via acoustic null steering. Piston 3 is rear-facing, least GRQ efficient in high frequencies, and only couples with the forward piston in low-mid frequencies to boost the sound-power level.

* The sound power levels in the spinorama heavily attenuates in high-frequency, and quickly rolls-off off-axis due to the increased directivity. In the mid-frequencies, sound power and early reflections increases while listening window power decreases due to the boosts in the forward and rear pistons. This results in a temporary dip to the SPDI.

In practice, this homogeneous piston array would be too directive in high-frequency for a center channel beamformer design. The transducer spacings are also too disparately angled to optimize for custom spatial windows used in specialized beamformer designs such as cross-talk cancellation, and ceiling/side wall surround beams. Many loudspeaker arrays also use heterogeneous pistons such as tweeters, midranges, and woofers that have varied directivity and disparate operating frequencies. Next, we consider the application of constrained GRQ beamformer designs.


### Channel and Frequency Weighted Directivity

We can regularize the constant GRQ constraint to penalize channel-dependent output power. Consider the following channel-dependent diagonal loading matrix $`\bf{\Sigma} = \textrm{diag} \left ( [\sigma_1, … , \sigma_N] \right )`$, and the channel output power diagonal matrix $`\bf{\Gamma} = \textrm{diag} \left ( [R_{11}, … , R_{nn}, ... R_{NN}] \right )`$  of the rejection covariance matrix $\bf{R}$. Loading the rejection covariance’s diagonal by $`\bf{R} + \bf{\Gamma} \bf{\Sigma}`$ penalizes the uncorrelated channel-dependent output power $\bf{\Gamma}$ per channel loadings $`\sigma_n`$ where $`0 \leq \sigma_n \leq \infty`$. Thus, the penalized rejection covariance matrix is ill-conditioned when only some channels are heavily loaded as is typical of loudspeaker channels with disparate operating frequency ranges. Fortunately, we can improve the matrix's conditioning when evaluating the GRQ under a change-of-variables. Substituting the penalized rejection covariance into the GRQ yields the following generalized Rayleigh penalized quotient (GRPQ)[^LUO_CDB] following

$$\textrm{GRQ}(\bf{A}, \bf{R} + \bf{\Gamma} \bf{\Sigma} , \bf{w}) =
\textrm{GRQ}( \bf{\Xi} \bf{A} \bf{\Xi}, \hat{\bf{R}}, \bf{y}), \qquad 
\bf{\Xi} = \left ( \bf{\Sigma} + \bf{I} \right )^{-\frac{1}{2}} = \textrm{diag}\left ( [\xi_1, … , \xi_N ] \right ), \quad
0 \leq \xi_n \leq 1,
$$

where $`\bf{y} = \bf{\Xi}^{-1} \bf{w}`$ is a change-of-variables, and the transformed rejection matrix $`\hat{\bf{R}} = \bf{\Gamma} + \bf{\Xi} (\bf{R} - \bf{\Gamma} ) \bf{\Xi}`$ attenuates the non-diagonal entries of $\bf{R}$ by $`\bf{\Xi}^2`$ and is therefore well-conditioned when computing the generalized eigenvector and eigenvalues. The transformed diagonal weighting matrix $\bf{\Xi}$ contain the operating channel weights $`\xi_n`$ and are bounded below unity where $`\xi_n =  1 \Longleftrightarrow \sigma_n = 0`$ indicates active (non-penalized) channel, and $`\xi_n =  0 \Longleftrightarrow \sigma_n = \infty`$  indicates inactive (penalized) channel. Thus, we can parameterize $`\xi_n`$ by the magnitude responses of smooth (non-ripple) filters such as Butterworth and Linkwitz-Riley[^LINKWITZ] crossovers. Lastly, the GRPQ can be shown to be bounded above by the GRQ following the relation

$$\textrm{GRQ}( \bf{\Xi} \bf{A} \bf{\Xi}, \hat{\bf{R}}, \bf{y}) = 
\textrm{GRQ}(\bf{A}, \bf{R}, \bf{w} ) \textrm{GRQ}(\bf{\Xi} \bf{R} \bf{\Xi}, \hat{\bf{R}}, \bf{y} ),
$$

where the right-most GRQ term is bounded below unity. A target GRQ must therefore reside below that of the GRPQ's extrema. Thus, it is useful to differentiate GRPQ from GRQ of the unweighted spatial covariances by referring to the latter by its alternative name given by generalized directivity index (GDI), which contains the directivity index[^BARANEK] derived as follows:

For spatial covariances matrices $\bf{A}(\omega)$, and $\bf{R}(\omega)$ across angular frequencies $\omega$, let GDI be the GRQ $`\tau(\omega) = \textrm{GRQ}(\bf{A}(\omega), \bf{R}(\omega), \bf{w})`$. The formulation subsumes the conventional directivity index definition, where $`\bf{A} = \bf{a}(\Omega) \bf{a}^H(\Omega)`$ is a rank-1 outer-product of the steering vector $\bf{a}(\Omega)$, and $\bf{R}$ is the covariance between far-field steering vectors sampled over a uniform density of spherical coordinates $\Omega$. Both discrete and continuous density functions for constructing $\bf{A}$ and $\bf{R}$ are supported in the design function [bf_grq_sh.m](../bf/bf_grq_sh.m), which follows the previous [SH-MVDR integration](../tut/tutorial_sh_bf_arr.md) tutorial. 

Let us now generate sample operating channel weights for the previous array with partial overlapping operating frequencies:
```
[wt, h_fig_wt] = ft_freq_wt_chs(freq, [4000 8000; 100, 10000; 0 1000], [24 0; 12 6; 24 24], 'enable_disp', true, 'disp_dB_lim', [-48, 6], 'soft_knee', true, 'norm_sum', false);
```
|   Description  |  Operating Frequency Weights| 
| :---: |  :---: | 
| Ch 1: Tweeter (above 4 kHz) <br> Ch 2: Tapered midrange (100 Hz to 4000 Hz) <br> Ch 3: Bass + lower midrange (below 1 kHz) |<img src="../figs/figs_bf/grq_max_wt/freq_wt.png" width="400"/> | 

The maximizer of the GRPQ and its GDI are computed, and its beam pattern is shown:

```    
varargin = {'A_pdf', A_pdf, 'R_pdf', R_pdf, 'enable_disp', true, 'disp_dB_lim', [-48, inf], 'disp_tau_dB_lim', [-36, inf], 'disp_DI_dB_lim', [-inf, inf], 'disp_EP_dB_lim', [-48, inf], 'disp_spin', true};
[~, f_fig_bf] = bf_grq_sh(C_arr, false, freq, db2pow(inf), mode_sol, varargin{:}, 'wt', wt);
```

|   GRQ Range  |  Beamformer Response| 
| :---: |  :---: | 
<img src="../figs/figs_bf/grq_max_wt/bf_sensitivity_tau_Inf_wt_tau_tgt.png" width="400"/>   | <img src="../figs/figs_bf/grq_max_wt/bf_sensitivity_tau_Inf_wt_bp.png" width="400"/>   |

where the frequency weighting introduces the gap between the maximum GDI and the GDI of GRPQ's maximizer. Moreover, the GDI of GRPQ’s minimizer can exceed that of the minimizer due to the asymmetry of the frequency weighting applied to $\bf{R}$ and $\bf{A}$. The maximizer’s beam pattern exhibits two distinct regions between the mid-to-high frequencies due to the stepper frequency weighting roll-off and higher natural directivity of the third channel.

| Beamformer Statistics |  Spinorama | 
| :---: |  :---: | 
| <img src="../figs/figs_bf/grq_max_wt/bf_sensitivity_tau_Inf_wt_bs.png" width="400"/>   | <img src="../figs/figs_bf/grq_max_wt/bf_sensitivity_tau_Inf_wt_spin.png" width="400"/>   | 

The resulting electrical power of the beamformer weights trend towards that of the input frequency weighting similar to a traditional crossover network. 

## Constant Directivity Beamformer Designs

In broadband beamforming designs, it is often useful to enforce a target directivity (CD) $\tau(\omega)$ across frequency so that its beam pattern is regular. This can be specified via two quadratic equality constraints of the RQ's numerator and denominator respectively at each frequency in beamformer optimization:

$$\tau = \textrm{RQ}(\bf{M}, \bf{x}) \quad \Rightarrow \quad \bf{y}^H \bf{\Lambda} \bf{y} = \tau, \quad \bf{y}^H \bf{y} = 1,$$

or alternatively as a single homogeneous quadratic constraint given by

$$\tau = \textrm{GRQ}(\bf{A}, \bf{R}, \bf{w}) \quad \Rightarrow \quad \bf{w}^H \bf{D} \bf{w} = 0,$$

where Hermitian matrix $`\bf{D} = \bf{A} - \tau \bf{R} \in \mathbb{C}^{N \times N}`$ is indefinite for $`\lambda_{min} < \tau < \lambda_{max}`$. 

In practice, the GDI targets $\tau(\omega)$ must fit between the generalized eigenvalue's extrema over frequency. A smooth $\tau(\omega)$ is a necessary but not sufficient condition for a smooth set of beamformer weights $\bf{w}(\omega)$. When $\tau$ is neither extrema, multiple $\bf{y}$ satisfy the GRQ constraints, which affords optionality in the choice of objective functions, and domain constraints in beamformer designs. Bounded and smoothed GDI targets can be fitted to a desired input GDI curve given the spatial covariance matrices $`\bf{A}(\omega)`$, $`\bf{R}(\omega)`$ over frequency in the function [bf_grq_tau_tgt](../bf/bf_grq_tau_tgt.m). Several sample GDI targets are shown for the loudspeaker piston array example:

| $\tau=6$ dB |  $\tau=0$ dB |  $\tau=-12$ dB | 
| :---: |  :---: |  :---: | 
| <img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_6_tau_tgt.png" width="400"/>   | <img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_0_tau_tgt.png" width="400"/>   | <img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_-12_tau_tgt.png" width="400"/>   |  

whereby the fitted GDI target curve minimizes a least-squares fitting error penalized by the squared finite difference approximation of the curve’s second-derivative. The latter penalty term can be specified via the `tau_lambda` options in the spherical harmonic GRQ design function [bf_grq_sh.m](../bf/bf_grq_sh.m). Next, we consider several beamformer design objectives subject to the constant GDI constraint.

### Minimum Spectral Norm

Beamformer weights $`\bf{w} = \bf{R}^{-\frac{H}{2}} \bf{V} \bf{x}`$, for the matrix of generalized eigenvectors $\bf{V}$ of $\bf{A}$, $\bf{R}$, satisfying constant GDI can be constructed by minimizing the norm of the spectral weights $\bf{x}$ in the function [bf_grq_spec.m](../bf/bf_grq_spec.m). For complex $\bf{x}$, the optimization is under-constrained as each $x_n \in \bf{x}$ can have arbitrary phase-rotation at the solution, which results in non-smooth beamformer responses over frequency. Constraining $\bf{x}$ to the non-negative reals fortunately gives a unique solution to the optimization problem following 

$$\bf{x}_{\ast} = \textrm{arg min}_{\bf{x}}  \textrm{ } f(\bf{\lambda}, \bf{x} ) \quad \textrm{s.t.}  \quad
\sum_{n=1}^N \lvert x_n \rvert^2 = 1, \quad
\sum_{n=1}^N \lambda_n \lvert x_n \rvert^2  = \tau, \quad 
\bf{x} \in \mathbb{R}^{N \times 1}_{\geq 0},
$$

which also linearizes the constant GRI constraints in terms of the squared modulus $\lvert x_n \rvert^2$. Therefore, standard linearly constrained linear and quadratic programming solvers can efficiently find the minimizer for linear and quadratic spectral objective functions $f( \bf{\lambda}, \bf{x} )$:

| Largest Modulus `max` |  Weighted Modulus `wabs`  | Inverted Weighted Modulus `iwabs` | Square Modulus `sq`  | Weighted Square Modulus `wsq`  | Inverted Weighted Squared Modulus `iwsq` | 
| :---: |  :---: | :---: |  :---: | :---: |  :---: | 
| $`\max_n \lvert x_n \rvert^2`$ | $`\sum_{n=1}^N \lambda_n   \lvert x_n \rvert^2`$ | $`\sum_{n=1}^N  \lambda_n^{-1}  \lvert x_n \rvert^2`$  | $`\sum_{n=1}^N  \lvert x_n \rvert^4`$ | $`\sum_{n=1}^N \lambda_n^2 \lvert x_n \rvert^4`$ | $`\sum_{n=1}^N \lambda_n^{-2} \lvert x_n \rvert^4`$ |

In practice, the solutions are not guaranteed to be smooth in frequency due to sensitivity to changes in the eigen-spectra and associated eigenvectors. The beamformer’s responses are not steered towards any target, especially when the GDI target $\tau$ is distant from both extrema.

Sample GRQ (row 1) and GRPQ channel weighted (row 2) target solutions for the largest modulus `max` objective from the function [plot_bf_grq_sh_example.m](../plt/plot_bf_grq_sh_example.m), and via `spec` mode in the design function [bf_grq_sh.m](../bf/bf_grq_sh.m):

```
plot_bf_grq_sh_example('spec', 6, false);
plot_bf_grq_sh_example('spec', 6, true); 	% Channel Weighted
plot_bf_grq_sh_example('spec', 0, false);
plot_bf_grq_sh_example('spec', 0, true); 	% Channel Weighted
```

|   Beamformer Response $\tau = 6$ (dB) |  Beamformer Statistics $\tau = 6$ (dB) | Beamformer Response  $\tau = 0$ (dB)  |  Beamformer Statistics $\tau = 0$ (dB)| 
| :---: |  :---: | :---: |  :---: | 
<img src="../figs/figs_bf/grq_mixed/bf_spec_tau_6_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_spec_tau_6_bs.png" width="420"/>   |<img src="../figs/figs_bf/grq_mixed/bf_spec_tau_0_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_spec_tau_0_bs.png" width="420"/>   |
<img src="../figs/figs_bf/grq_mixed/bf_spec_tau_6_wt_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_spec_tau_6_wt_bs.png" width="420"/>   |<img src="../figs/figs_bf/grq_mixed/bf_spec_tau_0_wt_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_spec_tau_0_wt_bs.png" width="420"/>   |



### Minimum Norm Projection

We can express beamformer weights by their nearest projection $`\bf{w} = \bf{a} + \bf{u}`$ of the steering vector $\bf{a}$ onto the GRQ constraints via unknown vector $\bf{u} \in \mathbb{C}^{N \times 1}$. Minimizing the norm of the projection vector $\bf{u}$ in the function [bf_grq_proj.m](../bf/bf_grq_proj.m) yields the following optimization problem:

$$\bf{u}_{\ast} =  \textrm{arg min}_{\bf{u}} \textrm{ } \bf{u}^H\bf{u}   \quad \textrm{s.t.}  \quad 
\bf{w}^H \bf{D} \bf{w} = 0,
$$

where the solution is the root of a secular equation of Lagrangian multipliers[^LUO_CDB] closest to $0$. If the steering vector $\bf{a}$ is regular or varies smoothly over frequency, then the solutions $\bf{w}_{\ast}$ yield a regular beam pattern response. Such is the case in loudspeaker arrays where steering vectors are smoother when closer to a transducer’s surface normal of acoustic propagation after normalizing for wavelength and transducer/baffle sizes. 

Sample GRQ (row 1) and GRPQ channel weighted (row 2) target solutions from the function [plot_bf_grq_sh_example.m](../plt/plot_bf_grq_sh_example.m), and via `proj` mode in the design function [bf_grq_sh.m](../bf/bf_grq_sh.m):

```
plot_bf_grq_sh_example('proj', 6, false);
plot_bf_grq_sh_example('proj', 6, true); 	% Channel Weighted
plot_bf_grq_sh_example('proj', 0, false);
plot_bf_grq_sh_example('proj', 0, true); 	% Channel Weighted
```

|   Beamformer Response $\tau = 6$ (dB) |  Beamformer Statistics $\tau = 6$ (dB) | Beamformer Response  $\tau = 0$ (dB)  |  Beamformer Statistics $\tau = 0$ (dB)| 
| :---: |  :---: | :---: |  :---: | 
<img src="../figs/figs_bf/grq_mixed/bf_proj_tau_6_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_proj_tau_6_bs.png" width="420"/>   |<img src="../figs/figs_bf/grq_mixed/bf_proj_tau_0_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_proj_tau_0_bs.png" width="420"/>   |
<img src="../figs/figs_bf/grq_mixed/bf_proj_tau_6_wt_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_proj_tau_6_wt_bs.png" width="420"/>   |<img src="../figs/figs_bf/grq_mixed/bf_proj_tau_0_wt_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_proj_tau_0_wt_bs.png" width="420"/>   |


### Maximum Sensitivity

In loudspeaker beamformer designs where electrical headroom is limited due to transducer and amplifier specifications, we can maximize the array’s sensitivity given by the on-axis response’s output / input power ratio. This generalizes conventional loudspeaker sensitivity[^TOOLE] which measures a single transducer’s anechoic on-axis sound-pressure-level at one meter under constant voltage root-mean-square input.

We formulate the optimization in terms of the linear distortionless response constraint w.r.t. the steering vector $\bf{a}$, and GRQ constraints. Minimizing the electrical power or squared norm of the beamformer weights $\bf{w}$ gives the maximum sensitivity constant directivity (MSCD) beamformer design in the function [bf_grq_mscd.m](../bf/bf_grq_mscd.m), which solves for the following optimization problem:

$$\bf{w}_{\ast} = \textrm{arg min}_{\bf{w}}  \textrm{ } \bf{w}^H \bf{w} 	\quad  \textrm{s.t.} 	\quad 
\bf{w}^H \bf{D} \bf{w} = 0, \quad \bf{w}^H \bf{a} = 1,
$$

where the beamformer’s on-axis response is distortionless and therefore has constant power over frequency, its beampattern has constant GDI, and the total electrical power required to achieve the response is minimized. The solution is similarly found via Lagrangian multipliers[^LUO_CDB] in the minimum norm projection problem.

Sample GRQ (row 1) and GRPQ channel weighted (row 2) target solutions from the function [plot_bf_grq_sh_example.m](../plt/plot_bf_grq_sh_example.m), and via `sensitivity` mode in the design function [bf_grq_sh.m](../bf/bf_grq_sh.m):

```
plot_bf_grq_sh_example('sensitivity', 6, false);
plot_bf_grq_sh_example('sensitivity', 6, true); 	% Channel Weighted
plot_bf_grq_sh_example('sensitivity', 0, false);
plot_bf_grq_sh_example('sensitivity', 0, true); 	% Channel Weighted
```

|   Beamformer Response $\tau = 6$ (dB) |  Beamformer Statistics $\tau = 6$ (dB) | Beamformer Response  $\tau = 0$ (dB)  |  Beamformer Statistics $\tau = 0$ (dB)| 
| :---: |  :---: | :---: |  :---: | 
<img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_6_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_6_bs.png" width="420"/>   |<img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_0_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_0_bs.png" width="420"/>   |
<img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_6_wt_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_6_wt_bs.png" width="420"/>   |<img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_0_wt_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_sensitivity_tau_0_wt_bs.png" width="420"/>   |


### Maximum Efficiency

Beam patterns belonging to loudspeaker array’s are less regular than that of a single transducer due to spacing-induced lobing. Array sensitivity is less predictive of its acoustic efficiency in room environments due to reflected acoustic paths and modal behavior. Instead, the output power averaged over [spatial windows](#loudspeaker-directivity-and-array-example) per input unit of power is a more representative objective of the array’s efficiency. Thus, we can specify the acoustic power over a spatial window density $`f_{E}(\bf{\Omega})`$ by constructing another spatial covariance matrix $`\bf{E} = \mathbb{E}_{\bf{\Omega} \sim f_E(\bf{\Omega}) } \left [ \bf{a}(\bf{\Omega})  \bf{a}^H(\bf{\Omega})  \right ]`$. We can therefore maximize the acoustic power subject to constant GDI, and constant input power via the maximum efficiency constant directivity (MECD) beamformer design in the function [bf_grq_mecd.m](../bf/bf_grq_mecd.m), which solves the following optimization problem:

$$\bf{w}_{\ast} = \textrm{arg max}_{\bf{w}} \textrm{ } \bf{w}^H \bf{E} \bf{w} 	\quad \textrm{s.t.} \quad 
\bf{w}^H \bf{D} \bf{w} = 0, \quad \bf{w}^H \bf{w} = 1,
$$

which additionally substitutes the linear distortionless constraint of the MSCD beamformer with a quadratic unity power constraint. This moves array steering into the choice of the window function in covariance matrix $\bf{E}$, and reduces the sensitivity of the array to steering-vector mismatch along single directions. The solution can be found via projected ascent[^LUO_CDB] of a sequence of minimum matrix norm projections onto the unit sphere, or via a tight semi-definite relaxation[^PATAKI] of the problem via the trace operator.


Sample GRQ (row 1) and GRPQ channel weighted (row 2) target solutions from the function [plot_bf_grq_sh_example.m](../plt/plot_bf_grq_sh_example.m), and via `efficiency` mode in the design function [bf_grq_sh.m](../bf/bf_grq_sh.m):

```
plot_bf_grq_sh_example('efficiency', 6, false);
plot_bf_grq_sh_example('efficiency', 6, true); 	% Channel Weighted
plot_bf_grq_sh_example('efficiency', 0, false);
plot_bf_grq_sh_example('efficiency', 0, true); 	% Channel Weighted
```

|   Beamformer Response $\tau = 6$ (dB) |  Beamformer Statistics $\tau = 6$ (dB) | Beamformer Response  $\tau = 0$ (dB)  |  Beamformer Statistics $\tau = 0$ (dB)| 
| :---: |  :---: | :---: |  :---: | 
<img src="../figs/figs_bf/grq_mixed/bf_efficiency_tau_6_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_efficiency_tau_6_bs.png" width="420"/>   |<img src="../figs/figs_bf/grq_mixed/bf_efficiency_tau_0_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_efficiency_tau_0_bs.png" width="420"/>   |
<img src="../figs/figs_bf/grq_mixed/bf_efficiency_tau_6_wt_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_efficiency_tau_6_wt_bs.png" width="420"/>   |<img src="../figs/figs_bf/grq_mixed/bf_efficiency_tau_0_wt_bp.png" width="600"/>   | <img src="../figs/figs_bf/grq_mixed/bf_efficiency_tau_0_wt_bs.png" width="420"/>   |


In the next tutorial, we will cover beamformer realization methods in terms of optimal filter networks.


[^CTA_2034_A]: Consumer Electronics Association (CEA). Standard Method of Measurement for In-Home Loudspeakers. ANSI/CEA-2034-A (ANSI/CTA-2034-A), Feb. 2015.

[^LUO_CDB]: Luo, Y. (2024, August). [Constant Directivity Loudspeaker Beamforming](https://ieeexplore.ieee.org/abstract/document/10715159). In 2024 32nd European Signal Processing Conference (EUSIPCO) (pp. 246-250). IEEE.

[^TOOLE]: Floyd E. Toole. Sound Reproduction: The Acoustics and Psychoacoustics of Loudspeakers and Rooms (3rd Edition, Routledge, 2017)

[^PATAKI]: G. Pataki, “On the rank of extreme matrices in semidefinite programs and the multiplicity of optimal eigenvalues,” Mathematics of Operations Research, vol. 23, 05 1998

[^BARANEK]: Beranek, L. L., & Mellow, T. J. (2012). Acoustics: Sound Fields and Transducers. Academic Press.

[^LINKWITZ]: Linkwitz, Siegfried H. (February 1976). "Active Crossover Networks for Noncoincident Drivers". Journal of the Audio Engineering Society. 24 (1): 2–8. Retrieved 2024-05-05.