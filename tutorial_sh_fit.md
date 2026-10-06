# Tutorial: Spherical Harmonic (SH) Basis Fitting

A common task in signal processing and machine learning is to fit a smooth function to data collected over spatial-temporal or spatial-frequency domain. For spatial domains over fixed distances, the spherical coordinates are specified in terms of co-latitude $\theta \in [0, \pi]$ and azimuth $\phi \in [0, 2 \pi)$ domains, whereby the [spherical harmonic](guide_func.md) (SH) bases represent the smooth and orthogonal functions that span the domain. In this tutorial, we cover SH fitting methods and applications introduced in our paper
> Luo, Y. (2021). [Spherical Harmonic Covariance and Magnitude Function Encodings for Beamformer Design](https://link.springer.com/article/10.1186/s13636-021-00230-7). EURASIP Journal on Audio, Speech, and Music Processing, 2021(1), 41.

as well sum-of-magnitude squared extensions introduced in 

> Luo, Y. (2026). [Spherical Harmonic Sliced Wasserstein Displacement Interpolation for Acoustic Source and Reflection Density Modeling](https://arxiv.org/abs/2609.22028). arXiv:2609.22028.

The high-level methods are as follows:
* [Complex Response Fit](#complex-response-fit)
  * [Truncated Singular Value Decomposition](#truncated-singular-value-decomposition-least-squares)
  * [Tikhonov Regularization](#tikhonov-regularized-least-squares)
  * [Radial Basis Functions](#radial-basis-function-max-marginal-likelihood)
* [Magnitude Response Fit](#magnitude-response-fit)
  * [Sum-of-Magnitude Squared](#sum-of-magnitude-squared-least-squares)
  * [Magnitude Squared](#magnitude-squared-least-squares)
  * [Mixtures-of-Magnitude Squared](#mixture-of-magnitude-squared-least-squares)
  * [Mixture-Power Squared](#mixture-power-least-squares)

## Complex Response Fit


The squared error objective is given by

$$\displaystyle
\begin{split} 
F_{LS}(\bf{C}) & = \sum_{n=1}^N \left \lvert  \bf{Y}(\theta_n, \phi_n) \bf{C} - x_n \right \rvert^2
%%%%%%%
 = \left ( \bf{Y} \bf{C} - \bf{x} \right )^H \left ( \bf{Y} \bf{C} - \bf{x} \right ), \\\\
%%%%%%
\bf{Y} & = [Y(\theta_1, \phi_1); … ; Y(\theta_N, \phi_N)] \in \mathbb{C}^{N \times N_C}, \quad 
%%%%%
\bf{C}  = [C_0^0, C_1^{-1}, C_1^0, C_1^1, … C_L^L ]^T \in \mathbb{C}^{N_C \times 1}, 
\end{split}
$$

where row-matrix $\bf{Y}$ are the SH evaluations over $N$ spherical coordinates $(\theta, \phi)$ for $N_C = (L + 1)^2$ of max-degree $L$, he SH expansion coefficients are contained in the column-vector $\bf{C}$, and the target responses in column-vector $\bf{x} = [x_1, … x_N]^T \in \mathbb{C}^{N \times 1}$.

### Truncated Singular Value Decomposition Least Squares

The singular value decomposition (SVD) is given by

$$\bf{Y} = \bf{U} \bf{\Sigma} \bf{V}^H, \quad 
\bf{U} \in \mathbb{C}^{N \times N_C}, \quad 
\bf{\Sigma}, \bf{V} \in \mathbb{C}^{N_C \times N_C}, $$ 

where $`\bf{U} = [\bf{u}_{1}, …, \bf{u}_{N_C}]`$ and $`\bf{V} = [\bf{v}_{1}, …, \bf{v}_{N_C}]`$ are orthonormal matrices, and $`\bf{\Sigma} = \bf{\Lambda}^{\frac{1}{2}}`$ is the diagonal matrix of the square-root of eigenvalues belonging to the normal matrix $`\bf{Y}^H \bf{Y} = \bf{V} \bf{\Lambda} \bf{V}^H`$ for $`\bf{\Lambda} = \textrm{diag}(\bf{\lambda})`$. The least squares solution $`\bf{C}_{\ast}`$ follows from inversion of the normal matrix given by 

$$\bf{C}_{\ast} = (\bf{Y}^H \bf{Y})^{-1} \bf{Y}^H \bf{x} = \bf{V} \bf{\Sigma}^{-1} \bf{U}^H \bf{x}, $$

which is sensitive to noise in $\bf{x}$ when projected onto singular vectors corresponding to small singular values. The condition number of $\bf{Y}$ increases when the sampling spherical coordinates $(\theta, \phi)$ are non-uniform. Therefore, one choice of regularization is to truncate the small singular values $\sqrt{\lambda} \leq \tau \max_{\lambda \in \bf{\lambda}} \sqrt{ \lambda }$ for some fraction $\tau$ of the max singular value by zeroing those reciprocal singular values in $\bf{\Sigma}^{-1}$. Other selection criteria includes the fraction of total variance, the bottom (smallest) $100 \tau$ percent, and the discrete Picard condition[^HANSEN_DPC] (DPC) where $`\lim_{i \rightarrow \infty} |\bf{u}_i^H x| / \sigma_i = 0`$ for descending singular values  $`\sigma_i`$ indexed $`i > k_c`$ after some crossover index $`k_c`$. This is implemented in the function `sh_fit_svd.m`.

Let us walk through an example in `plot_sh_fit_svd_example.m`:

* Generate a random SH field, sample the field over randomized spherical coordinates, and add Gaussian noise to observations:
  ```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Generate and sample from random field
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  rng(441);
  max_odr = 6;
  is_real = false;
  N_C = (max_odr + 1)^2;
  N_pts = (max_odr + 1)^2;
  
  %Reference field
  C_ref = sh_rand(max_odr, 1, is_real);
  
  %Sample at spherical coordinates
  %[theta, phi] = sh_grd_fib(N_pts);       %Fibonnaci
  [theta, phi] = sh_grd_rand_unis(N_pts);  %Random over sphere
  
  X = sh_dec(C_ref, theta, phi, is_real);
  %Add noise
  X = X + (randn(size(X)) + randn(size(X)) * 1i) * 5e-2;
  ```

* Fit SH coefficients to the randomized samples with different `mode` and $\tau$ specified by `trunc_frac`:

  ```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fit across various modes and truncation fractions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mode_list = {'max', 'totalvar', 'bottom'};
  trunc_frac_list = [0.01, 0.1, 0.2];
  
  N_trunc = numel(trunc_frac_list);
  N_modes = numel(mode_list);
  
  C_list = zeros([N_C, N_modes, N_trunc]);
  err_list = zeros(N_modes, N_trunc);
  trunc_percent_list = zeros(N_modes, N_trunc);
  
  for i = 1:N_modes
      for j = 1:N_trunc
          [C_list(:, i, j), ~, trunc_percent_list(i, j)] = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac_list(j), mode_list{i});
          err_list(i, j) = err_NSHMSQ(C_ref, C_list(:, i, j)); %Normalized error
      end
  end

  % Non-regularized least squares
  [C_ls]  = sh_fit_svd(X, theta, phi, max_odr, is_real, 0, 'max', 'enable_disp', true);
  
  % Picard cross criterion
  trunc_frac_pc = 0.5;
  [C_pc, ~, ~, h_pc]  = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac_pc, 'picard', 'enable_disp', true);
  h_picard = h_pc{2};
  ```

* Fit SH coefficients after truncating successively larger singular values:
  ```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fit after truncating smallest singular values
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  trunc_frac_all = (0:N_C) ./ N_C;
  N_frac_all = N_C + 1;
  err_frac_all = zeros(N_frac_all, 1);
  trunc_percent_frac_all = zeros(N_frac_all, 1);
  for i = 1:N_frac_all
      [C_i, ~, trunc_percent_frac_all(i)] = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac_all(i), 'bottom');
      err_frac_all(i) = err_NSHMSQ(C_ref, C_i); %Normalized error
  end
  ```
  
* Plot the reference field, non-regularized least squares, Picard crossover, and fitting error:

| Reference | Non-Regularized Least Squares |  Picard Crossover TSVD |
| :-: | :-: | :-: | 
| <img src="./figs/figs_sh/fit/svd_ref.png" width="400"/> | <img src="./figs/figs_sh/fit/svd_ls.png" width="400"/> |  <img src="./figs/figs_sh/fit/svd_pc.png" width="400"/> | 

where the non-regularized least squares yields a poor fit due to sensitivity to noise of the non-uniform sampling grid. We can compare the decay of the projection of the noisy observations onto the left singular vectors $|\bf{u}_i^H x|$ to the descending singular values $\sigma_i$ in the following Picard plot: 

 |  Picard Plot | Fitting Error |
  | :-: | :-: |
<img src="./figs/figs_sh/fit/svd_picard.png" width="400"/> |<img src="./figs/figs_sh/fit/svd_err.png" width="600"/> |

where the singular values decay faster than the magnitude of the projections. This amplifies the noise in $\bf{x}$ in the least squares solution $`\bf{C}_{\ast} = \sum_{i=1}^{N_C} \bf{v}_i \frac{\bf{u}_i^H x}{\sigma_i}`$ for small $\sigma_i$ in the tail. We can identify when the decays diverge via reverse integration and normalizing their cumulative distribution functions (CDFs) as to account for arbitrary scaling of $\bf{x}$. Sorting the singular values in ascending order and computing the normalized CDF (NCDF), we see that the curve crosses that of the projection’s NCDF after 19 singular values, and exceeds half the latter after 5 singular values. More generally, we can define a crossover index for the DPC that removes ascending singular values below index $k_c$ for some adjustable fraction $\tau$ given by

$$k_c = \textrm{arg min }  i \quad  \textrm{s.t.}  \quad \textrm{NCDF}(\sigma_i) \geq \tau  \textrm{ NCDF}(\bf{u}^H_i \bf{x}).$$

Choosing a small $\tau$ truncates smaller singular values that have decayed further below the observation’s projection onto the corresponding singular vectors. In fact, truncating just 2/49 of the smallest singular values yields the minimum normalized magnitude squared error (NMSQE) over the spherical coordinates given by

$$\displaystyle
\begin{split} 
\textrm{NMSQE}(\bf{C}_{ref}, \bf{C}_{\ast}) & =  \frac{ \int_{\theta=0}^{\pi} \int_{\phi = 0}^{2\pi} \left \lvert  Y(\theta, \phi) (\bf{C}_{ref} - \bf{C}_{\ast})  \right \rvert^2 \sin (\theta) d \theta d \phi }{\int_{\phi = 0}^{2\pi} \left \lvert  Y(\theta, \phi) \hat{\bf{C}}_{ref}  \right \rvert^2 \sin (\theta) d \theta d \phi } \\\\
%%%%%%
& = \frac{ (\bf{C}_{ref} - \bf{C})^H (\bf{C}_{ref} - \bf{C}) }{\hat{\bf{C}}_{ref}^H \hat{\bf{C}}_{ref}}, \quad 
%%%%%%%%%
\hat{\bf{C}}_{ref} = \bf{C}_{ref} - [C_0^0, 0, … , 0]^T,
\end{split}$$

where the denominator normalizes for the variance of the reference field. We can compare the Picard crossover solution to other truncation modes shown as follows:

  | Mode | $\tau = 0.01$ | $\tau = 0.1$ | $\tau = 0.2$ |
  | :---: |  :---: |  :---:  |  :---: |
  | Max |<img src="./figs/figs_sh/fit/svd_fit_1_1.png" width="400"/>|<img src="./figs/figs_sh/fit/svd_fit_1_2.png"  width="400"/>|<img src="./figs/figs_sh/fit/svd_fit_1_3.png"  width="400"/>|
  | Total Variation |<img src="./figs/figs_sh/fit/svd_fit_2_1.png" width="400"/>|<img src="./figs/figs_sh/fit/svd_fit_2_2.png"  width="400"/>|<img src="./figs/figs_sh/fit/svd_fit_2_3.png"  width="400"/>|
  | Bottom |<img src="./figs/figs_sh/fit/svd_fit_3_1.png" width="400"/>|<img src="./figs/figs_sh/fit/svd_fit_3_2.png"  width="400"/>|<img src="./figs/figs_sh/fit/svd_fit_3_3.png"  width="400"/>|

### Tikhonov Regularized Least Squares

We can regularize the least square solutions by introducing a quadratic penalty term to the SH coefficients $\bf{C}$ in the squared error objective given by

$$\displaystyle
\begin{split}
\bf{C}_{\ast} & = \textrm{arg min}_{\bf{C}}  \left \Vert  \bf{Y} \bf{C} - \bf{x} \right \Vert^2 + \alpha \bf{C}^H \bf{Q} \bf{C} \\
%%%%%
& = (\bf{Y}^H \bf{Y} + \lambda \bf{Q} )^{-1} \bf{Y}^H \bf{x}, \quad \alpha \geq 0,
\end{split}
$$

where for diagonal penalty matrix $\bf{Q} = \textrm{diag}(q)$, larger entrants in $`\bf{q} = [q_0^0, q_{1}^{-1}, q_{1}^{0}, q_{1}^{1}, … , q_{L_C}^{L_C}] \in \mathbb{R}^{N_C}_{\geq 0}`$ penalize corresponding coefficients in $\bf{C}$. In the case of $\bf{q} = \bf{1}$, the penalty is the squared Euclidean norm of $\bf{C}$. For processes where the gradient of field is constrained, the larger SH degrees $l$ are penalized. Several options include quadratic $l^2$ and quadratic + linear $`q_{l}^m = l (l+1)`$[^DURAIS_TR], which are proportional to the normalization terms in SH functions. 

For general matrix $\bf{Q}$, Tikhonov regularization relates to SVD when $\bf{Q}$ is diagonalized by the right singular vectors $\bf{V}$ of $\bf{Y}$ as follows:

$$\bf{Q} = \bf{V} \bf{D} \bf{V}^H, \quad \bf{D} = \textrm{diag}(\bf{d}) \quad \Rightarrow \quad 
\bf{C}_{\ast} = \bf{V} (\bf{\Sigma}^2 + \alpha \bf{D}  )^{-1}  \bf{\Sigma}  \bf{U}^H \bf{x},
$$

where the entrants of $\bf{d}$ can go to infinity for truncating singular values. The regularized least squares solution is given by 

$$\bf{C}_{\ast} = \sum_{i=1}^{N_C} \bf{v}_i \frac{\sigma_i}{\sigma_i^2 + d_i}  \bf{u}_i^H x,$$ 

where we wish to minimize the amount of regularization whilst satisfying DPC. Therefore, let us find $\bf{d}$ with minimum norm that satisfies several constraints in the following linear optimization problem:

$$\displaystyle
\begin{split}
\min_{\bf{d}} \sum_{i=1}^{k_{c}} d_{i} \quad  \textrm{s.t.} \quad
%%%%%%%%
d_{i}  \geq 0, \quad  & \textrm{Non-negative solution} \\
%%%%%%%%
|\bf{u}_{i}^H \bf{x}| \frac{\sigma_{i}}{\sigma_{i}^2 + d_{i}} \leq |\bf{u}_{i+1}^H \bf{x}| \frac{\sigma_{i+1}}{\sigma_{i+1}^2 + d_{i+1}}, \quad & \textrm{Monotonic Picard ratio, } i < k_c \\
%%%%%%%%
|\bf{u}_{i}^H \bf{x}| \frac{\sigma_i}{\sigma_{i}^2 + d_{i}} \leq \frac{|\bf{u}_{i+1}^H \bf{x}|}{\sigma_{i+1}}, \quad & \textrm{Monotonic Picard ratio, } i = k_c
\end{split}
$$

where $k_c$ is the Picard crossover index found for ascending singular values $\sigma_i \leq \sigma_{i+1}$. The monotonoic constraints can be converted into linear constraints via cross-multiplying the denominator terms. They ensure that the sequence of Picard ratio of projections onto descending singular vectors / regularized singular values do not increase. The Tikhonov regularized least squares methods are implemented in `sh_fit_tr.m`.

Let us consider the example in `plot_sh_fit_tr_example.m`:

* Generate a random SH field, sample the field over randomized spherical coordinates, and add Gaussian noise to observations:
  ```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Generate and sample from random field
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rng(443);
  
  max_odr = 6;
  is_real = false;
  N_C = (max_odr + 1)^2;
  N_pts = (max_odr + 1)^2;
  
  %Reference field
  C_ref = sh_rand(max_odr, 1, is_real);
  
  %Sample at spherical coordinates
  %[theta, phi] = sh_grd_fib(N_pts);       %Fibonnaci
  [theta, phi] = sh_grd_rand(N_pts);  %Random over sphere
  
  X = sh_dec(C_ref, theta, phi, is_real);
  %Add noise
  X = X + (randn(size(X)) + randn(size(X)) * 1i) * 5e-2;
  ```
  
* Fit non-regularized, Tikhonov regularized, and truncated SVD with Picard cross-over fraction of 0.5. Evaluate their NMSQE compared to reference field:
  ```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fit
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  C_ls = sh_fit_tr(X, theta, phi, max_odr, is_real, 0);
  C_tr_id = sh_fit_tr(X, theta, phi, max_odr, is_real, 0.1, 'identity');
  [C_tr_pc, ~, h_tr_pc_list] = sh_fit_tr(X, theta, phi, max_odr, is_real, 0.5, 'picard', 'enable_disp', true);
  C_svd_pc = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.5, 'picard');
  
  % Plot reference
  dB_lim = [-32, 24];
  h_ref = sh_plt(C_ref, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Reference', 'disp_theta_phi', [theta, phi], 'disp_theta_phi_markersize', 12);
  
  h_ls = sh_plt(C_ls, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Non-Regularized Least Squares');
  
  h_tr_id = sh_plt(C_tr_id, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'TR Least Squares L2 Norm');
  
  h_tr_pc = sh_plt(C_tr_pc, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'TR Least Squares Picard');
  
  h_svd_pc = sh_plt(C_svd_pc, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'TSVD Picard');
  
  % Compute errors
  err_ls = err_NSHMSQ(C_ref, C_ls)
  err_tr_id = err_NSHMSQ(C_ref, C_tr_id)
  err_tr_pc = err_NSHMSQ(C_ref, C_tr_pc)
  err_svd_pc = err_NSHMSQ(C_ref, C_svd_pc)
  ```

* Compare Tikhonov regularized, and truncated SVD with varying Picard crossover fractions:

  ```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute errors for varying crossover fractions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  N_frac = 10;
  frac_list = linspace(0, 1, N_frac);
  
  err_tr_pc_list = zeros(1, N_frac);
  err_svd_pc_list = zeros(1, N_frac);
  for n = 1:N_frac
      C_tr_pc_n = sh_fit_tr(X, theta, phi, max_odr, is_real, frac_list(n), 'picard');
      C_svd_pc_n = sh_fit_svd(X, theta, phi, max_odr, is_real, frac_list(n), 'picard');
  
      err_tr_pc_list(n) = err_NSHMSQ(C_ref, C_tr_pc_n);
      err_svd_pc_list(n) = err_NSHMSQ(C_ref, C_svd_pc_n);
  end
  ```

* Plot the fitted fields and errors:
  | Reference | Picard Tikhonov Regularized Least Squares | 
  | :-: | :-: | 
  | <img src="./figs/figs_sh/fit/tr_ref.png" width="400"/> | <img src="./figs/figs_sh/fit/tr_pc.png" width="400"/> |
     Picard Plot | Varying Picard Crossover |
   <img src="./figs/figs_sh/fit/tr_picard.png" width="400"/> |  <img src="./figs/figs_sh/fit/tr_tsvd_picard.png" width="400"/>

    where the Tikhonov regularization smoothly lifts the singular values below the Picard crossover index above the singular vector projection in the NCDF. For small Piard crossover fractions, too few tiny singular values are regularized or truncated, which results in high sensitivity to noise in the resulting fit. Tikhonov regularization outperforms the TSVD when more singular values are regularized. 


### Radial Basis Function Max Marginal Likelihood

SH model-order selection (max-degree $L$) remains a difficult choice for least squares fitting methods on irregular sampling grids. When the underlying observations are smooth in spherical coordinates, as is common to many physical processes, regularization techniques do not explicitly account for covariances at neighboring inputs. Non-parametric models such as radial basis functions (RBF), krigging, and Gaussian processes address the latter. Such methods specify an interpolating function (RBF network, Gaussian process posterior mean) given by

$$\displaystyle
\begin{split}
f(\theta, \phi) &  =  \tilde{\bf{k}}(\theta, \phi)  (\bf{K} + \sigma^2 \bf{I})^{-1} \bf{x}, \\
%%%%%%
K_{ij} & = K(\theta_i, \phi_i, \theta_j, \phi_j), \\
%%%%%%%
k_{i}(\theta, \phi) & = K(\theta, \phi, \theta_i, \phi_i), \quad
%%%%%
\tilde{k}_i(\theta, \phi)  = \bf{Y}(\theta, \phi) \bf{C}_i \approx k_i(\theta, \phi),
\end{split}$$

where the covariance functions $`K(\theta, \phi, \theta_i, \phi_i)`$ have SH expansion coefficients $`\bf{C}_n`$ in $`\tilde{\bf{k}}(\theta, \phi) = [\tilde{k}_1(\theta, \phi), … ,\tilde{k}_{N}(\theta, \phi)] = \bf{Y}(\theta, \phi) [\bf{C}_1, … , \bf{C}_N] \in \mathbb{R}^{1 \times N}`$, and the Gram matrix $\bf{K} = \mathbb{R}^{N \times N}$ remains exact to ensure positive definiteness. The SH expansion coefficients are therefore given by

$$f(\theta, \phi) = \bf{Y}(\theta, \phi) \bf{C}_{RBF}, \quad \bf{C}_{RBF} = [\bf{C}_1, … , \bf{C}_N] (\bf{K} + \sigma^2 \bf{I})^{-1} \bf{x} \in \mathbb{C}^{N_C \times N}.$$

For stationary covariance functions, the interpolating function are a sum of local density expansions centered on $`(\theta_1, \phi_1), … , (\theta_N, \phi_N)`$, and the input observations in $\bf{x}$ following the Representer’s theorem. Moreover, Mercer’s theorem expands the positive-definite kernel functions along an integrable basis, where it suffices to show that our covariance functions have SH expansions. For example, the class of Matérn covariance functions of chordal distances are stationary and can be expanded along SH bases via the Legendre polynomials[^LUO_COV]. The chordal distance is simply the Euclidean distances between points on a unit sphere given by

$$d(\theta, \phi, \theta', \phi') = 2 \sin \left ( \frac{ \left | \cos^{-1} (\bf{v}^T \bf{v}') \right |   }{2} \right ) = \left \lVert \bf{v} - \bf{v}'  \right \rVert_2, $$

where $\bf{v} \in \mathbb{R}^{3}$ is a unit vector.  The Matérn covariance functions $k(\theta, \phi, \theta’, \phi')$ of the chordal distances $d(\theta, \phi, \theta’, \phi')$ include the squared exponential and exponential, and are listed in terms of decreasing smoothness (differentiability w.r.t. $d$) as follows:


|Squared Exponential | Matérn $\nu=5/2$ | Matérn $\nu=3/2$ | Exponential |
| :-: | :-: | :-: | :-: | 
|$\exp \left ( - \frac{d^2} {2 \ell^2} \right )$ | $\left ( 1 +  \frac{\sqrt{5}  d}{\ell} + \frac{5 d^2}{3 \ell^2}   \right )  \exp \left ( - \frac{\sqrt{5} d}{\ell} \right )$ | $\left ( 1 +  \frac{\sqrt{3}  d}{\ell}  \right )  \exp \left ( - \frac{\sqrt{3} d}{\ell} \right )$ | $\exp \left ( \frac{d}{\ell} \right )$ |

Thus, RBF fitting follows from optimizing the RBF’s hyper-parameters such as its length-scale $\ell$ and the unknown noise-variance term $\sigma^2$. For Gaussian processes, hyper-parameter optimization follows from maximizing the following log-marginal likelihood (LMH) objective[^RASMUSSEN_GP]:

$$F_{LMH} = -\frac{1}{2} \textrm{trace} \left ( \bf{x}^H (\bf{K} + \sigma^2 \bf{I} )^{-1} \bf{x} \right ) - \frac{1}{2} \log \left \vert \bf{K} + \sigma^2 \bf{I} \right \vert - \frac{N}{2} \log(2\pi), $$

which is differentiable w.r.t.  $\sigma^2$ and $\ell$  within covariance functions $k(\theta, \phi, \theta’, \phi’)$. In practice, the optimization is tractable for small to medium sized number of observations ($N \leq 4096$), as the computational costs are bounded above by the $O(N^3)$ matrix-inverse operation, although sub-cubic and fast algorithms are possible for structured matrices and multipole expansions. Therefore, RBF fitting is suitable unto max-degree $L=63$ SH expansions. The RBF fitting methods are implemented in `sh_fit_rbf.m`.

Let us consider the example in `plot_sh_fit_rbf_example.m`:

* Sample from noisy random field, generate RBF fits with varying length-scales $\ell$ of increasing smoothness, plot and compute error: 
  ```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compare RBF fits for varying length-scales
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rng(43);
  
  max_odr = 4;
  max_odr_fit = 16;
  is_real = false;
  N_pts = 50;
  
  % Reference field
  C_ref = sh_rand(max_odr, 1, is_real);
  C_ref = sh_resize(C_ref, max_odr_fit);
  
  % Sample at spherical coordinates
  [theta, phi] = sh_grd_rand(N_pts);  %Random over sphere
  
  X = sh_dec(C_ref, theta, phi, is_real);
  % Add noise
  noise_std = 5e-1;
  X = X + (randn(size(X)) + randn(size(X)) * 1i) * noise_std;
  
  % Plot ref.
  dB_lim =  [-60, 20];
  h_ref_ell = sh_plt(C_ref, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_theta_phi', [theta, phi], 'disp_theta_phi_markersize', 12);
  
  % Fit RBFs for various length-scales
  ell_list = [0.2, 0.4, 0.8];
  N_ell = numel(ell_list);
  C_fit_list = zeros((max_odr_fit + 1)^2, N_ell);
  err_list = zeros(1, N_ell);
  h_rbf_ell_list = cell(1, N_ell);
  varargin = {'max_iter', 0, 'lb', [0, noise_std^2], 'ub', [1, noise_std^2]};
  mode = 'Mat32';
  for n = 1:N_ell
      C_fit_list(:, n) = sh_fit_rbf(X, theta, phi, max_odr_fit, is_real, mode, ell_list(n), noise_std^2, varargin{:}); 
      err_list(n) = err_NSHMSQ(C_ref, C_fit_list(:, n));
      h_tmp = sh_plt(C_fit_list(:, n), 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', [mode, ' length-scale = ', num2str(ell_list(n))]);
      h_rbf_ell_list{n} = h_tmp{1};
  end
  err_list
  ```
  
  | Reference | $\ell = 0.2$, NMSQE =  0.4509 | 
  | :-: | :-: | 
  |<img src="./figs/figs_sh/fit/rbf_ell_ref.png" width="400"/>|<img src="./figs/figs_sh/fit/rbf_ell_1.png" width="400"/>|
  | $\ell = 0.4$, NMSQE =  0.2618 | $\ell = 0.8$, NMSQE =  0.2747 |
  <img src="./figs/figs_sh/fit/rbf_ell_2.png" width="400"/>|<img src="./figs/figs_sh/fit/rbf_ell_3.png" width="400"/>|
  
  where the fitted field for $\ell = 0.4$ achieves the lowest error, being neither too smooth or coarse. Note that the noise-variance $\sigma^2$ of the measurements is incorporated into RBF’s Gram matrix. We can optimize $\ell$ by increasing the number of iterations following
  
  ```
  varargin = {'max_iter', 1000, 'lb', [0, noise_std^2], 'ub', [1, noise_std^2]};
  ```
  which yields an best-fit $\ell \approx 0.302$.

* Compare RBF to TSVD and Tikhonov regularization:  Sample from a noisy SH random field of max-degree $L=5$, fit higher max-degree $L=10$ SH bases by optimizing $\ell$, and compute the normalized error.

  ```
  varargin = {'max_odr', 5, 'max_odr_fit', 10, 'N_pts', 100};
  
  [err_svd, h_ref, h_svd]		= tst_sh_fit('TSVD', 'NSHMSQ', varargin{:}, 'svd_mode', 'picard', 'svd_trunc_frac', 0.5);
  [err_tr, h_ref, h_tr]		= tst_sh_fit('TRegu', 'NSHMSQ', varargin{:}, 'tr_mode', 'picard', 'tr_lambda', 0.5);
  
  [err_SqExp, ~, h_SqExp]		= tst_sh_fit('SqExp', 'NSHMSQ', varargin{:});
  [err_Mat52, ~, h_Mat52]		= tst_sh_fit('Mat52', 'NSHMSQ', varargin{:});
  [err_Mat32, ~, h_Mat32]		= tst_sh_fit('Mat32', 'NSHMSQ', varargin{:});
  [err_Exp, ~, h_Exp]		= tst_sh_fit('Exp', 'NSHMSQ', varargin{:});

  err_svd
  err_tr
  err_SqExp
  err_Mat52
  err_Mat32
  err_Exp
  ```

  | Reference | Picard Tikhonov Regularization: NMSQE = 0.3072 | Squared Exponential RBF: NMSQE = 0.1193 |
    | :-: | :-: | :-: |
    |<img src="./figs/figs_sh/fit/rbf_ref.png" width="400"/>|<img src="./figs/figs_sh/fit/rbf_tr.png" width="400"/>|<img src="./figs/figs_sh/fit/rbf_sqexp.png" width="400"/>|
  | Matérn $\nu=5/2$ RBF: NMSQE = 0.0787  | Matérn $\nu=3/2$ RBF: NMSQE = 0.0936 |  Exponential RBF: NMSQE = 0.1981 |
    |<img src="./figs/figs_sh/fit/rbf_mat52.png" width="400"/>|<img src="./figs/figs_sh/fit/rbf_mat32.png" width="400"/>|<img src="./figs/figs_sh/fit/rbf_exp.png" width="400"/>|

	The Picard Tikhonov regularized solution prevents some over-fitting to the higher-order SH bases but still produces a function with spurious high-frequency features. The Picard TSVD solution gives a similar solution and therefore omitted. The Matérn $\nu=5/2$ RBF yields both the maximum log-marginal likelihood as well as the lowest NMSQE. The squared exponential RBF solution is too smooth, and the exponential RBF too coarse.

## Magnitude Response Fit

Some problems require that the solution minimizes the differences in the magnitude responses or power responses over the spherical coordinates. For example, beam-former designs can solve for weights belonging to multiple array element's acoustic frequency responses such that the latter’s weighted summation achieves a target directivity pattern. The beam-former’s phase responses over the spherical coordinates are ignored in the objective function. Another application is density function modeling where we require SH expansions to be non-negative everywhere in spherical coordinates. A non-negative SH expansion $\bf{D}$ can be constructed from the sum-of-magnitude squares (SOMS) of $S$ expansions of SH expansion coefficients $\bf{C}_s \in \mathbb{C}^{N_C \times 1}$  given by

$$ \bf{D} = \sum_{s=1}^{S} \bf{C}_s \diamond \bf{\tilde{C}}_s \quad  \Rightarrow  \quad Y(\theta, \phi) \bf{D} \geq 0,$$

where $\diamond$ is the SH product operator, and $\tilde{\bf{C}}$ is the SH conjugate operator. Several least squares magnitude fitting methods are presented in the function `sh_fit_msq.m`. Let’s walk through an example in `plot_sh_fit_msq_example.m`:

* Generate a sample nonnegative sum-of-magnitude square function with $S=3$ components of max-order $L_C = 4$ component expansion. Sample the function over a uniform Fibonnaci grid, and add randomized Gaussian noise. Fit `sum-of-magnitude square (SOMS)`, `magnitude square (MS)`, `mix-of-magnitude (MOMS)`, and `mixture power (MP)` functions to the samples:

  ```
  num_func = 3;
  max_odr_half = 4;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Generate sample sum-of-magnitude square function
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rng(options.rseed);
  max_odr = 2 * max_odr_half;
  N_C = (max_odr + 1)^2;
  N_pts = N_C;
  
  is_real = false;
  
  % Generate reference field
  C_refs = sh_rand(max_odr_half, num_func, is_real);
  D_ref = sum(sh_msq(C_refs, is_real), 2);
  
  % Sample Fibonnaci uniform points over sphere
  [theta, phi] = sh_grd_fib(N_pts);  
  
  % Sample from field, add random noise
  X = real(sh_dec(D_ref, theta, phi, is_real));
  X = X + randn(N_pts, 1) * 5e-1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fit non-negative SH expansions to noisy observations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [D_SOMS, C_SOMS]    = sh_fit_msq(X, theta, phi, max_odr, is_real, 'SOMS');
  [D_MS, C_MS]        = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MS', 'C0', sh_rand(max_odr_half, 20, is_real) );
  
  B = sh_rand(max_odr_half, N_C, is_real);
  [D_MOMS, C_MOMS]    = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MOMS', 'B', B);
  [D_MP, C_MP]    = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MP', 'B', B, 'B0', rand(N_C, 20));
  
  err_SOMS    = err_NSHMSQ(D_ref, D_SOMS)
  err_MS      = err_NSHMSQ(D_ref, D_MS)
  err_MOMS    = err_NSHMSQ(D_ref, D_MOMS)
  err_MP      = err_NSHMSQ(D_ref, D_MP)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plotting
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  varargin = {'disp_phase', false, 'dB_lim', options.dB_lim };
  h_ref   = sh_plt(D_ref, 'mercator', is_real, 'title_name_override', 'Reference', varargin{:}, 'disp_theta_phi', [theta, phi], 'disp_theta_phi_markersize', 8);
  
  h_SOMS  = sh_plt(D_SOMS, 'mercator', is_real, 'title_name_override', ['Sum-of-Magnitude Square Fit NMSQE: ', num2str(err_SOMS)], varargin{:});
  h_MS    = sh_plt(D_MS, 'mercator', is_real, 'title_name_override', ['Magnitude Square Fit NMSQE: ', num2str(err_MS)], varargin{:});
  
  h_MOMS  = sh_plt(D_MOMS, 'mercator', is_real, 'title_name_override', ['Mix-of-Magnitude Square Fit NMSQE: ', num2str(err_MOMS)], varargin{:});
  h_MP    = sh_plt(D_MP, 'mercator', is_real, 'title_name_override', ['Mixture Power Fit NMSQE: ', num2str(err_MP)], varargin{:});
  ```
  
  |Reference| Sum-of-Magnitude Squared | Magnitude Squared | Mix-of-Magnitude Squared | Mixture Power |
  | :-: | :-: | :-: | :-: | :-: |
  | <img src="./figs/figs_sh/fit/msq_ref_4_3.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_soms_4_3.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_ms_4_3.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_moms_4_3.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_mp_4_3.png" width="350"/> |
  
  The SOMS and MOMS solutions converge to the reference solution. The MS and MP solutions converge to rank-1 approximations due to lack of modeling capacity.

* In the single component case where  $S = 1$, the various methods converge to the same solution:

  |Reference| Sum-of-Magnitude Squared | Magnitude Squared | Mix-of-Magnitude Squared | Mixture Power |
  | :-: | :-: | :-: | :-: | :-: |
  | <img src="./figs/figs_sh/fit/msq_ref_4_1.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_soms_4_1.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_ms_4_1.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_moms_4_1.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_mp_4_1.png" width="350"/> |


### Sum-of-Magnitude Squared Least Squares

The sum-of-magnitude squared objective is given by

$$\displaystyle
\begin{split} 
F_{SOMS}(\bf{C}) &= \sum_{n=1}^N \left \lvert  \sum_{i=1}^{N_C} \left \lvert  \bf{Y}(\theta_n, \phi_n) \bf{C}_i  \right \rvert^2 - x_n \right \rvert^2 \\ 
%%%%%%%%%%%%%%%
 &= \sum_{n=1}^N \left \lvert  \textrm{Trace} \left (\bf{Q}   \bf{Y}^H(\theta_n, \phi_n)   \bf{Y}(\theta_n, \phi_n)   \right ) - x_n \right \rvert^2 \\ 
 %%%%%%%%%%%%%%%
& =  \sum_{n=1}^N \left \lvert   \bf{Y}(\theta_n, \phi_n) \bf{D} - x_n \right \rvert^2, \quad 
%%%%%%%%
\bf{Q} = \sum_{i=1}^{N_C} C_i C_i^H, \quad
\bf{D} = \sum_{i=1}^{N_C} \bf{C}_i \diamond \bf{\tilde{C}}_i,
\end{split}
$$

where the minimizer can be found by solving for positive semi-definite (PSD) $\bf{Q} \succeq 0$ via semi-definite programming (SDP), and recovering $\bf{C}_i = \sqrt{\lambda_i} \bf{v}_i$ for eigenvalue $\lambda_i$ and eigenvector $\bf{v}_i$ of the matrix $\bf{Q}$ solution.

**Example:** Fit density functions to numerical estimates of an optimal transport Wasserstein interpolation path between two radial basis density functions in `tst_sh_pdf_transport.m`:

```
tst_sh_pdf_transport('SlicedWass', 'FwdLeft', 'SqExp', 'FwdRight', 'SqExp', [0.25, 0.5, 0.75], 'preset_title', 'Left to Right');
```

  |<img src="./figs/figs_sh/fit/FwdLeftSqExp_FwdRightSqExp_1.png" width="280"/> | <img src="./figs/figs_sh/fit/FwdLeftSqExp_FwdRightSqExp_2.png" width="250"/> |<img src="./figs/figs_sh/fit/FwdLeftSqExp_FwdRightSqExp_3.png" width="250"/> |<img src="./figs/figs_sh/fit/FwdLeftSqExp_FwdRightSqExp_4.png" width="250"/> |<img src="./figs/figs_sh/fit/FwdLeftSqExp_FwdRightSqExp_5.png" width="250"/> |
  | :-: | :-: | :-: | :-: | :-: |   

### Magnitude Squared Least Squares

In the instance where the SDP solution $`\bf{Q}_{\ast}`$ to the SOMS objective is tight ($`\bf{Q}_{\ast}`$ is rank-1), then the solution recovers a smooth square-root function approximation in $\bf{C}$. This is equivalent to the minimizer of the magnitude squared (MS) objective given by

$$F_{MS}(\bf{C}) = \sum_{n=1}^N \left \lvert  \bf{Y}(\theta_n, \phi_n) \bf{D} - x_n \right \rvert^2, \quad \bf{D} = \bf{C} \diamond \bf{\tilde{C}},$$

which has applications for spatial phase-retrieval via  $\angle \bf{Y}(\theta, \phi) \bf{C}$. In this regard, SOMS is a semi-definite relaxation (SDR) of the MS problem, which in practice may not be tight. 

**Example**: Acoustic measurements of a transducer may be difficult to time-sync. Recover the acoustic phase response of a random field from its power response (magnitude squared target) across spherical coordinates in the function `plot_sh_fit_msq_phase_recovery_example.m` for low-order expansions ($L_D \leq 5$):

```
plot_sh_fit_msq_phase_recovery_example(4);
plot_sh_fit_msq_phase_recovery_example(5);
```
  |Reference | Magnitude Square Fit |
  | :-: | :-: |
  | <img src="./figs/figs_sh/fit/msq_ms_phaserecover_ref_4.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_ms_phaserecover_fit_4.png" width="350"/> |
  |<img src="./figs/figs_sh/fit/msq_ms_phaserecover_ref_5.png" width="350"/> | <img src="./figs/figs_sh/fit/msq_ms_phaserecover_fit_5.png" width="350"/> | 



### Mixture-of-Magnitude Squared Least Squares

We can further constrain the expansions in SOMS to be a mixture of expansions $\bf{C}_i = \bf{B} \bf{w}_i$ from an input dictionary $\bf{B} \in \mathbb{C}^{N_C \times K}$ weighted by $\bf{w}_i \in \mathbb{C}^{K \times 1}$. This has the benefit of model-order reduction when $K << N_C$ (no dependence on spatial resolution), and the practical application of array designs where the dictionary elements or columns of $\bf{B}$ are empirical quantities such as frequency responses in the far-field. The mixture-of-magnitude squared (MOMS) objective is given by

$$\displaystyle
\begin{split} 
F_{MOMS}(\bf{W}) &= \sum_{n=1}^N \left \lvert  \sum_{i=1}^{K} \left \lvert  \bf{Y}(\theta_n, \phi_n) \bf{C}_i  \right \rvert^2 - x_n \right \rvert^2 \\ 
%%%%%%%%%%%%%%%
 &= \sum_{n=1}^N \left \lvert  \textrm{Trace} \left (\bf{Q}   \bf{B}^H \bf{Y}^H(\theta_n, \phi_n)   \bf{Y}(\theta_n, \phi_n) \bf{B}  \right ) - x_n \right \rvert^2 \\ 
%%%%%%
& =  \sum_{n=1}^N \left \lvert   \bf{Y}(\theta_n, \phi_n) \bf{D} - x_n \right \rvert^2, \quad
%%%%%%
\bf{Q} = \sum_{i=1}^{K} w_i w_i^H, \quad
\bf{D} = \sum_{i=1}^{K} \bf{C}_i \diamond \bf{\tilde{C}}_i, 
%%%%%%%%
\end{split}$$

where $\bf{W} = [\bf{w}_1, … \bf{w}_K] \in \mathbb{C}^{N_C \times K}$, and can be similarly found via SDP for PSD constrained matrix $\bf{Q} \succeq 0$ where $\bf{w}_i = \sqrt{\lambda_i} \bf{v}_i$ for eigenvalue $\lambda_i$ and eigenvector $\bf{v}_i$ of the matrix $`\bf{Q}_{\ast}`$ solution.  We note that each $`\bf{Y}(\theta, \phi) \bf{B}w_i`$ would be a separate array design instance, and that the summation of  their power responses fits $x_n$ only for arrays driven by uncorrelated input signals. 

**Example**: Fit a collection of beam-formers across a 4-element array whose total array gain power is constant (omni-directional) across spherical coordinates in the function `plot_sh_fit_msq_omni_array_example.m`:

```
plot_sh_fit_msq_omni_array_example

% Output
beamformer_coeffs =

  -0.0469 + 0.0463i   0.0008 - 0.0008i   0.0426 - 0.0000i  -0.0366 - 0.0000i
   0.0011 + 0.0000i   0.0660 + 0.0000i   0.0426 + 0.0000i   0.0366 + 0.0000i
   0.0469 - 0.0463i  -0.0008 + 0.0008i   0.0426 + 0.0000i  -0.0366 - 0.0000i
  -0.0011 + 0.0000i  -0.0660 + 0.0000i   0.0426 + 0.0000i   0.0366 + 0.0000i
```

| Piston 1 Response | Piston 2 Response | Piston 3 Response | Piston 4 Response|
  | :-: | :-: | :-: | :-: | 
  | <img src="./figs/figs_sh/fit/omni_array_dic_1.png" width="350"/> | <img src="./figs/figs_sh/fit/omni_array_dic_2.png" width="350"/> | <img src="./figs/figs_sh/fit/omni_array_dic_3.png" width="350"/> | <img src="./figs/figs_sh/fit/omni_array_dic_4.png" width="350"/> |
| Beam 1 Response | Beam 2 Response | Beam 3 Response | Beam 4 Response |
 | <img src="./figs/figs_sh/fit/omni_array_beam_1.png" width="350"/> | <img src="./figs/figs_sh/fit/omni_array_beam_2.png" width="350"/> | <img src="./figs/figs_sh/fit/omni_array_beam_3.png" width="350"/> | <img src="./figs/figs_sh/fit/omni_array_beam_4.png" width="350"/> |

| Fitted Total Power Response |
| :-: | 
|<img src="./figs/figs_sh/fit/omni_array_fit.png" width="350"/>|

### Mixture Power Least Squares

In the instance where the SDP solution $`\bf{Q}_{\ast}`$ to the MOMS objective is tight, then the solution $`\bf{w}_{\ast}`$ is the minimizer of the mixture power (MP) objective given by

$$F_{MP}(\bf{w}) = \sum_{n=1}^N \left \lvert \bf{Y}(\theta_n, \phi_n) \bf{D} - x_n \right \rvert^2, \quad 
%%%%%%%
\bf{D} = \bf{C} \diamond \bf{\tilde{C}}, \quad 
%%%%%%%
\bf{C} = \bf{B} \bf{w},
$$

which are simple mixture weights of the dictionary elements in $\bf{B}$. This is useful in array designs applications where a desired directivity pattern is specified through $\bf{x}_n = \left \lvert P(\theta_n, \phi_n) \right \rvert^2$ for pressure $P$, and fitted via weightings of the array element’s far-field transfer functions $\bf{Y}(\theta_n, \phi_n) \bf{B} \bf{w}$, and subject to additional constraints such as electrical power or white-noise gain $\bf{w}^H \bf{w} \leq \tau$. 

**Example**: Fit a 12-element array that targets a omni-directional + spatial null beam-pattern in `plot_sh_fit_msq_null_array_example.m`:

```
plot_sh_fit_msq_null_array_example

% Output
beamformer_coeffs =

  -0.2430 - 0.1759i
   0.1542 + 0.1424i
  -0.0941 - 0.1044i
   0.0703 + 0.0831i
  -0.0362 - 0.0599i
   0.0309 + 0.0423i
  -0.0481 - 0.0234i
   0.0865 + 0.0127i
  -0.1953 - 0.0010i
   0.2026 - 0.0168i
  -0.4028 - 0.1353i
   0.3106 + 0.1784i

```

| Target Power Response | Fitted Power Response | Fitted Frequency Response | 
| :-: | :-: |:-: |
| <img src="./figs/figs_sh/fit/null_array_tgt_pow.png" width=“400”/> | <img src="./figs/figs_sh/fit/null_array_fit_pow.png" width=“400”/> | <img src="./figs/figs_sh/fit/null_array_fit_resp.png" width=“400”/> 

[^HANSEN_DPC]: Hansen, P. C. (1990). The discrete Picard condition for discrete ill-posed problems. BIT Numerical Mathematics, 30(4), 658-672.

[^DURAIS_TR]: Duraiswami, R., Zotkin, D. N., & Gumerov, N. A. (2004, May). Interpolation and range extrapolation of HRTFs. In Proc. IEEE ICASSP (Vol. 4, pp. 45-48). Montreal, Canada.

[^LUO_COV]: Luo, Y. (2021). [Spherical Harmonic Covariance and Magnitude Function Encodings for Beamformer Design](https://link.springer.com/article/10.1186/s13636-021-00230-7). EURASIP Journal on Audio, Speech, and Music Processing, 2021(1), 41.

[^RASMUSSEN_GP]: Rasmussen, C. E. (2003). Gaussian processes in machine learning. In Summer school on machine learning (pp. 63-71). Berlin, Heidelberg: Springer Berlin Heidelberg.