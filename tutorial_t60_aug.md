# Tutorial: Spatial Room Impulse Response T60 Augmentation

In this tutorial, we cover sound decay time augmentation of room impulse responses (RIRs) from our paper 
>Yuancheng Luo, "Fast Time-Varying Exponentiated Convolution Methods for Generative Direction Dependent Reverberation", Proceedings of the 161th Audio Engineering Society Convention.

The high-level steps are as follows:
* [Sampling T60 from Gaussian Processes (GPs):](#sampling-t60-functions-from-gaussian-processes) Sample or design a sound decay time function $T_{60}(\omega, \theta, \phi)$ of angular frequency $\omega$, and spherical coordinates $(\theta, \phi)$ co-latitude and azimuth respectively
	* [GP prior mean specifications](#mean-function-prior-specifications)
    * [GP prior covariance specifications](#covariance-function-prior-specifications)
    * [Drawing T60 samples from GP prior and posterior distributions](#sampling-from-gp-prior-and-posterior-distributions)
    * [Optimizing GP covariance hyper parameters](#optimizing-gaussian-process-covariance-function-hyper-parameters)
* [Exponentiating FIR Optimization:](#fitting-fir-to-t60-functions) Fit short finite impulse response (FIR) exponentiating filters $\bf{g}$ to sampled $T_{60}(\omega, \theta, \phi)$ functions
* [Specifying or generating colorless RIR:](#generating-colorless-room-impulse-responses) Specify an input RIR or generate a colorless (constant T60) RIR $\bf{h}$
* [Augmenting RIR:](#applying-time-varying-exponentiated-convolution) Apply time-varying exponentiated convolution $\textbf{f} = f(\bf{g}, \bf{h})$

## Sampling T60 Functions from Gaussian Processes
We can sample smooth T60 functions $T_{60}(\omega, \theta, \phi)$ of frequency and spherical coordinates from Gaussian processes defined by prior mean function $\overline{T}_{60}(\omega)$, prior covariance function $k(\bf{x},\bf{x}')$, and observed log-T60 times $\bf{y} = [y_1, … , y_S]$ at $\bf{X} = \left \lbrace \bf{x}_1, …, \bf{x}_S \right \rbrace$ given by

$$y_n = \log T_{60}(\underline{\omega_n}, \underline{\theta_n}, \underline{\phi_n} ) - \log \overline{T}_{60}(\underline{\omega_n}), \quad  \bf{x}_n = (\underline{\omega_n}, \underline{\theta_n}, \underline{\phi_n}).
$$

The supported mean and covariance functions are specified as follows:

### Mean Function Prior Specifications

A simple T60 prior mean is the power law function given by

$$\overline{T}_{60}(\omega) =  \alpha  \max \left ( 1, \nu(\omega)  \right ) ^{-\beta}, \quad \nu(\omega) = \frac{\omega}{2 \pi}, $$

where $\alpha$ is the T60 (seconds) at 0 Hz, $\beta$ is the decay exponent, and $\nu(\omega)$ is ordinary frequency. Increasing the decay rate $\beta$ inflects the T60 curve towards 0 seconds for larger frequencies in the following figure: 

```
gp_plt('mu_pow');
```
<img src="./figs/figs_t60/mu_pow_1.png" alt="Power Function" width="600"/>

The power law function is monotonic decreasing for positive $\beta$ and whose derivative converges to $0$ for increasing frequency. We can introduce a 3rd shape parameter to control the tilt following an idealized low-pass filter response design given by

$$\overline{T}_{60}(\omega) =  \frac{\alpha}{1 + \left ( \frac{\nu(\omega)}{f_c}  \right )^{\beta}}, $$

where $f_c$ is the cross-over frequency relative to the flat T60 curve $(\beta = 0)$ where the absolute derivative is maximal, as shown in the following figures:

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

$$d(\theta, \phi, \theta', \phi') = 2 \sin \left ( \frac{ \left | \cos^{-1} (\bf{v}^T \bf{v}') \right |   }{2} \right ) = \left \| \bf{v} - \bf{v}'  \right \|_2 .$$

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

$$k(\bf{x},\bf{x}') = \sigma^2 \exp \left (  - \frac{d^2(\theta, \phi, \theta', \phi') }{2 \ell_c^2}  \right )  \exp \left (  - \frac{ \left | log(f) - log(f') \right | ^2 }{2 \ell_f^2}  \right ) ,$$

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

### Sampling from GP Prior and Posterior Distributions

We can sample T60 functions from either a GP prior or posterior distribution via the function `gp_t60_sample.m`. Let us walk through several cases from the function `plot_gp_T60_prior_post_example.m`.

* Drawing independent and identically distributed (IID) sample log-T60 functions from a GP prior distribution at a common spherical coordinate:
  ```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Setup GP prior mean and covariance
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  options_mu  = gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.25, 'mu_fc', 8000);
  options_cov = gp_cov_opts('cov_sigma', sqrt(2)/2, 'cov_gamma', 2/3, 'cov_ell', 343 * 1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Single point spherical coordinate evaluation grid, sampled 4 times, GP prior
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Fs = 48000; %Sample rate
  
  % Setup evaluation grid
  N_B = 128; % Number of angular frequencies
  N_E = 1;   % Number of spherical coordinates
  
  omega        = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
  [theta, phi] = sh_fib(N_E);
  
  % Sample 4 functions
  num_evals = 4;
  rng(21136 + 12); %Seed RNG
  [log_T60, h_fig_prior] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
          'options_mu', options_mu, 'options_cov', options_cov, ...
          'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1.75], 'disp_legend_num_cols', 1));
  ```
	<img src="./figs/figs_t60/sample_GP_prior.png" alt="Sample T60s drawn from GP prior" width="480"/>

	where by virtue of modeling log-T60 functions, the variance hyper parameter  $\sigma^2$ belonging to the covariance function is multiplicative w.r.t. the prior mean $\overline{T}_{60}(\omega)$, and the sampled functions are non-negative.

* Drawing a single log-T60 field at varying spherical coordinates on the horizontal plane from a GP prior distribution:
  ```
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Azimuth plane spherical coordinate evaluation grid, sampled once, GP prior
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup evaluation grid
    N_B = 64;
    omega   = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
    phi     = deg2rad([0:30:150])';
    theta   = pi/2 * ones(numel(phi), 1);
    
    % Sample 1 function
    rng(11136 + 13); 
    num_evals = 1;
    [log_T60, h_fig_prior_grid] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
            'options_mu', options_mu, 'options_cov', options_cov, ...
            'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 0.75], 'disp_legend_loc', 'northeast', 'disp_legend_num_cols',  2));
  ```
	<img src="./figs/figs_t60/sample_GP_prior_grid.png" alt="Sample T60s drawn from GP prior grid" width="480"/>

	where the T60s highly covary for lower frequencies by virtue of the non-stationary covariance function. 

* Specifying observed log-T60 via options in `gp_obs_opts.m` and drawing IID log-T60 functions from the GP posterior distribution:
	```
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Single point spherical coordinate evaluation grid, sampled 4 times, GP posterior
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Setup evaluation grid
  N_B = 128;
  N_E = 1;
  omega        =  2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
  [theta, phi] = sh_fib(N_E);
  
  % Observations
  N_S = 3; %3 observations in frequency
  omega_obs = logspace(log10(20), log10(Fs/2), N_S)' * 2 * pi;
  
  obs = gp_obs_opts('omega', omega_obs, ...
                      'theta', theta(1) * ones(N_S, 1), ...
                      'phi', phi(1) * ones(N_S, 1), ...
                      'T60', mu_pow(omega_obs, 2.5, 0.4), ...
                      'log_noise_std', 0.05 * ones(N_S, 1) ... %5 percent
                      );
  
  % Sample function
  rng(21136 + 18);
  num_evals = 4; %Sample 4 functions
  [log_T60, h_fig_post] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
              'options_mu', options_mu, 'options_cov', options_cov, ...
              'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1], 'disp_legend_loc', 'northeast', 'disp_legend_num_cols', 1));
    ```
	
	<img src="./figs/figs_t60/sample_GP_post.png" alt="Sample T60s drawn from GP posterior" width="480"/>

	where `obs` struct contains the three observed log-T60s specified at $[20, 692, 24000]$ Hz on a single spherical coordinate. The sampled log-T60 functions from the GP posterior distribution highly covary at the observed T60 frequencies due to the small log-noise standard deviation specification (5%).

* Drawing a single log-T60 field at varying spherical coordinates on the horizontal plane from a GP posterior distribution:
	```
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Azimuth plane spherical coordinate evaluation grid, sampled once, GP posterior
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Setup evaluation grid
  N_B = 64;
  omega   = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
  phi     = deg2rad([0:30:150])';
  theta   = pi/2 * ones(numel(phi), 1);
  
  %Observations
  N_S = 3;
  omega_obs = logspace(log10(20), log10(Fs/2), N_S)' * 2 * pi;
  
  obs = gp_obs_opts('omega', omega_obs, ...
                      'theta', deg2rad(90) * ones(N_S, 1), ...
                      'phi', deg2rad(0) * ones(N_S, 1), ...
                      'T60', mu_lpf(omega_obs, 1, 1, 4000), ...
                      'log_noise_std', 0.005 * ones(N_S, 1) ... %0.5 percent
                      );
  
  options_disp_gp = gp_disp_opts('disp_legend_samples', false, ...
      'disp_legend_mean', false, 'disp_legend_var', false, 'disp_legend_loc', 'northeast', ...
      'disp_legend_num_cols', 1, 'disp_ylim', [0, 1.1], ...
      'disp_sample_eval', true, 'disp_mean', false, 'disp_var', true, ...
      'disp_legend_compact', true, 'disp_position', [100, 100, 560 * 0.825, 480 * 0.666], ...
      'disp_legend_transparency', 0.75, 'disp_colororder', 'gem12', ...
      'disp_var_transparency', 0.025, ...
      'disp_sample_stride', 1);
  
  % Sample 1 function
  rng(11136 + 13); 
  num_evals = 1;
  [log_T60, h_fig_posterior_grid] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
          'options_mu', options_mu, 'options_cov', options_cov, ...
          'enable_disp', true, 'options_disp', options_disp_gp);
    ```
	<img src="./figs/figs_t60/sample_GP_post_grid.png" alt="Sample T60s drawn from GP posterior" width="480"/>

	where the sampled log-T60s revert back to the prior mean function for spherical coordinates far away from the observed coordinates at $(\theta = 90^{\circ}, \phi = 0^{\circ})$.

### Optimizing Gaussian Process Covariance Function Hyper Parameters

We can optimize the covariance function’s hyper parameters by maximizing the marginal data likelihood following the example in `plot_gp_optimize.m`. Let us draw some sample log-T60 functions from the following GP prior distribution:

```
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup GP prior mean and covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_mu  = gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.25, 'mu_fc', 8000);
options_cov = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', 0.5, 'cov_gamma', 0.75, 'cov_ell', 2 * 343);

Fs = 48000; %Sample rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuth plane spherical coordinate evaluation grid, sampled once, GP prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup evaluation grid
N_B     = 16;
N_E     = 24;
omega   = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
[theta, phi] = sh_fib(N_E);

options_disp = gp_disp_opts('disp_ylim', [0, 0.6], 'disp_legend_loc', 'southoutside', ...
    'disp_legend_num_cols',  4, 'disp_legend_compact', true, ...
    'disp_position', [100, 100, 800, 800], 'disp_legend_samples', false);

% Sample 1 field
rng(1213); 
num_evals = 1;

[log_T60, h_fig_prior_grid] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
        'options_mu', options_mu, 'options_cov', options_cov, ...
        'enable_disp', true, 'options_disp', options_disp); % [N_B x N_E x num_evals]
```

<img src="./figs/figs_t60/opt_GP_prior_grid.png" alt="Sample T60s drawn from GP posterior" width="480"/>

We now specify the drawn samples as the observed log-T60 and maximize the log-marginal likelihood w.r.t. the covariance hyper parameters as follows:

```
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega_grid = repmat(omega, [1, N_E]);
theta_grid = repmat(theta(:)', [N_B, 1]);
phi_grid = repmat(phi(:)', [N_B, 1]);
obs = gp_obs_opts('omega', omega_grid(:), 'theta', theta_grid(:), 'phi', phi_grid(:), 'T60', exp(log_T60(:)), ...
    'log_noise_std', 0.002 * ones(N_B * N_E, 1)); %0.2 percent 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit covariance hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_fmincon = optimoptions("fmincon", SpecifyObjectiveGradient=true, Display="iter", checkGradients=false, ...
        ScaleProblem=true, ...
        FunctionTolerance=1e-8, ConstraintTolerance=1e-8, OptimalityTolerance=1e-10, StepTolerance=1e-8, ...
        MaxIterations=1000, MaxFunctionEvaluations=10000);

options_cov_0 = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', ...
    'cov_sigma_lim', [0.2, 2], 'cov_gamma_lim', [0.4, 2], 'cov_ell_lim', [100, 1000]);

options_cov_fit = gp_t60_optimize(obs, 'options_mu', options_mu, 'options_cov', options_cov_0, ... 
    'options_fmincon', options_fmincon, 'num_start', 1);

% Compare covariance hyperparameters in structs
disp([newline, 'Prior cov_sigma: ', num2str(options_cov.cov_sigma)]);
disp(['Prior cov_ell: ', num2str(options_cov.cov_ell)]);
disp(['Prior cov_gamma: ', num2str(options_cov.cov_gamma)]);
```

The optimization converges and the fitted covariance hyper parameters are closer to that of the GP prior where the log-T60 functions were originally drawn as follows:

```
log-marginal likelihood: 1162.794
cov_sigma: 0.50844
cov_ell: 687.5028
cov_gamma: 0.74528

Prior cov_sigma: 0.5
Prior cov_ell: 686
Prior cov_gamma: 0.75
```

A sampled log-T60 field from the GP posterior distribution at twice the number of spherical coordinates and higher frequency resolution yields the following:
```
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample field from posterior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Draw sample log-T60 from fitted GP
options_disp_fitted = gp_disp_opts('disp_ylim', [0, 0.6], 'disp_legend_loc', 'southoutside', ...
    'disp_legend_num_cols', 4, 'disp_legend_compact', true, ...
    'disp_position', [100, 100, 800, 800], 'disp_legend_samples', false, ...
    'disp_legend_var', false, 'disp_legend_mean', false, ...
    'disp_marker_size', 4);

omega_fitted = 2 * pi * logspace(log10(20), log10(Fs/2), 100)';
[theta_fitted, phi_fitted] = sh_fib(2 * N_E);

[log_T60_fitted, h_fig_post_grid] = gp_t60_sample(omega_fitted, theta_fitted, phi_fitted, num_evals, obs, ...
    'options_mu', options_mu, 'options_cov', options_cov_fit, ...
    'enable_disp', true, 'options_disp', options_disp_fitted);
```
<img src="./figs/figs_t60/opt_GP_post_grid.png" alt="Sample T60s drawn from GP posterior" width="480"/>


## Fitting FIR to T60 Functions

## Generating Colorless Room Impulse Responses

## Applying Time-Varying Exponentiated Convolution

[^PACIOREK_NS]: Paciorek, C., & Schervish, M. (2003). Nonstationary covariance functions for Gaussian process regression. Advances in neural information processing systems, 16.

