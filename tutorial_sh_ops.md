# Tutorial: Spherical Harmonic (SH) Operators

In this tutorial, we cover sample usage of basic operators on the spherical harmonic expansions:

* [Decode](#decode)
* [Addition](#addition)
* [Multiplication](#multiplication)
* [Conjugation](#conjugation)
* [Magnitude Squared](#magnitude-squared)
* [Convolution](#convolution)
* [Normalization](#normalization)
* [Rotation](#rotation)
* [Reflection](#reflection)
* [Resample](#resample)
* [Resize](#resize)
* [Filter Time Domain](#filter-time-domain)
* [Filter Frequency Domain](#filter-frequency-domain)

Sample codes are contained in the function `plot_sh_ops_example.m`.

## Decode

Evaluate (decode) spherical harmonic expansion at select spherical coordinates. 

Sample usage `sh_dec.m`:
```
    is_real = false;
    max_odr = 6;    
    C = sh_enc_proj(max_odr, pi/2, 0, is_real);

    N_val = 128;
    phi = linspace(-pi, pi, N_val)';
    theta = pi/2 * ones(N_val, 1); 

    f = sh_dec(C, theta, phi, is_real);    
    dB_lim = [-20, 16];    
    h_1 = sh_plt(C(:, 1), 'mercator', is_real, 'title_name', 'Field C', ...
        'dB_lim', dB_lim, 'disp_theta_phi', [theta, phi], 'disp_theta_phi_markersize', 4);

    fontsize = 14;
    h_2 = figure;
    plot(rad2deg(phi), real(f), 'linewidth', 1.5);
    xlabel('Azimuth \phi', 'fontsize', fontsize); ylabel('Re[ Y(\theta, \phi) C ]', 'fontsize', fontsize);
    title('C Decoded', 'fontsize', fontsize + 1); set(gca, 'fontsize', fontsize - 1);
    grid on; axis tight;
```

|<img src="figs/figs_sh/ops/ops_dec_1.png"  width="450"> | <img src="figs/figs_sh/ops/ops_dec_2.png"  width="450"> |
| :-: | :-: | 

## Addition

Add multiple spherical harmonic expansions using Matlab’s `+` operator.

Sample usage:
```
    is_real = false;
    max_odr = 6;    
    C = sh_enc_proj(max_odr, [pi/2 pi/2], [0, pi/3], is_real);
    C_sum = C(:, 1) + C(:, 2);

    dB_lim = [-20, 16];    
    h_1 = sh_plt(C(:, 1), 'mercator', is_real, 'title_name', 'Field A', 'dB_lim', dB_lim);
    h_2 = sh_plt(C(:, 2), 'mercator', is_real, 'title_name', 'Field B', 'dB_lim', dB_lim);
    h_3 = sh_plt(C_sum, 'mercator', is_real, 'title_name', 'Field A + B', 'dB_lim', dB_lim);
```

| <img src="figs/figs_sh/ops/ops_add_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_add_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_add_3.png"  width="450"> |
| :-: | :-: | :-: | 


## Multiplication
Multiply spherical harmonic expansions via Clebsch-Gordan coefficients.

Sample usage `sh_mul.m`:
```
    is_real = false;
    max_odr = 6;    
    C = sh_enc_proj(max_odr, [pi/2 pi/2], [0, pi/3], is_real);
    C_mul = sh_mul(C(:, 1), C(:, 2));
    
    dB_lim = [-20, 16];
    h_1 = sh_plt(C(:, 1), 'mercator', is_real, 'title_name', 'Field A', 'dB_lim', dB_lim);
    h_2 = sh_plt(C(:, 2), 'mercator', is_real, 'title_name', 'Field B', 'dB_lim', dB_lim);
    h_3 = sh_plt(C_mul, 'mercator', is_real, 'title_name', 'Field A \times B', 'dB_lim', dB_lim);
```

| <img src="figs/figs_sh/ops/ops_mul_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_mul_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_mul_3.png"  width="450"> |
| :-: | :-: | :-: | 

## Conjugation
Conjugate spherical harmonic expansion evaluations.

Sample usage `sh_conj.m`:
```
    is_real = false;
    max_odr = 6;    
    C = sh_rand(max_odr, 1, is_real);
    C_conj = sh_conj(C);

    h_1 = sh_plt(C, 'mercator', is_real, 'title_name', 'Field A');
    h_2 = sh_plt(C_conj, 'mercator', is_real, 'title_name', 'Field conj(A)');
```

| <img src="figs/figs_sh/ops/ops_conj_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_conj_2.png"  width="450">  |
| :-: | :-: | 

## Magnitude Squared
Compute magnitude squared expansion of spherical harmonic expansion.

Sample usage `sh_msq.m`:
```
    rng(5324);
    is_real = false;
    max_odr = 6;    
    C = sh_rand(max_odr, 1, is_real);
    D = sh_msq(C, is_real);

    dB_lim = [-12, 36];
    h_1 = sh_plt(C, 'mercator', is_real, 'title_name', 'Field A', 'dB_lim',  dB_lim);
    h_2 = sh_plt(D, 'mercator', is_real, 'title_name', 'Field |A|^2', 'dB_lim', dB_lim);
```

| <img src="figs/figs_sh/ops/ops_msq_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_msq_2.png"  width="450">  |
| :-: | :-: | 

## Convolution
Convolve multiple spherical harmonic expansions.

Sample usage `sh_conv.m`:
```
    rng(5324);
    is_real = false;
    max_odr = 12;    
    C = sh_rand(max_odr, 1, is_real);
    D = 5 * sh_enc_rbf('SqExp', max_odr, 0, 0, 0.5, is_real);
    E = sh_conv(C, D);

    dB_lim = [-12, 24];
    h_1 = sh_plt(C, 'mercator', is_real, 'title_name', 'Field C', 'dB_lim',  dB_lim);
    h_2 = sh_plt(D, 'mercator', is_real, 'title_name', 'Field D', 'dB_lim', dB_lim);
    h_3 = sh_plt(E, 'mercator', is_real, 'title_name', 'Field Conv(C, D)', 'dB_lim', dB_lim);
```

| <img src="figs/figs_sh/ops/ops_conv_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_conv_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_conv_3.png"  width="450"> |
| :-: | :-: | :-: | 

## Normalization
Normalize N3D spherical harmonic expansion by its integral sum or power.

Sample usage `sh_nrm.m`:
```
    rng(5324);
    is_real = false;
    max_odr = 6;    
    C = sh_msq(sh_rand(max_odr, 1, is_real), is_real);
    D = sh_nrm(C, 'Sum');
    E = sh_nrm(C, 'Pow');

    h_1 = sh_plt(C, 'mercator', is_real, 'disp_phase', false, 'title_name', 'Field C', 'dB_lim',  [-20, 36]);
    h_2 = sh_plt(D, 'mercator', is_real, 'disp_phase', false, 'title_name', "Field norm(C, 'Sum')", 'dB_lim', [-40, 6]);
    h_3 = sh_plt(E, 'mercator', is_real, 'disp_phase', false, 'title_name', "Field norm(C, 'Pow')", 'dB_lim', [-40, 6]);
```

| <img src="figs/figs_sh/ops/ops_norm_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_norm_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_norm_3.png"  width="450"> |
| :-: | :-: | :-: | 

## Rotation

Rotate spherical harmonic expansion coefficients by zyx (yaw-pitch-roll), Tait-Bryan rotation, intrinsic.

Sample usage `sh_rot.m`:
```
    rng(1123);
    is_real = false;
    max_odr = 1;    
    C = sh_rand(max_odr, 1, is_real);
    C_z90 = sh_rot(C, deg2rad([90, 0, 0]), is_real);
    C_y90 = sh_rot(C, deg2rad([0, 90, 0]), is_real);
    C_x90 = sh_rot(C, deg2rad([0, 0, 90]), is_real);

    dB_lim = [-30, 6];
    h_1 = sh_plt(C, 'mercator', is_real, 'title_name', 'Field C', 'dB_lim',  dB_lim);
    h_2 = sh_plt(C_z90, 'mercator', is_real, 'title_name', 'Field Rotate(C, z=90)', 'dB_lim',  dB_lim);
    h_3 = sh_plt(C_y90, 'mercator', is_real, 'title_name', 'Field Rotate(C, y=90)', 'dB_lim',  dB_lim);
    h_4 = sh_plt(C_x90, 'mercator', is_real, 'title_name', 'Field Rotate(C, x=90)', 'dB_lim',  dB_lim);
```

| <img src="figs/figs_sh/ops/ops_rot_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_rot_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_rot_3.png"  width="450"> | <img src="figs/figs_sh/ops/ops_rot_4.png"  width="450"> |
| :-: | :-: | :-: | :-: | 

## Reflection

Reflect spherical harmonic expansion over a cardinal plane.

Sample usage `sh_refl.m`:
```
    rng(1953);
    is_real = false;
    max_odr = 1;    
    C = sh_rand(max_odr, 1, is_real);
    C_med = sh_refl(C, 'median', is_real);
    C_trans = sh_refl(C, 'transverse', is_real);
    C_front = sh_refl(C, 'frontal', is_real);

    dB_lim = [-30, 6];
    h_1 = sh_plt(C, 'mercator', is_real, 'title_name', 'Field C', 'dB_lim',  dB_lim);
    h_2 = sh_plt(C_med, 'mercator', is_real, 'title_name', 'Field Reflect(C, Median)', 'dB_lim',  dB_lim);
    h_3 = sh_plt(C_trans, 'mercator', is_real, 'title_name', 'Field Reflect(C, Transverse)', 'dB_lim',  dB_lim);
    h_4 = sh_plt(C_front, 'mercator', is_real, 'title_name', 'Field Reflect(C, Frontal)', 'dB_lim',  dB_lim);
```

| <img src="figs/figs_sh/ops/ops_refl_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_refl_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_refl_3.png"  width="450"> | <img src="figs/figs_sh/ops/ops_refl_4.png"  width="450"> |
| :-: | :-: | :-: | :-: | 


## Translation

Translate spherical harmonic expansion's origin to a new coordinate. Refit to higher max-order expansions.

Sample usage `sh_translate.m`:
```
    % Generate SH spherical piston normalized on-axis and translate the expansion
    Fs = 8000;
    is_real = false;
    N_taps = 128;
    w = linspace(0, 2 * pi, N_taps + 1); w = w(1:end-1);
    t = (0:(N_taps-1))/Fs;
    omega = w * Fs;
    freq = omega / (2 * pi);
    
    C = sh_enc_pist_sphere(30, pi/2, 0, 0.25, deg2rad(20), freq, 4, 'normalize_piston_axis', true);
    [D, h_1] = sh_translate(C, omega, is_real, [- 0.5/1000 * 343, 0, 0], 'enable_disp', true, 'disp_dB_lim', [-20, 6]); % 1/2 ms in -x direction
        
    C_td = ifft(C, [], 2, 'symmetric');
    h_2 = sh_plt(C_td, 'cardinal', is_real, 't', t, 'title_name', 'Original');
    D_td = ifft(D, [], 2, 'symmetric');
    h_3 = sh_plt(D_td, 'cardinal', is_real, 't', t, 'title_name', 'Translated');
```

| <img src="figs/figs_sh/ops/ops_trans_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_trans_2.png"  width="450"> |
| :-: | :-: |
|<img src="figs/figs_sh/ops/ops_trans_3.png"  width="450"> | <img src="figs/figs_sh/ops/ops_trans_4.png"  width="450"> |

## Resample

Regularized least squares fit FIR filters to SH coefficients per basis over frequency. Then compute frequency response of FIRs at target angular frequency for resampling.

Sample usage `ft_fit_sh.m` and `sh_filter_val.m`:
```
    % Resample log-frequency response of spherical piston to linear-frequency
    Fs = 48000;
    hz = logspace(log10(10), log10(Fs / 2), 256 + 1);
    w  = hz / Fs * (2 * pi);
    r = 3; % Evaluation distance
    radius = 0.25; % Spherical baffle radius
    C = sh_enc_pist_sphere(30, pi/2, 0, radius, deg2rad(20), hz, r);
    
    is_real = false;
    num_taps = 1024;
    h_fit_exp = ft_fit_sh(C, w, is_real, num_taps, 'exp', 'exp_mu', (r - radius) / 343 * Fs, 'exp_std', 100, ...
                      'enable_disp', false, 'disp_dB_lim', [-20, 60], 'disp_err_dB_lim', [-80, 10]);
    
    w_lin = linspace(0, pi, 1024);
    C_lin = sh_filter_val(h_fit_exp, w_lin);
    
    h_1 = sh_plt(C, 'horizontal', is_real,  'hz', w / (2 * pi) * Fs, 'dB_lim', [-20, 60], 'title_name', 'Original');
    h_2 = sh_plt(C_lin, 'horizontal', is_real,  'hz', w_lin / (2 * pi) * Fs, 'dB_lim', [-20, 60], 'title_name', 'Resampled');
```

| <img src="figs/figs_sh/ops/ops_resamp_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_resamp_2.png"  width="450">  |
| :-: | :-: | 

## Resize

Resize max order of spherical harmonic expansions via truncation or zero-padding.

Sample usage `sh_resize.m`:
```
    rng(534);
    P = 4;
    M = 1;
    is_real = false;
    C = sh_rand(P, M, is_real);
    D = sh_resize(C, 2); % Truncated 
    E = sh_resize(C, 6); % Zero-padded
    
    h_1 = sh_plt(C, 'mercator', is_real);
    h_2 = sh_plt(D, 'mercator', is_real);
    h_3 = sh_plt(E, 'mercator', is_real); % Match C
```

| <img src="figs/figs_sh/ops/ops_resize_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_resize_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_resize_3.png"  width="450"> |
| :-: | :-: | :-: | 

## Filter Time domain

Filter spherical harmonic expansion coefficients in the time-domain by numerator/denominator coefficients or second-order-sections.

Sample usage `sh_filter.m`:
```
    Fs    = 48000;
    M_w   = 256;
    freq  = logspace(log10(100), log10(Fs / 2), M_w);
    w     = freq / Fs * 2 * pi;
    
    fc = 2000;
    [b,a]     = butter(4,fc/(Fs/2));
    [z,p,k]   = butter(4,fc/(Fs/2));
    sos       = zp2sos(z,p,k);
    
    rng(234);
    max_odr = 3;
    num_func = 1024;
    
    is_real = true;
    C = sh_rand(max_odr, num_func, is_real, true);
    C_butter_ba  =  sh_filter(C, b, a);
    C_butter_sos =  sh_filter(C, sos);
    
    D     = sh_freqz(C, w);
    D_ba  = sh_freqz(C_butter_ba, w);
    D_sos = sh_freqz(C_butter_sos, w);
    
    dB_lim = [-10, 50];
    h_1 = sh_plt(D, 'horizontal', is_real, 'hz', freq, 'title_name', 'Real Random Field', 'dB_lim', dB_lim);
    h_2 = sh_plt(D_ba, 'horizontal', is_real, 'hz', freq, 'title_name', 'Real Random Field * LPF (b/a)', 'dB_lim', dB_lim);
    h_3 = sh_plt(D_sos, 'horizontal', is_real, 'hz', freq, 'title_name', 'Real Random Field * LPF SOS', 'dB_lim', dB_lim);
```

| <img src="figs/figs_sh/ops/ops_flt_t_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_flt_t_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_flt_t_3.png"  width="450"> |
| :-: | :-: | :-: | 

## Filter Frequency Domain

Filter spherical harmonic expansion coefficients in the frequency domain (Hardamard product) by numerator/denominator coefficients or second-order-sections.

Sample usage `sh_filter_freq.m`:
```
    Fs    = 48000;
    freq  = logspace(log10(100), log10(Fs / 2), 256);
    w     = freq / Fs * 2 * pi;
    
    fc = 2000;
    [b,a]     = butter(4,fc/(Fs/2));
    [z,p,k]   = butter(4,fc/(Fs/2));
    sos       = zp2sos(z,p,k);
    
    C = sh_enc_pist_sphere(40, pi/2, 0, 0.25, deg2rad(20), freq, 4);
    C_butter_ba  =  sh_filter_freq(C, w, b, a);
    C_butter_sos =  sh_filter_freq(C, w, sos);
    
    dB_lim = [-10, 70];
    h_1 = sh_plt(C, 'horizontal', false, 'hz', freq, 'title_name', 'Piston on Spherical Baffle', 'dB_lim', dB_lim);
    h_2 = sh_plt(C_butter_ba, 'horizontal', false, 'hz', freq, 'title_name', 'Piston on Spherical Baffle * LPF (b/a)', 'dB_lim', dB_lim);
    h_3 = sh_plt(C_butter_sos, 'horizontal', false, 'hz', freq, 'title_name', 'Piston on Spherical Baffle * LPF SOS', 'dB_lim', dB_lim);
```

| <img src="figs/figs_sh/ops/ops_flt_f_1.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_flt_f_2.png"  width="450">  |  <img src="figs/figs_sh/ops/ops_flt_f_3.png"  width="450"> |
| :-: | :-: | :-: | 