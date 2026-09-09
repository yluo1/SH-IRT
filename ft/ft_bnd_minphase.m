function [y, y_minphase, err, lambda] = ft_bnd_minphase(X_mag, num_taps, options)
%Fit FIR filter that minimizes the weighted least-squares error to the
%target minimum phase response of real-cepstrum
%subject to upper bound constraints on the frequency response's magnitude (unity)

%Solver: Second-order cone program

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X_mag:         [N x 1] target magnitude at DC to Fs uniform spaced frequency bins, positive valued, symmetric
%num_taps:      Number of taps in output FIR filter, must be <= N

%options:                Struct

%options.mode:           String, fitting method {'least_squares', 'minimax'}
%options.ub:             Scalar, magnitude upper bound, positive
%options.N_bins:         Number of frequency bins uniform spaced between DC, Nyquist 
%                        to constrain frequency response below unity

%options.wt:             [N x 1] frequency weights, [] for unity

%options.Fs:             Sampling rate
%options.enable_disp:    Logical, if true, plot FIR responses

%options.options_coneprog:  struct for coneprog options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%y:             [num_taps x 1] time-domain near minimum phase filter
%y_minphase:    [N x 1] time-domain minimum phase filter
%err:           [N x 1] frequency-response error  (minimum phase response - target response)
%lambda:        Augmented variable (objective)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Filter fit of varying number of taps

%X_mag_oneside = db2mag([0, -3, -0.1, -12])';
%X_mag = [X_mag_oneside; conj(flipud(X_mag_oneside(2:end-1)))];

%[y_ls, y_minphase]  = ft_bnd_minphase(X_mag, numel(X_mag), 'mode', 'least_squares', 'enable_disp', true);
%[y_ls]              = ft_bnd_minphase(X_mag, numel(X_mag)-1, 'mode', 'least_squares', 'enable_disp', true); %1 fewer tap
%[y_ls]              = ft_bnd_minphase(X_mag, numel(X_mag)-2, 'mode', 'least_squares', 'enable_disp', true); %2 fewer taps

%[y_minimax, y_minphase]     = ft_bnd_minphase(X_mag, numel(X_mag), 'mode', 'minimax',  'enable_disp', true);
%[y_minimax]                 = ft_bnd_minphase(X_mag, numel(X_mag)-1, 'mode', 'minimax',  'enable_disp', true); %1 fewer tap
%[y_minimax]                 = ft_bnd_minphase(X_mag, numel(X_mag)-2, 'mode', 'minimax',  'enable_disp', true); %2 fewer taps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Filter fit of varying frequency weighting
% 
% Fs = 48000;
% N_half = 64;
% rng(12);
% X_mag_oneside = db2mag(-rand(N_half, 1) * 20);
% X_mag = [X_mag_oneside; conj(flipud(X_mag_oneside(2:end-1)))];
% 
% hz = 0:(Fs / N):Fs; hz = hz(1:end-1)';
% wt_logQuad = ft_freq_wt(hz, 'logQuad');
% wt_logSqCov = ft_freq_wt(hz, 'logSqCov');
% wt_oct = ft_freq_wt(hz, 'octave', 'oct_freq', 100);
% wt_oct = wt_oct.^(0.25);
% 
% [y_ls, y_minphase]  = ft_bnd_minphase(X_mag, 64, 'mode', 'least_squares', 'enable_disp', true);
% [y_wls_oct]         = ft_bnd_minphase(X_mag, 64, 'mode', 'least_squares', 'wt', wt_oct, 'enable_disp', true);

% [y_minimax]         = ft_bnd_minphase(X_mag, 64, 'mode', 'minimax', 'enable_disp', true); 
% [y_minimax_oct]     = ft_bnd_minphase(X_mag, 64, 'mode', 'minimax', 'wt', wt_oct, 'enable_disp', true); 

arguments
    X_mag (:,1) double {mustBePositive} = [1];
    num_taps (1,1) double {mustBePositive, mustBeInteger}  = numel(X_mag);

    options.mode  {mustBeMember(options.mode, {'least_squares', 'minimax'})} =  'least_squares';

    options.ub (1,1) double {mustBePositive} = 1;
    options.N_bins (1,1) double {mustBePositive} = 512; 

    options.wt (:,1) double {mustBeNonnegative} = [];

    options.Fs (1,1) double {mustBePositive} = 48000;

    options.enable_disp (1,1) logical = false;

    options.options_coneprog = optimoptions("coneprog", MaxIterations=1000, ConstraintTolerance=1e-10, OptimalityTolerance=1e-8);
end

N = numel(X_mag);
if num_taps > N
    error('num_taps > N');
end
xhat = real(ifft(log(X_mag)));
odd = fix(rem(N,2));
wn = [1; 2 * ones([(N+odd)/2-1,1]) ; ones([1-rem(N,2),1]); zeros([(N+odd)/2-1,1])];
H_tgt = exp(fft(wn.*xhat)); %Minimum phase target

%Reference solution
y_minphase = real(ifft(H_tgt));

%Augmented variables
%x_aug = [x; lambda], lambda augmented error bound

% D = dftmtx(N);
% D = D(:, 1:num_taps);

% k_list = (0:(N-1))';
% n_list = 0:(num_taps-1);
% D_test = exp(-1i * 2*pi * k_list * n_list / N);
% norm(D - D_test)

k_list = (0:(N-1))';
n_list = 0:(num_taps-1);
D = exp(-1i * 2*pi * k_list * n_list / N); %[N x num_taps]
N_half = floor(N / 2) + 1; %DC to Nyquist

%Get frequency weights for error weighting
if isempty(options.wt)
    wt = ones(N, 1);
else
    wt = options.wt;
end

if strcmp(options.mode, 'minimax')

    %Minimize the maximum error: lambda is min(max(|err|))
    for n = 1:N_half %Single-sided spectrum

        Asoc_n = [[real(D(n, :)); imag(D(n, :))], zeros(2, 1)];
        bsoc_n = [real(H_tgt(n)); imag(H_tgt(n))];
        dsoc_n = [zeros(num_taps, 1); 1];
        gamma_n = 0;
       
        %Apply frequency weighting
        Asoc_n = wt(n) * Asoc_n;
        bsoc_n = wt(n) * bsoc_n;

        if n == 1
            socConstraints(1) = secondordercone(Asoc_n, bsoc_n, dsoc_n, gamma_n);
        else
            socConstraints(end+1) = secondordercone(Asoc_n, bsoc_n, dsoc_n, gamma_n);
        end
    end

else %'least_squares'

    %Minimize least square error: lambda is sqrt(sum(|err|.^2)))

    % Asoc = [[real(D); imag(D)], zeros(2 * N, 1)]; %Double spectrum
    % bsoc = [real(H_tgt); imag(H_tgt)];

    Asoc = [[real(D(1:N_half, :)); imag(D(1:N_half, :))], zeros(2 * N_half, 1)]; %Single-sided spectrum
    bsoc = [real(H_tgt(1:N_half)); imag(H_tgt(1:N_half))];
    dsoc = [zeros(num_taps, 1); 1];
    gamma = 0;

    %Apply frequency weighting
    W_mat = diag([wt(1:N_half); wt(1:N_half)]);
    Asoc = W_mat * Asoc;
    bsoc = W_mat * bsoc;

    socConstraints(1) = secondordercone(Asoc, bsoc, dsoc, gamma);
end

%Quadratic constraints
% D_FR = dftmtx(options.N_bins);
% D_FR = D_FR(:, 1:num_taps);
% N_bins_half = floor(options.N_bins / 2 ) + 1;

% k_list = (0:(options.N_bins-1))';
% n_list = 0:(num_taps-1);
% D_test = exp(-1i * 2*pi * k_list * n_list / options.N_bins);
% norm(D_FR - D_test)

k_list = (0:(options.N_bins-1))';
n_list = 0:(num_taps-1);
D_FR = exp(-1i * 2*pi * k_list * n_list / options.N_bins);
N_bins_half = floor(options.N_bins / 2 ) + 1;

for n = 1:N_bins_half
    Asoc_n = [[real(D_FR(n, :)); imag(D_FR(n, :))], zeros(2, 1)];
    bsoc_n = zeros(2, 1);
    dsoc_n = zeros(num_taps + 1, 1);
    gamma_n = -options.ub;
    socConstraints(end+1) = secondordercone(Asoc_n, bsoc_n, dsoc_n, gamma_n);
end

%Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

lb = [-inf(num_taps, 1); 0];
ub = [];

%Objective
f = [zeros(num_taps, 1); 1];

%Call solver
[x_aug] = coneprog(f, socConstraints, A, b, Aeq, beq, lb, ub, options.options_coneprog);
lambda = x_aug(end);
y = x_aug(1:end-1);

%Compute error
err = D(1:N_half, :) * y  - H_tgt(1:N_half);

%Display
if options.enable_disp

    freq_disp = logspace(log10(2), log10(options.Fs / 2), 512 );
    Y = freqz(y, 1, freq_disp, options.Fs);
    Y_minphase = freqz(y_minphase, 1, freq_disp, options.Fs);

    plot_cpx_resp([Y_minphase(:), Y(:)], freq_disp(:), 'resp_names', {'Minimum Phase', 'Near Minimum Phase'});

end