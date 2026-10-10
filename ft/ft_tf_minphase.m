function [b_mp, a_mp, b_ap, a_ap] = ft_tf_minphase(b, a, options)
% Decompose transfer function (numerator b / denominator a) into minimum phase
% and all-pass filter coefficients

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%options:               struct
%options.margin:        Scalar, projects any poles and zeros with modulus within 
%                       margin [1/margin, margin] onto nearest margin
%options.enable_disp:   Logical, if true, display frequency responses in fvtool
%options.Fs:            Scalar, sampling rate for display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%b_mp:  [1 x *] Numerator coefficients of minimum phase filter
%a_mp:  [1 x *] Denominator coefficients of minimum phase filter
%b_ap:  [1 x *] Numerator coefficients of all-pass filter
%a_ap:  [1 x *] Denominator coefficients of all-pass filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Minimum-phase allpass decomposition of IIR

%Fs = 48000;
%rng(5234);
%b = randn(1, 5);
%a = randn(1, 3);

%[b_mp, a_mp, b_ap, a_ap] = ft_tf_minphase(b, a, 'enable_disp', true);
%ft_tf_minphase(b, a, 'margin', 2, 'enable_disp', true); % With margin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Minimum-phase allpass decomposition of FIR

%Fs = 48000;
%rng(5234);
%b = randn(1, 10);
%a = 1;

%[b_mp, a_mp, b_ap, a_ap] = ft_tf_minphase(b, a, 'enable_disp', true);
%ft_tf_minphase(b, a, 'margin', 2, 'enable_disp', true); % With margin

arguments
    b (1,:) double = 1;
    a (1,:) double = 1;

    options.margin (1,1) double {mustBeNonnegative} = 1;
    options.enable_disp (1,1) logical = false;
    options.Fs (1,1) double {mustBeNonnegative} = 48000;
end

assert(options.margin >= 1, 'Require margin >= 1');

[z, p, k] = tf2zpk(b, a);

N_z = numel(z); % Number of zeros
N_p = numel(p); % Numoer of poles

% Flip poles and zeros outside of complex unit circle into unit-circle
z_mp = [];
p_mp = [];

z_ap = [];
p_ap = [];

for i = 1:N_z % Iterate over zeros
    
    z_i = z(i);
    d_i = abs(z_i);


    if d_i >= 1 / options.margin && d_i <= options.margin % Zero is within margin of exclusion
        if d_i > 1 % Outside
            z_i = z_i / d_i * options.margin;
        else % Inside or on unit circle
            z_i = z_i / d_i / options.margin;            
        end
    end

    if d_i >= 1   % Zero is outside or on unit circle

        z_conj_inv = 1 ./ conj(z_i);      
        
        z_mp    = [z_mp; z_conj_inv];
        z_ap    = [z_ap; z_i];
        p_ap    = [p_ap; z_conj_inv];

    else % Zero is inside unit circle

        z_mp    = [z_mp; z_i];

    end
end

for i = 1:N_p % Iterate over poles

    p_i = p(i);
    d_i = abs(p_i);
 
    if d_i >= 1 / options.margin && d_i <= options.margin % Pole is within margin of exclusion
        if d_i > 1 % Outside
            p_i = p_i / d_i * options.margin;
        else % Inside or on unit circle
            p_i = p_i / d_i / options.margin;            
        end
    end

    if abs(p_i) >= 1 % Pole is outside or on unit circle
        p_conj_inv = 1 ./ conj(p_i);  
        
        p_mp    = [p_mp; p_conj_inv];
        z_ap    = [z_ap; p_conj_inv];
        p_ap    = [p_ap; p_i];

    else % Pole is inside unit circle

        p_mp    = [p_mp; p_i];

    end
end

k_ap = real(prod(p_ap));
[b_mp, a_mp] = zp2tf(z_mp, p_mp, k / k_ap);
[b_ap, a_ap] = zp2tf(z_ap, p_ap, k_ap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.enable_disp && coder.target("Matlab")

    h_f_orig = fvtool(b, a, 'Fs', options.Fs);
    legend(h_f_orig, 'Original Transfer Function');

    h_f_mp_ap = fvtool(b_mp, a_mp, b_ap, a_ap, 'Fs', options.Fs);
    legend(h_f_mp_ap, 'Minimum Phase', 'All-Pass');

end