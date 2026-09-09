function [C, t] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, options)
%Generate all pair-wise combinations of spherical harmonic expansions
%between acoustic source and receiver in image-source-receiver room model.
%See paper reference:

%Luo, Y. and Kim, W., 2021, January. 
%Fast source-room-receiver acoustics modeling.
%In 2020 28th European Signal Processing Conference (EUSIPCO) (pp. 51-55). IEEE.

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%max_src_odr:   Max SH order of source expansions
%max_rec_odr:   Max SH order of receiver expansions
%is_real:       Logical, if true, C is real, otherwise, C is complex, non-negative

%T:             Scalar, time (sec)
%s:             [1 x 3] Source offset from origin (meters)
%r:             [1 x 3] Receiver offset from origin (Meters)
%l:             [1 x 3] Orthotope dimensions (meters)
%gamma_pos:     [L x 3] Filter (L-tap length) reflection IR for wall on +axis
%gamma_neg:     [L x 3] Filter (L-tap length) reflection IR for wall on -axis

%options,       struct
%options.Fs:                    Sampling rate
%options.kernel_sample_width:   Lanczos kernel half-window size
%options.jitter_coord_bnd:      [1 x 2] (min, max) bounds on jitter to image-source and receiver coordinates
%options.jitter_srand:          Random seed for jitter

%options.mode_enc:              String, encoding mode {'proj', 'proj_msq'}
%                                   'proj':         kernel projection onto SH
%                                   'proj_msq':     magnitude squared kernel projection onto SH
%                                   'proj_msq_mex': magnitude squared kernel projection onto SH, mex version
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:     [(max_src_odr + 1)^2 x (max_rec_odr + 1)^2 x M]
%t:     [1 x M] time (sec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate SH expansion of source-receiver expansions in shoebox room

% max_src_odr = 1;
% max_rec_odr = 5;
% is_real = true;
% T = 0.2;
% s = [1 2 1];
% r = [2 0.7 1];
% l = [5 6 3];
% gamma_pos = [0.8, 0.8, 0.8];
% gamma_neg = [0.8, 0.8, 0.8];

% [C, t] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, 'jitter_coord_bnd', [1 1] * 1e-3);
% 
% sh_plt(squeeze(C(1, :, :)), 'horizontal', is_real, 'disp_phase', false, 'dB_lim', [-120, -10], 't', t, 'disp_xaxis_ker_size', 1024);
% sh_plt(squeeze(C(1, :, ceil(0.15 * 48000)))', 'mercator', is_real);

arguments
    max_src_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 0;
    max_rec_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 3;
    is_real (1,1) logical = false;

    T (1,1) double {mustBeNonnegative} = 0.2;
    s (1,3) double = [1 2 1];
    r (1,3) double = [2 1 1];
    l (1,3) double {mustBePositive} = [5 6 3];
    
    gamma_pos (:,3) double = [0.8, 0.8, 0.8];
    gamma_neg (:,3) double = [0.8, 0.8, 0.8];


    options.Fs (1,1) double {mustBePositive} = 48000;
    options.kernel_sample_width  (1,1) double {mustBePositive, mustBeInteger}  = 10;
    options.jitter_coord_bnd  (1,2) double {mustBeNonnegative} = [0, 0];
    options.jitter_srand (1,1) double {mustBeInteger} = 6452;

    options.mode_enc (1,:) char {mustBeMember(options.mode_enc, {'proj', 'proj_msq', 'proj_msq_mex'} )} = 'proj';
    options.enable_mex (1,1) logical = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(all([(abs(s) < l/2),  (abs(r) < l/2), (l > 0)]), 'Invalid source, receiver coordinates, or room sizes');

assert(size(gamma_pos, 1) == size(gamma_neg, 1), 'gamma_pos, gamma_neg size mismatch');
L = size(gamma_pos, 1);

N = 3; %Number of dimensions (3)

c = 343; %Speed of sound
k = T * c;
q = k^2;
Ts = 1 / options.Fs;
c_Ts = c * Ts;

M = ceil(T * options.Fs);  %Number of taps

%h = zeros(M, 1);

N_src = (max_src_odr + 1)^2;
N_rec = (max_rec_odr + 1)^2;
C = complex(zeros([N_src, N_rec, M]));

%Compute bounds
ub = zeros(1, N);
lb = zeros(1, N);
x_list = cell(1, N);
for n = 1:N
    while u_vector(ub(n), s(n), r(n), l(n))^2 < q
        ub(n) = ub(n) + 1;
    end
    while u_vector(lb(n), s(n), r(n), l(n))^2 < q
        lb(n) = lb(n) - 1;
    end
    if coder.target('MATLAB')
        x_list{n} = lb(n):ub(n);
    else
        x_list{N-n+1} = lb(n):ub(n);
    end
end

%Compute lattice coordinates
if coder.target('MATLAB')
    V_set = table2array(combinations(x_list{:})); %[S x 3]
else %Equivalent order via ndgrid for codegen assuming reverse order x_list
    V_set_cell = cell(1, N);
    [V_set_cell{:}] = ndgrid(x_list{:});
    V_set = zeros(numel(V_set_cell{1}), N);
    for n = 1:N
        V_set(:, n) = V_set_cell{n}(:);
    end
    V_set = fliplr(V_set);
end

S = size(V_set, 1);

IS_set = u_vector(V_set, s, r, l);   %[S x 3] vectors from receiver to image-source
IR_set = u_vector_flip(V_set, s, r, l);   %[S x 3] vectors from source to image-receiver


if ~isequal(options.jitter_coord_bnd, [0, 0]) %Jitter the image-source coordinates
    rng(options.jitter_srand);

    rand_off_IS_LUT = unifrnd(min(options.jitter_coord_bnd), max(options.jitter_coord_bnd), size(IS_set));
    rand_off_IR_LUT = unifrnd(min(options.jitter_coord_bnd), max(options.jitter_coord_bnd), size(IR_set));

    rand_off_IS_LUT(ismember(V_set, zeros(1, N), 'rows'), :) = zeros(1, N); %Exclude original source
    rand_off_IR_LUT(ismember(V_set, zeros(1, N), 'rows'), :) = zeros(1, N); %Exclude original source

    IS_set = IS_set + rand_off_IS_LUT;   %Add random offset to image-source coordinates
    IR_set = IR_set + rand_off_IR_LUT;   %Add random offset to image-receiver coordinates
end
R_sq = sum(IS_set.^2, 2);            %[S x 1]
D_dist = sqrt(R_sq);

% %Sanity test
% R_sq2 = sum(IR_set.^2, 2);
% norm(R_sq - R_sq2)

mask_valid = (R_sq <= q); %[S x 1]

%Pre-compute reflection IR multiplicity,
%gamma_pos_mult{m, n} = self_conv(gamma_pos(:, n), m-1)
%gamma_neg_mult{m, n} = self_conv(gamma_neg(:, n), m-1)

%With cell array
% gamma_pos_mult = cell(max(max(abs(floor((V_set(mask_valid, :) + 1)/2) ))) + 1, N); %[* x N]
% gamma_neg_mult = cell(max(max(abs(floor((1 - V_set(mask_valid, n))/2) ))) + 1, N); %[* x N]
% for n = 1:N
%     gamma_pos_mult{1, n} = 1;
%     for m = 2:size(gamma_pos_mult, 1)
%         gamma_pos_mult{m, n} = conv(gamma_pos_mult{m-1, n}, gamma_pos(:, n)); 
%     end
%     gamma_neg_mult{1, n} = 1;
%     for m = 2:size(gamma_neg_mult, 1)
%         gamma_neg_mult{m, n} = conv(gamma_neg_mult{m-1, n}, gamma_neg(:, n));        
%     end
% end

%With array of structs for codegen
gamma_mult_struct.gamma = 0;
coder.varsize("gamma_mult_struct.gamma");
gamma_pos_mult = repmat(gamma_mult_struct, max(max(abs(floor((V_set(mask_valid, :) + 1)/2) ))) + 1, N);
gamma_neg_mult = repmat(gamma_mult_struct, max(max(abs(floor((1 - V_set(mask_valid, n))/2) ))) + 1, N);
for n = 1:N
    gamma_pos_mult(1, n).gamma = 1;
    for m = 2:size(gamma_pos_mult, 1)
        gamma_pos_mult(m, n).gamma = conv(gamma_pos_mult(m-1, n).gamma(:), gamma_pos(:, n)); 
    end
    gamma_neg_mult(1, n).gamma = 1;
    for m = 2:size(gamma_neg_mult, 1)
        gamma_neg_mult(m, n).gamma = conv(gamma_neg_mult(m-1, n).gamma(:), gamma_neg(:, n));        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute reflection IR per valid image-source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single-tap reflection filters
R_set = double(mask_valid);  %[S x 1]
for n = 1:N
    R_set(mask_valid) = R_set(mask_valid) .* ... 
        gamma_pos(n).^abs(floor((V_set(mask_valid, n) + 1)/2) )  .* ...
        gamma_neg(n).^abs(floor((1 - V_set(mask_valid, n))/2) );
end

%Multi-tap reflection filters
R_set_pos_idx = ones([S, N]); %[S x N]
R_set_neg_idx = ones([S, N]); %[S x N]
for n = 1:N
    R_set_pos_idx(mask_valid, n) = abs(floor((V_set(mask_valid, n) + 1)/2) ) + 1;
    R_set_neg_idx(mask_valid, n) = abs(floor((1 - V_set(mask_valid, n))/2) ) + 1;
end

%Valid image-sources 
idx_is = find(mask_valid); %Image-source indices [K x 1]

%Sample fractional indices of valid image-sources
idx_h_frac = D_dist(mask_valid) / c * options.Fs; %[K x 1]

%Compute spherical coordinates of incident angles of image-source - receiver
[az, elev] = cart2sph(IS_set(mask_valid, 1), IS_set(mask_valid, 2), IS_set(mask_valid, 3));
theta_IS = pi/2 - elev;
phi_IS = az;

%Compute spherical coordinates of incident angles of image-receiver - source
[az, elev] = cart2sph(IR_set(mask_valid, 1), IR_set(mask_valid, 2), IR_set(mask_valid, 3));
theta_IR = pi/2 - elev;
phi_IR = az;


%Sum of image-source, image-receiver, wall reflection contributions
a = options.kernel_sample_width; %kernel width

for i = 1:numel(idx_h_frac)
    
    idx = (ceil(idx_h_frac(i) - a) : floor(idx_h_frac(i) + a))'; %Write sample index
    idx(idx <= 0) = []; idx(idx >= M) = [];         %Trim at boundaries
    N_i = numel(idx);

    if N_i > 0 %Non-empty write indices
    
        if strcmp(options.mode_enc, 'proj')

            C_src_i = sh_enc_proj(max_src_odr, theta_IR(i), phi_IR(i), is_real);
            C_rec_i = sh_enc_proj(max_rec_odr, theta_IS(i), phi_IS(i), is_real);    

        elseif strcmp(options.mode_enc, 'proj_msq')

            C_src_i = sh_enc_proj_msq(max_src_odr, theta_IR(i), phi_IR(i), is_real);
            C_rec_i = sh_enc_proj_msq(max_rec_odr, theta_IS(i), phi_IS(i), is_real);

        elseif strcmp(options.mode_enc, 'proj_msq_mex')

            C_src_i = sh_enc_proj_msq_mex(max_src_odr, theta_IR(i), phi_IR(i), is_real, false);
            C_rec_i = sh_enc_proj_msq_mex(max_rec_odr, theta_IS(i), phi_IS(i), is_real, false);    

        else
            error('Unsuported mode');
        end
    
        if L == 1 %Single-tap reflection filters

            k_idx = c_Ts * idx; %Sample distance        
    
            %Iterate over sinc-kernel samples
            C_src_rec_i = R_set(idx_is(i)) * C_src_i(:) * C_rec_i(:).';
            for m = 1:N_i
                C(:, :, idx(m)) = C(:, :, idx(m)) + C_src_rec_i * lanczos_kernel(idx(m) - idx_h_frac(i), a) ./ (k_idx(m).^((N - 1)/2));
            end
    
        else      %Multi-tap reflection filters

            %Compute sinc-kernel delayed reflection filter
            r_i = lanczos_kernel(idx - idx_h_frac(i), a);
            for n = 1:N
                r_pos_n = R_set_pos_idx(idx_is(i), n);
                r_neg_n = R_set_neg_idx(idx_is(i), n);  

                %With cell
                % r_i = conv(r_i, gamma_pos_mult{r_pos_n, n});
                % r_i = conv(r_i, gamma_neg_mult{r_neg_n, n});

                %With array of structs for codegen
                r_i = conv(r_i, gamma_pos_mult(r_pos_n, n).gamma(:));
                r_i = conv(r_i, gamma_neg_mult(r_neg_n, n).gamma(:));
            end
            
            %Update idx
            idx = idx(1):(idx(1) + numel(r_i) - 1);
            idx(idx >= M) = [];
            N_i = numel(idx);
            k_idx = c_Ts * idx; %Sample distance 

            %Iterate over sinc-kernel samples * reflection filter
            C_src_rec_i = C_src_i(:) * C_rec_i(:).';
            for m = 1:N_i
                C(:, :, idx(m)) = C(:, :, idx(m)) + C_src_rec_i * r_i(m) ./ (k_idx(m).^((N - 1)/2));
            end
    
        end
    end

end

t = (0:(M-1)) / options.Fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = u_vector(v, s, r, l)
%Compute vector from receiver to image-source at lattice coordinate

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%v:         [N x M] or [M x N] Lattice coordinate for N dimensions
%s:         [N x 1] or [1 x N] Source coordinate
%r:         [N x 1] or [1 x N] Receiver coordinate
%l:         [N x 1] or [1 x N] Orthotope dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%u:         [N x M] or [1 x N] Vector from receiver to image-source

%u = l .* v + ((-1).^v) .* s - r;
u = bsxfun(@times, v, l) + bsxfun(@minus, bsxfun(@times, ((-1).^v), s), r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = u_vector_flip(v, s, r, l)
%Compute vector from source to image-receiver at lattice coordinate

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%v:         [N x M] or [M x N] Lattice coordinate for N dimensions
%s:         [N x 1] or [1 x N] Source coordinate
%r:         [N x 1] or [1 x N] Receiver coordinate
%l:         [N x 1] or [1 x N] Orthotope dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%u:         [N x M] or [1 x N] Vector from receiver to image-source

%u = l .* v + ((-1).^v) .* s - r;
v(rem(v, 2) == 0) = -v(rem(v, 2) == 0);

u = bsxfun(@times, v, l) + bsxfun(@minus, bsxfun(@times, ((-1).^v), r), s);
