function tst_rand_ortho_sh_pdf

rng(31641);
max_odr = 3;
N = (max_odr + 1)^2;

Q = zeros(N, N);
for l = 0:max_odr
    N_l = 2 * l + 1;
    [Q_l, R_l] = qr(randn(N_l));
    Q_l = Q_l * diag(sign(diag(R_l)));
    
    idx = l^2 + (1:(2*l+1));

%    Q(idx, idx) = Q_l;
    %Q(idx, idx) = (Q_l + Q_l') / 2; %Not orthogonal
 
    %Cpx->Re
    % m = (1:l)';
    % pos = [-1i * ones(l,1); sqrt(2)/2;  (-1).^m] / sqrt(2);    
    % neg = [1i * (-1).^flipud(m); sqrt(2)/2; ones(l,1)] / sqrt(2);
    % Q(idx, idx) = bsxfun(@times, Q_l, pos) + bsxfun(@times, flipud(Q_l), neg);

    %Re->cpx
    m = (1:l)';
    pos = [1i * ones(l,1); sqrt(2)/2;  (-1).^m] / sqrt(2);    
    neg = [ones(l,1); sqrt(2)/2; -1i * (-1).^m] / sqrt(2);
    Q(idx, idx) = bsxfun(@times, Q_l, pos) + bsxfun(@times, flipud(Q_l), neg);
end

Q(1,1) = 1; %Preserve positive

%R_blk = sh_rot_mat(max_odr, deg2rad([10 30 20]));

C = sh_enc_rbf('SqExp', max_odr, pi/2, 0, 1);
D = Q * C;

dB_lim = [-20, 0];

sh_plt(C, 'mercator', 'dB_lim', dB_lim);
sh_plt(D, 'mercator', 'dB_lim', dB_lim);