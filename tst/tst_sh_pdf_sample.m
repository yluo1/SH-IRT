function tst_sh_pdf_sample
%Test sh_pdf_sample.m

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1277);

max_odr = 6;
theta = pi/3;
phi = 0;
ell = 0.25;
is_real = false;
dB_lim = [-120, 0];

% C_pdf = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, theta, phi - pi/2, ell, is_real),  is_real), 'Sum') ...
%       + sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, theta + pi/2, phi + pi/4, ell, is_real), is_real), 'Sum');

C_pdf = sh_msq(sh_rand(floor(max_odr / 2), 1, is_real, true), is_real);

C_pdf = sh_nrm(C_pdf, 'Sum');

num_samples = 500;

rseed = 2341;

rng(rseed);
tic;
[theta_acceptRejectUniform, phi_acceptRejectUniform] = sh_pdf_sample(C_pdf, num_samples, 'acceptRejectUniform', is_real);
duration_acceptRejectUniform = toc

rng(rseed);
tic;
[theta_inverseTransform, phi_inverseTransform] = sh_pdf_sample(C_pdf, num_samples, 'inverseTransform', is_real);
duration_inverseTransform = toc

sh_plt(C_pdf, 'mercator', is_real, 'title_name', 'Accept Reject', 'disp_theta_phi', [theta_acceptRejectUniform, phi_acceptRejectUniform], 'disp_phase', false, 'dB_lim', dB_lim);
sh_plt(C_pdf, 'mercator', is_real, 'title_name', 'Inverse Transform', 'disp_theta_phi', [theta_inverseTransform, phi_inverseTransform], 'disp_phase', false, 'dB_lim', dB_lim);