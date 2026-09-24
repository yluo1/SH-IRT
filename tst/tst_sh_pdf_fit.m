function tst_sh_pdf_fit(mode)
%Test sh_pdf_fit.m

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Test pdf fit to various targets

%tst_sh_pdf_fit('two_pt')
%tst_sh_pdf_fit('rand')

if strcmp(mode, 'two_pt')

    max_odr = 8;
    N = (max_odr + 1)^2;
    is_real = false;
    rng(12463);
    dB_lim = [-60, 0];
    disp_phase = false;
    
    [theta, phi] = sh_fib(N);
    X = zeros(N, 1);  X(25) = 1; X(45) = 1;
   

elseif strcmp(mode, 'rand')

    max_odr = 8;
    N = (max_odr + 1)^2;
    is_real = false;
    rng(12463);
    dB_lim = [-60, 0];
    disp_phase = false;
    
    [theta, phi] = sh_fib(N);
    %X = zeros(N, 1);  X(25) = 1; X(45) = 1;
    C_pdf_ref = sh_nrm(sh_msq(sh_rand(floor(max_odr/2), 1, is_real), is_real), 'Sum');
    X = real(sh_dec(C_pdf_ref, theta, phi, is_real));  

    sh_plt(C_pdf_ref, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'PDF Ref.');
    
else
    error('Unsupported mode');
end

%Solve
[C_pdf_SqProjNNLS, err_SqProjNNLS]  = sh_pdf_fit(X, theta, phi, max_odr, is_real, 'SqProjNNLS');
[C_pdf_SqProjQP, err_SqProjQP]      = sh_pdf_fit(X, theta, phi, max_odr, is_real, 'SqProjQP');
[C_pdf_SqMagMS, err_SqMagMS]        = sh_pdf_fit(X, theta, phi, max_odr, is_real, 'SqMagMS', 'SqProjQP_C0', sh_rand(floor(max_odr/2), 50, true));
[C_pdf_SqMagSOMS, err_SqProjSOMS]   = sh_pdf_fit(X, theta, phi, max_odr, is_real, 'SqMagSOMS');

%Plot
sh_plt(C_pdf_SqProjNNLS, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'SqProjNNLS');
err_SqProjNNLS

sh_plt(C_pdf_SqProjQP, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'SqProjQP');
err_SqProjQP

sh_plt(C_pdf_SqMagMS, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'SqProjQP');
err_SqMagMS

sh_plt(C_pdf_SqMagSOMS, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'SqProjSOMS');
err_SqProjSOMS

