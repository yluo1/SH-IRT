function tst_sh_cdf_theta(mode)
%Test sh_cdf_theta.m

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:
%tst_sh_cdf_theta('monotonicity');

if strcmp(mode, 'monotonicity')
    rng(134);
    theta = linspace(0, pi, 100);
    max_odr = 3;
    is_real = false;

%    C_pdf = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real), is_real), 'Sum'); %Is PDF
    C_pdf = sh_nrm(sh_msq(sh_enc_proj(max_odr, pi, 0, is_real), is_real), 'Sum'); %Is PDF

%    C_pdf = sh_nrm(sh_rand(max_odr, 1, is_real), 'Sum'); %Not PDF
    %C_pdf

    D_pdf = sh_nrm(sh_msq(sh_enc_proj(max_odr, pi/2, 0, is_real), is_real), 'Sum'); %Is PDF

    sh_plt(C_pdf, 'mercator', is_real, 'title_name', 'Function f_C');
    sh_plt(D_pdf, 'mercator', is_real, 'title_name', 'Function f_D')

    CDF_C = sh_cdf_theta(C_pdf, theta);
    CDF_D = sh_cdf_theta(D_pdf, theta);

    %Plotting
    fontsize = 16;
    figure; 
    plot(theta, CDF_C, theta, CDF_D, 'linewidth', 1.5); grid on; axis tight; 
    h_lg = legend('$F_C(\theta)$', '$F_D(\theta)$', 'location', 'best', 'interpreter', 'latex');
    set(h_lg, 'fontsize', fontsize - 1);
    xlabel('\theta', 'fontsize', fontsize);
    ylabel('CDF(\theta)', 'fontsize', fontsize);
    set(gca, 'fontsize', fontsize - 1);
    title('Cumulative Distribution Function of Co-latitude');

end