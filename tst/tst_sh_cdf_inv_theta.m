function tst_sh_cdf_inv_theta
%Evaluate sh_cdf_inv_theta.m performance

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(134);
max_odr = 6;
is_real = false; 

C_pdf = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real), is_real), 'Sum'); %Is PDF
%C_pdf = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, pi/2, pi/3, 1, is_real), is_real), 'Sum'); %Is PDF

N_fwd = 100;
theta = linspace(0, pi, N_fwd);
N_inv = 11;
CDF_C = sh_cdf_theta(C_pdf, theta);
u = linspace(0, 1, N_inv)';

%Evaluations
[theta_CDF_inv, stats_bisection] = sh_cdf_inv_theta(C_pdf, u, 'mode', 'bisection');
disp_stats(stats_bisection);

[theta_CDF_inv_NRB, stats_NRB] = sh_cdf_inv_theta(C_pdf, u, 'mode', 'NewtonRaphsonBisection');
disp_stats(stats_NRB);


%Plotting
fontsize = 18;
h_cdf = figure;
h_cdf.Position = [100, 100, [600, 350] * 4/5];
plot(theta, CDF_C, '-', theta_CDF_inv, u, '*', 'linewidth', 3, 'markersize', 10); grid on; axis tight; 
h_lg = legend('$u = F_{\Theta}(\theta)$', '$\theta = F_{\Theta}^{-1}(u \sim \mathcal{U}(0,1))$', 'location', 'best', 'interpreter', 'latex');
set(h_lg, 'fontsize', fontsize );
xlabel('Co-latitude \theta', 'fontsize', fontsize);
ylabel('CDF(\theta)', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize - 1);
title('Cumulative Distribution Function', 'fontsize', fontsize + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display stats
function disp_stats(stats)
disp(['Mean +- Std. Number of iterations:       ', num2str(mean(stats.num_iter(:))), ' +- ', num2str(std(stats.num_iter(:)))  ]);
disp(['Mean +- Std. Number of function evals:   ', num2str(mean(stats.num_func_eval(:))), ' +- ', num2str(std(stats.num_func_eval(:)))  ]);
disp(['Mean +- Std. Number of derivative evals: ', num2str(mean(stats.num_deriv_eval(:))), ' +- ', num2str(std(stats.num_deriv_eval(:)))  ]);
disp(' ');
