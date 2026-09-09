function tst_sh_fib
%Evaluate condition number of SH evaluations over uniform points on sphere

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_list = [64, 256, 1024];
%N_list = ((1:6) + 1).^2;

for N = N_list
     [theta, phi] = sh_fib(N);
    
     max_odr = floor(sqrt(N) - 1);

     Y_cpx  = sh_val(max_odr, theta, phi, false); cond(Y_cpx)
     Y_real = sh_val(max_odr, theta, phi, true); cond(Y_real)

     sc_plt(theta, phi, ones(size(theta)));
end
