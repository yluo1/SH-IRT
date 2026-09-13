function pass = tst_gp_cov_check_psd(mode)
%Empirical checks of covariance matrix positive semi-definiteness 

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:      String, test mode {'rand', 'rand_dup_dir', 'const_dir'}
%                       'rand':         Randomize direction and frequency
%                       'rand_dup_dir': Random direction and frequency, duplicate direction once at half frequency
%                       'const_dir':    One direction, vary frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:      Test different direction x frequency combinations

%tst_gp_cov_check_psd('rand')
%tst_gp_cov_check_psd('rand_dup_dir')
%tst_gp_cov_check_psd('const_dir');

arguments
    mode (1,:) char {mustBeMember(mode, {'rand', 'rand_dup_dir', 'const_dir'})} = 'rand';
end

if strcmp(mode, 'rand')
    rng(1278);
 
%    N_S = 10;
    N_S = 50;
%    N_S = 100;
%    N_S = 200;
    
    [theta, phi] = sh_fib(N_S);
    omega = 2 * pi * (1000 * rand(N_S, 1) + 50);
    
    
    %Setup input set {X, y}
    freq_obs = max(1, omega / (2 * pi) ); %[N_S x 1]
    X = [2 * pi * freq_obs(:), theta(:), phi(:)];
        
    options_cov = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', sqrt(2)/2, 'cov_ell', 343, 'cov_gamma', 1);
    
    K = gp_cov(X, X, options_cov);
    min(eig(K))
    
    %Empirical test for PSD for varying ell, gamma
    N_ell = 100;
    N_gamma = 100;
    ell_list = linspace(1, 10000, N_ell);
    gamma_list = linspace(0.1, 10, N_gamma);
    
    pass = true;
    for i = 1:N_ell
        for j = 1:N_gamma
            options_cov = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', sqrt(2)/2, 'cov_ell', ell_list(i), 'cov_gamma', gamma_list(j));
            
            K = gp_cov(X, X, options_cov);
            
            %eig(K)
            min_eig = min(eig(K));
            if min_eig < -1e-8
                pass = false;
                min_eig
                ;
            end
        end
    end

elseif strcmp(mode, 'rand_dup_dir')
    rng(1278);
 
%    N_S = 10;
    N_S = 50;
%    N_S = 100;
    %N_S = 200;
    
    [theta, phi] = sh_fib(N_S);
    omega = 2 * pi * (1000 * rand(N_S, 1) + 50);
    
    %Setup input set {X, y}
    freq_obs = max(1, omega / (2 * pi) ); %[N_S x 1]
    X = [2 * pi * freq_obs(:), theta(:), phi(:)];
    
    X = [X; [X(:,1)/2, theta(:), phi(:)]];
    
    options_cov = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', sqrt(2)/2, 'cov_ell', 343, 'cov_gamma', 1);
    
    K = gp_cov(X, X, options_cov);
    min(eig(K))
    
    %Empirical test for PSD for varying ell, gamma
    N_ell = 100;
    N_gamma = 100;
    ell_list = linspace(1, 10000, N_ell);
    gamma_list = linspace(0.1, 10, N_gamma);
    
    pass = true;
    for i = 1:N_ell
        for j = 1:N_gamma
            options_cov = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', sqrt(2)/2, 'cov_ell', ell_list(i), 'cov_gamma', gamma_list(j));
            
            K = gp_cov(X, X, options_cov);
            
            %eig(K)
            min_eig = min(eig(K));
            if min_eig < -1e-8
                pass = false;
                min_eig
                ;
            end
        end
    end


elseif strcmp(mode, 'const_dir')

    N_omega = 100;    
    theta = pi/2;
    phi = 0;

    omega = 2 * pi * linspace(50, 1000, N_omega);
    
    %Setup input set {X, y}
    freq_obs = max(1, omega / (2 * pi) ); %[N_S x 1]
    X = [2 * pi * freq_obs(:), theta * ones(N_omega, 1), phi * ones(N_omega, 1)];
    
    options_cov = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', sqrt(2)/2, 'cov_ell', 343, 'cov_gamma', 1);
    
    K = gp_cov(X, X, options_cov);
    min(eig(K))
    
    %Empirical test for PSD for varying ell, gamma
    N_ell = 100;
    N_gamma = 100;
    ell_list = linspace(1, 10000, N_ell);
    gamma_list = linspace(0.1, 10, N_gamma);
    
    pass = true;
    for i = 1:N_ell
        for j = 1:N_gamma
            options_cov = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', sqrt(2)/2, 'cov_ell', ell_list(i), 'cov_gamma', gamma_list(j));
            
            K = gp_cov(X, X, options_cov);
            
            %eig(K)
            min_eig = min(eig(K));
            if min_eig < -1e-8
                pass = false;
                min_eig
                ;
            end
        end
    end


else
    error('Unknown mode');
end