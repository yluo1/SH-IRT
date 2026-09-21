function Y = legendre_poly(n, x, mode)
% Evaluate codegen supported Legendre polynomial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%n:     Polynomial order
%x:     [M x 1] Polynomial variable
%mode:  Compute method {'asc', 'bonnet'}
%           'asc':      Subset of associated Legendre polynomials (See asc_legendre.m)
%           'bonnet':   Bonnet recurrence relation

%Output
%Y:     [M x 1] Polynomial evaluation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('legendre_poly', '-o', 'sh/legendre_poly_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Plot lower-order Legendre polynomials

% M = 100;  % Number of sample points
% N = 5;    % Max order
% %mode = 'asc';
% mode = 'bonnet';
% x = linspace(-1, 1, M)';
% Y = zeros(M, N + 1);
% legend_str = {};
% for n = 0:N
%   Y(:, n + 1) = legendre_poly(n, x, mode);
%   legend_str{n + 1} = ['n = ', num2str(n)];
% end
% fontsize = 16;
% figure; plot(x, Y, 'linewidth', 1.5);
% xlabel('x', 'fontsize', fontsize); ylabel('P_n(x)', 'fontsize', fontsize);
% title('Legendre Polynomial', 'fontsize', fontsize + 1);
% h_lg = legend(legend_str, 'location', 'best', 'NumColumns', 2); set(h_lg, 'fontsize', fontsize - 1);
% set(gca, 'fontsize', fontsize - 1); grid on; axis tight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Compare high-order Legendre polynomials

% M = 100;  % Number of sample points
% N = 50;   % Max-order
% x = linspace(-1, 1, M)';
% Y_asc = zeros(M, N + 1);
% Y_bon = zeros(M, N + 1);
% for n = 0:N
%   Y_asc(:, n + 1) = legendre_poly(n, x, 'asc');
%   Y_bon(:, n + 1) = legendre_poly(n, x, 'bonnet');
% end

% err_asc_bon = norm(Y_asc - Y_bon)

arguments
    n (1,1) double = 0;
    x (:,1) double = 0;
    mode (1,:) char {mustBeMember(mode, {'asc', 'bonnet'})} = 'asc';
end

if strcmp(mode, 'asc')
    P = asc_legendre(n, x);
    Y = P(1,:).';

elseif strcmp(mode, 'bonnet')  % Recurrence relation
    % P_0(x) = 1
    % P_1(x) = x
    % P_n(x) = ((2 * n - 1) * x * P_{n - 1} - (n - 1) * P_{l - 2}) / n
    
    if n == 0
        Y = ones(size(x));
    elseif n == 1
        Y = x;
    else
        Y_prev_prev = ones(size(x));
        Y_prev = x;
        Y = zeros(size(x));
        for l = 2:n          
            Y = ((2 * l -1) * x .* Y_prev - (l - 1) * Y_prev_prev) / l;
            Y_prev_prev = Y_prev;
            Y_prev = Y;
        end
    end

else
    error('Unsupported mode');
end