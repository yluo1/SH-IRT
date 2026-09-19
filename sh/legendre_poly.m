function Y = legendre_poly(n, x, mode)
% Evaluate codegen supported Legendre polynomial from asociated Legendre polynomial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%n:     Polynomial order
%x:     [M x 1] Polynomial variable
%mode:  Compute method {'asc'}
%           'asc':      Subset of associated Legendre polynomials (See asc_legendre.m)

%Output
%Y:     [M x 1] Polynomial evaluation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Plot lower-order Legendre polynomials

% M = 100;  % Number of sample points
% N = 5;    % Max order
% x = linspace(-1, 1, M)';
% Y = zeros(M, N + 1);
% legend_str = {};
% for n = 0:N
%   Y(:, n + 1) = legendre_poly(n, x);
%   legend_str{n + 1} = ['n = ', num2str(n)];
% end
% fontsize = 16;
% figure; plot(x, Y, 'linewidth', 1.5);
% xlabel('x', 'fontsize', fontsize); ylabel('P_n(x)', 'fontsize', fontsize);
% title('Legendre Polynomial', 'fontsize', fontsize + 1);
% h_lg = legend(legend_str, 'location', 'best', 'NumColumns', 2); set(h_lg, 'fontsize', fontsize - 1);
% set(gca, 'fontsize', fontsize - 1); grid on; axis tight;

arguments
    n (1,1) double = 0;
    x (:,1) double = 0;
    mode (1,:) char {mustBeMember(mode, 'asc')} = 'asc';
end

if strcmp(mode, 'asc')
    P = asc_legendre(n, x);
    Y = P(1,:).';
else
    error('Unsupported mode');
end