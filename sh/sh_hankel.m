function y = sh_hankel(nu, z, kind)
%Evaluate spherical Hankel function of the first or second kind

%Author:    Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%nu:        [1 x N] Function orders
%z:         [1 x M] Wave numbers
%kind:      String, first or second kind {'first', 'second'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%y:         [N x M] Function evaluations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate sample evaluations

%nu = -5:5;
%z  = -10:10;
%y = sh_hankel(nu, z);

arguments
    nu   (1,:) double = 0;
    z    (1,:) double = 1;
    kind (1,:) char {mustBeMember(kind, {'first', 'second'})} = 'first';
end

N = numel(nu);
M = numel(z);

besselj_list = complex(zeros(N, M));
bessely_list = complex(zeros(N, M));

for m = 1:M
    besselj_list(:, m) = besselj(nu + 0.5, z(m));
    bessely_list(:, m) = bessely(nu + 0.5, z(m));
end

if strcmp(kind, 'first')
    y = bsxfun(@times, besselj_list + 1i * bessely_list, sqrt(pi ./ (2 * z(1:M)) ));
elseif strcmp(kind, 'second')
    y = bsxfun(@times, besselj_list - 1i * bessely_list, sqrt(pi ./ (2 * z(1:M)) ));
else
    error('Unsupported kind');
end

