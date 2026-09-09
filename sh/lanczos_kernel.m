function f = lanczos_kernel(x, a)
%f = sinc(x) * sinc(x/a) = a*sin(pi*x) .* sin(pi*x/a) ./ (pi^2 * x.^2)
%for |x| <= a, otherwise 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%x:     [M x N]
%a:     Scalar, non-negative

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:     [M x N]

arguments
    x (:,:) double = 0;
    a (1,1) double {mustBeNonnegative} = 10;    
end

f = a*sin(pi*x) .* sin(pi*x/a) ./ (pi^2 * x.^2);

f(abs(x) > a) = 0;
f(x == 0) = 1;
