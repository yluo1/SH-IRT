function [theta, phi, v] = sh_plat(name)
%Get spherical coordinates of platonic solid vertices

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%name:          String, solid name {'tetrahedron', 'octohedron', 'cube', 'icosohedron', 'dodecahedron'}
%               Number of vertices        4,            6,         8,         12,             20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:         [N x 1]  Co-latitude [0, pi]
%phi:           [N x 1]  Azimuth [0, 2 * pi)

%v:             [N x 3]  Cartesian coordinates, unit length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate platonic vertex coordinates,
%project onto first order SH bases and evaluate condition number

% max_odr = 1;
% name_list = {'tetrahedron', 'octohedron', 'cube', 'icosohedron', 'dodecahedron'};
% for n = 1:numel(name_list)
%     [theta, phi, v] = sh_plat(name_list{n});
%     Y_cpx = sh_val(max_odr, theta, phi, false); cond(Y_cpx)
%     Y_real = sh_val(max_odr, theta, phi, true); cond(Y_real)
% 
%     figure; scatter3(v(:, 1), v(:, 2), v(:, 3), 60, 'filled');
%     axis square;
% end

arguments
    name (1,:) char {mustBeMember(name, {'tetrahedron', 'octohedron', 'cube', 'icosohedron', 'dodecahedron'})} = 'tetrahedron';
end

if strcmp(name, 'tetrahedron')
    v = tetrahedron_vertices;
elseif strcmp(name, 'octohedron')
    v = octohedron_vertices;
elseif strcmp(name, 'cube')
    v = cube_vertices;
elseif strcmp(name, 'icosohedron')
    v = icosohedron_vertices;
elseif strcmp(name, 'dodecahedron')
    v = dodecahedron_vertices;
else
    error('Unsupported N');
end

v = bsxfun(@rdivide, v, sqrt(sum(v.^2, 2))); %Normalize
[phi, theta] = cart2sph(v(:, 1), v(:, 2), v(:, 3));
theta = pi/2 - theta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = tetrahedron_vertices
%Output
%v:     [4 x 3]

sqrt_8_9 = sqrt(8.0 / 9.0);
sqrt_2_9 = sqrt(2.0 / 9.0);
sqrt_2_3 = sqrt(2.0 / 3.0);

v = zeros(4, 3);

v(1, :) = [sqrt_8_9, 0, -1.0/3.0];
v(2, :) = [-sqrt_2_9, sqrt_2_3, -1.0/3.0];
v(3, :) = [-sqrt_2_9, -sqrt_2_3, -1.0/3.0];
v(4, :) = [0, 0, 1.0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = octohedron_vertices
%Output
%v:     [6 x 3]

v = [   1, 0, 0; ...
       -1, 0, 0; ...
        0, 1, 0; ...
        0,-1, 0; ...
        0, 0, 1; ... 
        0, 0,-1; ... 
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = cube_vertices
%Output
%v:     [8 x 3]

a = 1/sqrt(3);

v(1, :) = [a, a, a];
v(2, :) = [a, a, -a];
v(3, :) = [a, -a, a];
v(4, :) = [a, -a, -a];
v(5, :) = [-a, a, a];
v(6, :) = [-a, a, -a];
v(7, :) = [-a, -a, a];
v(8, :) = [-a, -a, -a];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = icosohedron_vertices
%Output
%v:     [12 x 3]

a = 2 / sqrt(10 + 2*sqrt(5));
b = (1+sqrt(5)) / sqrt(10 + 2 * sqrt(5));


v = zeros(12, 3);

v(1,:) = [0 a b];
v(2,:) = [0 a -b];
v(3,:) = [0 -a b];
v(4,:) = [0 -a -b];
v(5,:) = [a b 0];
v(6,:) = [a -b 0];
v(7,:) = [-a b 0];
v(8,:) = [-a -b 0];
v(9,:) = [b 0 a];
v(10,:) = [-b 0 a];
v(11,:) = [b 0 -a];
v(12,:) = [-b 0 -a];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = dodecahedron_vertices
%Output
%v:     [20 x 3]

a = 1 / sqrt(3);
b = (1+sqrt(5))/2 / (sqrt( ((1+sqrt(5))/2)^2 +  (2 / (1+sqrt(5)))^2 ));
c = 2/(1+sqrt(5)) / (sqrt( ((1+sqrt(5))/2)^2 +  (2 / (1+sqrt(5)))^2 ));


v = zeros(20, 3);

v(1,:) = [a, a, a];
v(2,:) = [a, a, -a];
v(3,:) = [a, -a, a];
v(4,:) = [a, -a, -a];
v(5,:) = [-a, a, a];
v(6,:) = [-a, a, -a];
v(7,:) = [-a, -a, a];
v(8,:) = [-a, -a, -a];

v(9,:) = [0, b, c];
v(10,:) = [0, b, -c];
v(11,:) = [0, -b, c];
v(12,:) = [0, -b, -c];

v(13,:) = [c, 0, b];
v(14,:) = [c, 0, -b];
v(15,:) = [-c, 0, b];
v(16,:) = [-c, 0, -b];

v(17,:) = [b, c, 0];
v(18,:) = [b, -c, 0];
v(19,:) = [-b, c, 0];
v(20,:) = [-b, -c, 0];