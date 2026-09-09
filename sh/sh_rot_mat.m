function [R_list_out, Rxyz] = sh_rot_mat(max_odr, theta_zyx, is_real, options)
%Compute spherical harmonic rotation matrices from 
%https://github.com/polarch/Spherical-Harmonic-Transform

%Archontis Politis, Microphone array processing for parametric spatial audio techniques, 2016
%Doctoral Dissertation, Department of Signal Processing and Acoustics, Aalto University, Finland

% which based its implementation on 

% Ivanic, J., Ruedenberg, K. (1996). Rotation Matrices for Real 
% Spherical Harmonics. Direct Determination by Recursion. The Journal 
% of Physical Chemistry, 100(15), 6342?6347.

% Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real 
% Spherical Harmonics. Direct Determination by Recursion Page: Additions 
% and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100.

% Modifications by Yuancheng Luo, 2026:
% * Combine into single script
% * Add block-matrix support via array of structs
% * Add intrinsic/extrinsic rotation options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%theta_zyx:     [1 x 3]  Radians, counter clock-wise rotation about [+x, +y, +z] axis
%max_odr:       Max order spherical harmonic
%is_real:       Logical, if true, evaluate real SH

%options:       struct
%options.convention:            String, rotation ordering  {'zyx', 'zxy', 'yxz', 'yzx', 'xyz', 'xzy'}  
%options.rot_intrinsic:         Logical, if false, reverse angle and convention for extrinsic order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%R_list_out:            [1 x max_odr + 1] cell 
%Rxyz:                  [3 x 3] Euler rotation matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute 3rd order SH rotation matrix

%max_odr = 3;
%theta_zyx = deg2rad([0, 0, 0]);
%is_real = true;
%R_list_out = sh_rot_mat(max_odr, theta_zyx, is_real);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    arguments
        max_odr (1,1) double {mustBeNonnegative} = 1;
        theta_zyx (1,3) double = [0 0 0];
        is_real (1,1) logical = false;  
    
        options.convention (1,3) char {mustBeMember(options.convention, {'zyx', 'zxy', 'yxz', 'yzx', 'xyz', 'xzy'})} = 'zyx';
        options.rot_intrinsic (1,1) logical = true;
    end
    
    R_struct.R = 0;
    coder.varsize("R_struct.R");
    R_list = repmat(R_struct, 1, max_odr + 1);
    for l = 1:(max_odr + 1)
        R_list(l).R = zeros([2 * (l-1) + 1, 2 * (l-1) + 1]);
    end                       
    
    if options.rot_intrinsic
        Rxyz = euler2rotationMatrix(theta_zyx(3), theta_zyx(2), theta_zyx(1), reverse(options.convention));
    else %Extrinsic
        Rxyz = euler2rotationMatrix(theta_zyx(1), theta_zyx(2), theta_zyx(3), options.convention);
    end
    
    % initialize zeroth and first band rotation matrices for recursion
    %	Rxyz = [Rxx Rxy Rxz
    %           Ryx Ryy Ryz
    %           Rzx Rzy Rzz]
    %
    % zeroth-band (l=0) is invariant to rotation
    R_list(1).R(1,1) = 1;
    
    if max_odr > 0
        % the first band (l=1) is directly related to the rotation matrix
        R_list(2).R(-1+2,-1+2) = Rxyz(2,2);       
        R_list(2).R(-1+2, 0+2) = Rxyz(2,3);
        R_list(2).R(-1+2, 1+2) = Rxyz(2,1);
        R_list(2).R( 0+2,-1+2) = Rxyz(3,2);
        R_list(2).R( 0+2, 0+2) = Rxyz(3,3);
        R_list(2).R( 0+2, 1+2) = Rxyz(3,1);
        R_list(2).R( 1+2,-1+2) = Rxyz(1,2);
        R_list(2).R( 1+2, 0+2) = Rxyz(1,3);
        R_list(2).R( 1+2, 1+2) = Rxyz(1,1);

        % compute rotation matrix of each subsequent band recursively
        for l = 2:max_odr                 
            for m=-l:l
                for n=-l:l
                    % compute u,v,w terms of Eq.8.1 (Table I)
                    d = (m==0); % the delta function d_m0
                    if abs(n)==l
                        denom = (2*l)*(2*l-1);
                    else
                        denom = (l*l-n*n);
                    end
                    u = sqrt((l*l-m*m)/denom);
                    v = sqrt((1+d)*(l+abs(m)-1)*(l+abs(m))/denom)*(1-2*d)*0.5;
                    w = sqrt((l-abs(m)-1)*(l-abs(m))/denom)*(1-d)*(-0.5);
        
                    % computes Eq.8.1
                    if u~=0, u = u * U(l, m, n); end
                    if v~=0, v = v * V(l, m, n); end
                    if w~=0, w = w * W(l, m, n); end
                    
                    R_list(l+1).R(m+l+1, n+l+1) = u + v + w;
                end
            end
        end
    end
    
    %Output rotation matrix
    R_struct_out.R = complex(0);
    coder.varsize("R_struct_out.R");
    R_list_out = repmat(R_struct_out, 1, max_odr + 1);
    for l = 1:(max_odr + 1)
        R_list_out(l).R = complex(R_list(l).R);
    end

    if ~is_real %Compute real to complex matrix conversion
    
        for l = 1:max_odr
            m = (1:l)';
            % form the diagonal
            diagT = [1i * ones(l,1); sqrt(2)/2; (-1).^m]/sqrt(2);
            % form the antidiagonal
            adiagT = [-1i*(-1).^m(end:-1:1); sqrt(2)/2; ones(l,1)]/sqrt(2);
            % form the transformation matrix for the specific band l
            T = diag(diagT) + fliplr(diag(adiagT));
    
            R_list_out(l+1).R = T.' * R_list_out(l+1).R * conj(T);
        end
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ret = U(l,m,n)
        % functions to compute terms U, V, W of Eq.8.1 (Table II)
        ret = P(0,l,m,n);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ret = V(l,m,n)
    
        if (m==0)
            p0 = P(1,l,1,n);
            p1 = P(-1,l,-1,n);
            ret = p0+p1;
        else
            if (m>0)
                dd = (m==1);
                p0 = P(1,l,m-1,n);
                p1 = P(-1,l,-m+1,n);        
                ret = p0*sqrt(1+dd) - p1*(1-dd);
            else
                dd = (m==-1);
                p0 = P(1,l,m+1,n);
                p1 = P(-1,l,-m-1,n);        
                ret = p0*(1-dd) + p1*sqrt(1+dd);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ret = W(l,m,n)
    
        if (m==0)
            error('should not be called')
        else
            if (m>0)
                p0 = P(1,l,m+1,n);
                p1 = P(-1,l,-m-1,n);        
                ret = p0 + p1;
            else
                p0 = P(1,l,m-1,n);
                p1 = P(-1,l,-m+1,n);        
                ret = p0 - p1;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function to compute term P of U,V,W (Table II)
    function ret = P(i,l,a,b)
        
        ri1 = R_list(2).R(i+2,1+2);
        rim1 = R_list(2).R(i+2,-1+2);
        ri0 = R_list(2).R(i+2,0+2);
        
        if (b==-l)
            ret = ri1 * R_list(l).R(a+l, 1) + rim1 * R_list(l).R(a+l, 2*l-1);
        else
            if (b==l)
                ret = ri1 * R_list(l).R(a+l, 2*l-1) - rim1 * R_list(l).R(a+l, 1);        
            else
                ret = ri0 * R_list(l).R(a+l, b+l);
            end
        end

    end


end    
               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = euler2rotationMatrix(alpha, beta, gamma, convention)

    %EULER2ROTATIONMATRIX Construct the rotation matrix from Euler angles (clockwise)
    %
    %   alpha:  first angle of rotation
    %   beta:   second angle of rotation
    %   gamma:  third angle of rotation
    %
    %   connvention: definition of the axes of rotation, e.g. for the y-convention
    %                this should be 'zyz', for the x-convnention 'zxz', and for
    %                the yaw-pitch-roll convention 'zyx', (clockwise)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Archontis Politis, 7/2/2015
    %   archontis.politis@aalto.fi
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Rx = @(theta) [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
    Ry = @(theta) [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
    Rz = @(theta) [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
    
    switch convention(1)
        case 'x'
            R1 = Rx(alpha);
        case 'y'
            R1 = Ry(alpha);
        case 'z'
            R1 = Rz(alpha);
    end
    
    switch convention(2)
        case 'x'
            R2 = Rx(beta);
        case 'y'
            R2 = Ry(beta);
        case 'z'
            R2 = Rz(beta);
    end
    
    switch convention(3)
        case 'x'
            R3 = Rx(gamma);
        case 'y'
            R3 = Ry(gamma);
        case 'z'
            R3 = Rz(gamma);
    end
    
    R = R3*R2*R1;
end