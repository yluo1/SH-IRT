function [ineqnonlin, eqnonlin, gradineqnonlin, gradeqnonlin] = sh_msq_pdf_nonlcon(C, relax_eq_constr)
%Equality constraints for unity integration of magnitude squared spherical harmonic
%expansion density functions

%int abs(Y(\omega) * C)^2 d\omega = 1 
%<=>
%C'*C = 1

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C: [(P + 1)^2 x 1]     Real spherical harmonic expansion coefficients
%relax_eq_constr:       Logical, if true, relax equality constraint C'*C = 1 
%                       to inequality constraints C'*C <= 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%ineqnonlin:        Scalar or empty
%eqnonlin:          Scalar or empty
%gradineqnonlin:    [(P + 1)^2 x 1] or empty
%gradeqnonlin:      [(P + 1)^2 x 1] or empty

arguments
     C (:,1) double = 0;     
     relax_eq_constr (1,1) logical = false;
end

if relax_eq_constr
    ineqnonlin = C'*C - 1;
    eqnonlin =  [];

    if nargout > 2
        gradineqnonlin = 2 * C;
        gradeqnonlin = [];
    end

else
    ineqnonlin = [];
    eqnonlin =  C'*C - 1;

    if nargout > 2
        gradineqnonlin = [];
        gradeqnonlin = 2 * C;
    end

end

