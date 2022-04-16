function f = init_calib(x, MOM, params,list, pol)
% initial guess for calibration


pg=exp(x(1));
Y=exp(x(2));

[~, ~, pe, ~, E]=aux_calib_Prod(MOM, pg, Y, C, params,list, pol);
f(1) = E*pe/Y-MOM.EpeY;
%- resource constraint
C   = Y-xn-xf-xg;
end



