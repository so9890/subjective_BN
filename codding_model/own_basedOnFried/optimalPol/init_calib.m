function f = init_calib(x, MOM, params,list, pol)
% initial guess for calibration


pg=exp(x(1));

[ ~, pe, ~, EY]=aux_calib_Prod(MOM, pg, params,list, pol);
f(1) = EY*pe-MOM.EpeY;

end



