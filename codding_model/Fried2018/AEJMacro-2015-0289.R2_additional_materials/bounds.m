%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: bounds.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Ensures parameter values stay within reasonable bounds for
% minDist (Nelder-Mead) routine. 
% Inputs: 
% Vector of calibrated paramters from minDist: pp
% Outputs:
% Vector of calibrated paramters for minDist
% Note: Calibration is robust to economicially reasonable changes in the bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newP = bounds(pp)

alphag = pp(1);  eta = pp(2); deltaye = pp(3); 
epsf = pp(4); deltaff = pp(5);  nu1 = pp(6); 


if eta > 0.98
    eta = 0.98;
elseif eta < .1;
    eta =.1;
end


if alphag > 0.99
    alphag = 0.99;
elseif alphag < .1;
    alphag =.1;
end

if epsf < 1.099
    epsf =1.099;
elseif epsf> 15
    epsf=15; 
end

if nu1 < .1 
    nu1 =1;
elseif nu1 > 1; 
    nu1 =1; 
end


if nu1 < .1 
    nu1 =.1;
elseif nu1 > 1; 
    nu1 = 1; 
end

if deltaff > .8
    deltaff =.8;
elseif deltaff < .2; 
    deltaff = .2;
end

pp(1) = alphag; 
pp(2) = eta;
pp(3) =deltaye;
pp(4) = epsf;
pp(5) = deltaff;
pp(6) =nu1;
 
newP = pp;