function allo_trans=trans_allo_out(indexx, solution, params, list)

% function to transform solution of fsolve or vpasolve with transformed
% model.

% input
% indexx: structure of logical vector ==1 if respective transformation has
%           to be performed
% solution: raw, untransformed solution from fsolve or vpasolve
% params: structure or double of parameter values



%-- transformation 

allo_trans=solution;
allo_trans(indexx.sqr)=solution(indexx.sqr).^2;
allo_trans(indexx.exp)=exp(solution(indexx.exp));
allo_trans(indexx.lab)=params(list.params=='upbarH')/(1+exp(solution(indexx.lab)));

end
