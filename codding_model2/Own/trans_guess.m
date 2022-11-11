function guess_trans=trans_guess(indexx, guess)

% function to transform initial values to pass to function with 
% transformed variables

%-- extract necessary params

guess_trans=guess;

% guess_trans(indexx.sqr)=sqrt(guess(indexx.sqr));
% guess_trans(indexx.exp)=log(guess(indexx.exp));
guess_trans(indexx.oneabove)=log(1-guess(indexx.oneabove));
end