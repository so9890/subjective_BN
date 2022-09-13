function guess_trans=trans_guess(indexx, guess, params, listt)

% function to transform initial values to pass to function with 
% transformed variables

%-- extract necessary params

guess_trans=guess;

guess_trans(indexx.sqr)=sqrt(guess(indexx.sqr));
guess_trans(indexx.exp)=log(guess(indexx.exp));
guess_trans(indexx.lab)=log((params(listt=='upbarH')-guess(indexx.lab)*0.999)./(guess(indexx.lab)*0.999));
guess_trans(indexx.oneab)=log((1-guess(indexx.oneab))./guess(indexx.oneab));


if isfield(indexx, 'BN')
    guess_trans(indexx.BN)=log((params(listt=='B')-guess(indexx.BN))./guess(indexx.BN));
end
if isfield(indexx, 'BNh')
    guess_trans(indexx.BNh)=log((params(listt=='Bh')-guess(indexx.BNh))./guess(indexx.BNh));
end
if isfield(indexx, 'BNl')
    guess_trans(indexx.BNl)=log((params(listt=='Bl')-guess(indexx.BNl))./guess(indexx.BNl));
end
end