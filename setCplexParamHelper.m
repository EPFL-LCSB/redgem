function [mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters)


if strcmp(CplexParameters,'CplexDefault')
    mipTolInt = 1e-5;
    scalPar = 0;
    feasTol = 1e-6;
    emphPar = 0;
elseif strcmp(CplexParameters,'LCSBDefault')
    % Empty will also just fetch the LCSB default
    mipTolInt = 1e-9;
    scalPar = -1;
    feasTol = 1e-9;
    emphPar = 1;
elseif strcmp(CplexParameters,'LCSBScaling')
    mipTolInt = 1e-9;
    scalPar = 0;
    feasTol = 1e-9;
    emphPar = 1;
elseif strcmp(CplexParameters,'redGEM_m7')
    mipTolInt = 1e-7;
    scalPar = -1;
    feasTol = 1e-7;
    emphPar = 1;
elseif strcmp(CplexParameters,'redGEM_m8')
    mipTolInt = 1e-8;
    scalPar = -1;
    feasTol = 1e-8;
    emphPar = 1;
else
    error('Wrong option!')
end


end