function[surrTensor] = computeTME(dataTensor,surrogate_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <TME>
% Copyright (C) 2017 Gamaleldin F. Elsayed and John P. Cunningham 
%       (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demonstration of how to use this code package 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rng('shuffle', 'twister') % randomize the seed

% surrogate_type = 'surrogate-N';
% this data is a preprocessed version of the data available online from 
% http://stat.columbia.edu/~cunningham/pdf/ChurchlandNature2012_code.zip


%% quantify primary features of the original data
[targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures(dataTensor);
%% sample many surrogates and build null distribution of summary statistics
numSurrogates = 100;
params = [];


if strcmp(surrogate_type, 'surrogate-T')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = [];
    params.margCov{3} = [];
    params.meanTensor = M.T;
elseif strcmp(surrogate_type, 'surrogate-TN')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = [];
    params.meanTensor = M.TN;
elseif strcmp(surrogate_type, 'surrogate-TNC')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = targetSigmaC;
    params.meanTensor = M.TNC; 
elseif strcmp(surrogate_type, 'surrogate-N')
    params.margCov{1} = [];
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = [];
    params.meanTensor = M.N;
elseif strcmp(surrogate_type, 'surrogate-NC')
    params.margCov{1} = [];
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = targetSigmaC;
    params.meanTensor = M.NC; 
else
    error('please specify a correct surrogate type') 
end

maxEntropy = fitMaxEntropy(params);             % fit the maximum entropy distribution

for i = 1:numSurrogates
%     fprintf('surrogate %d from %d \n', i, numSurrogates)
    [surrTensor] = sampleTME(maxEntropy);       % generate TME random surrogate data.
end


end