function scalar = getPertbScalar(freq,pDuration,trialInd,maxScalar)
% freq is in the unit of cycles per experiment


scalar = -maxScalar .* sin(2 .* pi .* freq ./ pDuration .* (trialInd-1));
end