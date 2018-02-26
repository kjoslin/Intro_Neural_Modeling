function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials)
dt = 1/1000; %time step in seconds
nBins = floor(tSim/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = 0:dt:tSim-dt;

% 
% spikeInd = find(spikeMat == 1);
% 
% tSpike = tVec(spikeInd);
% 
% s = diff(tSpike);
% 
% CV = std(tSpike)/mean(tSpike);




