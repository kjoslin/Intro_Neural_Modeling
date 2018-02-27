function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials)

% POISSONSPIKEGEN accepts a frequency (fr), length of simulation (tSim),
% and the number of trials (nTrials) and then generates a matrix with
% nTrials of poisson spike trains of the given frequency. It ouputs the
% spike train matrix (spikeMat) and the time vector (tVec) of the spike
% trains.
% Call format: [spikeMat, tVec] = PoissonSpikeGen(fr, tSim, nTrials)


dt = 1/1000; %time step in seconds
nBins = floor(tSim/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = 0:dt:tSim-dt;




