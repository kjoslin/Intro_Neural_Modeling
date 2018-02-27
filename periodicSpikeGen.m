function [spikeMat, tVec] = periodicSpikeGen(fr, tSim, nTrials)

% PERIODICSPIKEGEN accepts a frequency (fr), length of simulation (tSim),
% and the number of trials (nTrials) and then generates a matrix with
% nTrials of periodic spike trains of the given frequency. It ouputs the
% spike train matrix (spikeMat) and the time vector (tVec) of the spike
% trains.
% Call format: [spikeMat, tVec] = periodicSpikeGen(fr, tSim, nTrials)

dt = 0.001; %time step in seconds

tVec = 0:dt:tSim-dt;
spike = sin(2*pi*fr*tVec); % making a wave with fireing rate frequency

[pks,locs] = findpeaks(spike);

spikeVec = zeros(1,length(tVec));
spikeVec(locs) = 1;

spikeMat = zeros(nTrials,length(spikeVec));

for i = 1:nTrials
    spikeMat(i,:) = spikeVec;
end
    

