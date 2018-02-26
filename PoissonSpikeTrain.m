myFirstSpikeTrain = [0 1 0 0 0 1 0 0 1 0];

%% Single Poisson Spike Train Vector
fr = 100; % frequency in Hz
dt = 1/1000; % time step in seconds
nBins = 10; % 10 ms spike train
myPoissonSpikeTrain = rand(1,nBins) < fr*dt;

%% Poisson Spike Train Mulitple

fr = 100; % frequency in Hz
dt = 1/1000; % time step in seconds
nBins = 10; % 10 ms spike train
nTrials = 20; % # of simulations
spikeMat = rand(nTrials, nBins) < fr*dt;

%% 

[spikeMat, tVec] = poissonSpikeGen(30, 1, 20);
plotRaster(spikeMat, tVec*1000);
xlabel('Time(ms)');
ylabel('Trial Number');

%% POISSON SPIKE TRAIN

% simulate the baseline period
[spikeMat_base, tVec_base] = poissonSpikeGen(6, 0.5, 20);
tVec_base = (tVec_base - tVec_base(end))*1000 - 1;
% simulate the stimulus period
[spikeMat_stim, tVec_stim] = poissonSpikeGen(30, 1, 20);
tVec_stim = tVec_stim*1000;
% put the baseline and stimulus periods together
spikeMat = [spikeMat_base spikeMat_stim];
tVec = [tVec_base tVec_stim];
% plot the raster and mark stimulus onset
figure(1)
plotRaster(spikeMat, tVec);
hold all;
plot([0 0], [0 size(spikeMat, 1)+1]);
% label the axes
xlabel('Time (ms)');
ylabel('Trial number');


%PLOT AVERAGED (PSTH) DATA FOR ONE PARTICULAR 
figure(2)
TrialSum_vect = sum(spikeMat); %adds up spike trains across trials figure(2)
plot(tVec,TrialSum_vect,'.')
xlabel('time (ms)')
ylabel('total number of spikes')
figure(3)
plot(tVec,1000*TrialSum_vect/(20*dt),'.')
xlabel('time (ms)')
ylabel('Average firing rate (Hz)')

% CALCULATE FANO FACTOR
% FF should be equal to 1 under ideal simulation, with only 20
% trials we see some deviation from 1, but if you increase to
% 100 trials you get close to 1
NumCount_vec = sum(spikeMat,2); %This finds the number of spikes in each of the trials
FanoFactor = (std(NumCount_vec)^2)/mean(NumCount_vec) %the stdev squared divided by mean


%FIND INTER SPIKE INTERVAL (ISI) AND MAKE HISTOGRAM
tSpike1 = tVec(find(spikeMat(1,:)==1)) %this is finding what times the first trial has a spike (a spike occurs when spikeMat = 1)
isi_vec = diff(tSpike1) %interspike interval vector
figure(4)
hist(isi_vec,8) %Plotting a histogram of isi in 8 equal bins
xlabel('isi (ms)')
ylabel('Number of occurrences')

%FIND MEAN, STDEV, COEFFICIENT OF VARIANCE OF ISI

isi_mean = mean(isi_vec)
isi_std = std(isi_vec)
isi_CV = isi_std/isi_mean
isi_AvgFR = 1000*(1/isi_mean) %the average firing rate for the entire trial








