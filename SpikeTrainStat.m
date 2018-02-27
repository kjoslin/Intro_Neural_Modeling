function [FanoFactor, isi_mean_trial, CV_isi_trial, AvgFR_isi_trial] = ...
    SpikeTrainStat(spikeMat, tVec, trialNum)

% SPIKETRAINSTAT accepts a matrix of several spike train trials
% (spikeMat), the time vector for the trials (tVec), and a specified
% trial number (trialNum) and calculates the spike train statistics.
% It outputs the Fano Factor (FanoFactor)of the spike train matrix and 
% the the average ISI (isi_mean_trial), the CV of the ISI (CV_isi_trial),
% and the average firing rate (AvgFR_isi_trial). It also creates a
% histogram showing thedistribution of the ISI.
% Call format: [FanoFactor, isi_mean_trial, CV_isi_trial, AvgFR_isi_trial] = ...
%    SpikeTrainStat(spikeMat, tVec, trialNum)


% This converts the spikes into binary (1 for spike, 0 for not spike)
% if not already in binary
BinarySpikeMat = SpikeCurveToBinary(spikeMat,tVec);

figure(1)
plotRaster(BinarySpikeMat, tVec*1000);
xlabel('Time(ms)');
ylabel('Trial Number');
title('Raster Plot of the Individual Spikes in the Spike Train Trials');

tSpike1 = tVec(BinarySpikeMat(trialNum,:)==1);
isi_vec = diff(tSpike1); %interspike interval vector

string = sprintf('ISI Distribution for Trial %d', trialNum);

figure(2)
hist(isi_vec,8); %Plotting a histogram of isi in 8 equal bins
title(string);
xlabel('isi (ms)');
ylabel('Number of occurrences');

% CALCULATE FANO FACTOR
% FF should be equal to 1 under ideal simulation, with only 20
% trials we see some deviation from 1, but if you increase to
% 100 trials you get close to 1

fprintf('\nFano Factor for the Entire Simultation:\n')

NumCount_vec = sum(BinarySpikeMat,2); %This finds the number of spikes in each of the trials
FanoFactor = (std(NumCount_vec)^2)/mean(NumCount_vec) %the stdev squared divided by mean

%FIND MEAN, STDEV, COEFFICIENT OF VARIANCE OF ISI

fprintf('\nSpike Train Statistics for Trial %d:\n',trialNum)

isi_mean_trial = mean(isi_vec)
isi_std_trial = std(isi_vec);
CV_isi_trial = isi_std_trial/isi_mean_trial
AvgFR_isi_trial = 1000*(1/isi_mean_trial) %the average firing rate for the entire trial



