function SpikeTrainStat(spikeMat, tVec, trialNum)

tSpike1 = tVec(find(spikeMat(trialNum,:)==1));
isi_vec = diff(tSpike1); %interspike interval vector
figure(4)
hist(isi_vec,8); %Plotting a histogram of isi in 8 equal bins
title('ISI Distribution for Trial 1');
xlabel('isi (ms)');
ylabel('Number of occurrences');

% CALCULATE FANO FACTOR
% FF should be equal to 1 under ideal simulation, with only 20
% trials we see some deviation from 1, but if you increase to
% 100 trials you get close to 1
NumCount_vec = sum(spikeMat,2); %This finds the number of spikes in each of the trials
FanoFactor = (std(NumCount_vec)^2)/mean(NumCount_vec) %the stdev squared divided by mean


%FIND MEAN, STDEV, COEFFICIENT OF VARIANCE OF ISI

isi_mean_trial = mean(isi_vec)
isi_std_trial = std(isi_vec);
CV_isi_trial = isi_std_trial/isi_mean_trial
AvgFR_isi_trial = 1000*(1/isi_mean_trial) %the average firing rate for the entire trial



