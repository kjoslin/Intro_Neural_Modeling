function [spikeMat, tVec] = periodicSpikeGen(fr, tSim, nTrials)

dt = 0.01; %time step in seconds

tVec = 0:dt:tSim-dt;
spike = sin(2*pi*fr*tVec); % making a wave with fireing rate frequency

figure(1)
plot(tVec,spike,'g')

[pks,locs] = findpeaks(spike);

spikeVec = zeros(1,length(tVec));
spikeVec(locs) = 1;

spikeMat = zeros(nTrials,length(spikeVec));

for i = 1:nTrials
    spikeMat(i,:) = spikeVec;
end
    

% figure(2)
% plot(tVec,spikeVec,'g')
%  
%  %FIND INTER SPIKE INTERVAL (ISI) AND MAKE HISTOGRAM
% tSpike = tVec(locs); %this is finding what times the first trial has a spike (a spike occurs when spikeMat = 1)
% isi_vec = diff(tSpike); %interspike interval vector
% figure(3)
% hist(isi_vec,8) %Plotting a histogram of isi in 8 equal bins
% xlabel('isi')
% ylabel('Number of occurrences')
% 
% %FIND MEAN, STDEV, COEFFICIENT OF VARIANCE OF ISI
% 
% isi_mean = mean(isi_vec)
% isi_std = std(isi_vec);
% isi_CV = isi_std/isi_mean
% isi_AvgFR = 1000*(1/isi_mean) %the average firing rate for the entire trial
% 
% 

    
    




