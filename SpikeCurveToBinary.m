function BinarySpikeMat = SpikeCurveToBinary(spikeMat,tVec)
% SPIKECURVETOBINARY takes a spike train vector input (spikeMat) along 
% with its time vector (tVec) and converts the spike to a value of 1 and
% all other values to 0.
% Call Format: BinarySpikeMat = SpikeCurveToBinary(spikeMat,tVec)

if sum(sum(abs(spikeMat) > 1))
    BinarySpikeMat = zeros(size(spikeMat));
    for i = 1:size(spikeMat,1)
        [pks,locs] = findpeaks(double(spikeMat(i,:)));
        spikeVec = zeros(1,length(tVec));
        spikeVec(locs) = 1;
        BinarySpikeMat(i,:) = spikeVec;
    end
else
    BinarySpikeMat = spikeMat;
end

end