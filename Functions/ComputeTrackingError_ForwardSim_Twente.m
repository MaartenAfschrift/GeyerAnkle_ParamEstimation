function [TrackingError] = ComputeTrackingError_ForwardSim_Twente(Dat,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% first gait cycle with perturbations
i0 = 5;
if ~isempty(varargin)
    i0 = varargin{1};
end


ncycl = length(Dat);

% unperturbed gait cycles
if i0>1
    TrackingUnp = nan(i0-1,1);
    for i=1:i0-1
        Terror = Dat(i).Tankle-Dat(i).Tid_ankle;  
        TrackingUnp(i) = rms(Terror);
    end
else
    TrackingUnp = NaN;
end


% Perturbed gait cycles
TrackingP = nan(ncycl-i0+1,1);
for i=i0:ncycl
    Terror = Dat(i).Tankle-Dat(i).Tid_ankle;
    TrackingP(i) = rms(Terror);
end

   
TrackingError.UnpAll = nanmean(TrackingUnp);
TrackingError.UnpAll_gc = TrackingUnp;
TrackingError.PertAll = nanmean(TrackingP); 
TrackingError.PertAll_gc = TrackingP;
    


end

