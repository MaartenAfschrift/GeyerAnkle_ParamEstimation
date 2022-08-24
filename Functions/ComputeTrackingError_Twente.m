function [RMSE_Unp,RMSE_P,RMSE_Push,RMSE_Pull] = ComputeTrackingError_Twente(Dat,nUnp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

iUnp = 1 : Dat.iDatSet(nUnp+1);
iP = Dat.iDatSet(nUnp+1):Dat.iDatSet(end);

% get indices of forward and backward perturbations
iPush = find(Dat.Set.FileInfo>4);
iPull = find(Dat.Set.FileInfo<5 & Dat.Set.FileInfo>1);
IndexPush = [];
IndexPull = [];
for i = 1:length(iPush)
    IndexPush = [IndexPush Dat.iDatSet(iPush(i)):Dat.iDatSet(iPush(i)+1)-1];
end
for i = 1:length(iPull)
    IndexPull = [IndexPull Dat.iDatSet(iPull(i)):Dat.iDatSet(iPull(i)+1)-1];
end

if isfield(Dat.Set,'ScaleTroque') && Dat.Set.ScaleTroque
    Tid = Dat.Tid./(Dat.Set.BodyMass*Dat.Set.Height).*(70*1.75);
else
    Tid = Dat.Tid;
end

RMSE_Unp = nanmean(abs(Tid(iUnp)-Dat.Tmus(iUnp)));
RMSE_P = nanmean(abs(Tid(iP)-Dat.Tmus(iP)));
if ~isempty(IndexPush)
    RMSE_Push = nanmean(abs(Tid(IndexPush)-Dat.Tmus(IndexPush)));
else
    RMSE_Push = NaN;
end
if ~isempty(IndexPull)
    RMSE_Pull = nanmean(abs(Tid(IndexPull)-Dat.Tmus(IndexPull)));
else
    RMSE_Pull = NaN;
end


end

