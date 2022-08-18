function [RMSE_Unp,RMSE_P] = ComputeTrackingError_Twente(Dat,nUnp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

iUnp = 1 : Dat.iDatSet(nUnp+1);
iP = Dat.iDatSet(nUnp+1):Dat.iDatSet(end);

if isfield(Dat.Set,'ScaleTroque') && Dat.Set.ScaleTroque    
    Tid = Dat.Tid./(Dat.Set.BodyMass*Dat.Set.Height).*(70*1.75);
    RMSE_Unp = nanmean(abs(Tid(iUnp)-Dat.Tmus(iUnp)));
    RMSE_P = nanmean(abs(Tid(iP)-Dat.Tmus(iP)));
else
    RMSE_Unp = nanmean(abs(Dat.Tid(iUnp)-Dat.Tmus(iUnp)));
    RMSE_P = nanmean(abs(Dat.Tid(iP)-Dat.Tmus(iP)));
end


end

