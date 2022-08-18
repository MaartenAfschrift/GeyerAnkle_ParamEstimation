function [RefCycle] = GetRefCycle_Twente(AllRefData,headers,varargin)
%GetRefCycle_Twente Gets all the input information of the unperturbed walking cycling
% in the Vlutters 2018 dataset
%   Detailed explanation goes here

% mass
if ~isempty(varargin)
    Mass = varargin{1};
    if length(Mass)> 1
        error('Update function GetRefCycle_Twente, please update the inpute arguments');
    end
else
    Mass = 70;
end


% velocity bounds filter to exclude "bad" COM data (typically errors in labelling)
Bound_maxVel = 2;
Bound_minVel = -1;
if length(varargin)>1
    Bounds = varargin{2};
    Bound_maxVel = max(Bounds);
    Bound_minVel = min(Bounds);
end


ntr = length(AllRefData.trial);
COMRef = nan(100,ntr);
COMdRef = nan(100,ntr);
COMddRef = nan(100,ntr);
qaRef = nan(100,ntr);
TaRef = nan(100,ntr);

Bound_maxVel = 1.4;
Bound_minVel = -1;
for i=1:ntr
    if isempty(AllRefData.trial(i).NaNs)
        tSel = linspace(AllRefData.trial(i).tspan(1), AllRefData.trial(i).tspan(end),100);
        DatSel = ppval(AllRefData.trial(i).DataSpline,tSel)';
        COM = DatSel(:,strcmp(headers,'COM'));
        COMd = DatSel(:,strcmp(headers,'COMd'));
        qa = DatSel(:,strcmp(headers,'qa'));
        Ta = DatSel(:,strcmp(headers,'Ta'));
        Fx1 =DatSel(:,strcmp(headers,'GRFLy'));
        Fx2 =DatSel(:,strcmp(headers,'GRFRy'));
        COMdd = (Fx1+Fx2)./Mass;
        if any(COMd>Bound_maxVel | COMd<Bound_minVel)
            % do not use these signals
        else
            COMRef(:,i) = COM - COM(1,:);
            COMdRef(:,i) = COMd;            
        end
        qaRef(:,i) = qa;
        TaRef(:,i) = Ta;
        COMddRef(:,i) =  COMdd;        
    end
end

RefCycle.COMmean  = nanmean(COMRef,2);
RefCycle.COMdmean = nanmean(COMdRef,2);
RefCycle.COMddmean = nanmean(COMddRef,2);

RefCycle.COMstd   = nanstd(COMRef,[],2);
RefCycle.COMdstd  = nanstd(COMdRef,[],2);
RefCycle.COMddstd  = nanstd(COMddRef,[],2);


tInt = linspace(0,1,100);
RefCycle.Spline_COM         = spline(tInt,RefCycle.COMmean);
RefCycle.Spline_COMd        = spline(tInt,RefCycle.COMdmean);
RefCycle.Spline_COMdd        = spline(tInt,RefCycle.COMddmean);

RefCycle.Spline_COMstd      = spline(tInt,RefCycle.COMstd);
RefCycle.Spline_COMdstd     = spline(tInt,RefCycle.COMdstd);
RefCycle.Spline_COMddstd     = spline(tInt,RefCycle.COMddstd);


RefCycle.qa                 = qaRef;
RefCycle.Ta                 = TaRef;
RefCycle.COMRef             = COMRef;
RefCycle.COMdRef            = COMdRef;
RefCycle.COMddRef           = COMddRef;


end

