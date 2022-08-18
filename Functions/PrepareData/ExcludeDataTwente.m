function [BoolInclude] = ExcludeDataTwente(Dat,headers,Speed,Sid)
%ExcludeDataTwente Removes some data of the Vlutters 2018 dataset because of problems in labelling
% The approach we used here was to exclude all trials that could have problems in
%   - event detection (walking with feet on 1 plate)
%   - impossible COM velocities (due to missing markers typically)
%   Detailed explanation goes here
BoolInclude = true;


% Low pass filter (needed for COM velocity)
Freq = 1000;
tVect = Dat.tspan(1):1/Freq :Dat.tspan(2);
DatSpline = Dat.DataSpline;
DatInt = ppval(DatSpline,tVect)';
[a,b] = butter(4,6./(Freq *0.5),'low');
DatInt_filt = filtfilt(a,b,DatInt);

% get low pass filtered COM velocity
COMd = DatInt_filt(:,strcmp(headers,'COMd'));

% get ankle angle
qa = DatInt_filt(:,strcmp(headers,'qa'));

% get ankle moment
Ta = DatInt_filt(:,strcmp(headers,'Ta'));

% get stride time
dt_stride = Dat.tspan(2)-Dat.tspan(1);

% Rules - stride time
if strcmp(Speed,'Fast') && (dt_stride>1.3  || dt_stride<0.75)
    BoolInclude = false;
end
if strcmp(Speed,'Slow') && (dt_stride>1.8  || dt_stride<1.0)
    BoolInclude = false;
end

% Rules COM velocity
if strcmp(Speed,'Fast') && (max(COMd)>1.6 || min(COMd)<0.6)
    BoolInclude = false;
end
if strcmp(Speed,'Fast') && (max(COMd)>1.4 && min(COMd)<0.9)
    BoolInclude = false;
end

if strcmp(Speed,'Slow') && (max(COMd)>1 || min(COMd)<0.2)
    BoolInclude = false;
end

% Rules ankle angle
if min(qa(1:400))<-0.3 || max(qa)>0.45 || min(qa)<-0.7 || qa(1)<-0.4 || qa(1)>0.15
    BoolInclude = false;
end


if Sid == 8 && max(qa)> 0.3 && strcmp(Speed,'Fast')
    BoolInclude = false;
end

if Sid == 7 &&  min(qa)< -0.38 && strcmp(Speed,'Fast')
    BoolInclude = false;
end

if Sid == 7 &&  min(qa)< -0.38 && strcmp(Speed,'Fast')
    BoolInclude = false;
end
    
    
if Sid == 5 &&  min(qa)< -0.38 && strcmp(Speed,'Fast')
    BoolInclude = false;
end

    
if Sid == 4 &&  min(qa(320:350))< -0.3 && strcmp(Speed,'Fast')
    BoolInclude = false;
end


    
    


% Set of Rules:
%----------------
% related to stride time:
% fast walking 1s +/- 0.25
% slow walking 1.4 +/- 0.3

% related to COM velocity
% only low pass filtered velocites.
% speedBound = [0.4 0.9;
%    0.6 1.6];

% angle angle
% no angles below -0.3 rad during first 40% of gait cycle
% no angle above 0.3 rad in general
% no angles below -0.5 rad in general
% initial angle below -0.15 and above 0.1 is wrong (same for final angle)


end

