function [Set] = GetDefaultSettings(Set)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(Set,'BoolPerturb')
    Set.BoolPerturb = 0;
end

% additional reflex based on COM feedback
if ~isfield(Set,'COMfb')
    Set.COMfb = 0;
end


% COM feedback based on deviation from reference trajectory
if ~isfield(Set,'PercGaitCycle')
    Set.PercGaitCycle = 1;
end

% bands on standard deviation when using Set.PercGaitCycle
if ~isfield(Set,'COM_std')
    Set.COM_std = false;
end


% Boolean conver tot radans
if ~isfield(Set,'BoolRad')
    Set.BoolRad = 1;
end

% treshold for event detection
if ~isfield(Set,'FyTreshold')
    Set.FyTreshold = 60;
end

% convert to radians
if ~isfield(Set,'InDeg')
    Set.InDeg = 1;
end

% mesh frequency
if ~isfield(Set,'MeshFreq')
    Set.MeshFreq = 150;
end

% muscle names
if ~isfield(Set,'MuscleNames')
    Set.MuscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};
end

% inequality constraint on feedback gains
if ~isfield(Set,'FB_ineq')
    % Bounds on inequality constraint feedback law
    Set.FB_ineq = 10^-4;
end

% scale factor for feedback gains (might improve numerical condition Jacobian)
if ~isfield(Set,'ScaleGains')
    Set.ScaleGains = 0.1;
end

% weights in objective function
if ~isfield(Set,'Weight')
   Set.Weight.Unp = 1;
   Set.Weight.Pert = 1;
   Set.Weight.PertPost = 1;
else
    if ~isfield(Set.Weight,'Unp')
        Set.Weight.Unp = 1;
    end
    if ~isfield(Set.Weight,'Pert')
        Set.Weight.Pert = 1;
    end
    if ~isfield(Set.Weight,'PertPost')
        Set.Weight.PertPost = 1;
    end
end

% inequality constraint on feedback law
if ~isfield(Set,'FB_ineq')
    Set.FB_ineq = 10^-4;
end

% plot Temp results
if ~isfield(Set,'BoolPlot')
    Set.BoolPlot = 0;
end

% minimal e0 reflex gain
if ~isfield(Set,'e0_min')
    Set.e0_min = 0.005;
end

if ~isfield(Set,'ReflexMinBound')
    Set.ReflexMinBound = 0.1;
end

if ~isfield(Set,'ReflexMaxBound')
    Set.ReflexMaxBound = 10;
end

if ~isfield(Set,'diaryName')
    Set.diaryName = fullfile(pwd,'EstimationInfo.txt');
end

% COM velocity feedback only
if ~isfield(Set,'COMfb_velOnly')
    Set.COMfb_velOnly = 0;
end

% time delay for COM feedback
if ~isfield(Set,'COM_TimeDelay')
    % default time delay COM feebdack in ms
    Set.COM_TimeDelay = 0.05;
end

if ~isfield(Set,'BoolGradientDS')
    Set.BoolGradientDS = true;
end

% smooth transition between Single stance and DS
if ~isfield(Set,'G_gradient')
    Set.G_gradient = -10;
end

% Normalised tendon stiffness
if ~isfield(Set,'ATendon')
    Set.ATendon = [20; 35];
end

% Select Cases for gain optimization
if ~isfield(Set,'GainsEstimation')
    % default implementation geyer model
    Set.GainsEstimation = 'GeyerDefault';

    % other options are:
    %   (1) Geyer_COM_COMd: Full Geyer model with COM position and velocity feedback
    %   (2) Geyer_COMd: Full Geyer model with only COM velocity feedback
    %   (3) GeyerDefault_Limited1: Default geyer model with only the following optimization variabls:
    %       - Sol.G, Tib.G, Tib.GSol
    %   (4) Geyer_COMd_Limited2  COMd feedback geyer model with only opt vars:
    %       - Sol.G, Tib.G, Tib.GSol, Sol.COMd, Tib.COMd
    %   (5) Geyer_COMd_GRF COMd feedback with soleus force scaled by the
    %       vertical ground reaction force
    %   (6) Geyer_Default_GRF Default Geyer model with GRFy scaling for the
    %       soleus force feedback
    %   (7) Geyer_COMd_GRF_lim: Same as Geyer_COMd_GRF but without Sol-e0
    %       as optimization parameter (result from identification analysis)
end

% default values used for Soleus feedback when these are not optimization variables
if ~isfield(Set,'Sol')
    Set.Sol.e0 = 0.02;
    Set.Sol.G = 1;
    Set.Sol.COM = 0;
    Set.Sol.COMd = 0;
end

% default values used for Tibialis feedback when these are not optimization variables
if ~isfield(Set,'Tib')
    Set.Tib.e0 = 0.02;
    Set.Tib.G = 0.1;
    Set.Tib.loff = 1;
    Set.Tib.GSol = 0;
    Set.Tib.COM = 0;
    Set.Tib.COMd = 0;
end

% set the intial guess of the feedback gains
if ~isfield(Set,'ReflexGuess')
    Set.ReflexGuess = [];
end

% set the intial guess of the feedback gains
if ~isfield(Set,'G_min')
    Set.G_min = 0;
end

% plot forward simulation with initial guess (in shooting optimization)
if ~isfield(Set,'PlotGuess')
    Set.PlotGuess = false;
end

% test covariance optimization variables (as in Vlugt2010)
if ~isfield(Set,'Vlugt2010')
    Set.Vlugt2010       = false;
end

% Plot figure with covariance
if ~isfield(Set,'PlotVlugt')
    Set.PlotVlugt       = false;
end

% optimize initial state in shooting appraoch
if ~isfield(Set,'OptInitState')
    Set.OptInitState= false;
end

% plot solution on top of guess
if ~isfield(Set,'PlotSolution')
    Set.PlotSolution = false;
end

% constanst STD band around unperturbed COM position and velocity
if ~isfield(Set,'ConstantSTD')
    Set.ConstantSTD = false;
end
if  Set.ConstantSTD
    if ~isfield(Set,'COMref_std')
        Set.COMref_std = 0.01; % 0.01 m band around average
    end
    if ~isfield(Set,'COMdref_std')
        Set.COMdref_std = 0.01; % 0.01 m/s band around average
    end
end

% minmize muscle activations
if ~isfield(Set,'MinA')
    Set.MinA = false;
end
% minmize muscle activations
if ~isfield(Set,'MinA_Weight')
    Set.MinA_Weight = 1;
end

if ~isfield(Set,'Height')
    Set.Height = 1.75;
end

% scales all ID torques to a reference person of 1.75m and 70kg
if ~isfield(Set,'ScaleTroque')
    Set.ScaleTroque = true;
end
    

end

