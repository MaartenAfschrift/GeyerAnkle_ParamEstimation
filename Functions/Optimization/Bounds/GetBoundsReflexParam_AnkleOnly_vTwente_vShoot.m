function [BoundsO] = GetBoundsReflexParam_AnkleOnly_vTwente_vShoot(Set,e0_min,G_min)
%GetBoundsReflexParam_AnkleOnly_vTwente_vShoot Bounds on reflex parameters
%   Detailed explanation goes here

% % specify the bounds of the default reflex gains
% Bounds = [e0_min 0.1;       % Sol.e0    (1)
%     G_min       1.4;        % Sol.G
%     e0_min      0.1;        % Gas.e0
%     G_min       1.4;        % Gas.G
%     e0_min      0.1;        % Tib.e0    (5)
%     G_min       0.5;        % Tib.G
%     0.8         1.2;        % Tib.loff
%     G_min       0.6;        % Tib.GSol
%     e0_min      0.1;        % Vas.e0
%     G_min       2;          % Vas.G     (10)
%     qk_max-0.2  qk_max;     % Vas.KneeOff
%     G_min       10;         % Vas.kOff
%     e0_min      0.1;        % Ham.e0
%     G_min       1;          % Ham.G
%     0.6         1;          % Ham.loff  (15)
%     e0_min      0.1;        % Glut.e0
%     G_min       1;          % Glut.G
%     G_min       0.2;        % Glut.kpTrunk
%     e0_min      0.1;        % Psoas.e0
%     G_min       1;          % Psoas.G   (20)
%     0.7         1.3;        % Psoas.loff
%     G_min       1;          % Psoas.Gham
%     G_min       10;         % Trunk.kp
%     min(qTrunk) max(qTrunk);% Trunk.qref
%     G_min       1.5;        % Trunk.kd  (25)
%     min(qTrunk) max(qTrunk);% Trunk.qrefToeOff
%     G_min       0.2;        % Weight.kbw
%     G_min       0.3;        % Weight.S
%     G_min       5;          % Weight.klean
%     G_min       0.02;       % Glut.HipStim (30)
%     G_min       0.02;       % Psoas.HipStim
% ];

% specify the bounds of the default reflex gains
% Bounds = [e0_min 0.1;       % Sol.e0    (1)
%     0.7         2;          % Sol.G
%     e0_min      0.1;        % Tib.e0    (5)
%     G_min       2;          % Tib.G
%     0.6         1.2;        % Tib.loff
%     G_min       2;        % Tib.GSol
% ];


Bounds.Sol.e0   = [e0_min   0.2];
Bounds.Sol.G    = [0.7      1.7];
Bounds.Tib.e0   = [e0_min   0.4];
Bounds.Tib.G    = [G_min    4];
Bounds.Tib.loff = [0.6      1.2];
Bounds.Tib.GSol = [G_min    1.5];
Bounds.Sol.COM  = [-2       2];
Bounds.Sol.COMd = [-2       2];
Bounds.Sol.COMdd = [-2       2];
Bounds.Tib.COM  = [-2       2];
Bounds.Tib.COMd = [-2       2];
Bounds.Tib.COMdd = [-2       2];


if strcmp(Set.GainsEstimation,'GeyerDefault') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
    BoundsO = [Bounds.Sol.e0;
        Bounds.Sol.G;
        Bounds.Tib.e0;
        Bounds.Tib.G;
        Bounds.Tib.loff;
        Bounds.Tib.GSol];

elseif strcmp(Set.GainsEstimation,'Geyer_COM_COMd')
    BoundsO = [Bounds.Sol.e0;
        Bounds.Sol.G;
        Bounds.Tib.e0;
        Bounds.Tib.G;
        Bounds.Tib.loff;
        Bounds.Tib.GSol;
        Bounds.Sol.COM;
        Bounds.Sol.COMd;
        Bounds.Tib.COM;
        Bounds.Tib.COMd];
elseif strcmp(Set.GainsEstimation,'Geyer_COMd') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF')
    BoundsO = [Bounds.Sol.e0;
        Bounds.Sol.G;
        Bounds.Tib.e0;
        Bounds.Tib.G;
        Bounds.Tib.loff;
        Bounds.Tib.GSol;
        Bounds.Sol.COMd;
        Bounds.Tib.COMd];
elseif strcmp(Set.GainsEstimation,'GeyerDefault_Limited1')
    BoundsO = [Bounds.Sol.G;
        Bounds.Tib.G;
        Bounds.Tib.GSol];
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2')
    BoundsO = [Bounds.Sol.G;
        Bounds.Tib.G;
        Bounds.Tib.GSol;
        Bounds.Sol.COMd;
        Bounds.Tib.COMd];
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
    BoundsO = [Bounds.Sol.G;
        Bounds.Tib.e0;
        Bounds.Tib.G;
        Bounds.Tib.GSol;
        Bounds.Sol.COMd;
        Bounds.Tib.COMd];
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited3')
    BoundsO = [Bounds.Sol.G;
        Bounds.Tib.e0;
        Bounds.Tib.G;
        Bounds.Tib.loff;
        Bounds.Tib.GSol;
        Bounds.Sol.COMd;
        Bounds.Tib.COMd];
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
    BoundsO = [Bounds.Sol.G;
        Bounds.Tib.e0;
        Bounds.Tib.G;
        Bounds.Tib.loff;
        Bounds.Tib.GSol;
        Bounds.Sol.COMd;
        Bounds.Tib.COMd];
elseif  strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
    BoundsO = [Bounds.Sol.e0;
        Bounds.Sol.G;
        Bounds.Tib.e0;
        Bounds.Tib.G;
        Bounds.Tib.loff;
        Bounds.Tib.GSol;
        Bounds.Sol.COM;
        Bounds.Tib.COM;
        Bounds.Sol.COMd;
        Bounds.Tib.COMd;
        Bounds.Sol.COMdd;
        Bounds.Tib.COMdd];
elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
    BoundsO = [Bounds.Sol.e0;
        Bounds.Sol.G;
        Bounds.Tib.e0;
        Bounds.Tib.G;
        Bounds.Tib.loff;
        Bounds.Tib.GSol;
        Bounds.Sol.COMd;
        Bounds.Tib.COMd;
        Bounds.Sol.COMdd;
        Bounds.Tib.COMdd];
end
end

