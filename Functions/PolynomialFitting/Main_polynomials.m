% This function generates polynomials to approximate muscle-tendon lengths
% and moment arms. The code is from Wouter Aerts and is adapted to be
% used with CasADi.
%
% Author: Antoine Falisse
%
% Datum: 03/04/2018
%
clear all
close all
clc

%% Filenames

MainPath = 'C:\Users\u0088756\Documents\FWO\Software\GitProjects\Geyer_InverseID';
ModelPath = fullfile(MainPath,'Data','PredSim_Data_Antoine','subject1_scaled_MuscleAdj_Antoine.osim');
output_path = fullfile(MainPath,'Results','PolyInfo','MuscleAdj_Antoine');

%% Run the muscle analysis on the dummy motion

% do a muscle analysis
output_path_MA = fullfile(output_path,'MA');
if ~isfolder(output_path_MA); mkdir(output_path_MA); end
% OpenSim_Muscle_Analysis(fullfile(pwd,'dummy_motion_gait1018.mot'),ModelPath,output_path_MA,[0 50])

%% Import data
% We optimize the polynomial coefficients based on muscle-tendon lengths and
% moment arms generated with OpenSim for a wide range of joint positions. To
% this aim, we first create a dummy motion and then perform a muscle
% analysis in OpenSim. Here, we saved the needed data in the structure
% MuscleData_subject1.mat. We provide code for easy reproduction. The user
% should run a muscle analysis based on the setup file provided in
% MuscleAnalysis/SetupMA_subject1.xml
% Extract time and angles from dummy motion
pathmain = pwd;
name_dummymotion = 'dummy_motion_gait1018.mot';
path_dummymotion = pwd;
path_resultsMA = output_path_MA;
dummy_motion = importdata(fullfile(path_dummymotion,name_dummymotion));
% 3 dofs: hip flex r, knee flex r, ankle r
q = dummy_motion.data(:,5:7).*(pi/180);
% Generate random numbers between -1000 and 1000 (°/s)
load(fullfile(path_dummymotion,'dummy_qdot_1018.mat'));
qdot = dummy_qdot(:,:);
% Import data
% lMT
lMT = importdata([path_resultsMA,'/dummy_motion_gait1018_MuscleAnalysis_Length.sto']);
% hip flexion r
MA.hip.flex = importdata([path_resultsMA,...
    '/dummy_motion_gait1018_MuscleAnalysis_MomentArm_hip_flexion_r.sto']);
% knee flexion r
MA.knee.flex = importdata([path_resultsMA,...
    '/dummy_motion_gait1018_MuscleAnalysis_MomentArm_knee_angle_r.sto']);
% ankle r
MA.ankle.flex = importdata([path_resultsMA,...
    '/dummy_motion_gait1018_MuscleAnalysis_MomentArm_ankle_angle_r.sto']);
% Organize MuscleData
MuscleData.dof_names = dummy_motion.colheaders(5:7);
muscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};
MuscleData.muscle_names = muscleNames;
for m = 1:length(muscleNames)
    % lMT
    MuscleData.lMT(:,m) = lMT.data(:,strcmp(lMT.colheaders,muscleNames{m}));
    % hip flexion r
    MuscleData.dM(:,m,1) = ...
        MA.hip.flex.data(:,strcmp(lMT.colheaders,muscleNames{m}));
    % knee flexion r
    MuscleData.dM(:,m,2) = ...
        MA.knee.flex.data(:,strcmp(lMT.colheaders,muscleNames{m}));
    % ankle flexion r
    MuscleData.dM(:,m,3) = ...
        MA.ankle.flex.data(:,strcmp(lMT.colheaders,muscleNames{m}));
end
MuscleData.q = q;
MuscleData.qdot = qdot;

%% Optimize polynomial coefficients
[muscle_spanning_joint_INFO,MuscleInfo] = PolynomialFit(MuscleData);
save(fullfile(output_path,'MuscleData_CE.mat'),'MuscleData');
save(fullfile(output_path,'muscle_spanning_joint_INFO_CE.mat'),'muscle_spanning_joint_INFO');
save(fullfile(output_path,'MuscleInfo_CE.mat'),'MuscleInfo');

%% Create CasADi functions
import casadi.*
% Order mobilities: hip_flex, knee_angle, ankle_angle,
NMuscle = length(MuscleInfo.muscle);
q_leg = 3;
qin     = SX.sym('qin',1,q_leg);
qdotin  = SX.sym('qdotin',1,q_leg);
lMT     = SX(NMuscle,1);
vMT     = SX(NMuscle,1);
dM      = SX(NMuscle,q_leg);
for i=1:NMuscle
    index_dof_crossing  = find(muscle_spanning_joint_INFO(i,:)==1);
    order               = MuscleInfo.muscle(i).order;
    [mat,diff_mat_q]    = n_art_mat_3_cas_SX(qin(1,index_dof_crossing),order);
    lMT(i,1)            = mat*MuscleInfo.muscle(i).coeff;
    vMT(i,1)            = 0;
    dM(i,1:q_leg) = 0;
    nr_dof_crossing     = length(index_dof_crossing);
    for dof_nr = 1:nr_dof_crossing
        dM(i,index_dof_crossing(dof_nr)) = ...
            (-(diff_mat_q(:,dof_nr)))'*MuscleInfo.muscle(i).coeff;
        vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr)) * ...
            qdotin(1,index_dof_crossing(dof_nr)));
    end
end
f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});

%% Check polynomial fit
lMT_out_r = zeros(size(MuscleData.q,1),NMuscle);
vMT_out_r = zeros(size(MuscleData.q,1),NMuscle);
dM_out_r = zeros(size(MuscleData.q,1),NMuscle,q_leg);
for i = 1:size(MuscleData.q,1)
    [out1_r,out2_r,out3_r] = ...
        f_lMT_vMT_dM(MuscleData.q(i,:),MuscleData.qdot(i,:));
    lMT_out_r(i,:) = full(out1_r);
    vMT_out_r(i,:) = full(out2_r);
    dM_out_r(i,:,1) = full(out3_r(:,1));
    dM_out_r(i,:,2) = full(out3_r(:,2));
    dM_out_r(i,:,3) = full(out3_r(:,3));
end

%% lMT
% right
figure()
subplot(3,3,1)
scatter3(MuscleData.q(:,1),MuscleData.q(:,2),lMT_out_r(:,1)); hold on;
scatter3(MuscleData.q(:,1),MuscleData.q(:,2),MuscleData.lMT(:,1));
xlabel('q hip');
ylabel('q knee');
title('HAM');
subplot(3,3,2)
scatter(MuscleData.q(:,2),lMT_out_r(:,2)); hold on;
scatter(MuscleData.q(:,2),MuscleData.lMT(:,2));
xlabel('q knee');
title('BFSH');
subplot(3,3,3)
scatter(MuscleData.q(:,1),lMT_out_r(:,3)); hold on;
scatter(MuscleData.q(:,1),MuscleData.lMT(:,3));
xlabel('q hip');
title('GLU MAX');
subplot(3,3,4)
scatter(MuscleData.q(:,1),lMT_out_r(:,4)); hold on;
scatter(MuscleData.q(:,1),MuscleData.lMT(:,4));
xlabel('q hip');
title('ILIO-PSOAS');
subplot(3,3,5)
scatter3(MuscleData.q(:,1),MuscleData.q(:,2),lMT_out_r(:,5)); hold on;
scatter3(MuscleData.q(:,1),MuscleData.q(:,2),MuscleData.lMT(:,5));
xlabel('q hip');
ylabel('q knee');
title('RF');
subplot(3,3,6)
scatter(MuscleData.q(:,2),lMT_out_r(:,6)); hold on;
scatter(MuscleData.q(:,2),MuscleData.lMT(:,6));
xlabel('q knee');
title('VAS');
subplot(3,3,7)
scatter3(MuscleData.q(:,2),MuscleData.q(:,3),lMT_out_r(:,7)); hold on;
scatter3(MuscleData.q(:,2),MuscleData.q(:,3),MuscleData.lMT(:,7));
xlabel('q knee');
ylabel('q ankle');
title('GAS');
subplot(3,3,8)
scatter(MuscleData.q(:,3),lMT_out_r(:,8)); hold on;
scatter(MuscleData.q(:,3),MuscleData.lMT(:,8));
xlabel('q ankle');
title('SOL');
subplot(3,3,9)
scatter(MuscleData.q(:,3),lMT_out_r(:,9)); hold on;
scatter(MuscleData.q(:,3),MuscleData.lMT(:,9));
xlabel('q ankle');
title('TA');
legend('Polynomial','Model');
suptitle('lMT right');

%% Assert results
for i = 1:NMuscle
    assertLMT.all(:,i) = abs(lMT_out_r(:,i) - MuscleData.lMT(:,i));
    assertdM.hip.flex(:,i) = abs(dM_out_r(:,i,1) - MuscleData.dM(:,i,1));
    assertdM.knee(:,i) = abs(dM_out_r(:,i,2) - MuscleData.dM(:,i,2));
    assertdM.ankle(:,i) = abs(dM_out_r(:,i,3) - MuscleData.dM(:,i,3));
end
assertLMTmax_r = max(max(assertLMT.all));
assertdM.hip.flexmax = max(max(assertdM.hip.flex));
assertdM.kneemax = max(max(assertdM.knee));
assertdM.anklemax = max(max(assertdM.ankle));

%% save polynomial information for individual muscles


% we want to save the information from the soleus and tibialis anterior
% here since we want to use these muscles in the control of an ankle-foot
% exoskeleton.


iSol = find(strcmp(muscleNames,'soleus_r'));
SoleusPoly.coeff = MuscleInfo.muscle(iSol).coeff;
SoleusPoly.order = MuscleInfo.muscle(iSol).order;

iTib = find(strcmp(muscleNames,'tib_ant_r'));
TibPoly.coeff = MuscleInfo.muscle(iTib).coeff;
TibPoly.order = MuscleInfo.muscle(iTib).order;


save(fullfile(output_path,'AnklePoly.mat'),'SoleusPoly','TibPoly');









