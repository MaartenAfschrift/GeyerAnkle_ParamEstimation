%% Test evaluate polynomials

clear all; close all; clc;

MainPath = 'C:\Users\u0088756\Documents\FWO\Software\GitProjects\Geyer_InverseID';
ModelPath = fullfile(MainPath,'Data','PredSim_Data_Antoine','subject1_scaled_MuscleAdj_Antoine.osim');
output_path = fullfile(MainPath,'Results','PolyInfo','MuscleAdj_Antoine');

load(fullfile(output_path,'AnklePoly.mat'),'SoleusPoly','TibPoly');


q = 0.1;
qd = 0;
i =1;
ndof = 1;
[mat,diff_mat_q]    = n_art_mat_3_eval(q,TibPoly.order,ndof);
% [mat,diff_mat_q]    = n_art_mat_3(q,TibPoly.order);
lMT(i)              = mat*TibPoly.coeff;
dM(i)               =(-(diff_mat_q))*TibPoly.coeff;
vMT(i)              = -dM(i) *qd;


