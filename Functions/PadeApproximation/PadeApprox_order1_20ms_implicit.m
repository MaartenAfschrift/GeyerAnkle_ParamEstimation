function [u,err] = PadeApprox_order1_20ms_implicit(e,y,dydt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% implicit dynamics
err = -100.*y + 16.*e - dydt;

% output 
u = 12.5.*y-e;


end

