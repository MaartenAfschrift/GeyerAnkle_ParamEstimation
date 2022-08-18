function [u,dydt] = PadeApprox_order1_20ms(e,y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dydt = -100.*y + 16.*e;
u = 12.5.*y-e;


end

