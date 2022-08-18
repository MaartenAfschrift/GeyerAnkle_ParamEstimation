function [u,dydt] = PadeApprox_order1_100ms(e,y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dydt = -20.*y + 8.*e;
u = 5.*y-e;


end

