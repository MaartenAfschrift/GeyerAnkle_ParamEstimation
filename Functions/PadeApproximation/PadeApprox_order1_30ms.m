function [u,dydt] = PadeApprox_order1_30ms(e,y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dydt = -66.67.*y + 16.*e;
u = 8.33.*y-e;


end

