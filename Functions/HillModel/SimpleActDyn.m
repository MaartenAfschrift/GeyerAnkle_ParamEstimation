function [dadt] = SimpleActDyn(e,a,tau)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dadt = (e-a)./tau;
end

