function [gains] = UnpackAnkleReflex(R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% vector with reflexes
gains = [R.Sol.e0 R.Sol.G R.Tib.e0 R.Tib.G R.Tib.loff R.Tib.GSol R.Sol.COMd R.Tib.COMd];


end

