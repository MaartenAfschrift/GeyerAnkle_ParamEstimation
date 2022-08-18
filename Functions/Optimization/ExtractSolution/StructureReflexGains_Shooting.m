function [Reflex] =StructureReflexGains_Shooting(Set,OptReflex,ConstReflex)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



if strcmp(Set.GainsEstimation,'GeyerDefault') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
    Sol.e0      = OptReflex(1);
    Sol.G       = OptReflex(2);
    Tib.e0      = OptReflex(3);
    Tib.G       = OptReflex(4);
    Tib.loff    = OptReflex(5);
    Tib.GSol    = OptReflex(6);
    Sol.COMd    = 0;
    Tib.COMd    = 0;
elseif strcmp(Set.GainsEstimation,'Geyer_COM_COMd')
    Sol.e0      = OptReflex(1);
    Sol.G       = OptReflex(2);
    Tib.e0      = OptReflex(3);
    Tib.G       = OptReflex(4);
    Tib.loff    = OptReflex(5);
    Tib.GSol    = OptReflex(6);
    Sol.COM     = OptReflex(7);
    Sol.COMd    = OptReflex(8);
    Tib.COM     = OptReflex(9);
    Tib.COMd    = OptReflex(10);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF')
    Sol.e0      = OptReflex(1);
    Sol.G       = OptReflex(2);
    Tib.e0      = OptReflex(3);
    Tib.G       = OptReflex(4);
    Tib.loff    = OptReflex(5);
    Tib.GSol    = OptReflex(6);
    Sol.COMd    = OptReflex(7);
    Tib.COMd    = OptReflex(8);
elseif strcmp(Set.GainsEstimation,'GeyerDefault_Limited1')
    Sol.e0      = ConstReflex(1);
    Sol.G       = OptReflex(1);
    Tib.e0      = ConstReflex(2);
    Tib.G       = OptReflex(2);
    Tib.loff    = ConstReflex(3);
    Tib.GSol    = OptReflex(3);
    Sol.COMd    = 0;
    Tib.COMd    = 0;
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2')
    Sol.e0      = ConstReflex(1);
    Sol.G       = OptReflex(1);
    Tib.e0      = ConstReflex(2);
    Tib.G       = OptReflex(2);
    Tib.loff    = ConstReflex(3);
    Tib.GSol    = OptReflex(3);
    Sol.COMd    = OptReflex(4);
    Tib.COMd    = OptReflex(5);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
    Sol.e0      = ConstReflex(1);
    Sol.G       = OptReflex(1);
    Tib.e0      = OptReflex(2);
    Tib.G       = OptReflex(3);
    Tib.loff    = ConstReflex(2);
    Tib.GSol    = OptReflex(4);
    Sol.COMd    = OptReflex(5);
    Tib.COMd    = OptReflex(6);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
    Sol.e0      = ConstReflex(1);
    Sol.G       = OptReflex(1);
    Tib.e0      = OptReflex(2);
    Tib.G       = OptReflex(3);
    Tib.loff    = OptReflex(4);
    Tib.GSol    = OptReflex(5);
    Sol.COMd    = OptReflex(6);
    Tib.COMd    = OptReflex(7);
elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
    Sol.e0      = OptReflex(1);
    Sol.G       = OptReflex(2);
    Tib.e0      = OptReflex(3);
    Tib.G       = OptReflex(4);
    Tib.loff    = OptReflex(5);
    Tib.GSol    = OptReflex(6);
    Sol.COM     = OptReflex(7);
    Tib.COM     = OptReflex(8);
    Sol.COMd    = OptReflex(9);
    Tib.COMd    = OptReflex(10);
    Sol.COMdd   = OptReflex(11);
    Tib.COMdd   = OptReflex(12);
elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
    Sol.e0      = OptReflex(1);
    Sol.G       = OptReflex(2);
    Tib.e0      = OptReflex(3);
    Tib.G       = OptReflex(4);
    Tib.loff    = OptReflex(5);
    Tib.GSol    = OptReflex(6);
    Sol.COMd    = OptReflex(7);
    Tib.COMd    = OptReflex(8);
    Sol.COMdd   = OptReflex(9);
    Tib.COMdd   = OptReflex(10);
end

Reflex.Sol = Sol;
Reflex.Tib = Tib;
    
 
 
end

