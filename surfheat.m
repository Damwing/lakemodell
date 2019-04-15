function [Hsurf] = surfheat(Tw,Ta,Wsp,RH,P,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURFHEAT.M:
% Compute the surface heat flux Hsurf [W/m^2] at a given time.
% INPUTS:
% - Tw: lake temperature [°C]
% - Ta: air temperature [°C]
% - Wsp: wind speed [m/s]
% - RH: relative humidity [%]
% - P: atmospheric pressure [Pa]
% - C: couldiness [-]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Parameters
% Parameters:
AL=0.03; % [-]
Ew=0.972; % [-], emissivity of water
sigma=5.67*10^(-8); % [W/m^2/K^4], Boltzmann constant 
Cpa=1005; % [J/kg/K], specific heat of air
Lv=2470*10^3; % [J/kg]
a=1.09;
e_sat=@(T) 611.2*exp(17.62*T./(243.12+T)); % [Pa] with T in [°C]
P=P/100; % [hPa]


%% 1. Compute the heat fluxes 
ea=e_sat(Ta).*RH./100; % [Pa]

% a. Sensible heat
f=4.8+1.98*Wsp+0.28*(Tw-Ta); % [W/m^2/mbar]
Hsens=-1*(-Cpa*P/(0.622*Lv).*f.*(Tw-Ta)); % [W/m^2], P in [hPa]

% b. Latent heat
Hlat=-1*(-f.*(e_sat(Tw)-ea)/100); % [W/m^2]

% c. Shortwave radiation: not a surface flux

% d. Longwave radiation
Ea=a*(1+0.17*C^2)*1.24*(ea/100./(Ta+273.15)).^(1/7); % [-]
HlwIN=-1*((1-AL)*sigma*Ea.*(Ta+273.15).^4); % [W/m^2]
HlwOUT=(Ew*sigma*(Tw+273.15).^4);

% e. Total
Hsurf=HlwOUT+HlwIN+Hsens+Hlat;
end

