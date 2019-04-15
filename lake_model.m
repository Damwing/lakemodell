close all % Close all windows
clear % Clear all variables
clc % Clear the command line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. General constants:
g=9.81; % [m/s^2], gravity
Kcst=0.4; % [-], von Kármán constant

% 2. Lake geometry:
A0=7.8*10^6; % [m^2], lake surface area
hmax=8; % [m], mean depth
V=A0*hmax; % [m^3], lake volume

% 3. Meteorological data:
meteodata = xlsread('Meteo_data_NASA.xlsx');
Ta=meteodata(:,8); % [°C], air temperature
RH=meteodata(:,6); % [%], relative humidity
Wsp=meteodata(:,9); % [m/s], wind speed
R=meteodata(:,5)*10^6/(24*3600); % [W/m^2], shortwave incoming radiation
P=meteodata(:,7)*1000; % [Pa], atmospheric pressure

% 4. Atmospheric variables:
rho_air=1.2; % [kg/m^3], air density
CD=0.001; % [-], drag coefficient

% 5. Heat fluxes variables:
albedo=0.1; % [-], albedo at the lake surface
% Shortwave radiation passing through the lake surface:
Hsw0=(1-albedo)*R; % [W/m^2]
C=0; % cloudiness
Cpw=4200; % [J/kg/K], specific heat of water
Dth=1.4*10^(-7); % [m^2/s], molecular heat diffusivity

% 6. Water properties:
Twinter=4; % [°C], lake temperature at ice-on/ice-off
nu=10^(-6); % [m^2/s], kinematic viscosity
% Thermal expansivity of water as a function of temperature:
alphaT=@(T)10^(-6)*(-65.4891+17.12544*T-0.178155*T^2); % [/K]
% Water density as a function of temperature:
rho_T=@(T) 999.84298 + (65.4891*T-8.56272*T^2+0.059385*T^3)*10^(-3); % [kg/m^3]
% Water density of as a function of temperature and solids concentration:
rho_TCss=@(T,Css)rho_T(T)+Css*10^(-3)*(1-rho_T(T)/rho_p); % [kg/m^3]

% 7. Particles properties:
rho_p=2650; % [kg/m^3], density of particles 
Dp=1.6*10^(-6); % [m], mean diameter of particles 
Vs=@(rho)g*(rho_p-rho)/rho*Dp^2/(18*nu); % [m/s], settling velocity 
C_FFT=10^4; % [mg/L], solids concentration of FFT below the hypolimnion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THERMOCLINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=(1:365)'; % [DOY], time
% Onset of stratification (rise of the thermocline):
stratif_start=datenum(2015,5,4)-datenum(2015,01,01); % [DOY]
% Maximum stratification (shallowest thermocline):
stratif_max=datenum(2015,06,08)-datenum(2015,01,01);
% End of the stratification (fall turnover):
stratif_end=datenum(2015,9,3)-datenum(2015,01,01);
hepi_min=5; % [m], shallowest thermocline
% Depth of the thermocline over time:
h_epi=seasonal_therm(t,[stratif_start stratif_max stratif_end],...
    [hepi_min hmax]); % [m]
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Heat fluxes
H_epi=zeros(length(t)-1,1); % [W/m^2], total heat flux in the epilimnion
H_hypo=zeros(length(t)-1,1); % [W/m^2], total heat flux in the hypolimnion
Hsurf=zeros(length(t)-1,1); % [W/m^2], surface heat flux
% Shortwave radiation absorbed in the epilimnion:
Hsw_epi=zeros(length(t)-1,1); % [W/m^2]
% Shortwave radiation absorbed in the hypolimnion:
Hsw_hypo=zeros(length(t)-1,1); % [W/m^2]
% Vertical diffusive flux between the epilimnion and the hypolimnion:
H_diff=zeros(length(t)-1,1); % [W/m^2]

% 2. Turbulence variables
B0=zeros(length(t)-1,1); % [W/kg], destabilizing buoyancy flux
ustar=zeros(length(t)-1,1); % [m/s], friction velocity
epsilon=zeros(length(t)-1,1); % [W/kg], rate of turbulent dissipation 
Kz=zeros(length(t)-1,1); % [m^2/s], vertical turbulent diffusivity
N2=zeros(length(t)-1,1); % [s^(-2)], buoyancy frequency

% 3. Temperature:
T_epi=zeros(size(t)); T_epi(1)=Twinter; % [°C], epilimnion 
T_hypo=zeros(size(t)); T_hypo(1)=Twinter; % [°C], hypolimnion 

% 4. Solids concentration:
Css_epi=zeros(size(t)); Css_epi(1)=100; % [mg/L], epilimnion
Css_hypo=zeros(size(t)); Css_hypo(1)=200; % [mg/L], hypolimnion

% 5. Density:
rho_epi=zeros(size(t)); % [kg/m^3], epilimnion
rho_hypo=zeros(size(t));% [kg/m^3], hypolimnion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITERATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(t)-1
    if h_epi(k)<hmax % Two boxes
        
        % Densities:
        rho_epi(k)=rho_TCss(T_epi(k),Css_epi(k)); % [kg/m^3]
        rho_hypo(k)=rho_TCss(T_hypo(k),Css_hypo(k)); % [kg/m^3]
        N2(k)=2*g/(rho_epi(k)+rho_hypo(k))*...
            (rho_hypo(k)-rho_epi(k))/(0.5*hmax); % [s^(-2)]
        
        % Surface heat flux:
        Hsurf(k)=-surfheat(T_epi(k),Ta(k),Wsp(k),RH(k),P(k),C); % [W/m^2]

        % Other heat fluxes (to fill):
        
        
        % Mass fluxes (to fill):
        
        
        % New temperatures (to fill):
        
        
        % New concentrations (to fill):

%%%%%%%%%%%%%%%%%%%%%%%%%%        
    elseif T_epi(k)>Twinter % One box but no ice
        if T_epi(k)~=T_hypo(k)
            T_hypo(k)=T_epi(k);
        end
        rho_epi(k)=rho_TCss(T_epi(k),Css_epi(k)); % [kg/m^3]
        rho_hypo(k)=rho_epi(k); % [kg/m^3]
        
        % Surface heat flux:
        Hsurf(k)=-surfheat(T_epi(k),Ta(k),Wsp(k),RH(k),P(k),C); % [W/m^2]
        
        % Other heat fluxes, mass fluxes, new temperatures and
        % concentrations (to fill):
        
%%%%%%%%%%%%%%%%%%%%%%%%%%         
    else % Ice (to fill):
        
    end
end

