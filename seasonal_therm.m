function [h] = seasonal_therm(t,tlim,hlim)
%SEASONAL_THERM Calculate the seasonal thermocline for each time value in t.
% INPUTS : - t (n*1 numerical array): time values;
%           - tlim (3*1 numerical array): values t0 (time of ice-off),
% t1 (time of highest thermocline), t2 (time of ice-on);
%           - hlim (2*1 numerical array): minimum thermocline and lake
%           depth.
%   OUTPUTS: - h (n*1 numerical array): depth of seasonal thermocline

if length(tlim)~=3 || length(hlim)~=2
    error('Wrong size of input arguments')
else
    tlim=sort(tlim); hlim=sort(hlim);
end

t0=tlim(1); t1=tlim(2); t2=tlim(3);
hmin=hlim(1); D=hlim(2);

h=ones(size(t))*D;


% Rise of the thermocline
ar=(D-hmin)/(t1-t0)^2;
br=-2*ar*t1;
cr=hmin+ar*t1^2;
tr=t(t>=t0 & t<=t1);
h(t>=t0 & t<=t1)=ar*tr.^2+br*tr+cr;



% Deepening of the thermocline
ad=(D-hmin)/(t2-t1)^2;
bd=-2*ad*t1;
cd=D+ad*(2*t1*t2-t2^2);
td=t(t>=t1 & t<=t2);
h(t>=t1 & t<=t2)=ad*td.^2+bd*td+cd;



end

