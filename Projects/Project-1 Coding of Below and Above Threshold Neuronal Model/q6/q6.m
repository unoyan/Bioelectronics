close all
clear all
clc
global Vr
%% FIXED PARAMETERS =======================================================
dt = 1e-5;         % 0.01 ms : time increment
VR = -105e-3;       % resting voltage (V)
Vr = VR*1e3;          % resting voltage (mV)
VNa = 52.4e-3;     % reversal voltage for Na+ (V)
VK = -72.1e-3;     % reversal voltage for K+ (V)
VL = -49.2e-3;     % reversal voltage for leakage Nernst potential.(V)     
Vthr = VR + 15e-3 ;         % Vthreshold=15;
Cm = 1e-6;         % membrane capacitance/area  (F.m^-2)
R = 10;
tmin = 0;          % starting time
tmax = 10e-3;      % finishing time (s)  default  5e-3

gKmax = 36e-3;     % K+ conductance (ohm^-1.cm^-2)
gNamax = 120e-3;   % Na+ conductance (ohm^-1.cm.-2)
gLmax = 0.3e-3;    % max leakage conductance (ohm-1.cm-2) 

Jext_max = 0;  % max current density for ext stimulus (A.cm^-2) -1e-4; 

ts = 0;            % stimulus ON  
tf = 0.15e-3;       % stimulus OFF 

sf = 1e3;          % scale factor for consersion  v to mV and s to ms
T = 6.3;           % temperature (deg C) default 18.5

fs = 14;           %fontsize

%% SETUP ==================================================================
t = tmin:dt:tmax;
num = length(t);
num1 = min(find(t > ts));       % index for stimulus ON
num2 = min(find(t > tf));       % index for stimulus OFF

Jext = zeros(num,1);       % external current density (A.cm^-2)
JNa  = zeros(num,1);       % Na+ current density (A.cm^-2)
JK   = zeros(num,1);       % K+  current density (A.cm^-2)
JL   = zeros(num,1);       % leakage current density (A.cm^-2)
Jm   = zeros(num,1);       % membrane current (A.cm^-2)
V    = zeros(num,1);       % membrane potential (V)
gNa  = zeros(num,1);       % Na+ conductance
gK   = zeros(num,1);       % K+ conductance
gL   = ones(num,1);        % gL conductance
n    = zeros(num,1);       % K+ gate parameter
m    = zeros(num,1);       % Na+ gate parameter
h    = zeros(num,1);       % Na+ gate parameter

V(1) = VR;                   % initial value for membrane potential

Jext(num1:num2) = Jext_max;  % external stimulus current
%Jext(num3:num4) = Jext_max;  % external stimulus current

% Initial Values
n(1) = 0.2803;
m(1) = 0.0393;
h(1) = 0.6798;

gK(1)  = gKmax * n(1)^4;
gNa(1) = gNamax * m(1)^3 * h(1);
gL = gLmax .* gL;

JK(1)  = gK(1)  * (V(1) - VK);
JNa(1) = gNa(1) * (V(1) - VNa);
JL(1)  = gL(1) * (V(1) - VL); %VR - 10.6e-3
Jm(1)  = JNa(1) + JK(1) + JL(1);

V(1) = VR + (dt/Cm) * (-JK(1) - JNa(1) - JL(1) + Jext(1));

for cc = 1 : num-1
    
[ An Am Ah ] = alpha(V(cc)*1000, T);
[ Bn Bm Bh ] = beta(V(cc)*1000, T);
An = sf * An;   Am = sf * Am;   Ah = sf * Ah;  
Bn = sf * Bn;   Bm = sf * Bm;   Bh = sf * Bh; 

n(cc+1) = n(cc) + dt * (An *(1-n(cc)) - Bn * n(cc)); 
m(cc+1) = m(cc) + dt * (Am *(1-m(cc)) - Bm * m(cc)); 
h(cc+1) = h(cc) + dt * (Ah *(1-h(cc)) - Bh * h(cc)); 

gK(cc+1) = n(cc+1)^4 * gKmax;
gNa(cc+1) = m(cc+1)^3 * h(cc+1) * gNamax;

JK(cc+1)  = gK(cc+1)  * (V(cc) - VK);
JNa(cc+1) = gNa(cc+1) * (V(cc) - VNa);
JL(cc+1)  = gL(cc+1) * (V(cc) - VL);
Jm(cc+1)  = JNa(cc+1) + JK(cc+1) + JL(cc+1);

V(cc+1) = V(cc) + (dt/Cm) * (-JK(cc+1) - JNa(cc+1) - JL(cc+1) + Jext(cc+1));

end

figure(1)     % current --------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.02 0.65 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'current densities   (mA.cm ^{-2})';
%title_main = 'Charging a Capacitor';
%tt1 = 'time constant   \tau  =';
%tt2 = num2str(1e3*tau,3);   tt3 = '  ms';
%tt = [tt1 tt2 tt3];

x = t.*sf;   y = Jext.*sf;
plot(x,y,'linewidth',2);   %  Current - ext
%text(1.5,50,tt);
xlabel(title_x); ylabel(title_y);
%title(title_main);

hold on
x = t.*sf;   y = JNa.*sf;
plot(x,y,'r','linewidth',2);   %  Current - Na+

x = t.*sf;   y = JK.*sf;
plot(x,y,'m','linewidth',2);   %  Current - K+

x = t.*sf;   y = JL.*sf;
plot(x,y,'c','linewidth',2);   %  Current - leakage

x = t.*sf;   y = Jm.*sf;
plot(x,y,'k','linewidth',2);   %  Current - K+

%h_L = legend('J_{ext}','J_{Na}','J_K','J_L','J_m');

figure(2)     % voltage --------------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.3 0.65 0.25 0.25]);
title_x = 'time  t   (ms)';   title_y = 'membrane voltage  V (mV)';

x = t.*sf;   y = V.*sf;
plot(x,y,'k','linewidth',2);   % membrane voltage
xlabel(title_x); ylabel(title_y);
grid on
hold on 

figure(3)     % conductancies  --------------------------------------------
set(gcf,'units','normalized');
set(gcf,'position',[0.58 0.65 0.25 0.25]);

title_x = 'time  t   (ms)';   title_y = 'conductance  g  ( mmho.cm^{-2})';

x = t.*sf;   y = gNa.*sf;
plot(x,y,'r','linewidth',2);   % conductance  Na+
hold on
x = t.*sf;   y = gK.*sf;
plot(x,y,'m','linewidth',2);   % conductance  K+

xlabel(title_x); ylabel(title_y);
grid on
legend('g_{Na}','g_K')


figure(4)
set(gcf,'units','normalized');
set(gcf,'position',[0.1 0.1 0.8 0.8]);

subplot(3,1,1)
set(gca,'fontsize',fs);
title_x = 'time  t   (ms)';   title_y = 'J   (mA.cm ^{-2})';

x = t.*sf;   y = Jext.*sf;
plot(x,y,'linewidth',2);   %  Current - ext
%text(1.5,50,tt);
xlabel(title_x); ylabel(title_y);
%title(title_main);

hold on
x = t.*sf;   y = JNa.*sf;
plot(x,y,'r','linewidth',2);   %  Current - Na+

x = t.*sf;   y = JK.*sf;
plot(x,y,'m','linewidth',2);   %  Current - K+

x = t.*sf;   y = JL.*sf;
plot(x,y,'c','linewidth',2);   %  Current - leakage

x = t.*sf;   y = Jm.*sf;
plot(x,y,'k','linewidth',2);   %  Current - K+

h_L = legend('J_{ext}','J_{Na}','J_K','J_L','J_m');

subplot(3,1,2)
set(gca,'fontsize',fs);
title_x = 'time  t   (ms)';   title_y = '  g  ( mmho.cm^{-2})';

x = t.*sf;   y = gNa.*sf;
plot(x,y,'r','linewidth',2);   % conductance  Na+
hold on
x = t.*sf;   y = gK.*sf;
plot(x,y,'m','linewidth',2);   % conductance  K+

xlabel(title_x); ylabel(title_y);
grid on
legend('g_{Na}','g_K')

subplot(3,1,3)
set(gca,'fontsize',fs);
title_x = 'time  t   (ms)';   title_y = ' V_m (mV)';

x = t.*sf;   y = V.*sf;
plot(x,y,'linewidth',2);   % membrane voltage
xlabel(title_x); ylabel(title_y);
grid on
%hold on 
figure()
plot(t,n,'linewidth',2)
set(gca,'fontsize',fs);
%xlim([0 20])
hold on
plot(t,m,'g','linewidth',2)
%xlim([0 20])
hold on
plot(t,h,'r','linewidth',2)
%xlim([0 20])
legend('n','m','h')
xlabel('time(s)')
ylabel('gates')
grid on
%[N out1] = spike_times(V,Vthr);
%disp('Time interval from start of stimulus to the peak of subsequent action potential length for stimulus current with amplitude 200:')
%disp(t(out1));

