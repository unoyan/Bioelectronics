clear
close all
clc
%% Constants
global V_rest
global T
T=6.3;
V_rest = -60 ; % resting membrane voltage (mV)

% Nernst voltages
E_K = -72.1; % potassium Nernst voltage (mV)
E_Na = 52.4; % sodium Nernst voltage (mV)
E_L = -49.2; % leak Nernst voltage (mV)
C_m = 1; % membrane capacitance (uF/cm^2)

% ion conductivities
G_Na_max = 120; % max sodium conductance (mS/cm^2)
G_K_max = 36; % max potassium conductance (mS/cm^2)

% leakage conductivity
G_L = 0.3; % leakage conductance (mS/cm^2)

d = 300*1e-4; % diameter of cylindric fiber (cm)
%ri = 30 / (pi*d*d);
%re = 20 / (3*pi*d*d);
ri = 0.1061*1e5;
re = 0.2357*1e4;
cons = 1/((2*pi*d)*(ri+re));
c_m = C_m*2*pi*d;
%% Simulation Properties



delta_x = 0.05; % cm
x_min = 0; %  cm
x_max = 4; %  cm
positions = x_min:delta_x:x_max;
n_positions = length(positions);

delta_t = 0.01;
%delta_t = ((2*pi*d)*(ri+re)*C_m*delta_x^2)*0.5; %2.5461
t_end = 125;
timesteps = 0:delta_t:t_end;
timesteps1 = 0:delta_t:t_end+delta_t;
n_timesteps = length(timesteps);

%initial values
n0 = 0.2803;
m0 = 0.0393;
h0 = 0.6798;
V0 = V_rest;

mesh = delta_t/(ri*C_m*delta_x^2);


%% Simulation

% variables
alpha_n = zeros(n_positions, n_timesteps);
beta_n = zeros(n_positions, n_timesteps);
alpha_m = zeros(n_positions, n_timesteps);
beta_m = zeros(n_positions, n_timesteps);
alpha_h = zeros(n_positions, n_timesteps);
beta_h = zeros(n_positions, n_timesteps);

g_k = zeros(n_positions, n_timesteps);
g_na = zeros(n_positions, n_timesteps);

n = zeros(n_positions, n_timesteps);
n_dot = zeros(n_positions, n_timesteps);
m = zeros(n_positions, n_timesteps);
m_dot = zeros(n_positions, n_timesteps);
h = zeros(n_positions, n_timesteps);
h_dot = zeros(n_positions, n_timesteps);
G_K = zeros(n_positions, n_timesteps);
G_Na = zeros(n_positions, n_timesteps);
I_K = zeros(n_positions, n_timesteps);
I_Na = zeros(n_positions, n_timesteps);
I_L = zeros(n_positions, n_timesteps);
V_m_dot = zeros(n_positions, n_timesteps);
E_l = E_L + zeros(n_positions, n_timesteps);
I_m = zeros(n_positions, n_timesteps);
i_p = zeros(n_positions, n_timesteps);
Is = zeros(n_positions, n_timesteps);
Iion = zeros(n_positions, n_timesteps);
Iamp = 24;
stim_duration = 10;
t1s = 545;
t1f = t1s+stim_duration;
t2s = t1f+200;
t2f = t2s+stim_duration;
%11172
Is(2,t1s:t1f)=Iamp;
Is(n_positions-1,t1s:t1f)=-Iamp;
% Is(100:200,length(x)-1)=-100;
Is(2,t2s:t2f)=Iamp;
%11427
Is(n_positions-1,t2s:t2f)=-Iamp;
%Is(2,3100:3200) = 200; 
%Is(3100:3200,length(x)-1)=-200;

% initialize
V_m = V_rest + zeros(n_positions, n_timesteps);
n(:, 1) = n0;
m(:, 1) = m0;
h(:, 1) = h0;
t = 0;
% iterate
for x_index=2:n_positions-1
    for t_index=1:n_timesteps
        %calculate ion gates
        I_m(x_index,t_index) = cons*(((V_m(x_index-1,t_index) -2*V_m(x_index,t_index)+V_m(x_index+1,t_index))/(delta_x^2))-re*i_p(x_index,t_index)) + Is(x_index, t_index);
        g_k(x_index,t_index) = G_K_max*n(x_index,t_index)^4;
        g_na(x_index,t_index) = G_Na_max*m(x_index,t_index)^3*h(x_index,t_index);
        I_K(x_index, t_index) = (V_m(x_index, t_index)-E_K) * G_K_max * n(x_index, t_index)^4;
        I_Na(x_index, t_index) = (V_m(x_index, t_index)-E_Na) * G_Na_max * m(x_index, t_index)^3 * h(x_index, t_index);
        
        E_l(x_index,t_index) = ((g_k(x_index,t_index)+g_na(x_index,t_index) + G_L)*V_rest -(g_k(x_index,t_index)*E_K)+ (g_na(x_index,t_index)*E_Na))/G_L;
        Iion(x_index,t_index) = g_k(x_index,t_index)*(V_m(x_index,t_index)-E_K) + g_na(x_index,t_index)*(V_m(x_index,t_index)-E_Na) + G_L*(V_m(x_index,t_index)-E_l(x_index,t_index));
        I_L(x_index, t_index) = (V_m(x_index, t_index)- E_l(x_index,t_index)) * G_L;
        G_K(x_index,t_index) = G_K_max*n(x_index,t_index)^4;
        G_Na(x_index,t_index) = G_Na_max*m(x_index,t_index)^3*h(x_index,t_index);
        
        % calcultate membrane voltage change
        V_m_dot(x_index, t_index) =  (delta_t/C_m) * (I_m(x_index,t_index)- Iion(x_index,t_index));
        V_m(x_index, t_index+1) = V_m(x_index, t_index) +  V_m_dot(x_index, t_index);
           
        [alpha_n(x_index, t_index), beta_n(x_index, t_index), alpha_m(x_index, t_index), beta_m(x_index, t_index), alpha_h(x_index, t_index), beta_h(x_index, t_index)] = calc_gates(V_m(x_index, t_index));
        n_dot(x_index, t_index) = alpha_n(x_index, t_index) * (1 - n(x_index, t_index)) - beta_n(x_index, t_index) * n(x_index, t_index);
        m_dot(x_index, t_index) = alpha_m(x_index, t_index) * (1 - m(x_index, t_index)) - beta_m(x_index, t_index) * m(x_index, t_index);
        h_dot(x_index, t_index) = alpha_h(x_index, t_index) * (1 - h(x_index, t_index)) - beta_h(x_index, t_index) * h(x_index, t_index);

        n(x_index, t_index+1) = n(x_index, t_index) + n_dot(x_index, t_index) * delta_t;
        m(x_index, t_index+1) = m(x_index, t_index) + m_dot(x_index, t_index) * delta_t;
        h(x_index, t_index+1) = h(x_index, t_index) + h_dot(x_index, t_index) * delta_t;
    
    end
%     plot(timesteps1, V_m(x_index,:))
%     hold on
%     imagesc(V_m)
%     hold on
%     pause(0.1)
end
   
[N out1] = spike_times(V_m(min(find(positions>0.05*x_max)),:),-45);            
%num = min(find(timesteps())); 
num=out1;
%% Plot Results

figure()
subplot(5,1,1)
grid on
plot(timesteps1,V_m(min(find(positions>0.05*x_max)),:),'b', 'Linewidth',2)
xlabel('Time (ms)')
ylabel(' (mV)')
title('Transmembrane voltage vs time at start of fiber')
subplot(5,1,2)
grid on
plot(timesteps1,V_m(min(find(positions>0.95*x_max)),:),'g','Linewidth',2)
xlabel('Time (ms)')
ylabel(' (mV)')
title('Transmembrane voltage vs time at end of fiber')
subplot(5,1,3)
grid on
plot(timesteps1,n(min(find(positions>0.05*x_max)),:), 'Linewidth',2)
hold on
plot(timesteps1, m(min(find(positions>0.05*x_max)),:), 'g', 'Linewidth',2)
hold on
plot(timesteps1,h(min(find(positions>0.05*x_max)),:),'r', 'Linewidth',2)
hold off
legend('n','m','h')
xlabel('Time (ms)')
ylabel('n m h')
title('Gating variables')
% subplot(4,1,4)
% grid on
% plot(timesteps,n(min(find(positions>0.75*x_max)),:), 'Linewidth',2)
% hold on
% plot(timesteps, m(min(find(positions>0.75*x_max)),:), 'g', 'Linewidth',2)
% hold on
% plot(timesteps,h(min(find(positions>0.75*x_max)),:),'r', 'Linewidth',2)
% hold off
% legend(['n','m','h'])
% xlabel('Time (ms)')
% ylabel('Gating variables')
subplot(5,1,4)
grid on
plot(timesteps,I_K(min(find(positions>0.05*x_max)),:),'c', 'Linewidth',2)
hold on
plot(timesteps, I_Na(min(find(positions>0.05*x_max)),:),'r', 'Linewidth',2)
hold on
plot(timesteps,I_L(min(find(positions>0.05*x_max)),:),'g', 'Linewidth',2)
hold off
legend('I_K','I_{Na}','I_L')
xlabel('Time (ms)')
ylabel('(\mu A/cm^2)')
title('Current density')
subplot(5,1,5)
grid on
plot(timesteps,Is(2,:),'b', 'Linewidth',2)
hold on
plot(timesteps,Is(n_positions-1,:),'r', 'Linewidth',2)
hold off
legend('Is(start)','Is(end)')
xlabel('Time (ms)')
ylabel(' (\mu A/cm^2)')
title('Stimulus Current density')
figure()

subplot(3,1,1)

grid on
semilogy(timesteps,G_K(min(find(positions>0.05*x_max)),:),'b', 'Linewidth',2)
hold on
semilogy(timesteps, G_Na(min(find(positions>0.05*x_max)),:),'r', 'Linewidth',2)
hold off
title_x = 'time  t   (ms)';   
title_y = '  g  ( mmho.cm^{-2})';
xlabel(title_x); ylabel(title_y);
legend('g_K', 'g_{Na}')
title('Membrane Conductances at left end of fiber')
subplot(3,1,2)
grid on
semilogy(timesteps,G_K(min(find(positions>0.5*x_max)),:),'b', 'Linewidth',2)
hold on
semilogy(timesteps, G_Na(min(find(positions>0.5*x_max)),:),'r', 'Linewidth',2)
hold off
title_x = 'time  t   (ms)';   
title_y = '  g  ( mmho.cm^{-2})';
xlabel(title_x); ylabel(title_y);
legend('g_K', 'g_{Na}')
title('Membrane Conductances at center of fiber')
subplot(3,1,3)
grid on
semilogy(timesteps,G_K(min(find(positions>0.95*x_max)),:),'b', 'Linewidth',2)
hold on
semilogy(timesteps, G_Na(min(find(positions>0.95*x_max)),:),'r', 'Linewidth',2)
hold off
title_x = 'time  t   (ms)';   
title_y = '  g  ( mmho.cm^{-2})';
xlabel(title_x); ylabel(title_y);
legend('g_K', 'g_{Na}')
title('Membrane Conductances at right end of fiber')
% figure()
% grid on
% plot(positions,V_m(:,num), 'Linewidth',2)
% xlabel('Distance (cm)')
% ylabel('Transmembrane voltage (mV)')
% figure()
% grid on
% plot(positions,n(:,num), 'Linewidth',2)
% hold on
% plot(positions, m(:,num), 'g', 'Linewidth',2)
% hold on
% plot(positions,h(:,num),'r', 'Linewidth',2)
% hold off
% legend(['n','m','h'])
% xlabel('Distance (cm)')
% ylabel('Gating variables')
% figure()
% grid on
% plot(positions,I_K(:,num), 'Linewidth',2)
% hold on
% plot(positions, I_Na(:,num),'r', 'Linewidth',2)
% hold on
% plot(positions,I_L(:,num),'g', 'Linewidth',2)
% hold off
% legend(['I_K','I_Na','I_L'])
% xlabel('Distance (cm)')
% ylabel('Current density (\mu A/cm^2)')

