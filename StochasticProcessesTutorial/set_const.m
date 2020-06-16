%Fundamental LIF Cell Properties
global capacitance;
capacitance = 100*10^-12
global g_leak;
g_leak = 10 * 10^-9
global e_leak;
e_leak = -65 * 10^-3
global v_threshold;
v_threshold = -55 *10^-3  
global v_reset;
v_reset = -60 * 10^-3     
global t_refractory;
t_refractory = 2 * 10^-3

%Properties for Exp. IF Model
global delta;
delta = 2 * 10^-3;

%Properties for including Synaptic Conductances
global e_exc;
e_exc = 0 * 10^-3
global e_inh;
e_inh = -80 * 10^-3

%Values for practical running of simulation
global timestep;
timestep = 0.1 * 10^-3

%Stimulus timing properties.
global stim_start;
stim_start = 0.5
global stim_end;
stim_end = stim_start+3;
global duration;
duration = stim_end+0.5;