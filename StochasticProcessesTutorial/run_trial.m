
function [spikecounts,spikerates] = run_trial(mu,sigma,c,g_exc1,g_inh1,g_exc2,g_inh2,g_excall,g_inhall)
%
% computes spike counts and spike rates of two correlated cells for a single trial run
% mu and sigma are mean and variance of the input current
% c is the correlation level (netween 0 and 1)
% the g_e's and g_i's are excitatory and inhibitory currents for cell 1, 2 and correlated input (*all)
% 


spikecounts = [0;0];
spikerates = [0;0];
    
[v_m,spiketimes] = corr_leaky_if(mu,sigma,c,g_exc1,g_inh1,g_excall,g_inhall);
spikecounts(1) = length(find(spiketimes>0.5&spiketimes<3.5));
spikerates(1) = spikecounts(1)/3.0;
clear spiketimes; clear v_m;

[v_m,spiketimes] = corr_leaky_if(mu,sigma,c,g_exc2,g_inh2,g_excall,g_inhall);
spikecounts(2) = length(find(spiketimes>0.5&spiketimes<3.5));
spikerates(2) = spikecounts(2)/3.0;
clear spiketimes; clear v_m;

