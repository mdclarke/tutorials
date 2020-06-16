
function [g_exc1,g_inh1,g_exc2,g_inh2,g_excall,g_inhall] = make_trials(ntrials,syn_exc,syn_inh)
%
% make_trials generates conductances for a simulation of paired LIF cells with correlated input
% ntrials is the number of trials you wish to simulate
% syn_exc and syn_inh are synaptic excitatory and inhibitory conductnace templates
% you can obtain these using gen_syn_template.m
%

global timestep
global duration

% for cell 1
g_exc1 = zeros(ntrials,(duration/timestep));
g_inh1 = zeros(ntrials,(duration/timestep));

% for cell 2
g_exc2 = zeros(ntrials,(duration/timestep));
g_inh2 = zeros(ntrials,(duration/timestep));

% common input
g_excall = zeros(ntrials,(duration/timestep));
g_inhall = zeros(ntrials,(duration/timestep));


for i = 1:ntrials

    % for cell 1
    pe1 = gen_poisson_spikes(7000,4.0); 
    pi1 = gen_poisson_spikes(3000,4.0); 
    t1 = conv(pe1,syn_exc); t1i = conv(pi1,syn_inh);
    t1(length(pe1)+1:end) = []; t1i(length(pi1)+1:end) = [];
    g_exc1(i,:) = t1'; 
    g_inh1(i,:) = t1i';

    % for cell 2
    pe2 = gen_poisson_spikes(7000,4.0); 
    pi2 = gen_poisson_spikes(3000,4.0);         
    t2 = conv(pe2,syn_exc); t2i = conv(pi2,syn_inh);                  
    t2(length(pe2)+1:end) = []; t2i(length(pi2)+1:end) = [];              
    g_exc2(i,:) = t2';
    g_inh2(i,:) = t2i';

    % for common input
    pec = gen_poisson_spikes(7000,4.0);
    pic = gen_poisson_spikes(3000,4.0);     
    tc = conv(pec,syn_exc); tci = conv(pic,syn_inh);              
    tc(length(pec)+1:end) = [];	tci(length(pic)+1:end) = [];          
    g_excall(i,:) = tc';
    g_inhall(i,:) = tci';

end
