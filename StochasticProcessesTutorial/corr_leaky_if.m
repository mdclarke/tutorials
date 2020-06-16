function [v_m, spiketimes] = corr_leaky_if(mu,sigma,cor,g_exc,g_inh,g_excall,g_inhall)
%
% corr_leaky_if returns the vm and spikes for an LIF cell with correlated and uncorrelated (independent) inputs
% sigma sets the variance of the input current
% cor is the onput correlation level (between 0 and 1)
% the g_e's and g_i's are excitatory and inhibitory conductances for independent and correlated input respectively
%
                                          
%Bring in the required parameters from the workspace
global capacitance;
global g_leak;
global e_leak;
global v_threshold;
global v_reset;
global t_refractory;
global timestep;
global e_exc;
global e_inh;

%Set up a vector to store the membrane potential
%It needs to be one longer than the vector for the current input,
%and start at the resting membrane potential.

i_inj = step_current(mu);
v_m = zeros(1,(length(i_inj)+1));
i_tot = zeros(1,(length(i_inj)+1));
v_m(1) = e_leak;

%Just a variable to count the number of spikes - starts at zero.
num_spikes = 0;

%The loop that works its way through the simulation.
index = 2;
while index<=length(v_m)
    %The numerical scheme itself
    
    % i_common is the shared current input
    i_common = sqrt(0.010)*sigma*sqrt(cor)*(g_excall(index-1)*(e_exc - v_m(index-1)) + g_inhall(index-1)*(e_inh - v_m(index-1)));

    % i_ind is the independent input for this cell
    i_ind = sqrt(0.010)*sigma*sqrt(1-cor)*(g_exc(index-1)*(e_exc - v_m(index-1)) + g_inh(index-1)*(e_inh - v_m(index-1)));

    v_m(index) = v_m(index-1) + timestep*(i_inj(index-1) + i_ind + i_common + g_leak*(e_leak - v_m(index-1)))/capacitance;

    %Test for threshold-crossing, then register spike and reset v_m, with refractory period.
    if(v_m(index)>v_threshold)
            num_spikes = num_spikes+1;
            spiketimes(num_spikes) = index*timestep;
            v_m(index) = v_threshold;
            v_m((index+1):(index+(t_refractory/timestep))) = v_reset;
            index = index+(t_refractory/timestep);
    end
    index = index+1;
end

%We want to return a list of spike times.
%In case there were no spikes, we set this to NaN
if(num_spikes==0)
    spiketimes = [];
end

%Print the number of spikes, for convenience
%%num_spikes

%In case it's convenient to plot out the current and voltage directly
%set the figure_plot_on value to 1.
%This will probably become annoying though in later simulations.
figure_plot_on = 0;
if(figure_plot_on == 1)
    t_ref_i = (0:timestep:((length(i_inj)-1)*timestep));
    t_ref_v = (0:timestep:((length(v_m)-1)*timestep));
    figure;
    subplot(2,1,1);
    plot(t_ref_i,i_inj);
    subplot(2,1,2);
    plot(t_ref_v,v_m);
end
