function synaptic_template = gen_syn_template(tau_short,tau_long,syn_amp)

global timestep

x_ref = 0:timestep:(6*tau_long);
synaptic_template  = exp(-x_ref/tau_long) - exp(-x_ref/tau_short);
synaptic_template = syn_amp*synaptic_template/max(synaptic_template);