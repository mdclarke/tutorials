function i_inj = step_current(curr_amp)

global timestep
global stim_start
global stim_end
global duration

i_inj = zeros(1,(duration/timestep));
i_inj((stim_start/timestep):(stim_end/timestep)) = curr_amp * 10^-12;
