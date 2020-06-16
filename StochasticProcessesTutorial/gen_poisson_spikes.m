function poisson_train = gen_poisson_spikes(rate,duration)

global timestep

%Generate poisson process.
%Generates series of interspike intervals; advances over time until
%end duration is reached.

poisson_train = zeros(1,(duration/timestep));
if(rate>0)
    spike_time = - log(1-rand(1))/rate;
    index = 1;
    while(spike_time < duration)
        poisson_train(ceil(spike_time/timestep)) = poisson_train(ceil(spike_time/timestep)) +1;
        spike_time = spike_time -log(1-rand(1))/rate;
    end
end
