
function [r_sc, v, p_sc] = run_ntrials(ntrials,mu,sigma,c,g_exc1,g_inh1,g_exc2,g_inh2,g_excall,g_inhall)
%
% run_ntrials runs many trials of a paired LIF simualtion with correlated input
% and returns the spike count correlation r_sc(mu,sigma) with p-value p_sc(mu,sigma) 
% and mean output spike rate for the pair v(mu,sigma)
% ntrials is the number of trials to run
% mu and sigma are mean and variance of the input current to the cells
% these can be vectors as well, then the output variables will be matrices with entries for each (mu,sigma) pair
% c is the correlation level of the input (between 0 and 1)
% the g_e's and g_i's are matrices of excitatory and inhibitory conductances for cell 1, 2 and correlated input (*all)
% NOTE: there should be *ntrials* number of conductances (can use make_trials to generate these) 
%

 r_sc = zeros(length(mu),length(sigma));
 p_sc = zeros(length(mu),length(sigma));
 v = zeros(length(mu),length(sigma));

% need to call run_trial ntrials times for each pair of mu, sigma values 

 for l = 1:length(mu)

    l                     % keeping track of progress
 
   for k = 1:length(sigma)  

    k              % keeping track of progress
    spikecounts = zeros(2,ntrials);
    spikerates = zeros(2,ntrials);

      for j = 1:ntrials
    
       [spikecounts(:,j),spikerates(:,j)] = run_trial(mu(l),sigma(k),c,g_exc1(j,:),g_inh1(j,:),g_exc2(j,:),g_inh2(j,:),g_excall(j,:),g_inhall(j,:));
    
      end
     
    % calculate pearson's correlation coefficient 
    [r_mat p_mat] = corr(spikecounts');
    r_sc(l,k) = r_mat(1,2);
    % p_sc returns the p-value
    p_sc(l,k) = p_mat(1,2);

    % calculate geometric mean spike rate
    v(l,k) = sqrt(mean(spikerates(1,:))*mean(spikerates(2,:)));
  
    clear spikecounts; clear spikerates;

   end

 end
