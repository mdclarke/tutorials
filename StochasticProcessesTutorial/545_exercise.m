
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stochasticProcessesTutorial.m
%
% Michael N. Shadlen
% Neubeh545  February, 2007

% Version history
% From gamSpikeTrainTutorial (MNS for CSH course 1998)
% update by G. Horwitz (2000)
% update for NB545 2003 by MNS. Many new bits
% update 2005 by MNS. Markov chain section added.
% update Feb 2007. Mainly explanations 
% update Feb 2009. Minor stuff.
% update Jan 2012. Text, function handle, other minor annoyances
%
% Dependencies: 
% 	plot1ras: written by Mike Shadlen
%   plot2ras: written by Greg Horwitz
%   many functions in stats toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Introduction
%
% This tutorial serves as an introduction to the study of stochastic
% processes and their application to spike trains in neurons. It introduces
% a few facts about renewal processes. We introduce the topic with an
% example of a very irregular sequence of spikes and learn about the basic
% statistical description of intervals and counts.  We then apply these
% statistics to the important Poisson point process and extend this to
% gamma processes. Through this series of exercises, we get a sense of
% "disorderliness". We introduce two ideas to  additional concepts to
% describe irregularity: the hazard function and entropy. In the final
% section of the tutorial, we introduce Markov processes and chains.
%
% Before you jump into stochastic processes, make sure that you understand
% what a random variable is, what is meant by its expectation (or mean),
% its variance and its standard deviation. Do you know what a probability
% density function is? Do you know how to use a PDF to calculate
% expectation and variance? If not, take a moment to run through the
% probabilityBasicsTutorial.m
%
% The most critical
% concept that you need to know is how to calculate the expectation of
% some function of a random variable, say g(x), where x is the random
% variable, using the probability density function, f(x). This is covered
% in probabilityBasicsTutorial.m and in the Berg chapter, which covers many
% other essential concepts.


%% Preliminary matlab stuff.
% Make sure you've got plot1ras.m and plot2ras.m on your path. The course
% website also provides other functions that you may need if you do not
% have the matlab stats package on your computer. Type 'help stats' in the
% command line to see if you've got it. 
clear all
rand('state',sum(100*clock))
format short g


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I. Descriptive statistics.
%
% What is a point process?   It is a list of events that are
% indistinguishable from one another, like spikes or radioactive decay, or
% the onset of an alarm,  or a grain of silver on a piece of film. We can
% think of each *event* as having a time or a location.   That means we can
% describe the process by a list of times, sometimes called arrival times.
% We typically characterize a point process by a list of times. In some
% very important cases, we might choose instead to characterize the
% intervals between these times, and/or the number of events in some epoch.
% These latter two descriptions are especially useful when the random
% process that gives rise to the events can be described in some way that
% is not changing as a function of time. An important case is known as a
% renewal process. 
 
% What is a renewal process? It is a point process in which every event
% starts the process all over again. For example, suppose the events under
% consideration are the failure of your hard drive, that is, the times of
% these failures.  Every time your hard drive fails, we replace it with a
% new one. Whatever process describes the time to failure, it just plays
% itself out over and over again. So one way to characterize a renewal
% process is this: a point process whose intervals are drawn randomly from
% a common distribution.  Each interval is independent from the others.  So
% we say that the intervals are iid (independent and identically
% distributed).  They are described by a random variable. Note that saying
% that the intervals are iid does not mean that there is not some
% relationship between pairs of events. For example, it might be the case
% that your hard drive has some likelihood to fail about 1 hour into use,
% but if not, then it might have some random lifespan of months to years.
% The distribution of intervals would be some complicated thing with a mode
% near 1 hour and some second mode in the months to years range. So long as
% we think that every time a hard drive starts out new, it's time to
% failure is described by the same interval distribution, we have a renewal
% process on our hands.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 1. A renewal that looks like a spike train
%
% Here is an example of a renewal process.
% We're going to draw interevent intervals from a lognormal
% distribution.  A lognormal distribution look like this:
dx = .001;
x = [dx:dx:200];
u = 2.5;  % play with this if you want, but I chose 2.5 for a reason
s = 1;
y = lognpdf(x,u,s);
figure(1), clf
plot(x,y);
title('The Lognormal Distribution');
xlabel('Interval duration (ms)');
ylabel('Probability');

% Notice that the lognormal density has non-zero probability for only
% positive interval durations.  This makes sense because the time
% between two events can't possibly be negative -- lengths (whether in 
% time or space) are always non-negative.

% (As aside: if X has a lognormal distribution then log(X) has a normal
% (Gaussian) distribution.  That's where the name, lognormal, comes from.)
% You can appreciate this by changing the abscissa to a log scale 
set(gca,'XScale','log')

% Let's imagine that a neuron fires spikes according to a renewal process. 
% Further, let's imagine that the interspike intervals from this neuron 
% have a lognormal distribution.  We can simulate spikes trains by drawing
% a bunch of interspike intervals at random from the lognormal distribution, 
% and calling each event a spike.  We'll repeat this 100 times to simulate
% 100 trials.

ntrials = 100
n_intervals = 300
intervals = lognrnd(u,s,n_intervals,ntrials);

% Here's a histogram of the interspike intervals:

figure(1), clf
x = [0:500]';
[n,bins] = hist(intervals(:),x);
n = n/(sum(n));
% bar(x(1:end-1),n(1:end-1));
bar(x,n);
set(gca,'XLim',[min(x) max(x)]);
xlabel('Interval')
ylabel('Frequency')
% The relatively high number of interspike intervals that are very large is
% due to the fact that the  lognormal distribution has a very long tail --
% a characteristic  that we'll come back to soon.)

% As expected, the histogram looks a lot like the lognormal 
% distribution that the interspike intervals were drawn from.  
% Let's superimpose the lognormal density to verify the similarity.

figure(1), hold on;
plot(x,lognpdf(x,u,s),'r-','LineWidth',2);
hold off;

% We can use the 'cumsum' (cumulative sum) function to add up the
% interspike intervals and thus find the arrival time of each spike.

arrivals = cumsum(intervals);
tmax = 1000
arrivals(arrivals > tmax) = nan;
figure(1), clf
for i=1:ntrials
	plot1ras(arrivals(:,i),i); hold on;
end
hold off;
set(gca,'YLim',[.5 ntrials+.5],'Box','off')
xlabel('Time (ms)');
ylabel('Trial number');

% Notice that the arrivals look pretty irregular.
% Let's listen to this irregular point process by playing the events 
% as clicks in the computer speaker.  This is similar to what an
% electrophysiologist hears over the audio monitor when recording 
% from a cell in vivo.  We'll play the spikes from the first five trials
% and wait 1 second between trials.
% This may not work (it fails on my computer)
for i = 1:5
    % sound(hist(arrivals(:,i), [0:max(arrivals(:,i))]),1000)
    sound(hist(arrivals(:,i), [0:max(arrivals(:,i))]),1000)
	pause(1)
end

% Now that we've generated some simulated spike trains we can do 
% some statistical analyses of our simulated data.
% First we'll calculate the mean and standard deviation of the 
% interspike intervals.

mean_int = mean(intervals(:))
% which implies that the rate (in events per second) is
1000 ./ mean_int
var_int = var(intervals(:))

% The ratio of standard deviation to mean is called the  'coefficient of
% variation', or CV for short.  In the right context, it is the reciprocal
% of the 'signal to noise ratio'. Can you think of some situations where
% signal to noise ratio is the reciprocal of CV?

CV_int = sqrt(var_int)/mean_int

% Now we compute the theoretical value of the CV (with enough simulated
% spikes the empirical CV will tend to this value).

[m,v] = lognstat(u,s)
sqrt(v)/m

% Keep this value in the back of your mind. In a moment, we'll compare it
% to the CV from a very special  interval distribution: the exponential
% distribution. 

% Neurophysiologists often calculate statistics from spike count data. In
% this next bit we're going to calculate spike count statististics from our
% simulated spike trains. The expected spike count in an epoch is simply
% the duration of the interval divided by the expected interval. Of course,
% we don't expect to get the same count every time we observe the point
% process. Our interest in counts is in the variability in the number of
% events that we actually observe in an epoch.

% We generated 100 spikes for each of the trials in our simulation.  
% Because the spike times are random, this means that each trial
% ends at a different time.  We will now calculate spike count statistics
% during the first 500 ms of each trial.  This is a brief enough epoch
% that even the shortest trial exceeds it.

epochDur = 500;
if (min(max(arrivals)) < epochDur);
	error('Not enough spikes!');
end
counts = sum(arrivals < epochDur);

% Here is the frequency histogram of spike counts.

figure(1)
hist(counts,[min(counts):2:max(counts)]);
xlabel('Number of spikes in epoch')
ylabel('Number of trials');

% We can now calculate the mean and variance of the distribution 
% of spike counts (just as we previously did for the distribution
% of interspike intervals).
mean_counts = mean(counts)
% ... which implies the rate (in events/sec) is
mean_counts / (epochDur/1000)
% The variance of the counts is 
var_counts = var(counts)

% For any renewal process there's a remarkable asymptotic relationship
% between the mean and variance of the interspike interval distribution
% and the mean and variance of the distribution of spike counts.
% First, the mean of the spike count is (asymptotically) equal to the 
% duration of the counting window divided by the mean of the interval
% distribution.

epochDur/mean_int
mean_counts

% Second, the variance of the spike count is (asymptotically) equal
% to the duration of the counting window times the variance of the 
% interval distribution, divided by the mean of the interval distribution
% cubed.  (Note: this converges a lot more slowly than the formula for the
% mean, above).

epochDur*var_int/mean_int^3
var_counts

% Another way of saying this is that, for a renewal process, the square
% of the CV of the interval distribution is asymptotically equal to 
% the Fano factor (variance to mean ratio) of the counts.

CV_int^2
var_counts/mean_counts

% Before we leave renewal processes in general it's worth pointing 
% that as the length of time over which the process is obsevered becomes
% very long (with respect the mean interarrival time) something cool 
% happens to the distribution of counts within that window: it becomes
% more and more Gaussian.  This is a result of the Central Limit theorem.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Example 2.  A Poisson process. 
%
% In this exercise we will use the same simulation technique to 
% generate spike trains from a few different hypothetical neurons.
% The first is a Poisson neuron that fires spikes according to 
% a Poisson process -- a special type of renewal process.  
%
% A Poisson process is a renewal process whose inter-event times
% have an exponential distribution.  To simulate the spiking of
% a Poisson neuron, we are going to draw random numbers from 
% an exponential distribution and place spikes according to these
% random numbers.
%

ntrials = 100;
spike_rate = 50;                 % spikes/sec
meanint = 1000 / spike_rate;     % mean interspike interval (ms)
intervals = exprnd(meanint,100,ntrials);

% Here's the histogram of the interspike intervals that we've 
% drawn.  To no one's great surprise, it looks pretty exponential.
figure(1), clf
[n,bins] = hist(intervals(:),x);
n = n/(sum(n)*.1);
bar(x(1:end-1),n(1:end-1));
set(gca,'XLim',[min(x) max(x)]);
xlabel('Interval duration (msec)')
ylabel('Probability')

% Now we add up the interspike intervals to find the arrival time
% of each spike.  We're also going to play the spikes over the 
% speaker and draw the raster to the screen.
arrivals = cumsum(intervals);
tmax = 1000;
arrivals(arrivals>tmax) = nan;
figure(1),clf
for i=1:ntrials
	plot1ras(arrivals(:,i),i); hold on;
end
hold off;
xlabel('Time (ms)');
ylabel('Trial number');
% Listen to the 1st five trials
% Sorry, the sound is not working on my mac. You may need to skip these lines.
for i = 1:5
    sound(hist(arrivals(:,i), [0:max(arrivals(:,i))]),1000)
	pause(1)
end

% Notice that the spike rasters look qualitatively similar 
% to the neuron whose interspike intervals had a lognormal
% distribution.  The spike trains sound about the same over the 
% speaker.  This is because our eyes and ears aren't so good
% at discriminating between these two kinds of random processes.
% The statistics of the Poisson spike trains are rather different from 
% those we simulated at the begining of this tutorial.
% To illustrate this point, we'll now look at the statistics of the 
% spike counts across trials.

epochDur = 500;
if (min(max(arrivals)) < epochDur);
	error('Not enough spikes!');
end
poisscounts = sum(arrivals < epochDur);

% Here is the frequency histogram of spike counts.
figure(2),clf
bin0 = min([counts(:); poisscounts(:)]);
binLast = max([counts(:); poisscounts(:)]);
binwidth = 2;
subplot(2,1,1)
[n,x] = hist(poisscounts,[bin0:binwidth:binLast]);
bar(x,n/sum(n));
ylabel('Proportion of trials');
xlabel('Number of spikes in epoch')
title('Poisson neuron','FontSize',14);

% For comparison, we'll also display the frequency histogram of
% spike counts for the neuron we simulated earlier -- the one
% whose interspike intervals come from the lognormal distribution.

subplot(2,1,2)
[n,x] = hist(counts,[bin0:binwidth:binLast]);
bar(x,n/sum(n));
xlabel('Number of spikes in epoch')
ylabel('Proportion of trials');
title('non-Poisson neuron','FontSize',14);

% Notice that the number of spikes fired by the Poisson neuron, although
% random from trial to trial, is more consistent that the number of spikes
% fired by the neuron from the begining of the tutorial. To appreciate
% this, compare the width's of the distributions.  This is a consequence of
% the fact that the  coefficient of variation of the lognormal distribution
% (with the  parameters that we used) is bigger than the coefficient of 
% variation of the exponential distribution.
%
% OK, so what's the big deal about the Poisson process? 
% What's so great about having exponentially distributed 
% interevent times?  
%
% The Poisson process has a number of interesting qualities
% two of which we'll mention now.  The first is that 
% the number of events (spikes in our case) arriving during
% an epoch of fixed duration has a distribution that can be
% derived analytically: the probability that the neuron fires 
% 'k' spikes during an interval of length 't' is:

%                    (lambda*t)^k e^(-lambda*t)
% Prob (X(t) = k) = --------------------------
%                                k!

% It can be viewed as a limiting case of sequence of independent events
% whos probability of success is controlled only by knowledge of overall
% rate of occurrence. Suppose I tell you that the overall rate is 10 per
% sec and ask you what the probability is that you will see an event in any
% one msec bin? Your answer had better be 1 in 100. If every msec is
% independent of every other one, then
% The expected number of events is
N = 1000    % possible msec's to observe an event.
p = 10/N
q = 1-p
% and there are 
p*N
% The actual number is random of course. Its distribution is the binomial
% distribution
figure(3), clf
stem([0:30],binopdf([0:30],N,p))
xlabel('Count')
ylabel('Probability')
% Why did I use the stem function for the plot? 
% Here's the mean and variance
[binoMean, binoVar] = binostat(N,p)
% The Poisson distribution of counts arises when we allow the "bin" size to
% become infinitessimally small and we therefore make N infinitely large.
% As N gets big, p gets small and the limit leads us to an exponential. See
% Berg. 
%
% To illustrate the Poisson distribution, we'll superimpose the
% analytically derived probabilities over our histogram.
figure(2)
subplot(2,1,1);
hold on;
x = [bin0:binwidth:binLast];
plot(x,binwidth*poisspdf(x,spike_rate*(epochDur/1000)),'m*');

% Notice that this theoretical distribution matches the spike counts
% from the Poisson neuron pretty well, but it doesn't match the
% spikes counts from the other neuron at all.  Sadly, the spike
% counts from the non-Poisson neuron don't have a distribution with
% a nice analytical solution, at least as far as I know.

% You may have noticed that that superimposed Poisson distribution 
% looks almost like a normal distribution.  The fact is that the 
% Poisson distribution converges to a normal distribution as the
% duration of the counting window (or the rate of events) goes to 
% infinity.  So does the distribution of spike coutns from the other
% neuron.  As was mentioned briefly earlier, the distribution of
% counts from any renewal process converges to a normal distribution 
% as the number of counts becomes large.

%% Distribution of intervals in a Poission process
% The second remarkable property of the Poisson process pertains to its
% distribution of intervals: the exponential distribution. This
% distribution is special because it conforms to what we really mean by a
% process in which the exact time of the next event is unpredictable as
% time passes. Consider an event (spike) that happened to occur at t = 57
% ms, and suppose the rate is 50 events per second. If we consider the
% probability that an event is about to occur -- that is, we think about
% waiting for this "next" event -- in the next moment, this probability is
% the same in the 58th ms and in the 90th etc., so long as we happen to be
% waiting. In other words, we learn nothing about the probability of an
% event in the next moment simply by virtue of the fact that time has
% passed and we're still waiting. Because of this property, the Poisson process
% is sometimes called "memoryless".

% To illustrate this we'll first plot the exponential 
% density function.
figure(1); clf;
t = [0:.1:200];
y = exppdf(t,meanint);
plot(t,y);
xlabel('Interval duration');
ylabel('Probability');

%% Hazard function
% Looking at the probability distribution you might think that  if a spike
% just occurred another one is probably going to happen pretty soon.  After
% all, the greatest probability density is for very short intervals (at the
% left of the graph).  The fact is that the probability of a spike
% happening doesn't  depend at all on when the last spike happened.  They
% are  completely independent, for the Poisson process.
%
% To see this let's consider a conditional probability: the  probability of
% firing a spike at time 't' given that we haven't fired a spike from time
% 0 to 't'.  Assume that a spike was fired at time 0 exactly.  This first
% probability is just the  exponential probability density (it is the
% probability  of observing an interspike interval 't' ms long).  The
% second probability is one minus the cumulative exponential distribution
% from 0 to 't' (the probability that a spike hasn't occurred between time
% 0 and 't').   To find the conditional probability, we divide one by the
% other.

hazard = exppdf(t,meanint)./(1-expcdf(t,meanint));  % note the PDF in numerator and 1-CDF in denominator
clf;
plot(t, hazard)
xlabel('Time (ms)')
ylabel('Prob spike | Prev spike at t=0')
set(gca,'Ylim',[0 .1]);

% The so-called "hazard function" is totally flat.  This means that 
% the probability of firing a spike in the next small time increment
% is the same irrespective of how much time has gone by since the 
% previous spike.  The product of the hazard function at time 't' and 
% a small increment, 'dt', is the probability that an event occurs
% sometime between 't' and 't+dt'.

% Let's compare that to the lognormal distribution (in red)
hazlogn = lognpdf(t,u,s) ./ (1 - logncdf(t,u,s));
hold on
plot(t,hazlogn,'r')
% This hazard tells us that as time passes from the last spike there is a
% favored time of the next spike (~12 msec), but once this passes, the
% likelihood of a spike in the next moment declines, and the longer the
% wait, the less likely a spike will come.

% The hazard rate is really a conditional probability: the probability that
% an event occurs at t, given that it has not occurred yet.
%   h(t) = p(t | "not yet")
%
% Recall the rule for conditional probability
% p(t,y) = p(t|y)p(y)
% Here y stands for "not yet". So
%
%    h(t) = p(t,y) / p(y)
%
% Let's unpack this. p("not yet") is the probability that the event has not
% occured at time t. It's a function of time that's related to the
% cumulative probability of seeing an event by time t. If f(t) is the
% probability density of observing intervals (e.g., the exponential
% distribution) then "not yet" is 1-F(t), where F(t) is the cumulative
% probability of observing an interval less than t. Now, here's the
% potentially confusing part. What is p(t,y)? It's the probability that an
% event occurs at t AND it has not occurred yet. Well, that's exacly what
% an interval distribution describes. So p(x,t) is just f(t). Substituting
% f(t) and F(t) into the expresssion above, we get
%
%    h(t) = f(t) / [1 - F(t)]


% Homework question 1.
% A sprinter is waiting for the 'go' signal to start the race. The 'go' is
% the second of two pips (i.e., clicks). Suppose the interval between the
% pips is 1 sec, on average, but is random such that any time from 0.5 to
% 1.5 sec is equally probable (this is called a uniform distribution).
% Sketch the shape of the hazard function. You don't have to derive it (but
% it's easy if you try). If the sprinter has the knowledge you have, do you
% think it could affect her start time? How could we tell?

% Homework question 2.
% What is the general relationship between the probability function (or probability
% density), f(t), and its associated cumulative probaility function, F(t)? 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 3.  Gamma-distributed interevent intervals. 
%
% In this section we're going to simulate a renewal process by a different
% method.  Imagine that you have a neuron that integrates excitatory
% post-synaptic potentials that arrive as a Poisson process.  Once "n" of
% these EPSPs have been integrated (i.e., accumulated or counted), the
% neuron fires a spike and resets its membrane voltage to its resting
% state.  This is an example of a renewal process: interevent intervals are
% independent and come from a single underlying distribution.  We haven't
% explictly stated the parametric form of the interevent interval
% distribution, but that's OK.  It'll turn out that the interevent interval
% distribution is the gamma distribution.  Lucky us!

% To start, let's say that the neuron has to receive
% exactly 5 EPSPs to fire a spike.  We'll further assume that EPSPs arrive
% as a Poisson process and we'll adjust the rate of these Poisson processes
% so that our model neuron will fire, on average, 50 spikes per second (if
% should be clear that if we don't adjust the rate of the EPSP processes
% our model neuron will fire progressively more slowly as we require that
% it integrate more and more EPSPs to  fire a spike).

nstepstothresh = 5;
spike_rate = 50;
mean_int = 1000/spike_rate;
intervals = exprnd(meanint/nstepstothresh,nstepstothresh,1000);
intervals = sum(intervals,1);

% Let's take a look at the histogram of interspike intervals and the
% spike rasters from the first 500 ms of the simulation.
figure(2), clf
plot2ras(2,intervals,nstepstothresh);

% Notice that we don't get many very short interspike intervals.
% This is because our neuron has to wait for the 5th EPSP to fire
% a spike.  In principle, 5 EPSPs could arrive immediately after our
% neuron fires a spike, in which case it would fire a second spike,
% but this is exceedingly unlikely to happen.

% For comparison, let's simulate spikes from a neuron that only needs
% to integrate 1 EPSP before it fires a spike.  It should come as 
% no surprise that this neuron fires as a Poisson process -- it waits 
% for an EPSP (which arrives as a Poisson process by assumption) and
% then fires a spike.

nstepstothresh = 1
intervals = exprnd(meanint/nstepstothresh,nstepstothresh,1000);
intervals = sum(intervals,1);
plot2ras(1,intervals,nstepstothresh);

% By now the fact that the interspike interval histogram looks 
% exponential should come as no surprise.

% Notice that we occasionally get very short (and very 
% long) interspike intervals.  It's a lot more likely for pairs
% of Poisson events to occur nearly simultaneously than for a
% packet of five Poisson events to occur nearly simultaneously.
% Likewise, pairs of Poisson events may have a substantial
% delay between them, whereas the time between the fifth
% Poisson event and the 10th (assuming that the events are coming
% 5x as quickly) is unlikely to be as long.

% Finally, let's look at a neuron that needs to integrate 20 Poisson
% EPSPs before firing a spike.

nstepstothresh = 20
intervals = exprnd(meanint/nstepstothresh,nstepstothresh,1000);
intervals = sum(intervals,1);
plot2ras(3,intervals,nstepstothresh);

% Now compare vertically across the panels.  As we require that our neuron
% integrate progressively greater numbers of EPSPs, the interspike interval
% histogram gets progressively tighter (less variable) and the spike trains
% become progressively more regular looking.
% The CV of the interspike interval distribution decreases as we
% increase the number of EPSPs to count.  As was shown above, this has 
% the consequence of decreasing the mean to variance ratio of the spike
% counts.  We fixed the mean spike count in this simulation, so we see
% a decrease in the spike count variance.

% This can be a problem for the simple integrate and fire model of 
% neuronal firing.  It predicts much less variable responses than one
% actually observes in in vivo cortical recordings.  There are ways
% of getting around this, however.  For instance, the fact that real
% neurons receive both excitation *and inhibition* changes the situation
% dramatically.

% It's interesting to note that the interspike interval distributions
% for the integration model has a gamma density.  I'll superimpose some
% gamma density functions on the interspike interval histograms to try to 
% convince you of this.

subplot(3,2,1); hold on;
plot([0:60],gampdf([0:60],1,meanint),'g-');
subplot(3,2,3); hold on;
plot([0:60],gampdf([0:60],5,meanint/5),'g-');
subplot(3,2,5); hold on;
plot([0:60],gampdf([0:60],20,meanint/20),'g-');

% This is a nice result because it means that we could have simulated
% this process simply by drawing interspike intervals from a gamma 
% distribution. 

% Now, let's look at the CV and hazard functions for these gamma
% distributions.
%
% The 1st order gamma is the exponential. We know how that's going to come
% out.
clear ax
figure(3), clf
k = 1;
ax(k) = subplot(3,2,k); k = k+1;
plot(t,gampdf(t,1,meanint))
[m,v] = gamstat(1,meanint);
cv = sqrt(v)/m;
s = sprintf('CV = %.2f', cv);
text(50,.07,s)
ax(k) = subplot(3,2,k); k = k+1;
haz = gampdf(t,1,meanint) ./ (1 - gamcdf(t,1,meanint));
plot(t,haz)
ax(k) = subplot(3,2,k); k = k+1;
plot(t,gampdf(t,5,meanint/5))
% next 4 lines calculate & display the CV
[m,v] = gamstat(5,meanint/5);
cv = sqrt(v)/m;
s = sprintf('CV = %.2f', cv);
text(50,.07,s)
ax(k) = subplot(3,2,k); k = k+1;
haz = gampdf(t,5,meanint/5) ./ (1 - gamcdf(t,5,meanint/5));
plot(t,haz)
ax(k) = subplot(3,2,k); k = k+1;
plot(t,gampdf(t,20,meanint/20))
[m,v] = gamstat(20,meanint/20);
cv = sqrt(v)/m;
s = sprintf('CV = %.2f', cv);
text(50,.07,s)
ax(k) = subplot(3,2,k); k = k+1;
haz = gampdf(t,20,meanint/20) ./ (1 - gamcdf(t,20,meanint/20));
plot(t,haz)
set(ax,'XLim',[0 100])
set(ax(1:2:end),'YLim',[0 .1])
set(ax(2:2:end),'YLim',[0 .8])
set(gcf,'CurrentAxes',ax(1))
title('Interval distribution (PDF)')
set(gcf,'CurrentAxes',ax(2))
title('Hazard rate')

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II. Entropy. 
%
% We have talked about the disorderliness of the spike trains -- some more
% disorderly than others -- in terms of the CV of the interval
% distribution. We can think of the Poisson point process (PPP) as the
% benchmark.  For the PPP, the intervals were described by an exponential
% distribution which has CV=1. Recall that our lognormal distribution was
% more disorderly in the sense that its CV > 1. The gamma distributions
% were less disorderly as reflected by their CV < 1. On the other hand, the
% hazard functions seemed to imply that the PPP was the most random in the
% sense that no one interval was preferred and nothing was learned about
% the possibility of an event in the next moment simply by virtue of the
% fact that it had not happened yet. That is a unique property of the
% exponential distribution.

% Here is another way to think about disorderliness. It is the entropy. The
% entropy is defined as the average value of -log base 2 of the probabilities.
% This is sum(-log2(p) .* p). 
fEnt = @(pdf) nansum(-log2(pdf).*pdf)  
% fEnt is a function handle. I use nansum instead of sum only because
% matlab does not know that the when p=0, p*log(p) should be 0. It returns
% NaN (not a number); nansum just ignores these terms.
%
% The entropy quantity has units 'bits' and is related
% intuitively to information in the following sense. Suppose you have a
% fair coin. The coin can take values 0 or 1. The uncertainty in the coin
% is 
fEnt([.5 .5])
% which is 1 bit. If I tell you that the coin came up heads, we have
% removed all uncertainty. That is, we have learned 1 bit of information.
% That is, we have learned 1 bit of information.
%
% Suppose the coin started out unfairly weighted toward heads with 3:1
% odds. Then the uncertainty is only
fEnt([.75 .25])

% If you learn that the coin came up either heads or tails you have learned
% only .81 bits. It may seem surprising at first that you learn just as
% much from a head or a tail. But remember, the entropy is an expectation.
% Although we may  be surprised by an unlikley event, it happens less
% often. 
% 
% A nice feature about this entropy measure is that it is additive. Suppose
% there are 8 possible outcomes, like the sequences of heads and tails from
% 3 flips of a coin. Then there are 3 bits of uncertainty
p  = repmat(1/8, 1, 8)
fEnt(p)
% And suppose I tell you that the first flip came out heads. That
% eliminates 1/2 of the possibilities. 
p(1:4) = 0
% The new amount of entropy is now 2 bits
fEnt(p)
nansum(-log2(p).*p)
% (Note: I use nansum because matlab does not realize that the limit as p->0 of p*log(p) equals 0.)

% Suppose I tell you that at least one of the tosses came up heads. How
% much uncertainty reduction did you achieve? Convince yourself that 
% the answer is approximately 0.19 bits.
p  = repmat(1/8, 1, 8); p(1) = 0; p = p/sum(p)
uncertaintyReduction = 3 - fEnt(p)


% Let's return to our stochastic point processes. In a sense, the entropy
% is a measure of the disorderliness of the random variable. Let's compare
% entropy vals for the distributions we considered earlier
dt = .01;
t = [dt:dt:1000];   % define 
% Exponential
p = exppdf(t,50);
p = p/sum(p);
entExp = fEnt(p)  
% the value you get depends on time grid (i.e., 1/dt). You'll get an extra
% bit if you double the sample points (dt = .005). We'll keep dt the same
% for our comparisons.

% Fifth order Gamma
p = gampdf(t,5,50/5);
entGam5 = fEnt(p/sum(p))

% 20th order Gamma
p = gampdf(t,20,50/20);
entGam20 = fEnt(p/sum(p))

% lognormal. This is the most interesting case! Remember the CV suggests
% that this intervals are more disorderly, but the hazard says that there
% is more predictability. The entropy measure says that it is in fact less
% disorderly than the exponential!
p = lognpdf(t,3.412,1);
entLogn = fEnt(p/sum(p))

% Display comparison. The exponential wins!
fprintf(1,'Dist:\tExp\tGam5\tGam20\tLogNormal\n');
fprintf(1,'     \t~~~\t~~~~\t~~~~~\t~~~~~~~~~\n');
fprintf(1,'     \t%.2f\t%.2f\t%.2f\t%.2f\n',entExp,entGam5,entGam20,entLogn);


% Homework question 3 (optional). Write a differential equation to capture
% the idea that the hazard function equals a constant, k. Let F be the
% cumuluative distribution function. Then F' = dy/dt is the PDF for the
% intervals, what we usually denote with lower case f(t). Show that F' must
% be an exponential distribution.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part III: Correlation in LIF neurons 
%   (Exercises created by Pamela Baker)
%
%
% In this section, you will be looking at how correlation of input to a pair of 
% leaky integrate and fire cells is translated into correlation of the output 
% spike counts. You will explore this by reproducing some results from the paper:
%
% De la Rocha J., Doiron B., Shea-Brown E., Kresimir J. and Reyes A. (2007) 
% Correlation between neural spike trains increases with firing rate. Nature 
% 448:802-806.
%
% The main results you will try to replicate are: 1) that spike count correlation
% r_SC is lower than the input correlation c (Figure 2A inset); and 2) r_SC 
% depends on output firing rate (Figure 2C). To do this, you need to simulate pairs 
% of LIF neurons that are driven with correlated inputs. 
%
% In the de la Rocha et al. paper, they injected currents to drive a pair of leaky 
% integrate and fire neurons. The currents they injected, I_i (for each neuron 
% i=1,2), have the form:
%
% I_i = mu + sigma * sqrt(tau_m) * ( sqrt(1-c)*I_ind(t) + sqrt(c)*I_common(t) )
%
% where c is the level of correlation of the input, I_ind is independent current to 
% each cell and I_common is shared input between the two simulated neurons. 
% Parameters mu and sigma set the mean and variance of the total input current 
% respectively, and are varied to obtain different output firing rates.
% The main difference between the simulations you will perform and the paper is that 
% instead of creating white noise currents for I_ind and I_common, the LIF cell you 
% will simulate, using corr_leaky_if.m, takes as input excitatory and inhibitory 
% conductances that are computed from randomly generated Poisson spike trains. This 
% means you will need to use different values of mu and sigma than the ones de la 
% Rocha et al. use in their paper (since their mu and sigma have are in mV). If you 
% prefer, you can also reproduce the paper results directly by creating white noise 
% currents as they have done, and then using their parameter values. 
%
% 1) First run set_const.m to set constant parameters for the whole simulation.
% Then use gen_syn_template.m to generate synaptic inputs for the simulation:
%

set_const
syn_exc = gen_syn_template(0.0001,0.004,(0.02*g_leak));
syn_inh = gen_syn_template(0.0001,0.008,(0.06*g_leak));

%
% 2) Tuning the model -- Use make_trials.m to generate inputs for ntrials, and 
% corr_leaky_if.m to run a simulation a single cell, test a range of different mu 
% and sigma values until you have found appropriate ranges of these parameters to 
% produce a mean output spike rate over the range that de la Rocha et al. tested, 
% around 1 - 150 spikes/s (Note: in the Figure 2c inset, they show the ranges of mu 
% and sigma they used to get the two curves. For the blue curve they set mu to a low
% value and vary sigma from low to high values. For the red curve, they set sigma 
% low and vary mu. Since you will reproduce these curves, try varying mu and sigma 
% in the same way to find suitable ranges for these parameters). For this tuning 
% stage of the model, you can test just one or two trials per set of 
% (mu, sigma) values. You can set c to any value between 0 and 1 for these runs.
%

c = 0.5;
mu = 5;
sigma = 20;
ntrials = 1;
[g_exc1,g_inh1,g_exc2,g_inh2,g_excall,g_inhall] = make_trials(ntrials,syn_exc,syn_inh);
[v_m, spiketimes] = corr_leaky_if(mu, sigma, c, g_exc1, g_inh1, g_excall, g_inhall);
spikecounts = length(find(spiketimes>0.5&spiketimes<3.5));
spikerates = spikecounts/3.0

%
% 3) Running a simulation of paired cells -- run_ntrials.m runs multiple 
% trials for two paired cells and calculates the resulting spike count correlation 
% r_SC and p-value for the correlation (using the MATLAB corr function) as well as 
% the mean output spike rate, v. This function takes as input the number of trials 
% you wish to run, as well as mu values (can take a vector for several mu values at 
% once), sigma (again, will take as input single or multiple values), input 
% correlation value c and the input conductances you can generate with make_trials 
% (make sure you use the same number of trials when you run make_trials and 
% run_ntrials). Note: for the matrices returned by run_ntrials, e.g. r_sc(1,1) will 
% be the value for the first mu and first sigma values in the vectors you pass in, 
% and so on.
%

mu = [5 10];
sigma = [10 20];
c = 0.4;
ntrials = 200;
[g_exc1,g_inh1,g_exc2,g_inh2,g_excall,g_inhall] = make_trials(ntrials,syn_exc,syn_inh);
[r_sc, v, p_sc] = run_ntrials(ntrials,mu,sigma,c,g_exc1,g_inh1,g_exc2,g_inh2,g_excall,g_inhall);


%
% Homework 4A: Plot the spike count correlation r_SC for input c values between 0 
% and 1 (as in the inset to Figure 2A from the paper). Use a value for mu and sigma 
% at the low end of the ranges you chose. Run enough trials to get a significant 
% correlation value (keeping in mind for low c values, this may not happen, so at 
% most try maybe 200 trials). Hand in this plot. What does this plot look like if 
% you use high mu and sigma values and why?
%
% Homework 4B: Generate plots of r_SC vs. output spike rate v (Figure 2C in de la 
% Rocha et al.). For this you will need to run the paired simulation run_ntrials at 
% 1) a low mu value with a range of sigmas and 2) a low sigma value with a range of 
% mus. Do this for a couple of different c values (note: experimentally observed 
% values fall between 0-0.6). Hand in these plots. Also make up some plots where 
% you use slightly higher fixed sigma or mu values while varying the other parameter
%  - what happens? Is the dependence of r_SC on v a robust effect?
%
