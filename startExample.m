%% Generalized linear model fit to calcium signals
% includes example test case
%
% Alex Kwan, 2018
%
% AK wrote the code based on the procedures described in
% Pinto and Dan, Neuron 2015

clear all;

% set up figure plotting
setup_figprop;

%% Generate mock data

% --- create mock behavior time stamps
numTrial = 100; %number of trials
trialDur = 10;  %each trial is 10-second long

dt = 0.01;
t = [0:dt:numTrial*trialDur]';

stimTime = [trialDur+1:trialDur:t(end)]'; %stimulus appears 1 s after start of trial (except the first trial)
choiceTime = stimTime + 0.5 + 2*rand(size(stimTime));      %animal responds with variable delay

% --- create mock calcium signals
dff=zeros(length(t),1);

% each stimulus would elicit a synaptic-potential-like dF/F
t_kernel=[-5:dt:5]';

tau = 0.5; ampl = 1;
kernel_stim=ampl.*t_kernel./tau.*exp(-t_kernel/tau);
kernel_stim(1:ceil(end/2))=0;
for j=1:length(stimTime)
    dff=dff+[zeros(round((stimTime(j)-t(1)+t_kernel(1))/dt),1); kernel_stim; zeros(length(dff)-length(zeros(1,round((stimTime(j)-t(1)+t_kernel(1))/dt)))-length(kernel_stim),1)];
end

% each choice would elicit a synaptic-potential-like dF/F
% (at half the initial amplitude, and twice the time constant)
tau = 1; ampl = -0.5;
kernel_choice=ampl.*t_kernel./tau.*exp(-t_kernel/tau);
kernel_choice(1:ceil(end/2))=0;
for j=1:length(choiceTime)
    dff=dff+[zeros(round((choiceTime(j)-t(1)+t_kernel(1))/dt),1); kernel_choice; zeros(length(dff)-length(zeros(1,round((choiceTime(j)-t(1)+t_kernel(1))/dt)))-length(kernel_choice),1)];
end

% add noise
dff = dff + 0.05*randn(size(dff));

% --- plot the mock data, first 100 seconds
figure;
subplot(4,1,1); hold on;
for j=1:length(stimTime)
    h1=plot(stimTime(j)*[1 1],[0 1],'k');
end
for j=1:length(choiceTime)
    h2=plot(choiceTime(j)*[1 1],[0 1],'r');
end
ylabel('Stim and Choice');
xlim([t(1) 100]);
legend([h1;h2],[{'Stim'};{'Choice'}]);

subplot(4,1,2); hold on;
plot(t,dff,'k-','LineWidth',2);
xlabel('Time (s) [first 100 s]');
ylabel('dF/F');
xlim([t(1) 100]);

print(gcf,'-djpeg','mockdata-raw');    %jpeg format

%% Compare ideal kernel with the trial-averaged signals

% Trial averaging at the event times
window=[-3 3];
[sigbyTrial_stim, tbyTrial_stim] = align_signal(t,dff,stimTime,window);
[sigbyTrial_choice, tbyTrial_choice] = align_signal(t,dff,choiceTime,window);

figure;
subplot(2,1,1); hold on;
plot(tbyTrial_stim,nanmean(sigbyTrial_stim,2),'r');
plot(t_kernel,kernel_stim,'k--');
xlabel('Time from stim (s)');
legend('Estimate using trial-averaging','Kernel used to create the data');
xlim(window); ylim([-0.3 0.6]);
subplot(2,1,2); hold on;
plot(tbyTrial_choice,nanmean(sigbyTrial_choice,2),'b');
plot(t_kernel,kernel_choice,'k--');
xlabel('Time from choice (s)');
legend('Estimate using trial-averaging','Kernel used to create the data');
xlim(window); ylim([-0.3 0.6]);

print(gcf,'-djpeg','mockdata-trialavg');    %jpeg format

%% Generalized linear regression

params=[];

% --- define the predictors

% the task events that will be used to construct GLM to predict neural activity
%predictor 1: stim
params.eventTime{1} = stimTime;
params.eventLabel{1} = 'Stimulus';
%predictor 2: choice
params.eventTime{2} = choiceTime;
params.eventLabel{2} = 'Choice';

%how far do the predictors go backward and forward in time?
params.window = window;

% --- define the time periods during which the predictors will be used to
% predict the calcium signals (for real data, will ignore periods when the
% animal was inactive and not performing any actions)

%only perform the analysis between intervals define by these times
params.trigTime(:,1) = stimTime - 3;
params.trigTime(:,2) = stimTime + 3;

% --- fit GLM and plot
[glreg, glreg_fit]=gl_regr( dff, t, params.trigTime, params.eventTime, params );

%plot_gl_regr_fit(glreg_fit,glreg,params);

plot_gl_regr_beta(glreg,params);

% --- also plot the kernel on the same plots, for comparison
subplot(2,1,1); hold on;
plot(t_kernel,kernel_stim,'k--');
legend('Estimate using GLM','Kernel used to create the data');
xlim(window); ylim([-0.3 0.6]);
subplot(2,1,2); hold on;
plot(t_kernel,kernel_choice,'k--');
legend('Estimate using GLM','Kernel used to create the data');
xlim(window); ylim([-0.3 0.6]);

print(gcf,'-djpeg','mockdata-glmfit');    %jpeg format