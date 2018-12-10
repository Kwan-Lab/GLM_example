function [ output, output_fit ] = gl_regr ( signal, t, trigTime, eventTime, params)
% % linear_regr %
%PURPOSE:   Generalized linear regression
%           see Pinto and Dan, Neuron 2015
%AUTHORS:   AC Kwan 171101
%
%INPUT ARGUMENTS
%   signal:         time-series signal (e.g., dF/F)
%   t:              time points corresponding to the signal
%   trigTime:       bounds of periods that will be used to train/test GLM
%   eventTime:      times of events that are predictive of the signal
%   params.window:  temporal window for predictors
%
%OUTPUT ARGUMENTS
%   output:         structure containing regression analysis results
%   output_fit:     structure containing specifically info regarding signal
%                   fit for this cell (large in size, and not useful to 
%                   summarize for population)
%
% To plot the output, use plot_gl_regr().
%
% Note -- when to use RIDGE REGRESSION: when there are many predictors, and
% the predictors are correlated. Problem with standard linear regression
% is that it is based on least-square errors, which gives an unbiased
% estimate of the parameters (beta's). However if predictors are corrected,
% then variability in the extracted beta's could be large. Instead, with
% ridge regression, estimation of the predictors are biased and typically
% under-estimated, but the problem is more well-conditioned and the
% extracted parameters are more consistent, yielding "better" prediction.
% (https://tamino.wordpress.com/2011/02/12/ridge-regression/)

%% the predictors

%the 'basis set' of the events
for k = 1:numel(eventTime)  %for each event type
    basisEvent{k} = zeros(size(t));  %default value for each time bin is 0
    
    basisEvent{k}(1) = 0;
    for j = 2:numel(t)   %if event occurs in that time bin, puts a 1 there
        basisEvent{k}(j) = sum(eventTime{k}>t(j-1) & eventTime{k}<=t(j));
    end
end

%to construct the full predictor set, we will translate basis set by going
%backward and forward in time (to ask if there is an event occurred prior
%to or after a particular image frame); total number of predictors to
%regress against is #event * #time_shift (e.g., if events are bout onset
%and bout offset, and time shift is +/- 2 seconds at 30Hz, then there are
%242 predictors)
frameRate = 1/mode(diff(t));
idxWinStart=round(frameRate * params.window(1));     %how much to translate, in units of image frames
idxWinLength=round(frameRate * (params.window(end) - params.window(1)));
idx=[idxWinStart:1:idxWinStart+idxWinLength];
t_Window = idx/frameRate;

X = []; %matrix of all the predictors
mask = [];
for k = 1:numel(eventTime)  %for each event type
    predictor = zeros(numel(t),numel(idx));
    
    for j=1:numel(idx)
        if idx(j)<0       %neural activity precedes event
            predictor(:,j)=[basisEvent{k}(-idx(j)+1:end); zeros(-idx(j),1)];
        elseif idx(j)>=0  %neural activity follows event
            predictor(:,j)=[zeros(idx(j),1); basisEvent{k}(1:end-idx(j))];
        end
    end
    
    X = [X predictor];
    
    %later, we will get a 1-dimensional arrays for the beta for each
    %predictor. We will want to separate the predictors by the event type
    %let's say idx is 3 element long, and there are 2 event types;
    %mask{1} = [1 1 1 0 0 0]; mask{2} = [0 0 0 1 1 1]
    mask{k} = [false((k-1)*numel(idx),1); true(numel(idx),1); false((numel(eventTime)-k)*numel(idx),1)];
end

%% training and testing sets
% use 20% data for testing; remainder 80%, use 80% for training, 20% for validation

numTrial = size(trigTime,1);

%which trials are picked for test/train
numTest = round(0.2*numTrial);
drawNum = randsample(numTrial,numTrial,'false'); %each time draw another set without replacement
testTrials = drawNum(1:numTest);
trainTrials = drawNum(numTest+1:end);

%take out the first and last events, so don't have to deal with boundary
%issues when identifying imaging frames
testTrials = testTrials(testTrials~=1 & testTrials~=numTrial);
trainTrials = trainTrials(trainTrials~=1 & trainTrials~=numTrial);

%what imaging frames do those test/train/validation trials correspond to
testIdx=false(size(signal,1),1);
for j=1:numel(testTrials)
    idx1=sum(trigTime(testTrials(j),1)>=t);
    idx2=sum(trigTime(testTrials(j),2)>=t);
    testIdx(idx1:idx2)=true;
end
trainIdx=false(size(signal,1),1);
for j=1:numel(trainTrials)
    idx1=sum(trigTime(trainTrials(j),1)>=t);
    idx2=sum(trigTime(trainTrials(j),2)>=t);
    trainIdx(idx1:idx2)=true;
end

%% run the regression

%identify the best regularization parameter lambda, using training/validation set
testLambda=10.^(-(5:-0.5:1));
testCC=[]; sqerr=[];

numRepeat = 5;

for jj = 1:numRepeat  %do 5 repeats for each lambda value

    %among the train trials, use 80% to train and 20% to test the lambda
    numTrainTrial = numel(trainTrials);
    numVal = round(0.2*numTrainTrial);
    
    drawNum = randsample(numTrainTrial,numTrainTrial,'false'); %each time draw another set without replacement

    trainvalTrials = trainTrials(drawNum(1:numVal));
    traintrainTrials = trainTrials(drawNum(numVal+1:end));

    traintrainIdx=false(1,size(signal,1));
    for j=1:numel(traintrainTrials)
        idx1=sum(trigTime(traintrainTrials(j),1)>=t);
        idx2=sum(trigTime(traintrainTrials(j),2)>=t);
        traintrainIdx(idx1:idx2)=true;
    end
    trainvalIdx=false(1,size(signal,1));
    for j=1:numel(trainvalTrials)
        idx1=sum(trigTime(trainvalTrials(j),1)>=t);
        idx2=sum(trigTime(trainvalTrials(j),2)>=t);
        trainvalIdx(idx1:idx2)=true;
    end
    
    for ll=1:numel(testLambda) %test a range of regularization parameter lambda
        %L2-penalized linear least square regressoin (ridge regression)
        MDl=fitrlinear(X(traintrainIdx,:)',signal(traintrainIdx)','ObservationsIn','columns','Learner','leastsquares','Solver','sgd','Regularization','ridge','Lambda',testLambda(ll),'FitBias',true);
        yfit=predict(MDl,X(trainvalIdx,:));
        
        %computer corr. coefficient between predicted and measured activity
        tempCC=corrcoef(signal(trainvalIdx),yfit);
        testCC(ll,jj)=tempCC(1,2);
        
        %computer squared error between predicted and measured activity
        %note -- in my hands, this seems to identify the same index more consistently than CC
%        sqerr(ll,jj) = sum((yfit - signal(trainvalIdx)).^2);
    end
    
end

%identify the fit that has the highest mean correlation coefficient
[~,idx] = max(mean(testCC,2));

%identify the fit that has the smallest mean least-square error
%[~,idx]=min(mean(sqerr,2));

%using that specific lambda parameter, fit again using training +
%validation set
MDl=fitrlinear(X(trainIdx,:)',signal(trainIdx)','ObservationsIn','columns','Learner','leastsquares','Solver','sgd','Regularization','ridge','Lambda',testLambda(idx),'FitBias',true);

%evaluate fit using test set
ytest_fit=predict(MDl,X(testIdx,:));

output_fit.ytest = nan(size(signal));
output_fit.ytest(testIdx) = signal(testIdx);  %the test signal, measured
output_fit.ytest_fit = nan(size(signal));
output_fit.ytest_fit(testIdx) = ytest_fit;    %the test signal, predicted

%compute correlation coefficient between predicted and measured activity for test set
tempCC=corrcoef(signal(testIdx),ytest_fit);
output.ytest_CC=tempCC(1,2);

%the GLM coefficients
output.lambda = testLambda(idx);

for k = 1:numel(eventTime)  %for each event type
    output.beta{k} = MDl.Beta(mask{k});
end
output.t_beta = t_Window';  %the time associated with the beta's
output.beta_label = params.eventLabel;

end
