%doing a specific survival analysis in ALS

dataTable=readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/ALS_all_Clinical_data.xlsx');


ptID=categorical(dataTable.INDDID);
visitDate=dataTable.VisitDate;
Weight=dataTable.Weight;
FVCSeated=dataTable.FVCSeated;
CNSLS=dataTable.CNSLS;
FRS=dataTable.FRSTotal;
Co2=dataTable.TranscutaneousCO2;
snp=dataTable.rs199347;
mutationSum=dataTable.Mutation_Summary;
ageAtOnset=dataTable.AgeatOnset;
onsetDate=dataTable.ALSSymptomOnsetDate;
biPapDate=dataTable.BiPAPDate;
pegDate=dataTable.PEGDate;
trachDate=dataTable.TracheostomyDate;
deathDate=dataTable.DOD;
ageAtDeath=dataTable.AgeatDeath;
bioSex=dataTable.Sex;
testDate=dataTable.VisitDate;
birthDay=dataTable.DOB;
initDiff=testDate-birthDay;

hours2use=hours(initDiff); 

ageAtTest=hours2use/(24*365.25);




emptySNP=cellfun(@isempty,snp);

%% cox survival analysis

% Get unique patient IDs
uniquePatients = unique(dataTable.INDDID);

sumMat=[];

% Loop through each patient
for gg = 1:length(uniquePatients)
    patientID = uniquePatients(gg);
    
    % Extract all rows for this patient
    patientData = dataTable(dataTable.INDDID == patientID, :);
    
    % extract SNP
    snpStatus = unique(patientData.rs199347);
    
if isempty(snpStatus{1})
        snpStatAtPlay=nan;
    elseif snpStatus{1}=='CC'
    snpStatAtPlay=0;
    elseif snpStatus{1}=='CT'
    snpStatAtPlay=1;
    elseif snpStatus{1}=='TT'
    snpStatAtPlay=2;
end


    % Take the first available row for date values (assuming they are the same)
    trachDate = unique(patientData.TracheostomyDate);
    pegDate = unique(patientData.PEGDate);
    biPapDate = unique(patientData.BiPAPDate);
    deathDate = unique(patientData.DOD);
    onsetDate=unique(patientData.ALSSymptomOnsetDate);
    diagDate=unique(patientData.DiagnosisDate);
    ageAtOnset= unique(patientData.AgeatOnset);

if isnan(ageAtOnset(1))
    ageAtOnset=nan;
end


    if ~isnat(onsetDate(1))
        date2use=onsetDate;
    elseif isnat(onsetDate(1)) && ~isnat(diagDate(1))
    date2use=diagDate;
    else
        date2use=NaT;
    end



    % Find the earliest date among the available ones
    allDates = [trachDate;  deathDate];
    validDates = allDates(~isnat(allDates)); % Remove NaT (missing values)
    if isempty(validDates)
        earliestEndDate = NaT;
    else
        earliestEndDate = min(validDates);
    end


    timeToEvent= earliestEndDate-date2use;
%convert to months 
hours2use=hours(timeToEvent)/24; weeks=hours2use/7; timeToEvent=weeks/4;

    
    % Store results in summary matrix
sumMat=[sumMat; [patientID,snpStatAtPlay,timeToEvent,ageAtOnset] ];






end


remValz=isnan( sumMat(:,2) + sumMat(:,3) ) ;

scoredValz=zscore(sumMat(~remValz,3));

snpAtPlay=sumMat(~remValz,2);
eventAtPlay=sumMat(~remValz,3);

scoredValz=zscore(sumMat(~remValz,3));


[b,logL,H,stats] = coxphfit(sumMat(scoredValz>1 & ~remValz,2 ) ,sumMat(scoredValz>1  & ~remValz )  );


%% running each SNP separately 


snpGroups = unique(sumMat(~remValz,2)); % Extract SNP categories

figure;
hold on;
colors = {'c', 'm', 'g'}; % Assign colors to SNP groups

for gg = 1:length(snpGroups)
    snpValue = snpGroups(gg);
    
    % Filter data for the current SNP group

    idx = sumMat(:,2) == snpValue  & ~remValz;
    survTimes = sumMat(idx,3); % Survival times
    scoredAtPlay=zscore(survTimes);
    eventOccurred = ones(size(survTimes  )); % Assume all events are observed (if censoring not available)

    % Fit Cox model only for this SNP group
    [~, ~, H, ~] = coxphfit(ones(size(survTimes )), survTimes , 'Censoring', 1-eventOccurred); 

    % Plot survival curve
    stairs(H(:,1), exp(-H(:,2)), 'Color', colors{gg}, 'LineWidth', 2);
end

xlabel('Time (Months)','FontSize',20);
ylabel('Survival Probability','FontSize',20);
title(['Survival Probability Stratified by rs199347 status','p=', num2str(stats.p)],'FontSize',20);
legend({'CC','CT','TT'});
grid on;
hold off;
xlim([0,200])

%%

scale_param = mean(H(:,1)) / log(2); % Approximate median survival
B = 1.5; % Shape parameter 

% Define time points
xx = linspace(0, max(H(:,1)), 100);

% Compute Weibull survival function
weibull_survival = exp(-(xx / scale_param).^B); 

% Plot Cox vs Weibull
figure;
hold on;
stairs(H(:,1), exp(-H(:,2)), 'LineWidth', 2, 'DisplayName', 'Cox Estimated Survivor Function');
plot(xx, weibull_survival, 'r', 'LineWidth', 2, 'DisplayName', 'Weibull Survivor Function');

xlim([0, max(H(:,1))]);
xlabel('Time');
ylabel('Survival Probability');
title('Cox vs Weibull Survivor Function');
legend;
grid on;
hold off;
xlim([0,200])


%% ok cool now onto modelling decline



ptID=(dataTable.INDDID);
visitDate=dataTable.VisitDate;
Weight=dataTable.Weight;
FVCSeated=dataTable.FVCSeated;
CNSLS=dataTable.CNSLS;
FRS=dataTable.FRSTotal;
ageAtDeath=dataTable.AgeatDeath;
bioSex=dataTable.Sex;
snp=dataTable.rs199347;
CNSLS=dataTable.CNSLS;
Co2=dataTable.TranscutaneousCO2;
geneStat=dataTable.Mutation_Summary;
c9Pos=contains(geneStat,'C9orf72');
testDate=dataTable.VisitDate;
birthDay=dataTable.DOB;
initDiff=testDate-birthDay;
emptySNP=cellfun(@isempty,snp);
trachDate=dataTable.TracheostomyDate;




hours2use=hours(initDiff); 

ageAtTest=hours2use/(24*365.25);


% model will ultimately be scores ~ 1/ID + time*gene + sex + ageAtDeath
% +startScore


%1 calculate time passed and get a Mat with the start scores for our model
uniquePt=unique(ptID);
scoreToTest=FRS;

startScoreMat=[];
timePassedMat=[];
numVisits=[];
numVisitInd=[];
for dd=1:length(uniquePt)
    datesAtPlay=visitDate(ptID==uniquePt(dd));
    scoresAtPlay=scoreToTest(ptID==uniquePt(dd));
        if isempty(scoresAtPlay) || isscalar(scoresAtPlay) 
            scoresAtPlay=nan;
        end
     [startDt,DtLoc]=min(datesAtPlay) ;  %finding earliest date (sometimes I think things come out of order)
     startScore=scoresAtPlay(DtLoc);
    timePassed=datesAtPlay-startDt; %order doesn't matter here
    %converting to months
    hours2use=hours(timePassed)/24; weeks=hours2use/7; timePassed=weeks/4;
    if timePassed==0
    timePassed=nan;
    end

    %generating our output Mats

    timePassedMat=[timePassedMat;timePassed];
    scoreHold=ones(length(timePassed),1); scoreHold=scoreHold* startScore;
    startScoreMat=[startScoreMat;scoreHold];
   visitHold=ones(length(timePassed),1); visitHold=visitHold* length(scoresAtPlay);
    numVisits=[numVisits;visitHold];
    numVisitInd=[numVisitInd,length(scoresAtPlay)];
end


filter2use= ptID~=112690 & ptID~=104346 &ptID~=110052 %~isnan(ageAtDeath) %~isnan(ageAtDeath) | ~isnat(trachDate);


% Create table
modelData = table(ptID(filter2use), timePassedMat(filter2use), bioSex(filter2use), ageAtTest(filter2use),   snp(filter2use), scoreToTest(filter2use), startScoreMat(filter2use), ageAtDeath(filter2use), ...
    'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest',  'SNP', 'Scores', 'startScore', 'ageAtDeath'});

% Convert categorical variables to factors
modelData.Sex = categorical(modelData.Sex);
modelData.ptID = categorical(string(modelData.ptID));
modelData.SNP=categorical(modelData.SNP); %I'm not entirely sure if this is correct but I'm going categorical here because I'm not convinced a relationship would be purely linear



% Fit GLME 
glme = fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtTest + SNP * timePassed + (1|ptID)', ...
    'Distribution', 'Normal', 'Link', 'Identity');

glmeAlt= fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtTest + (1|ptID)', ...
    'Distribution', 'Normal', 'Link', 'Identity');

compare(glmeAlt, glme) % cool, significant in terms of model comp as well

disp(glme);



% glme = fitglme(modelData, ...
%     'ageAtDeath ~ 1 + Sex + startScore  + SNP * timePassed + (1|ptID)', ...
%     'Distribution', 'Normal', 'Link', 'Identity');
% 

%%

filter2use= isnan(ageAtDeath); %~isnan(ageAtDeath) | ~isnat(trachDate);


% Create table
modelDataTest = table(ptID(filter2use), timePassedMat(filter2use), bioSex(filter2use), ageAtTest(filter2use),   snp(filter2use), scoreToTest(filter2use), startScoreMat(filter2use), ageAtDeath(filter2use), numVisits(filter2use), trachDate(filter2use), ...
    'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest',  'SNP', 'Scores', 'startScore', 'ageAtDeath','numVisits', 'TrachDate'});

%% ok coolio but let's be very hardcore w/ a permutation test on this glme

permNumber=1000;

pValTrue=glme.Coefficients(9,6);
pValTrue=pValTrue.pValue;
falseScore= scoreToTest(  (randperm(length(scoreToTest))   ) );    

permMat=nan(1,permNumber);
    for dd=1:permNumber
falseScore= scoreToTest(  (randperm(length(scoreToTest))   ) );    
modelData.Scores=falseScore;
glmeFalse = fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtDeath + SNP * timePassed + (1|ptID)', ...
    'Distribution', 'Normal', 'Link', 'Identity');
falsePval=glmeFalse.Coefficients(9,6);
permMat(dd)=falsePval.pValue;
    end

permP=sum(permMat< pValTrue)/permNumber;


%% plotting

%checking out model fit

% Get residuals from the GLME model
residuals2use = residuals(glme, 'ResidualType', 'Pearson'); % Pearson residuals

% Plot residuals against fitted values
fittedVals = fitted(glme);

figure;
scatter(fittedVals, residuals2use, 'o');
xlabel('Fitted Values','FontSize',25);
ylabel('Residuals','FontSize',25);
title('Residual Plot','FontSize',25);
refline(0,0); % Add ref line

% looks good to me! There is slight heteroscedasticity with lower variance
% at higher values but this makes intuitive sense to me since there is more
% intrinsic variability when people are getting sicker... 


% plotting the output of the likelihood ratio test w/ log likelihoods

% Calculate explained variance (using log-likelihood)
logLikelihoodFull = glme.LogLikelihood; 
logLikelihoodAlt = glmeAlt.LogLikelihood;   
% Plotting the bar graph
figure;
bar([logLikelihoodFull, logLikelihoodAlt]);
set(gca, 'XTickLabel', {'Full Model', 'Alt Model'},'FontSize',15);
ylabel('Explained Variance (Log-Likelihood)', 'FontSize',20);
title('Model Comparison:  Log-Likelihood w/ and w/o SNP interaction term', 'FontSize',15);
grid on;



%plotting fits for decline

% Get fitted values from the GLME
fittedVals = fitted(glme);

% Get unique SNP groups
uniqueSNPs = unique(modelData.SNP);
uniqueSNPs=uniqueSNPs(1:3);
% Define time bins
timeBins = 0:20:120;

% Initialize arrays for mean fitted values and standard errors
meanFittedOverTime = nan(length(timeBins)-1, length(uniqueSNPs));
stdErrorOverTime = nan(length(timeBins)-1, length(uniqueSNPs));

% Calculate mean fitted scores and standard error in each time bin
for i = 1:length(uniqueSNPs)
    for j = 1:length(timeBins)-1
        idx = modelData.SNP == uniqueSNPs(i) & modelData.timePassed >= timeBins(j) & modelData.timePassed < timeBins(j+1);
        
        if any(idx)
            meanFittedOverTime(j, i) = nanmean(fittedVals(idx));
            stdErrorOverTime(j, i) = nanstd(fittedVals(idx)) / sqrt(sum(idx)); % Standard error
        end
    end
end

% Plot the fitted values with error bars
figure;
hold on;
for i = 1:length(uniqueSNPs)
    errorbar(timeBins(1:end-1) + 5, meanFittedOverTime(:, i), stdErrorOverTime(:, i), '-o', 'LineWidth', 2, 'DisplayName', char(uniqueSNPs(i)));
end
xlabel('Time Passed (Months)','FontSize',20);
ylabel('Fitted Weight','FontSize',20);
title('Fitted Weight Over Time by Genotype','FontSize',20);
legend('Location', 'NorthEast','FontSize',30);
grid on;
hold off;

%xlim([0,120])




%% kaplan meir


% Extract SNP groups
snpGroups = unique(sumMat(~remValz,2)); % Extract SNP categories

figure;
hold on;
colors = {'c', 'm', 'g'}; % Colors for SNP groups
legendLabels = {}; % Store legend labels

for gg = 1:length(snpGroups)
    snpValue = snpGroups(gg);
    
    % Filter data for the current SNP group
    idx = sumMat(:,2) == snpValue & ~remValz;
    survTimes = sumMat(idx,3); % Survival times
    eventOccurred = ones(size(survTimes)); % Assume no censoring for now

    % Kaplan-Meier Estimator
    [f, x] = ecdf(survTimes, 'Function', 'survivor'); % Kaplan-Meier survival function

    % Plot survival curve
    plot(x, f, 'Color', colors{gg}, 'LineWidth', 2);
    legendLabels{end+1} = ['SNP ', num2str(snpValue)];
end

xlabel('Time (Months)','FontSize',20);
ylabel('Survival Probability','FontSize',20);
title('Kaplan-Meier Survival Curves by SNP Status','FontSize',20);
legend(legendLabels);
grid on;
hold off;
xlim([0,200]);
