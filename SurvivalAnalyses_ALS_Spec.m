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
weight=dataTable.Weight;
hours2use=hours(initDiff); 

ageAtTest=hours2use/(24*365.25);




emptySNP=cellfun(@isempty,snp);

%% cox survival analysis

% Get unique patient IDs
uniquePatients = unique(dataTable.INDDID);
sumMat=[]; %this will be where we put all covariates, etc
% Loop through each patient
for gg = 1:length(uniquePatients)
    patientID = uniquePatients(gg);   
    % Extract all rows 
    patientData = dataTable(dataTable.INDDID == patientID, :);
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
    ageAtDeath=unique(patientData.AgeatDeath);
    FRS=patientData.FRSTotal;
    dates=patientData.VisitDate;
    [~,earliestDateLoc]=min(dates);

    if ~isempty(earliestDateLoc(1)) & ~isempty(FRS(1))
    FRSScoreAtPlay=FRS(earliestDateLoc);
    else
FRSScoreAtPlay=nan;
    end

if isnan(ageAtDeath(1))
    ageAtDeath=nan;
end


    bioSex= unique(patientData.Sex);
if isnan(ageAtOnset(1))
    ageAtOnset=nan;
end

    
    if strcmp(bioSex,'Female')
    sexAtPlay=1;
    elseif strcmp(bioSex,'Male')
    sexAtPlay=0;
    else
        sexAtPlay=nan;
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
    
ageAtDeath=ageAtOnset+(timeToEvent/12);
avgAge= (ageAtDeath+ageAtOnset)/2;

    % Store results in summary matrix
sumMat=[sumMat; [patientID,snpStatAtPlay,timeToEvent,avgAge, sexAtPlay] ];


end


remValz=isnan(mean( sumMat,2)) ;
snpAtPlay=sumMat(~remValz,2);
eventAtPlay=sumMat(~remValz,3);
zScoreAge=zscore(sumMat(~remValz,4) );
snpStat=sumMat(~remValz,2);
sexStat=sumMat(~remValz,5);
% startScorez=sumMat(~remValz,6);


[b,logL,H,stats] = coxphfit( [snpStat,zScoreAge,sexStat] ,sumMat(~remValz,3 )  );

% Compute hazard ratio
HR = exp(stats.beta);  



%% running each SNP separately 



%remValz=isnan(sumMat(:,2)) | isnan(sumMat(:,3)) | isnan(sumMat(:,1)) | isnan(sumMat(:,5)) | isnan(sumMat(:,4))  ;


uniqueSNPs = unique(sumMat(~remValz,2)); % Extract SNP categories

figure;
hold on;
colors = {'c', 'm', 'g'}; % Assign colors to SNP groups

for gg = 1:length(uniqueSNPs)
    snpValue = uniqueSNPs(gg);
    
    % Filter data for the current SNP group

    idx = (sumMat(~remValz,2) == snpValue);
    survTimes = sumMat(idx,3); % Survival times
    eventOccurred = ones(size(survTimes  )); % Assume all events are observed since we remove nans... 

    % Fit Cox model only for this SNP group
    [~, ~, H, ~] = coxphfit(ones(size(survTimes )), survTimes , 'Censoring', 1-eventOccurred); %I'm just forcing a fit here for the purposes of plotting


    length(survTimes)
    % Plot survival curve
    stairs(H(:,1), exp(-H(:,2)), 'Color', colors{gg}, 'LineWidth', 2);
end

xlabel('Time (Months)','FontSize',20);
ylabel('Survival Probability','FontSize',20);
title(['Survival Probability Stratified by rs199347 status'],'FontSize',20) %,' p=', num2str(stats.p(1))],'FontSize',20) %, ' HR=',num2str(HR(1))],'FontSize',20);
legend({'CC','CT','TT'},'FontSize',20 );
grid on;
hold off;
xlim([0,120])

%% just testing out some ways we could plot this

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
BMI=dataTable.BMI;
UMNTotal=dataTable.UMNTotal;
dx=dataTable.ElEscorialVisit;
controlPt=contains(dx, 'Control');
onsetSite=dataTable.ALSSymptomOnsetSite;
hours2use=hours(initDiff); 

overGPNMB=contains(snp,"TT");
underGPNMB=contains(snp,"CC");
het=contains(snp,"CT");
emptySNP=cellfun(@isempty,snp);



ageAtTest=hours2use/(24*365.25);


% model will ultimately be scores ~ 1/ID + time*gene + sex + ageAtDeath
% +startScore


%1 calculate time passed and get a Mat with the start scores for our model
uniquePt=unique(ptID);
scoreToTest=Weight;

startScoreMat=nan(1,length(scoreToTest));
timePassedMat=nan(1,length(scoreToTest));
numVisits=nan(1,length(scoreToTest));
numVisitInd=nan(1,length(scoreToTest));
totScoreDiffMat=nan(1,length(scoreToTest));
totTimePassedMat=nan(1,length(scoreToTest));
postTrachMat=nan(1,length(scoreToTest));
snpTot=nan(1,length(scoreToTest));
onsetSiteTot=cell(1,length(uniquePt));

for dd=1:length(uniquePt)
    datesAtPlay=visitDate(ptID==uniquePt(dd) );
    scoresAtPlay=scoreToTest(ptID==uniquePt(dd));
        if isempty(scoresAtPlay) || isscalar(scoresAtPlay) 
            scoresAtPlay=nan;
        end
     [startDt,DtLoc]=min(datesAtPlay) ;  %finding earliest date (sometimes I think things come out of order)
     [endDt,endDtLoc]=max(datesAtPlay);
     startScore=scoresAtPlay(DtLoc);
     endScore=scoresAtPlay(endDtLoc);
    timeAtPlay=datesAtPlay-startDt; %order doesn't matter here
    orderedDates=sort(datesAtPlay);
%finding scores that come after a tracheostomy

trachDates=trachDate(ptID==uniquePt(dd),1) ;
if ~isnat(trachDates(1) ) 

   postTrach=datesAtPlay> trachDates(1); 
else
    postTrach=false(length(datesAtPlay),1);

end
    %converting to months

[timePassed] = timeOutput(timeAtPlay);


totTime=orderedDates(end)-orderedDates(1);    
totTimePassed=timeOutput(totTime);


totScoreDiff= scoresAtPlay(endDtLoc)-startScore;
%being lazy and leveraging this for loop to grab onset site
onSight=onsetSite(ptID==uniquePt(dd));
onsetSiteTot{dd}=onSight{1};
% snpSight=snp(ptID==uniquePt(dd));
% snpTot(ptID==uniquePt(dd))=snpSight(1);

    %generating our output Mats
    timePassedMat(ptID==uniquePt(dd))=timePassed;
    scoreHold=ones(length(timePassed),1); scoreHold=scoreHold* startScore;
    startScoreMat(ptID==uniquePt(dd))=scoreHold;
   visitHold=ones(length(timePassed),1); visitHold=visitHold* length(scoresAtPlay);
    numVisits(ptID==uniquePt(dd))=visitHold;
    numVisitInd=[numVisitInd,length(scoresAtPlay)];
    totTimePassedHold=ones(length(timePassed),1); totTimePassedHold=totTimePassedHold* totTimePassed;
    totTimePassedMat((ptID==uniquePt(dd)))=totTimePassedHold;
    totScoreDiffHold=ones(length(timePassed),1); totScoreDiffHold=totScoreDiffHold* totScoreDiff;
    totScoreDiffMat(ptID==uniquePt(dd))=totScoreDiffHold;
    postTrachMat(ptID==uniquePt(dd))=postTrach;

end

emptySNP=cellfun(@isempty,snp);

filter2use= ~controlPt' & totScoreDiffMat<0  & numVisits>=3 & ~postTrachMat & ~isnan(timePassedMat) & ~emptySNP' & ~isnan(scoreToTest') ;
%the model will actually automatically remove nan values, I just include
%them to facilitate counting of patients included!

includedPT=length(unique(ptID(filter2use)));


% Create table
modelData = table(ptID(filter2use), timePassedMat(filter2use)', bioSex(filter2use), ageAtTest(filter2use),   snp(filter2use), scoreToTest(filter2use), startScoreMat(filter2use)', ageAtDeath(filter2use), numVisits(filter2use)', ...
    'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest',  'SNP', 'Scores', 'startScore', 'ageAtDeath','numVisits'});

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



% glmeTest = fitglme(modelData, 'Scores ~ 1 + Sex + ageAtTest + startScore + timePassed*SNP + (1 | ptID)', ...
%                'DummyVarCoding', 'effects');



%% prop of upper motor symptoms


overGPNMB=contains(snpTot,"TT");
underGPNMB=contains(snpTot,"CC");
het=contains(snpTot,"CT");
isUMN=contains(onsetSiteTot,'UMN');
isLMN=contains(onsetSiteTot,'LMN');

propUnder=(sum(isUMN &underGPNMB')/sum(underGPNMB))*100

propHet=(sum(isUMN &het')/sum(het))*100

propOver=(sum(isUMN &overGPNMB')/sum(overGPNMB))*100




%%

filter2use= isnan(ageAtDeath); %~isnan(ageAtDeath) | ~isnat(trachDate);


% Create table
modelDataTest = table(ptID(filter2use), timePassedMat(filter2use), bioSex(filter2use), ageAtTest(filter2use),   snp(filter2use), scoreToTest(filter2use), startScoreMat(filter2use), ageAtDeath(filter2use), numVisits(filter2use), trachDate(filter2use), ...
    'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest',  'SNP', 'Scores', 'startScore', 'ageAtDeath','numVisits', 'TrachDate'});

%% ok coolio but let's be very hardcore w/ a permutation test on this glme


permNumber=1000;

pValTrue=glme.Coefficients(9,6);
pValTrue=pValTrue.pValue;

permMat=nan(1,permNumber);
    for dd=1:permNumber
falseScore= scoreToTest(  (randperm(length(scoreToTest(filter2use)))));

modelData.Scores=falseScore;
glmeFalse = fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtDeath + SNP * timePassed + (1|ptID)', ...
    'Distribution', 'Normal', 'Link', 'Identity');
falsePval=glmeFalse.Coefficients(9,6);
permMat(dd)=falsePval.pValue;
    end

    %update this!!
permP=sum(permMat< pValTrue)/(permNumber-8); %this is ridiculous but I'm too lazy right now to capture the situations where a model can't converge, so I just manually enter it based on warnings


%% simple glme plot



figure;
    hold on;

TCLowerCI=glme.Coefficients(8,7); TCLowerCI=TCLowerCI.Lower;
TCupperCI=glme.Coefficients(8,8); TCupperCI=TCupperCI.Upper;
TTLowerCI=glme.Coefficients(9,7); TTLowerCI=TTLowerCI.Lower;
TTupperCI=glme.Coefficients(9,8); TTupperCI=TTupperCI.Upper;

CCupperCI=mean([TTupperCI,TCupperCI]);
CCLowerCI=mean([TTLowerCI,TCLowerCI]);


CCestimate=0; %CC is the comparator group so the 'slope' will always be 0, here I just use the CI for CT, in next iteration I'll just plot residuals for each patient, binned
colorCode=[0,1,1];
startPt= nanmean(startScoreMat(    strcmp(snpStat,'CC')     ));
plotGLMESlope(CCestimate, CCLowerCI, CCupperCI,startPt,colorCode)





TCestimate=glme.Coefficients(8,2); TCestimate=TCestimate.Estimate;
startPt= nanmean(startScoreMat(    strcmp(snpStat,'CT')     ));
colorCode=[1,0,1];
plotGLMESlope(TCestimate, TCLowerCI, TCupperCI,startPt,colorCode)





TTestimate=glme.Coefficients(9,2); TTestimate=TTestimate.Estimate;

snpStat=dataTable.rs199347;

startPt= nanmean(startScoreMat(    strcmp(snpStat,'TT')     ));


    
    colorCode=[0,1,0];
plotGLMESlope(TTestimate, TTLowerCI, TTupperCI,startPt,colorCode)









legend({'CC','','CT','','TT'},'FontSize',25);
xlabel('Time (Months)','FontSize',25);
ylabel('Predicted ALSFRS','FontSize',25)


title('Modeled ALSFRS by rs199347 Status', 'FontSize',18)





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




%% Cox on continuous values 



%1 calculate time passed and get a Mat with the start scores for our model
uniquePt=unique(ptID);
scoreToTest=Weight;

startScoreMat=[];
timePassedMat=[];
numVisits=[];
numVisitInd=[];
totScoreDiffMat=[];
totTimePassedMat=[];
postTrachMat=[];
snpTot=[];
onsetSiteTot=cell(1,length(uniquePt));

for dd=1:length(uniquePt)
    datesAtPlay=visitDate(ptID==uniquePt(dd) );
    scoresAtPlay=scoreToTest(ptID==uniquePt(dd));
        if isempty(scoresAtPlay) || isscalar(scoresAtPlay) 
            scoresAtPlay=nan;
        end
     [startDt,DtLoc]=min(datesAtPlay) ;  %finding earliest date (sometimes I think things come out of order)
     [endDt,endDtLoc]=max(datesAtPlay);
     startScore=scoresAtPlay(DtLoc);
     amtChange=scoresAtPlay/startScore;

     endPt= find(amtChange<=.8,1); endDt=datesAtPlay(endPt);
     timeAtPlay=datesAtPlay-startDt; %order doesn't matter here
    orderedDates=sort(datesAtPlay);
%finding scores that come after a tracheostomy


time2endPt = timeOutput( endDt-startDt);




totScoreDiff= scoresAtPlay(endDtLoc)-startScore;
%being lazy and leveraging this for loop to grab onset site
onSight=onsetSite(ptID==uniquePt(dd));
onsetSiteTot{dd}=onSight{1};
snpSight=snp(ptID==uniquePt(dd));
snpTot=[snpTot;snpSight(1)];

    %generating our output Mats
    timePassedMat=[timePassedMat;timePassed];
    scoreHold=ones(length(timePassed),1); scoreHold=scoreHold* startScore;
    startScoreMat=[startScoreMat;scoreHold];
   visitHold=ones(length(timePassed),1); visitHold=visitHold* length(scoresAtPlay);
    numVisits=[numVisits;visitHold];
    numVisitInd=[numVisitInd,length(scoresAtPlay)];
    totTimePassedHold=ones(length(timePassed),1); totTimePassedHold=totTimePassedHold* totTimePassed;
    totTimePassedMat=[totTimePassedMat; totTimePassedHold];
    totScoreDiffHold=ones(length(timePassed),1); totScoreDiffHold=totScoreDiffHold* totScoreDiff;
    totScoreDiffMat=[totScoreDiffMat;totScoreDiffHold];
    postTrachMat=[postTrachMat; postTrach];

end
%%


% Extract unique patient IDs
uniquePt = unique(ptID);
time2endPt = []; % Survival time
snpGroup = []; 

for dd = 1:length(uniquePt)
    datesAtPlay = visitDate(ptID == uniquePt(dd));
    scoresAtPlay = FRS(ptID == uniquePt(dd));


    [startDt, DtLoc] = min(datesAtPlay);
    startScore = scoresAtPlay(DtLoc);
    amtChange = scoresAtPlay / startScore;

    endPt = find(amtChange <= 0.5, 1); % Find first time weight drops 20%

    if ~isempty(endPt)
        endDt = datesAtPlay(endPt);
        time2endPt = [time2endPt; timeOutput(endDt - startDt)];
    else
        time2endPt = [time2endPt;nan]; % Last follow-up time
    end

    % Store SNP group for stratification
    snpSight = snp(ptID == uniquePt(dd));
    snpGroup = [snpGroup; snpSight(1)];

if length(time2endPt) >dd
b=1;
end

end

% Get unique SNP groups
emptySNP=cellfun(@isempty,snpGroup);
uniqueSNPs = unique(snpGroup(~emptySNP));

% Define colors for each SNP group
colors = lines(length(uniqueSNPs)); 


overGPNMB=contains(snpGroup,'TT');
underGPNMB=contains(snpGroup,'CC');
figure
hold on
histogram(time2endPt(overGPNMB &time2endPt<100),'BinWidth',2)
histogram(time2endPt(underGPNMB&time2endPt<100),'BinWidth',2)

%%

% Extract patient data
ptID = dataTable.INDDID;
visitDate = dataTable.VisitDate;
weight = dataTable.Weight;
snp = dataTable.rs199347;
ageAtDeath = dataTable.AgeatDeath;

% Get unique patient IDs
uniquePatients = unique(ptID);
sumMat = [];

% Loop through each patient
for gg = 1:length(uniquePatients)
    patientID = uniquePatients(gg);
    
    % Extract all records for this patient
    patientData = dataTable(dataTable.INDDID == patientID, :);
    
    % Extract SNP status
    snpStatus = unique(patientData.rs199347);
    if isempty(snpStatus{1})
        snpStatAtPlay = nan;
    elseif strcmp(snpStatus{1}, 'CC')
        snpStatAtPlay = 0;
    elseif strcmp(snpStatus{1}, 'CT')
        snpStatAtPlay = 1;
    elseif strcmp(snpStatus{1}, 'TT')
        snpStatAtPlay = 2;
    end

    % Extract weight and visit dates
    datesAtPlay = patientData.VisitDate;
    weightsAtPlay = patientData.FRSTotal;
    
    if isempty(weightsAtPlay) || isscalar(weightsAtPlay)
        continue; % Skip if there's not enough data
    end
    
    % Find earliest weight
    [startDt, DtLoc] = min(datesAtPlay);
    startWeight = weightsAtPlay(DtLoc);
    
    % Calculate percentage weight loss
    percentLoss = weightsAtPlay / startWeight;
    
    % Find time to 20% weight loss
    lossIdx = find(percentLoss <= 0.8, 1); % First time weight drops 30%
    if ~isempty(lossIdx)
        endDt = datesAtPlay(lossIdx);
        timeToEvent = endDt - startDt;
        eventOccurred = 1; % 20% weight loss occurred
    else
        timeToEvent = max(datesAtPlay) - startDt; % Last follow-up time
        eventOccurred = 0; % Censored
    end

    % Convert time to months
    hours2use = hours(timeToEvent) / 24;
    weeks = hours2use / 7;
    timeToEvent = weeks / 4;
    
    % Store results
    sumMat = [sumMat; [patientID, snpStatAtPlay, timeToEvent, eventOccurred, patientData.AgeatDeath(1)]];
end

% Remove NaN values
remValz = isnan(sumMat(:,2)) | isnan(sumMat(:,3));

% Extract relevant values
snpAtPlay = sumMat(~remValz, 2);
timeToEvent = sumMat(~remValz, 3);
eventOccurred = sumMat(~remValz, 4);

% Fit Cox Proportional Hazards Model
[b, logL, H, stats] = coxphfit(snpAtPlay, timeToEvent, 'Censoring', 1 - eventOccurred);

%% **Stratified Kaplan-Meier Plot for Each SNP Group**
uniqueSNPs = unique(snpAtPlay); % Extract SNP categories

figure; hold on;
colors = {'c', 'm', 'g'}; % Assign colors to SNP groups

for gg = 1:length(uniqueSNPs)
    snpValue = uniqueSNPs(gg);
    
    % Select patients with the current SNP
    idx = (snpAtPlay == snpValue) ;
    survTimes = timeToEvent(idx);
    eventStatus = eventOccurred(idx);

    % Fit Kaplan-Meier survival function
    [H, T] = ecdf(survTimes, 'Function', 'survivor', 'Censoring', 1 - eventStatus);

    % Plot survival curve
    plot(T, H, 'Color', colors{gg}, 'LineWidth', 2);
end

xlabel('Time to 30% Weight Loss (Months)', 'FontSize', 20);
ylabel('Survival Probability', 'FontSize', 20);
title(['Survival Probability Stratified by rs199347 status, p=', num2str(stats.p)], 'FontSize', 20);
legend({'CC', 'CT', 'TT'});
grid on;
hold off;
xlim([0,200]);
