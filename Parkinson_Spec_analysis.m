%Sigh, eventually I'll make all of these one pipeline, but the way I've had
%to output INDD datasets makes this easier for me to think about since I have to do lots of fussy things for each patient group... 

%This contains all the broad data levadopa doses, etc
refDataTable=readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/PD_Broad_Overview_data_including_SNP.xlsx');

score2test=1; % 1 is for DRS, 2 is for UPDRSII

if score2test==1
    % DRS scores
    testData=readtable("/Volumes/PC60/InqueryDatasets/FinalizedSets/PDC_Psychometrics_all.xlsx");
    scoreAtPlay=testData.DRSTotalAge;
elseif score2test==2
    %UPDRSII scores
    testData= readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/UPDRS_Total.xlsx');
    scoreAtPlay=testData.P3Total;
end




%cogScore=dataTable.MoCATotal;
%cogScore=dataTable.DRSTotalAge;
if score2test==1
testDate=testData.TestDate;
elseif score2test==2
testDate=testData.VisitDate;
end




%% so annoying, ammending one table to another because of INDD output


% Create a mapping from refDataTable to testData based on INDDID
[isMatched, idx] = ismember(testData.INDDID, refDataTable.INDDID);

% Extract valid indices
validIdx = idx(isMatched);

% Initialize new variables with appropriate data types
testData.MotorDx1 = cell(height(testData), 1);   
testData.rs199347 = cell(height(testData), 1);  
testData.Education = nan(height(testData), 1);  
testData.Sex = cell(height(testData), 1);      
testData.AgeatDeath = nan(height(testData), 1);  
testData.MotorAgeOnset = nan(height(testData), 1);  

% Assign only for matched patients
testData.MotorDx1(isMatched) = refDataTable.MotorDx1(validIdx);
testData.rs199347(isMatched) = refDataTable.rs199347(validIdx);
testData.Education(isMatched) = refDataTable.Education(validIdx);
testData.Sex(isMatched) = refDataTable.Sex(validIdx);
testData.AgeatDeath(isMatched) = refDataTable.AgeatDeath(validIdx);
testData.MotorAgeOnset(isMatched) = refDataTable.MotorAgeOnset(validIdx);


% ok now for the harder task of matching levadopa doses



% Initialize the new variable in testData
testData.LevodopaDose = nan(height(testData), 1);
% Loop through each testData row to find the most recent VisitDate
for dd = 1:height(testData)
    patientID = testData.INDDID(dd);
    testDate2use = testDate(dd);
    % Get all visit records 
    patientVisits = refDataTable(refDataTable.INDDID == patientID, :);
    % Find the most recent VisitDate before or equal to TestDate so even if
    % a date is closer to the later levodopa dose we assume earlier was
    % being given
    validVisits = patientVisits(patientVisits.VisitDate <= testDate2use, :);
    if ~isempty(validVisits)  % If at least one valid visit is found
        [~, latestIdx] = max(validVisits.VisitDate); % Find the latest visit
        testData.LevodopaDose(dd) = validVisits.Levodopa(latestIdx);
    end
end



%% designing table and running glme
%1 calculate time passed and get a Mat with the start scores for our model
singleID=unique(testData.INDDID);

snpArray=nan(1,height(testData));

startScoreMat=nan(height(testData),1);
timePassedMat= nan(height(testData),1);
idMat=nan(height(testData),1);
snpMat=nan(height(testData),1);
numVisits=nan(height(testData),1);

for dd= 1:length(singleID)
    datesAtPlay=testDate(testData.INDDID==singleID(dd));
    scoresAtPlay=scoreAtPlay(testData.INDDID==singleID(dd));
        if isempty(scoresAtPlay) || isscalar(scoresAtPlay) %if nothing was collected or only one val was collected
            scoresAtPlay=nan;
            numVisitsAtPlay=nan;
            snpArray=nan;
        end

     [startDt,DtLoc]=min(datesAtPlay) ;  %finding earliest date (sometimes I think things come out of order)
     startScore=scoresAtPlay(DtLoc);
     numVisitsAtPlay=length(scoresAtPlay);
    time2use=datesAtPlay-startDt; %order doesn't matter here
    %converting to months
    [timePassed] = timeOutput(time2use);

    %generating our output Mats
    timePassedMat(testData.INDDID==singleID(dd))=timePassed;


   startScoreHold=ones(length(timePassed),1); startScoreHold=startScoreHold* startScore;

    startScoreMat(testData.INDDID==singleID(dd))= startScoreHold; 
    visitHold=ones(length(startScoreHold),1); numVisitsAtPlay=visitHold*length(startScoreHold);
    numVisits(testData.INDDID==singleID(dd))=numVisitsAtPlay;



snpArray=testData.rs199347(testData.INDDID==singleID(dd));



if isempty(snpArray(1))
    snpMat(testData.INDDID==singleID(dd))=nan;
elseif ~isempty(snpArray(1))
snpMat(testData.INDDID==singleID(dd))=nan;
end



end




emptySex=cellfun(@isempty, testData.Sex); %idk why but one patient doesn't have this recorded. 
% Ensure MotorDx1 is a cell array
    % Replace empty cells with an empty string
    emptyCells = cellfun(@isempty, testData.MotorDx1); 
    testData.MotorDx1(emptyCells) = {''}; % Replace empty entries with ''



if score2test==1 % 1 is for DRS, 2 is for UPDRSII
  filter2use = contains(string(testData.MotorDx1), 'Parkinson') & ~emptySex & numVisits>=4 & startScoreMat>4 ;

% Create table
modelData = table(testData.INDDID(filter2use), timePassedMat(filter2use), testData.Sex(filter2use), testData.AgeatTest(filter2use), testData.rs199347(filter2use), scoreAtPlay(filter2use), startScoreMat(filter2use), ...
    'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest', 'SNP', 'Scores', 'startScore'});

% Convert categorical variables to factors
modelData.Sex = categorical(modelData.Sex);
modelData.ptID = categorical(string(modelData.ptID));
modelData.SNP=categorical(modelData.SNP); %I'm not entirely sure if this is correct but I'm going categorical here because I'm not convinced a relationship would be purely linear
% Fit GLME 
glme = fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtTest + SNP * timePassed + (1|ptID)', ...
    'Distribution','normal','Link','identity')

else %if doing the motor score include levodopa dose as a cofactor 

    filter2use = contains(string(testData.MotorDx1), 'Parkinson') & ~emptySex & numVisits>4  ;


% Create table
modelData = table(testData.INDDID(filter2use), timePassedMat(filter2use), testData.Sex(filter2use), testData.AgeatTest(filter2use), testData.rs199347(filter2use), scoreAtPlay(filter2use), startScoreMat(filter2use), testData.LevodopaDose(filter2use), ...
    'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest', 'SNP', 'Scores', 'startScore', 'Levodopa'});

% Convert categorical variables to factors
modelData.Sex = categorical(modelData.Sex);
modelData.ptID = categorical(string(modelData.ptID));
modelData.SNP=categorical(modelData.SNP); %I'm not entirely sure if this is correct but I'm going categorical here because I'm not convinced a relationship would be purely linear
% Fit GLME 
glme = fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtTest + Levodopa + SNP * timePassed + (1|ptID)', ...
    'Distribution','normal','Link','identity')

disp(glme);



end










glmeAlt= fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtTest + (1|ptID)', ...
    'Distribution', 'Normal', 'Link', 'Identity');

compare(glmeAlt, glme) % cool, significant in terms of model comp as well

disp(glme);

%% plotting glme output and eventually running survival analysis


TTLowerCI=glme.Coefficients(10,7); TTLowerCI=TTLowerCI.Lower;
TTupperCI=glme.Coefficients(10,8); TTupperCI=TTupperCI.Upper;
TTestimate=glme.Coefficients(10,2); TTestimate=TTestimate.Estimate;

snpStat=testData.rs199347;

startPt= nanmean(startScoreMat(    strcmp(snpStat,'TT') & filter2use    ));


    figure;
    hold on;
    colorCode=[.2,.6,.8];
plotGLMESlope(TTestimate, TTLowerCI, TTupperCI,startPt,colorCode)


TCLowerCI=glme.Coefficients(9,7); TCLowerCI=TCLowerCI.Lower;
TCupperCI=glme.Coefficients(9,8); TCupperCI=TCupperCI.Upper;
TCestimate=glme.Coefficients(9,2); TCestimate=TCestimate.Estimate;
startPt= nanmean(startScoreMat(    strcmp(snpStat,'CT') & filter2use      ));
colorCode=[.9,.4,.3];
plotGLMESlope(TCestimate, TCLowerCI, TCupperCI,startPt,colorCode)
legend({'TT','','CT'});
xlabel('Time (Months)');
ylabel('Estimated Age Adjusted DRS')