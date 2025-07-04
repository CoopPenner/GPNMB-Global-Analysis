%basic check of rs199347 frequency and cognitive decline

addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');


%% reading in data
dataTable=readtable('/Users/pennerc/Downloads/InqueryDatasets/FinalizedSets/ADC_MMSE_TotalPatients_withUpdate.xlsx');

%dataTable=readtable('/Users/pennerc/Downloads/InqueryDatasets/AllPatients_extraSNP_withMMSE.xlsx');




%  dataTable=readtable('/Volumes/PC60/InqueryDatasets/AllPatients_extraSNP_DRS.xlsx');
%  dataTable=readtable('/Volumes/PC60/InqueryDatasets/allPatients_extraSNP_BNT.xlsx');
%  dataTable=readtable('/Volumes/PC60/InqueryDatasets/AllPatients_WithMoca_11_24.xlsx');

cogScore=dataTable.MMSETotal;
testDate=dataTable.TestDate;
dx=dataTable.ADCCDx1;
snp=dataTable.rs199347;
globalID=dataTable.INDDID;
Edu=dataTable.Education;
bioSex=dataTable.Sex;
ageATTest=dataTable.AgeatTest;
ageAtDeath=dataTable.AgeatDeath;
apoeHaplo=dataTable.APOE;


%% starting with data visualizations, lets plot cognitive slope across each of our snps by tertile


snpGroups = {'TT', 'CT', 'CC'};
colors = lines(3); % for SNPs
dxGroups = {'Normal', 'MCI', 'Alzheimer'}; % for labeling
dxMarkers = ["normal", "mci -memory", "ad"]; % strings to match

for dxType = 1:3 % we iterate through each diagnosis subtype
    figure;
    hold on;
    allSlopes = {};

    for jj = 1:length(snpGroups)
        snpAtPlay = snpGroups{jj};
        snpGroup = strcmp(snp, snpAtPlay);
        groupIDs = unique(globalID(snpGroup));
        slopes = [];

        for ii = 1:length(groupIDs) % now all of our patients
            IDAtPlay = groupIDs(ii);
            ptIdx = globalID == IDAtPlay & snpGroup;

            % Get diagnosis for this patient
            dxString = lower(string(dx(find(ptIdx, 1, 'first'))));

            % Check if diagnosis matches current group this is probably not
            % the most efficient way to implement this...
            if ~contains(dxString, 'mci') && ~contains(dxString, '- memory')
                continue;
            end

            if sum(ptIdx) > 1 % if there are at least 3 visits
              
                 datesAtPlay = datenum(testDate(ptIdx)); 
                scoresAtPlay = cogScore(ptIdx);
                falseIdx = abs(scoresAtPlay) > 30; %sometimes there are recorded values of 1000

                datesAtPlay = (datesAtPlay - min(datesAtPlay)) / 365.25; %converting to years (accoutning for leap years : )

                % Fit linear slope
                p = polyfit(datesAtPlay(~falseIdx), scoresAtPlay(~falseIdx), 1);
                if abs(p(1)) > 10 | length(datesAtPlay)<3
                    p(1) = nan;
                end

                slopes(end+1) = p(1);
            end
        end

        allSlopes{jj} = slopes;
% now plotting
        scatter(repelem(jj, length(slopes)), slopes, 50, 'filled', ...
            'MarkerFaceColor', colors(jj,:), 'MarkerEdgeColor', 'k');

        % Display mean slope
        disp([dxGroups{dxType} ', SNP ' snpAtPlay, ...
            ' mean slope: ', num2str(nanmean(slopes))])
    end

    xlim([0.5, 3.5]);
    xticks(1:3);
    xticklabels(snpGroups);
    ylabel('Slope of MMSE decline (points/year)');
    title(['Cognitive decline in ', dxGroups{dxType}, ' group']);
    grid on;
end

%% running my GLME
% Preprocess: Calculate time passed and baseline score per subject
singleID = unique(globalID);
scoreToTest = cogScore;

% Initialize matrices
startScoreMat = nan(length(globalID),1);
timePassedMat = nan(length(globalID),1);

for dd = 1:length(singleID)
    ptAtPlay = singleID(dd);
    idxAtPlay = globalID == ptAtPlay;

    datesAtPlay = testDate(idxAtPlay);
    scoresAtPlay = scoreToTest(idxAtPlay);

    if isempty(scoresAtPlay) || length(scoresAtPlay)<3
        continue;
    end

    % Find earliest test date and corresponding MMSE
    [startDt, DtLoc] = min(datesAtPlay);
    startScore = scoresAtPlay(DtLoc);

    % Calculate time passed in months
    timeDelta = datesAtPlay - startDt;
timePassed = days(timeDelta) / 30.44;  % average time passed in  months

    % Remove unusable cases )
    if sum(timePassed == 0) > 1 && length(timePassed) < 3
        timePassed(:) = nan;
    end

    % Assign values to  matrices
    timePassedMat(idxAtPlay) = timePassed;
    startScoreMat(idxAtPlay) = startScore;

end

% Filter for analysis
emptySNP = cellfun(@isempty, snp);
filter2use = contains(dx,'AD')&...%contains(dx, 'MCI') & contains(dx, ' - memory')& ... % only mci patients
             scoreToTest <= 30 & ... % I filter for this in the above code but it made me paranoid so I included it here too lol
             startScoreMat >= 10 & ... % patients that begin at the cutoff for "moderate dementia" are excluded
             ~isnan(scoreToTest) & ... % obviously we need scores
             ~emptySNP; % and geneotypes

% make model  table
modelData = table( ...
    categorical(string(globalID(filter2use))), ...
    timePassedMat(filter2use), ...
    categorical(bioSex(filter2use)), ...
    ageATTest(filter2use), ...
    categorical(snp(filter2use)), ...
    scoreToTest(filter2use), ...
    startScoreMat(filter2use), ...
    'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest', 'SNP', 'Scores', 'startScore'} ...
);

  modelData.SNP = reordercats(modelData.SNP, {'TT', 'CT', 'CC'});
% modelData.SNP = reordercats(modelData.SNP, {'AA', 'AC', 'CC'});

% Fit GLME  and alt model
glme = fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtTest + SNP * timePassed + (1|ptID)', ...
    'Distribution','normal','Link','identity');

glmeAlt = fitglme(modelData, ...
    'Scores ~ 1 + Sex + startScore + ageAtTest  + (1|ptID)', ...
    'Distribution','normal','Link','identity');

disp(glme);
compare(glmeAlt, glme);

%% simple glme plot


CCLowerCI=glme.Coefficients(9,7); CCLowerCI=CCLowerCI.Lower;
CCupperCI=glme.Coefficients(9,8); CCupperCI=CCupperCI.Upper;
CCEstimate=glme.Coefficients(9,2); CCEstimate=CCEstimate.Estimate;

snpStat=snp;

startPt= nanmean(startScoreMat(    strcmp(snpStat,'CC')     ));


    figure;
    hold on;
    colorCode=[.2,.6,.8];
plotGLMESlope(CCEstimate, CCLowerCI, CCupperCI,startPt,colorCode)


TCLowerCI=glme.Coefficients(8,7); TCLowerCI=TCLowerCI.Lower;
TCupperCI=glme.Coefficients(8,8); TCupperCI=TCupperCI.Upper;
TCestimate=glme.Coefficients(8,2); TCestimate=TCestimate.Estimate;
startPt= nanmean(startScoreMat(    strcmp(snpStat,'CT')     ));
colorCode=[.9,.4,.3];
plotGLMESlope(TCestimate, TCLowerCI, TCupperCI,startPt,colorCode)




TTEstimate=0; %TT is the comparator group so the 'slope' will always be 0, here I just use the CI for CT, in next iteration I'll just plot residuals for each patient, binned
colorCode=[0,1,1];
startPt= nanmean(startScoreMat(    strcmp(snpStat,'TT')     ));
plotGLMESlope(TTEstimate, CCLowerCI, CCupperCI,startPt,colorCode)

legend({'CC','','CT','','TT'},'FontSize',20);
xlabel('Time (Months)','FontSize',20);
ylabel('Modeled MMSE','FontSize',20)

title('Modeled MMSE Slope as a function of rs199347 status in MCI patients','FontSize',20)


%% plotting residuals etc

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
        idxAtPlay = modelData.SNP == uniqueSNPs(i) & modelData.timePassed >= timeBins(j) & modelData.timePassed < timeBins(j+1);
        
        if any(idxAtPlay)
            meanFittedOverTime(j, i) = nanmean(fittedVals(idxAtPlay));
            stdErrorOverTime(j, i) = nanstd(fittedVals(idxAtPlay)) / sqrt(sum(idxAtPlay)); % Standard error
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
ylabel('Fitted MMSE Scores','FontSize',20);
title('Fitted MMSE Over Time by Genotype','FontSize',20);
legend('Location', 'NorthEast','FontSize',30);
grid on;
hold off;











%% code cemetery ( I include this here for people who are really curious, this is all first pass pre cleanup read on if you dare)

% %% designing table and running glme
% %1 calculate time passed and get a Mat with the start scores for our model
% singleID=unique(globalID);
% scoreToTest=cogScore;
% 
% startScoreMat=nan(length(globalID),1);
% timePassedMat= nan(length(globalID),1);
% timePassedMatTest=[];
% idMat=[];
% for dd= 1:length(singleID)
%     datesAtPlay=testDate(globalID==singleID(dd));
%     scoresAtPlay=scoreToTest(globalID==singleID(dd));
%         if isempty(scoresAtPlay) || isscalar(scoresAtPlay) 
%             scoresAtPlay=nan;
%         end
%      [startDt,DtLoc]=min(datesAtPlay) ;  %finding earliest date (sometimes I think things come out of order)
%      startScore=scoresAtPlay(DtLoc);
%     timePassed=datesAtPlay-startDt; %order doesn't matter here
%     %converting to months
%     hours2use=hours(timePassed)/24; weeks=hours2use/7; timePassed=weeks/4;
% 
% if sum(timePassed==0)>1 & length(timePassed)<4
% timePassed=timePassed*nan;
% end
% 
%     %generating our output Mats
%     %timePassedMat=[timePassedMat;timePassed];
% timePassedMat(globalID==singleID(dd))=timePassed;
% 
% timePassedMatTest=[timePassedMatTest;timePassed]; %this was a really gnarly bug, the INDD dataset isn't output in numerical order!!! Luckily we caught it
% startScoreMat(globalID== singleID(dd))= startScore;
% 
% 
%     if length(timePassed)~= length(scoresAtPlay)
%     disp(['Patient: ', num2str(singleID(dd)), ' - Length timePassed: ', num2str(length(timePassed)), ' - Length scoresAtPlay: ', num2str(length(scoresAtPlay))]);
%     end
% end
% 
% 
% emptySNP=cellfun(@isempty,snp);
% 
% filter2use= contains(dx, 'Alzheimer')  & scoreToTest<=30 & startScoreMat>= 14   & ~isnan(scoreToTest)   & ~emptySNP; % define patient to use
% 
% %currently using the cutoff for moderate alzheimers 
% 
% 
% % Create table
% modelData = table(globalID(filter2use), timePassedMat(filter2use), bioSex(filter2use), ageATTest(filter2use), snp(filter2use), scoreToTest(filter2use), startScoreMat(filter2use), ...
%     'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest', 'SNP', 'Scores', 'startScore'});
% 
% % Convert categorical variables to factors
% modelData.Sex = categorical(modelData.Sex);
% modelData.ptID = categorical(string(modelData.ptID));
% modelData.SNP=categorical(modelData.SNP); %I'm not entirely sure if this is correct but I'm going categorical here because I'm not convinced a relationship would be purely linear
% 
% 
% 
% % Fit GLME 
% glme = fitglme(modelData, ...
%     'Scores ~ 1 + Sex + startScore + ageAtTest + SNP * timePassed + (1|ptID)', ...
%     'Distribution','normal','Link','identity')
% 
% 
% 
% 
% glmeAlt= fitglme(modelData, ...
%     'Scores ~ 1 + Sex + startScore + ageAtTest + (1|ptID)', ...
%     'Distribution', 'Normal', 'Link', 'Identity');
% 
% compare(glmeAlt, glme) % cool, significant in terms of model comp as well
% 
% disp(glme);
% 


% 
% 
% 
% 
% 
% 
% 
% 
% %cogScore=dataTable.MoCATotal;
% cogScore=dataTable.MMSETotal;
% %cogScore=dataTable.DRSTotalAge;
% testDate=dataTable.TestDate;
% dx=dataTable.GlobalDx;
% snp=dataTable.rs199347;
% globalID=dataTable.INDDID;
% Edu=dataTable.Education;
% bioSex=dataTable.Sex;
% ageATTest=dataTable.AgeatTest;
% ageAtDeath=dataTable.AgeatDeath;
% 
% 
% 
% slopeCollect=nan(1, length(singleID));
% startScore=nan(1, length(singleID));
% endScore=nan(1,length(singleID));
% year2DecCollect=nan(1, length(singleID));
% diagCollect=cell(1, length(singleID));
% snpCollect=cell(1, length(singleID));
% ageAtTestCollect=nan(1, length(singleID));
% edCollect=nan(1, length(singleID));
% sexCollect=cell(1, length(singleID));
% emptyIndex=false(1, length(singleID));
% IDCollect=startScore;
% 
% 
% 
% for tt= 1:length(singleID) %iterating through each patient in the dataset
% 
% relDates=testDate(   globalID==singleID(tt)    ); %all the test dates
% relScores=cogScore(globalID==singleID(tt) ); %all the scores
% %sigh these are sometimes out of order  
% [testYear,testMonth]=ymd(relDates);        
% pureTime=datenum(relDates); %convert it to pure time to sort
% [sorted,n]=sort(pureTime); %sorting on absolute value of time
% testYear=testYear(n); testMonth=testMonth(n); relScores=relScores(n); %then just outputting the reorder
% 
% testYear=testYear-testYear(1); %everything relative to start date;
% testYear=testYear*12;
% 
% testMonth=testMonth-testMonth(1);
% 
% allTime=testYear+testMonth;
% 
%         allTime((isnan(relScores)))=[];
%         relScores(isnan(relScores))=[]; 
%         [uniqueA, i,j] = unique(allTime,'first');
%         idxToRemove = find(not(ismember(1:numel(allTime),i))  | relScores'>30 )   ;
%         allTime(idxToRemove)=[]; relScores(idxToRemove)=[]; %sometimes for some reason there are duplicate times
%         allTime=allTime/12; %everything contextualized in years, since that is the standard by which slopes are presented
% 
%         if ~isempty(relScores) & length(relScores)>=3  
%         scoreVec=diff(relScores);  sigChange=find( (relScores-relScores(1) ) <=-2 ,1); %here we identify when scores sink two points below first measure.      
%         P=polyfit(allTime,relScores,1);
%         %P= (relScores(end)-relScores(1))/max(allTime); 
%         slopeCollect(tt)=P(1);
% if P(1)<-30
%     b=1;
% end
% 
% 
%         diagCollect(tt)=unique(dx(globalID==singleID(tt)));
%         snpCollect(tt)=unique(snp(globalID==singleID(tt)));
% 
% allAge=ageATTest(globalID==singleID(tt)   ) ;  ageAtTestCollect(tt)= mean(allAge, "omitnan")  ; %I take the age at first test
% allEd=Edu(globalID==singleID(tt)); edCollect(tt)=allEd(1); % I take the starting educational attainment
% sexCollect{tt}=unique(bioSex(globalID==singleID(tt))); %no categorization available in database for intersex individuals
% IDCollect(tt)= singleID(tt);
% 
%         % ageAtTest(tt)=ageAtTest(globalID==globalID(tt)) ;
% 
%             if ~isempty(sigChange)      
%         year2DecCollect(tt)= allTime(sigChange)  ;
%             else
%          year2DecCollect(tt)= nan ;
%             end
%         startScore(tt)=relScores(1); %collect starting point
%         endScore(tt)=relScores(end);   
%         else
%                 slopeCollect(tt)=nan;
%                 IDCollect(tt)= nan;
%                 ageAtTestCollect(tt)=nan;
%                 edCollect(tt)=nan;
%                 diagCollect{tt}='nan';
%                 snpCollect{tt}='nan';
%                 sexCollect{tt}='nan';
%                 endScore(tt)=nan;   
%                 emptyIndex(tt)=true;
% 
%         end
% 
% end
% 
% 
% 
% %% Plotting some basic clinical metrics
% 
% cellfun(@isempty, diagCollect)
% %setting up diff disease subgroups
% Parkinson= (contains(diagCollect, 'Parkinson'));
% Alzheimer= (contains(diagCollect, 'Alzheimer')   )  ;
% ALS= (contains(diagCollect, 'Amyotrophic'));
% DemLewy= (contains(diagCollect, 'Dementia with Lewy Bodies')); 
% MCI= (contains(diagCollect, 'Mild cognitive impairment')); 
% corticoBasal= (contains(diagCollect, 'Corticobasal syndrome')); 
% bvFTD= (contains(diagCollect, 'bvFTD-FTLD')); 
% PPA= (contains(diagCollect, 'PPA')); 
% supraNuc= (contains(diagCollect, 'Progressive supranuclear palsy')); 
% neuroPanel= Alzheimer + DemLewy  +Parkinson+ corticoBasal +bvFTD ;
% ParkinsonianPanel= Parkinson +corticoBasal +supraNuc +DemLewy;
% ofInterest=Parkinson+ALS+Alzheimer+bvFTD;
% 
% ParkinsonianDem=corticoBasal+DemLewy;
% other= ~neuroPanel;
% HC=(contains(diagCollect, 'Normal')); 
% 
% 
% overGPNMB= contains(snpCollect,'TT'); %the major allele
% het=contains(snpCollect, 'CT');
% underGPNMB=contains(snpCollect,'CC'); %the minor allele
% 
% 
% % figure
% % histogram(slopeCollect((overGPNMB) ),'BinWidth',1, 'FaceColor','b')
% % hold on
% % histogram(slopeCollect((underGPNMB)),'BinWidth',1, 'FaceColor','r')
% % histogram(slopeCollect((het)),'BinWidth',1, 'FaceColor','g')
% % xlabel('MOCA slope')
% 
% 
% pt2Plot=[ Parkinson  ];
% val2Plot=slopeCollect;
% SNP=snpCollect;
% remVals= isnan(val2Plot) | cellfun(@isempty,SNP)  ;
% 
% figure
% b=bar(1,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot )))  ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.5;
% b.CData(1,:) = [.8 .2 .5]; 
% ylabel('MOCA Slope')
% a=gca; a.XTickLabel=[];
% hold on
% 
% b=bar(3,nanmean(val2Plot(~remVals & het & pt2Plot )))  ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.5;
% b.CData(1,:) = [0 0.7 .25];
% 
% 
% b=bar(5,nanmean(val2Plot(~remVals & underGPNMB & pt2Plot )))  ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.5;
% b.CData(1,:) = [0.3 0.1 .6];
% 
% hold on
% a=scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot))+.5, val2Plot(~remVals & overGPNMB & pt2Plot), 'Marker', 'o' );
% a.CData=[.8 .2 .5]; 
% b=scatter(rand(1, sum(~remVals & het & pt2Plot))+2.5, val2Plot(~remVals & het & pt2Plot), 'Marker', 'o' );
% b.CData(1,:) = [0 0.7 .25];
% 
% c=scatter(rand(1, sum(~remVals & underGPNMB & pt2Plot))+4.5, val2Plot(~remVals & underGPNMB & pt2Plot), 'Marker', 'o' );
% c.CData(1,:) = [0.3 0.1 .6];
% 
% legend({'AA (over Production)', 'GC','GG (under Production)'})
% 
% 
% anovan(val2Plot(~remVals), {underGPNMB(~remVals)})
% 
% 
% 
% %% ok let's fit a glme to this function
% 
% % fixed effects of age sex SNP and educational attainment... consider using
% % INDDID as a random effect
% 
% %first let's hop in by formatting our table, there is actually very few
% %datapoints, at least for moca so let's avoid confusion
% 
% 
% remVals= isnan(val2Plot) | cellfun(@isempty,SNP) | startScore<14 | ~Parkinson | abs(val2Plot) >30   ;
% 
% 
% sexVar=categorical(string(sexCollect(~remVals)))';
% startVar=startScore(~remVals)';
% endVar=endScore(~remVals)';
% slopeVar=slopeCollect(~remVals)';
% % transformedSlope=slopeVar- (min(slopeVar)) +1;
% % transformedSlope=boxcox(transformedSlope);
% 
% snpVar=categorical(snpCollect(~remVals))';
% 
% ageVar=ageAtTestCollect(~remVals)';
% edVar=edCollect(~remVals)';
% diagVar=categorical(diagCollect(~remVals))';
% IDVar=(IDCollect(~remVals)  )';
% 
% varNames=["Sex","cogSlope", "rs199347","ageAtTest","edLevel",'dx','ID','startScore','endScore'];
% 
% glmeTable = table(sexVar ,slopeVar,snpVar,ageVar,edVar,diagVar,IDVar, startVar, endVar, 'Variablenames',varNames);
% 
% glmeTable.rs199347 = reordercats(glmeTable.rs199347, {'TT', 'CC', 'CT'});
% 
% 
% 
% glme = fitglme(glmeTable,...
% 		'cogSlope ~ 1  + rs199347 + Sex + startScore + ageAtTest',...
% 		'Distribution','normal','Link','identity','FitMethod','Laplace',...
% 		'DummyVarCoding','reference')
% 
% 
% glmeHeldOut=fitglme(glmeTable,...
% 		'cogSlope ~ 1 + Sex + startScore + (1|ageAtTest)   ',...
% 		'Distribution','normal','Link','identity','FitMethod','Laplace',...
% 		'DummyVarCoding','reference');
% 
% 
% compare(glme,glmeHeldOut)
% 
% 
% mdl = fitlm(glmeTable,'cogSlope')
% 
% 
% residuals = glme.residuals;
% figure
% histogram(residuals, 20); % Residual distribution
% plotResiduals(glme, 'fitted');
% 
% %nice, some interesting stuff, let's just quickly write this out so we can
% %hang onto it
% 
% 
% 
% 
% 
% 
% 
% 
% %% reorienting data correctly
% 
% origData2use=origData(:,12:end);
% origNames=origData2use.Properties.VariableNames;