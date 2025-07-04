function [qVals,pValVec] = quickAnovaForGWASAnalysis(usableGenes,easyAccessData,filter2use,genePres,origNames,dx,age,sexAtPlay,INDDID,normFactor)



pValVec=nan(1,sum(genePres));

        for tt=1:sum(genePres)
        
             origLoc = (strcmp(origNames, usableGenes{tt}));
             val2use=table2array(easyAccessData(:,origLoc));
             val2use(val2use<0)=0; %because of log2 normalization some things can be 0, I'll just consider this below the limit of detection.

             if sum(val2use>0)==0
                 pValVec(tt)=1;
                 continue
             end

             


% make model  table
modelData = table( ...
    categorical(string(INDDID(filter2use))), ...
    categorical((sexAtPlay(filter2use))), ...
    categorical(dx(filter2use)), ...
    age(filter2use), ...
    val2use(filter2use),...
    normFactor(filter2use)',...
    'VariableNames', {'ID', 'Sex', 'dx', 'ageAtTest', 'ProtLevel','NormFactor'} ...
);

% Fit GLME  and alt model
glme = fitglme(modelData, ...
    'ProtLevel ~ 1 + Sex + ageAtTest + dx  + NormFactor + (1|ID)', ...
    'Distribution','normal','Link','identity');
 disp(glme)



 pValVec(tt)=glme.Coefficients.pValue(find(contains(glme.Coefficients.Name, 'dx'))); 



             
            %[pValVec(tt),~,~]=anovan(val2use( ~filter2use) , {dx(~filter2use)},'display','off');


        end

        qVals = bh_fdr(pValVec);


% 
% 
% 
% %% Initialize data
% 
% addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/')
% 
% % Load reference data with levodopa doses and SNP info
% refDataTable = readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/PD_Broad_Overview_data_including_SNP.xlsx'); %sigh have to ammend one data structure to another cause of INDD crashes
% 
% 
% %  !!USER DECISION!!
% score2test = 1; % 1 = DRS, 2 = UPDRSII
% 
% 
% switch score2test
%     case 1
%         testData = readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/upDALL_CognitiveTesting.xlsx');
%         testScore = testData.DRSTotalAge;
%         testDate = testData.TestDate;
%     case 2
%         testData = readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/upDALL_MotorTesting.xlsx');
%         testScore = testData.P3Total;
%         testDate = testData.VisitDate;
% end
% 
% % Match reference info to testData
% [isMatched, idxAtPlay] = ismember(testData.INDDID, refDataTable.INDDID); %matching test data to data containing levodopa
% validIdx = idxAtPlay(isMatched);
% 
% % Initialize new variables
% fieldsToAssign = {'MotorDx1', 'rs199347', 'Sex'};
% for dd = fieldsToAssign
%     testData.(dd{1}) = repmat({''}, height(testData), 1);
%     testData.(dd{1})(isMatched) = refDataTable.(dd{1})(validIdx);
% end
% testData.Education = nan(height(testData),1);
% testData.AgeatDeath = nan(height(testData),1);
% testData.MotorAgeOnset = nan(height(testData),1);
% testData.Education(isMatched) = refDataTable.Education(validIdx);
% testData.AgeatDeath(isMatched) = refDataTable.AgeatDeath(validIdx);
% testData.MotorAgeOnset(isMatched) = refDataTable.MotorAgeOnset(validIdx);
% 
% % Ok now for the hard part, matching levodopa doses this will only be part
% % of UPDRSII model
% testData.LevodopaDose = nan(height(testData),1);
% for ff = 1:height(testData)
%     ptAtPlay = testData.INDDID(ff);
%     dateAtPlay = testDate(ff);
%     vRows = refDataTable(refDataTable.INDDID == ptAtPlay & refDataTable.VisitDate <= dateAtPlay, :); % we use nearest prior date
%     if ~isempty(vRows)
%         [~, lastVisitIdx] = max(vRows.VisitDate);
%         testData.LevodopaDose(ff) = vRows.Levodopa(lastVisitIdx);
%     end
% end
% 
% %% We need to extract and calculate a few things for our model (time passed, first score num visits)
% singleID = unique(testData.INDDID); %all patients
% startScoreMat = nan(height(testData),1); 
% timePassedMat = nan(height(testData),1);
% numVisits = nan(height(testData),1);
% 
% for ff = 1:length(singleID)
%     ptAtPlay = singleID(ff);
%     idxAtPlay = testData.INDDID == ptAtPlay;
%     datesAtPlay = testDate(idxAtPlay); %dates of test
%     scoreAtPlay = testScore(idxAtPlay); %scores of test
% 
%     if isempty(scoreAtPlay) || isscalar(scoreAtPlay) % if the score is singular (one visit also filtered later) or nanned out
%         continue;
%     end
% 
%     [startDt, loc] = min(datesAtPlay); %first visit
%     startScore = scoreAtPlay(loc); 
%     timePassed = days(datesAtPlay - startDt) / 30.44; %time Passed since first visit in months
% 
%     startScoreMat(idxAtPlay) = repmat(startScore, sum(idxAtPlay), 1); %this will be the same obviously
%     timePassedMat(idxAtPlay) = timePassed;  %loading
%     numVisits(idxAtPlay) = repmat(length(scoreAtPlay), sum(idxAtPlay), 1); %loading
% end
% 
% % generating meaningful filters (nanned values actually get automatically
% % removed by glmefit but best to be intentional
% emptySex = cellfun(@isempty, testData.Sex); %no bio sex
% testData.MotorDx1(cellfun(@isempty, testData.MotorDx1)) = {''}; %no diagnosis
% 
% if score2test == 1 %DRS
%     filter2use = contains(testData.MotorDx1, 'Parkinson') & ~emptySex & numVisits >= 3 & startScoreMat > 6; %not demented at baseline at least 3 visits only PD patients
% else %UPDRSIII
%     filter2use = contains(testData.MotorDx1, 'Parkinson') & ~emptySex & numVisits >= 3;
% end
% 
% % Build model table
% modelData = table( ...
%     categorical(string(testData.INDDID(filter2use))), ...
%     timePassedMat(filter2use), ...
%     categorical(testData.Sex(filter2use)), ...
%     testData.AgeatTest(filter2use), ...
%     categorical(testData.rs199347(filter2use)), ...
%     testScore(filter2use), ...
%     startScoreMat(filter2use), ...
%     testData.Education(filter2use),...
%     testData.GlobalAgeOnset(filter2use),...
%     'VariableNames', {'ptID', 'timePassed', 'Sex', 'ageAtTest', 'SNP', 'Scores', 'startScore','Edu','ageOnset'} ...
% );
% 
% if score2test == 2 %add levodopa for the UPDRSII cases
%     modelData.Levodopa = testData.LevodopaDose(filter2use);
% end
% 
% modelData.SNP = reordercats(modelData.SNP, {'TT', 'CT', 'CC'});
% 
% 
% % Fit GLME model
% if score2test == 1
%     formula = 'Scores ~ 1 + Sex + startScore + ageAtTest + Edu + ageOnset + SNP * timePassed + (1|ptID)';
% else
%     formula = 'Scores ~ 1 + Sex + startScore + ageAtTest + Levodopa + ageOnset + SNP * timePassed + (1|ptID)';
% end
% glme = fitglme(modelData, formula, 'Distribution','normal','Link','identity');
% 
% disp(glme);
% % 
% % % Alternative model
% % glmeAlt = fitglme(modelData, 'Scores ~ 1 + Sex + startScore + ageAtTest + (1|ptID)', 'Distribution','normal','Link','identity');
% % compare(glmeAlt, glme);
% % 
% % 
% % 
% % %% plotting glme output and eventually running survival analysis
% % 
% % 
% % CCLowerCI=glme.Coefficients(11,7); CCLowerCI=CCLowerCI.Lower;
% % CCupperCI=glme.Coefficients(11,8); CCupperCI=CCupperCI.Upper;
% % CCestimate=glme.Coefficients(11,2); CCestimate=CCestimate.Estimate;
% % 
% % snpStat=testData.rs199347;
% % 
% % startPt= nanmean(startScoreMat(    strcmp(snpStat,'CC') & filter2use    ));
% % 
% % 
% %     figure;
% %     hold on;
% %     colorCode=[.2,.6,.8];
% % plotGLMESlope(CCestimate, CCLowerCI, CCupperCI,startPt,colorCode)
% % 
% % 
% % TCLowerCI=glme.Coefficients(10,7); TCLowerCI=TCLowerCI.Lower;
% % TCupperCI=glme.Coefficients(10,8); TCupperCI=TCupperCI.Upper;
% % TCestimate=glme.Coefficients(10,2); TCestimate=TCestimate.Estimate;
% % startPt= nanmean(startScoreMat(    strcmp(snpStat,'CT') & filter2use      ));
% % colorCode=[.9,.4,.3];
% % plotGLMESlope(TCestimate, TCLowerCI, TCupperCI,startPt,colorCode)
% % 
% % 
% % TCLowerCI=glme.Coefficients(9,7); TCLowerCI=TCLowerCI.Lower;
% % TCupperCI=glme.Coefficients(9,8); TCupperCI=TCupperCI.Upper;
% % TCestimate=glme.Coefficients(9,2); TCestimate=TCestimate.Estimate;
% % 
% % 
% % TTUpperCI=mean([CCupperCI,CCupperCI]);
% % TTLowerCI=mean([CCLowerCI,CCLowerCI]);
% % 
% % 
% % TTestimate=0; %TT is the comparator group so the 'slope' will always be 0, here I just use the CI for CT, in next iteration I'll just plot residuals for each patient, binned
% % colorCode=[0,1,1];
% % startPt= nanmean(startScoreMat(    strcmp(snpStat,'TT')     ));
% % plotGLMESlope(TTestimate, TTLowerCI, TTUpperCI,startPt,colorCode)
% % 
% % 
% % 
% % legend({'CC','','CT', '', 'TT'},'FontSize',20);
% % xlabel('Time (Months)','FontSize',20);
% % ylabel('Estimated Age Adjusted DRS','FontSize',20)
% % 
% % title('Modeled DRS Slope as a function of rs199347 status in PD patients','FontSize',20)
% % 
% % 
% % %% now plotting for visualization of slopes 
% % 
% % 
% % figure;
% % hold on;
% % allSlopes = {};
% % globalID = testData.INDDID;
% % cogScore = testScore;
% % 
% % for jj = 1:length(snpGroups)
% %     snpAtPlay = snpGroups{jj};
% %     snpGroup = strcmp(testData.rs199347, snpAtPlay);
% %     groupIDs = unique(globalID(snpGroup));
% %     slopes = [];
% % 
% %     for ii = 1:length(groupIDs)
% %         IDAtPlay = groupIDs(ii);
% %         ptIdx = globalID == IDAtPlay & snpGroup;
% % 
% %         if sum(ptIdx) > 1 %if there's data to plot
% %             datesAtPlay = datenum(testDate(ptIdx));
% %             scoresAtPlay = cogScore(ptIdx);
% %             badVal = abs(scoresAtPlay) > 30; %sometimes incorrect data added to the mix
% % 
% %             datesAtPlay = (datesAtPlay - min(datesAtPlay)) / 30.44; % in months
% %             p = polyfit(datesAtPlay(~badVal), scoresAtPlay(~badVal), 1);
% % 
% %             if  length(datesAtPlay)<3 | abs(p(1))>.4
% %                 p(1)=nan;
% %                 % testData.rs199347(find(globalID==IDAtPlay))
% %                 % testDate(find(globalID==IDAtPlay))
% %                 % testScore(find(globalID==IDAtPlay))
% %                 % slopeDope=p(1)
% %             end
% %             slopes(end+1) = p(1);
% %         end
% %     end
% % 
% %     allSlopes{jj} = slopes;
% %     scatter(repelem(jj, length(slopes)), slopes, 50, 'filled', ...
% %         'MarkerFaceColor', colors(jj,:), 'MarkerEdgeColor', 'k');
% % 
% %     disp(['SNP ', snpAtPlay, ' mean slope: ', num2str(nanmean(slopes))]);
% % end
% % 
% % xlim([0.5, 3.5]);
% % xticks(1:3);
% % xticklabels(snpGroups);
% % ylabel('Slope of decline (points/year)');
% % title('Longitudinal cognitive slope by SNP group');
% % grid on;
