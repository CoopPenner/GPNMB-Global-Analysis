%quick count for RO1


dataTable=readtable('/Users/pennerc/Downloads/CSFTAUADRCPatientList.xlsx');
%% First getting patients with csf/plasma samples
ptID=dataTable.INDDID;

dx=dataTable.ADCCDx1;

uniquePt=unique(ptID);

csfMeasure=dataTable.CSF_Abeta42_pg_mL_;
plasmaMeasure=dataTable.Plasma_Abeta42_pg_mL_;

ptCollector=false(1,length(uniquePt));
ADCollect=false(1, length(uniquePt));
MCICollect=false(1, length(uniquePt));

ADIdx= (contains(dx,'AD')) ; 


MCIidx=  (contains(dx,'MCI')) &...
                        (~contains(dx, 'non memory')) &...
                        (~contains(dx, 'non-memory')) &...
                        (~contains(dx,'Possible'));




   for dd=1:length(uniquePt)
        ptAtPlay=uniquePt(dd);

        csfMeasurezMCI=csfMeasure(ptID==ptAtPlay & MCIidx);
        csfMeasurezAD=csfMeasure(ptID==ptAtPlay & ADIdx);

        plasmaMeasurezMCI=plasmaMeasure(ptID==ptAtPlay&MCIidx);
        plasmaMeasurezAD=plasmaMeasure(ptID==ptAtPlay &ADIdx);
     

            if (sum(~isnan(unique(plasmaMeasurezMCI) )) )>0 || sum(~isnan( unique (csfMeasurezMCI)) )>0 

            MCICollect(dd)=true;
            end


           if (sum(~isnan(unique(plasmaMeasurezAD))) )>0 || (sum(~isnan(unique(csfMeasurezAD))) )>0 
    
                ADCollect(dd)=true;
    
           end

            if  sum( ~isnan(plasmaMeasure(ptID==ptAtPlay )  ))>0 || sum( ~isnan(csfMeasure(ptID==ptAtPlay)  )) >0
            ptCollector(dd)=true;
            end



  end
    

%% then calculating MCI/AD populations 

dataTable=readtable('/Users/pennerc/Downloads/InqueryDatasets/FinalizedSets/ADC_MMSE_TotalPatients.xlsx');

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
onsetAge=dataTable.ADCAgeOnset;
onsetDate=dataTable.ADCYearOnset;

%% starting with data visualizations, lets plot cognitive slope across each of our snps by tertile


snpGroups = {'TT', 'CT', 'CC'};
colors = lines(3); % for SNPs
dxGroups = {'Normal', 'MCI', 'Alzheimer'}; % for labeling
dxMarkers = ["normal", "mci", "ad"]; % strings to match This is obvi catching a lot of varieties of MCI, I just do this for an initial check

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
            if ~contains(dxString, dxMarkers(dxType))
                continue;
            end

            if sum(ptIdx) > 2 % if there are at least 3 visits
              
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


singleID = unique(globalID);
scoreToTest = cogScore;


%% now outputting everything for the R01

snpGroups = {'TT', 'CT', 'CC'};
colors = lines(3); % for SNPs
dxGroups = {'Normal', 'MCI', 'Alzheimer'}; % for labeling
dxMarkers = ["normal", "mci", "ad"]; % strings to match This is obvi catching a lot of varieties of MCI, I just do this for an initial check

  
    allSlopes = {};
    startScorez={};
startScoreMat=[];
onsetDelayMat=[];
dxMat={};
snpMat={};
  groupIDs=unique(globalID);

        for ii = 1:length(groupIDs) % now all of our patients

            IDAtPlay = groupIDs(ii);
            ptIdx = globalID == IDAtPlay ;
            % Get diagnosis for this patient


                dxString = unique(dx(ptIdx));


                snpAtPlay= unique(snp(ptIdx));

                onsetAgeAtPlay=unique( onsetAge(ptIdx));

                    if ~isnan(unique(onsetDate(ptIdx)  ))
                        onsetDateAtPlay=datenum(  datetime(   {['2-July-',num2str(unique(onsetDate(ptIdx)  )  )]}   )  ); %I take the middle date of a diagnosis year
                    else
                        onsetDateAtPlay=nan;
                    end

                datesAtPlay = datenum(   testDate(ptIdx)); 
                scoresAtPlay = cogScore(ptIdx);
                falseIdx = abs(scoresAtPlay) > 30; %sometimes there are recorded values of 1000

                onsetDelay= (min(datesAtPlay)- onsetDateAtPlay)/365.25;

                datesAtPlay = (datesAtPlay - min(datesAtPlay)) / 365.25; %converting to years (accoutning for leap years : )

                [~,firstDate]= min(datesAtPlay);
                startScoreAtPlay= scoresAtPlay(firstDate);

                % % Fit linear slope
                % p = polyfit(datesAtPlay(~falseIdx), scoresAtPlay(~falseIdx), 1);
                % if abs(p(1)) > 10 | length(datesAtPlay)<3
                %     p(1) = nan;
                % end

                % slopes(end+1) = p(1);
                startScoreMat(end+1)=startScoreAtPlay;
                if onsetDelay>0
                onsetDelayMat(end+1)=onsetDelay;
                else
                onsetDelayMat(end+1)=nan;
                end
                dxMat(end+1)=dxString;

                snpMat(end+1)=snpAtPlay;

        end



singleID = unique(globalID);
scoreToTest = cogScore;

snpPres=~cellfun(@isempty,snpMat);


MCI=contains(dxMat,'MCI') & contains(dxMat,'memory') & ~contains(dxMat,'non-memory');

Alzheimer=contains(dxMat,'AD') ;

nansum((MCI  &   startScoreMat>18))


%% Initialize matrices
startScoreMat = nan(length(globalID),1);
timePassedMat = nan(length(globalID),1);

for dd = 1:length(singleID)
    ptAtPlay = singleID(dd);
    idxAtPlay = globalID == ptAtPlay;

    datesAtPlay = testDate(idxAtPlay);
    scoresAtPlay = scoreToTest(idxAtPlay);

    if isempty(scoresAtPlay) || isscalar(scoresAtPlay)
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



