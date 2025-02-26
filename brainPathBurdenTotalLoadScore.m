function [] = brainPathBurdenTotalLoadScore(pt2use, pathTable, pathType,  brainAreas,  diseaseName,pathCut,brainRegionCut   )

% Initialize variables
globalDx = pathTable.GlobalDx;
pathID = pathTable.INDDID;
ageAtDeath = pathTable.AgeatDeath;
ageAtOnset = pathTable.GlobalAgeOnset; 
disDur = ageAtDeath - ageAtOnset;
bioSex = pathTable.Sex;
snpStat = pathTable.rs199347;

overGPNMB = contains(snpStat,'TT'); 
het = contains(snpStat, 'CT');
underGPNMB = contains(snpStat,'CC'); 
emptySNP = cellfun(@isempty, snpStat);

% Create empty arrays to store all data points across regions
allPathScores = [];
allRegionLabels = [];
allSNP = [];
allAgeAtDeath = [];
allSex = [];
allDiseaseDur = [];
allPtID = [];


%outputting pathlogical load for each of our species in a given patient
%group we will take all brain areas in which at least 10% of patients have
%a 'rare' amount of pathology.

%initializing all variables
globalDx=pathTable.GlobalDx;
pathID=pathTable.INDDID;
%initializing covariates for glme
ageAtDeath=pathTable.AgeatDeath;
ageAtOnset=pathTable.GlobalAgeOnset; disDur=ageAtDeath-ageAtOnset;
bioSex=pathTable.Sex;

%organizing snp status readout
snpStat=pathTable.rs199347;
overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele
emptySNP=cellfun(@isempty,snpStat);





numSamp=nan(1,length(brainAreas));
avgBurd=nan(1,length(brainAreas));
percentPres=nan(1,length(brainAreas));
totalPathScoresUnder=[];
totalPathScoresHet=[];
totalPathScoresOver=[];
totalPathScores=[];


for dd=1:length(brainAreas)
    brainAreaAtPlay=brainAreas{dd};
[asynCont,tauCont,aBetaCont,TDPCont,neuronDropCont,gliosisCont] = pathScoreGenerate(brainAreaAtPlay,pathTable);

    if pathType==1
        val2Test=asynCont;
        pathName='Asyn';
    elseif pathType==2
        val2Test=aBetaCont;
            pathName='aBeta';
    elseif pathType==3
        val2Test=tauCont;
        pathName='Tau';
    elseif  pathType==4
        val2Test=TDPCont;
        pathName='TDP43';
    elseif pathType==5
        val2Test=neuronDropCont;
        pathName='NeuronLoss';
    elseif pathType==6
        val2Test=gliosisCont;
        pathName='Gliosis';
    end


keepVals=   pt2use' & ~emptySNP';
numSamp(dd)= sum( keepVals & ~isnan(val2Test) ); % outputting the number of path samples in a given brain area
avgBurd(dd)=nanmean(val2Test(keepVals)); %outputting mean path burden
percentPres(dd)= sum(val2Test(keepVals)>=pathCut)/sum(keepVals); % percentage of cases that have at least rare path burden and at least 50 cases
totalPathScoresUnder(dd,:)=val2Test(underGPNMB & pt2use);
totalPathScoresHet(dd,:)=val2Test(het & pt2use);
totalPathScoresOver(dd,:)=val2Test(overGPNMB & pt2use);
totalPathScores(dd,:)=val2Test(keepVals);

end


% avgBurd(percentPres< brainRegionCut |numSamp<50  )=[];
% 
% [highBurd,i]= maxk(avgBurd, 15);


%ok, now we have the brain areas to assess and the path values, let's just
%aggregate them and evaluate. 

 totalPathScoresUnder(numSamp<50 | percentPres<brainRegionCut,:)=[];
 totalPathScoresOver(numSamp<50 | percentPres<brainRegionCut,:)=[];
 totalPathScoresHet(numSamp<50 | percentPres<brainRegionCut,:)=[];
 totalPathScores(numSamp<50 | percentPres<brainRegionCut,:)=[];


underAvg=mean(totalPathScoresUnder,'omitnan'); %this is burden across all included areas for each snp carrier
hetAvg=mean(totalPathScoresHet,'omitnan'); %this is burden across all included areas for each snp carrier
overAvg=mean(totalPathScoresOver,'omitnan'); %this is burden across all included areas for each snp carrier
totAvg=mean(totalPathScores, 'omitnan');

%% plotting

%first setting up glme 


% Create table
modelData = table(pathID(keepVals), ageAtDeath(keepVals), bioSex(keepVals), disDur(keepVals), snpStat(keepVals), totAvg',  ...
    'VariableNames', {'ptID', 'ageAtDeath', 'Sex', 'diseaseDuration', 'SNP', 'pathScore'});

% Convert categorical variables to factors
modelData.Sex = categorical(modelData.Sex);
modelData.ptID = categorical(string(modelData.ptID));
modelData.SNP=categorical(modelData.SNP); %I'm not entirely sure if this is correct but I'm going categorical here because I'm not convinced a relationship would be purely linear



% Fit GLME 
glme = fitglme(modelData, ...
    'pathScore ~ 1 + Sex + ageAtDeath + diseaseDuration + SNP  + (1|ptID)', ...
    'Distribution', 'Normal', 'Link', 'Identity');


disp(glme)


figure


b=bar(1,nanmean(overAvg )) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(hetAvg )) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(underAvg )) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

hold on
a=scatter(rand(1, length(overAvg)  )+.5, overAvg, 'Marker', 'o' );
a.CData=[.8 .2 .5]; 

b=scatter(rand(1, length(hetAvg)  )+2.5, hetAvg, 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, length(underAvg)  )+4.5, underAvg, 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

legend({'AA (over Production)', 'GA','GG (under Production)'}, 'FontSize', 15)



title([diseaseName, ' Patients ', pathName, ' Burden', ' Across All brain regions with Pathology' ], 'FontSize',15)

ylabel([pathName, ' pathological ratings'] ,'FontSize', 15)




% Loop through each brain region and path type
for dd = 1:length(brainAreas)
    brainAreaAtPlay = brainAreas{dd};
    [asynCont, tauCont, aBetaCont, TDPCont, neuronDropCont, gliosisCont] = pathScoreGenerate(brainAreaAtPlay, pathTable);

        switch pathType
            case 1, val2Test = asynCont; pathName = 'Asyn';
            case 2, val2Test = aBetaCont; pathName = 'aBeta';
            case 3, val2Test = tauCont; pathName = 'Tau';
            case 4, val2Test = TDPCont; pathName = 'TDP43';
            case 5, val2Test = neuronDropCont; pathName = 'NeuronLoss';
            case 6, val2Test = gliosisCont; pathName = 'Gliosis';
         end

        keepVals = pt2use' & ~emptySNP' ;

 numSamp= sum( keepVals & ~isnan(val2Test) ); % outputting the number of path samples in a given brain area
 avgBurd=nanmean(val2Test(keepVals)); %outputting mean path burden

if numSamp <50 || avgBurd< pathCut
    continue
end

try
        % Append data from each region and path type
        allPathScores = [allPathScores; val2Test(keepVals)'];
        allRegionLabels = [allRegionLabels; repmat({brainAreaAtPlay}, sum(keepVals), 1)];  % Brain area cell array
        allSNP = [allSNP; snpStat(keepVals)];
        allAgeAtDeath = [allAgeAtDeath; ageAtDeath(keepVals)];
        allSex = [allSex; bioSex(keepVals)];
        allDiseaseDur = [allDiseaseDur; disDur(keepVals)];
        allPtID = [allPtID; pathID(keepVals)]; 





catch
b=1;

end
end

% Create a table with all collected data
modelData = table(allPtID, allAgeAtDeath, allSex, allDiseaseDur, allSNP, allRegionLabels, allPathScores, ...
    'VariableNames', {'ptID', 'ageAtDeath', 'Sex', 'diseaseDuration', 'SNP', 'BrainRegion', 'pathScore'});

% Convert categorical variables
modelData.Sex = categorical(modelData.Sex);
modelData.ptID = categorical(string(modelData.ptID));
modelData.SNP = categorical(modelData.SNP);
modelData.BrainRegion = categorical(modelData.BrainRegion);

% Fit GLME with BrainRegion as a cofactor
glme = fitglme(modelData, ...
    'pathScore ~ 1 + Sex + ageAtDeath + diseaseDuration + SNP + BrainRegion + (1|ptID)', ...
    'Distribution', 'Normal', 'Link', 'Identity');

disp(glme);






