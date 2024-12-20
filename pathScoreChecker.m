%% step 1, parse out all possible combinations to facilitate group creation

addpath('/Users/pennerc/Documents/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis')
%reading in path scores and other requisite info
pathTable=readtable('/Volumes/PC60/InqueryDatasets/allPatients_extraSNP_pathInfo.xlsx');
allNames=pathTable.Properties.VariableNames;
globalDx=pathTable.GlobalDx;
snpStat=pathTable.rs199347;
pathID=pathTable.INDDID;

allNames=allNames(1:270); %removing supplementary SNP's and whatnot
allNames=allNames(2:end);

%outputting all brain areas
tauContain=find(contains(allNames, 'Tau'));
brainAreas=cell(1,length(tauContain));
    for dd=1:length(tauContain)
    brainAreas{dd}= allNames{tauContain(dd)}(1:end-3);
    end
%removing NeocorticalT because I don't really understand what that means
brainAreas(strcmp('NeocorticalT',brainAreas))=[];


% reading in cog scores
cogData=readtable('/Volumes/PC60/InqueryDatasets/allCogScoreSlopesMMSE.csv');
cogID=cogData.ID;
cogSlope=cogData.cogSlope;
endScorez=cogData.endScore;


%it's a little strange to initialize things here since I pass the dataTable
%object to each of the downstream functions, but I feel it helps with
%readability and fast iteration 

Parkinson= (contains(globalDx, 'Parkinson'));
Alzheimer= (contains(globalDx, 'Alzheimer')   )  ;
ALS= (contains(globalDx, 'Amyotrophic'));
DemLewy= (contains(globalDx, 'Dementia with Lewy Bodies')); 
MCI= (contains(globalDx, 'Mild cognitive impairment')); 
corticoBasal= (contains(globalDx, 'Corticobasal syndrome')); 
bvFTD= (contains(globalDx, 'bvFTD-FTLD')); 
PPA= (contains(globalDx, 'PPA')); 
supraNuc= (contains(globalDx, 'Progressive supranuclear palsy')); 
neuroPanel= Alzheimer + DemLewy  +Parkinson+ corticoBasal + ALS;
ParkinsonianPanel= Parkinson +corticoBasal +supraNuc +DemLewy;
ofInterest=ALS+corticoBasal;
ParkinsonianDem=corticoBasal+DemLewy;
other= ~neuroPanel;
HC=(contains(globalDx, 'Normal')); 


overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele


%% analysis 1 is there a relationship between SNP status and path burden
% we will look at all path and copath separately, with an emphasis on aSYN

brainAreaAtPlay='Amyg';
pt2use=Alzheimer ;
pathType= 1; %1 is Asyn, 2 is aBeta, 3 is Tau, 4 is TDP 5 is Neuron Loss 6 is Gliosis
disName='Alzheimer"s';
plotType='scatBar'; %currently just two options Stacked bar graph ('stackBar')
%and a scatter bar graph (scatBar);
snpStatPathBurdenTestr(pt2use,brainAreaAtPlay, pathTable,pathType,disName,plotType)






%% analysis 3 is there a relationship between path score and cognition and
% is this in some way mediated by SNP status... 
brainAreaAtPlay='Amyg';
pt2use=Alzheimer;
pathType=4; %1 is Asyn, 2 is aBeta, 3 is Tau, 4 is TDP
snpStatPathCogTestr(pt2use,brainAreaAtPlay, pathTable, cogData,pathType   )


%% analysis 4 looking specifically at the apparent inverse relationship between TDP burden and aSyn in alz
%additionally adding in functionality for neuronal dropout and whatnot for
%all path scores

brainAreaAtPlay='Amyg';
pt2use=neuroPanel;
figure
pathCogCompTestr(pt2use,brainAreaAtPlay, pathTable, cogData   )



%% analysis 2 previously I have shown
% that rs199347 status predicts cog decline. is this mediated in anyway by
% concomitant pathology?
secondaryCut=4;
val2Plot=cogScoreAtPlay;
SNP=snpStat;
remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)   ;

figure
b=bar(1,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel('end score')
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Plot(~remVals & het & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Plot(~remVals & underGPNMB & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

hold on
scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot & asynCont'<secondaryCut))+.5, val2Plot(~remVals & overGPNMB & pt2Plot & asynCont'<secondaryCut), 'Marker', 'o', 'MarkerEdgeColor', [.8 .2 .5] );
scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot & asynCont'>= secondaryCut ))+.5, val2Plot(~remVals & overGPNMB & pt2Plot & asynCont'>= secondaryCut ), 'Marker', 'o','MarkerFaceColor','k', 'MarkerEdgeColor', [.8 .2 .5]  );


scatter(rand(1, sum(~remVals & het & pt2Plot & asynCont'<secondaryCut))+2.5, val2Plot(~remVals & het & pt2Plot & asynCont'<secondaryCut), 'Marker', 'o', 'MarkerEdgeColor', [0 0.7 .25] );
scatter(rand(1, sum(~remVals & het & pt2Plot & asynCont'>= secondaryCut ))+2.5, val2Plot(~remVals & het & pt2Plot & asynCont'>= secondaryCut), 'Marker', 'o','MarkerFaceColor','k', 'MarkerEdgeColor', [0 0.7 .25] );


scatter(rand(1, sum(~remVals & underGPNMB & pt2Plot & asynCont'<secondaryCut))+4.5, val2Plot(~remVals & underGPNMB & pt2Plot & asynCont'<secondaryCut), 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.1 .6]);
scatter(rand(1, sum(~remVals & underGPNMB & pt2Plot & asynCont'>= secondaryCut ))+4.5, val2Plot(~remVals & underGPNMB & pt2Plot & asynCont'>= secondaryCut), 'Marker', 'o','MarkerFaceColor','k', 'MarkerEdgeColor', [0.3 0.1 .6] );





legend({'AA (over Production)', 'GC','GG (under Production)'})




%% analysis 3 is there a relationship between path score and cognition and
% is this in some way mediated by SNP status... 
brainAreaAtPlay='Ang';
pt2use=Parkinson;
snpStatPathCogTestr(pt2use,brainAreaAtPlay, pathTable, cogData   )
