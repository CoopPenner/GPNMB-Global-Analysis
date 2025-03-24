%% step 1, parse out all possible combinations to facilitate group creation

addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis')
%reading in path scores and other requisite info
pathTable=readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/allPathData_Cerad_Ammended.xlsx');
allNames=pathTable.Properties.VariableNames;
globalDx=pathTable.GlobalDx;
pdDX=pathTable.PDCDx;
CNDRDX=pathTable.CNDRDx;
alsDX=pathTable.ALSDx;
adDX=pathTable.ADCDx;
snpStat=pathTable.rs199347;
pathID=pathTable.INDDID;
allNames=allNames(1:270); %removing supplementary SNP's and whatnot
allNames=allNames(2:end);

over60=pathTable.AgeatDeath>60;
under60=pathTable.AgeatDeath<60;
emptySnp=cellfun(@isempty, snpStat);


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

basicData=readtable('/Volumes/PC60/InqueryDatasets/DetailedBasicInfo_allPatients.xlsx');


spreadTable=readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/allPathSpreadData.xlsx');



%it's a little strange to initialize things here since I pass the dataTable
%object to each of the downstream functions, but I feel it helps with
%readability and fast iteration 

Parkinson= (contains(globalDx, 'Parkinson') | contains(pdDX, 'Parkinson') | contains(CNDRDX, 'Parkinson')    );
Alzheimer= (contains(globalDx, 'Alzheimer')  );  
ALS= contains(globalDx, 'Amyotrophic'); %| contains(alsDX, 'Amyotrophic') | contains(CNDRDX, 'Amyotrophic')    );
DemLewy= (contains(globalDx, 'Dementia with Lewy Bodies')| contains(pdDX, 'Dementia with Lewy Bodies') | contains(CNDRDX, 'Dementia with Lewy Bodies')    );
MCI= (contains(globalDx, 'Mild cognitive impairment')); 
corticoBasal= contains(globalDx, 'Corticobasal syndrome'); %| contains(pdDX, 'Corticobasal syndrome') | contains(CNDRDX, 'Corticobasal syndrome')    );
bvFTD= (contains(globalDx, 'bvFTD-FTLD')); 
PPA= (contains(globalDx, 'PPA')); 
supraNuc= (contains(globalDx, 'Progressive supranuclear palsy')| contains(pdDX, 'Progressive supranuclear palsy') | contains(CNDRDX, 'Progressive supranuclear palsy')    );
neuroPanel= Alzheimer | DemLewy  | Parkinson| corticoBasal ;
ParkinsonianPanel= Parkinson |corticoBasal |supraNuc |DemLewy;
ofInterest=ALS|Parkinson |Alzheimer;
ParkinsonianDem=corticoBasal|DemLewy;
other= ~neuroPanel;
HC=(contains(globalDx, 'Normal')); 
allDis=~HC;

overGPNMB = contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele

emptySNP= cellfun(@isempty, snpStat);


sum(~emptySNP & ALS)
%% analysis 1 is there a relationship between SNP status and path burden
% we will look at all path and copath separately, with an emphasis on aSYN

brainAreaAtPlay='Amyg';
pt2use=DemLewy ;
pathType= 4; %1 is Asyn, 2 is aBeta, 3 is Tau, 4 is TDP 5 is Neuron Loss 6 is Gliosis
disName='ALS';
plotType='both'; %currently just two options Stacked bar graph ('stackBar')
%and a scatter bar graph (scatBar); both has both as subplots
snpStatPathBurdenTestr(pt2use,brainAreaAtPlay, pathTable,pathType,disName,plotType, basicData)




cortical1={'MC'}; cortical2={'Cing','Neocortical'}; cortical3={'SN','GP','CP'}; cortical4={'CS','Amyg','DG','EC'};


%% analysis 3 is there a relationship between path score and cognition and
% is this in some way mediated by SNP status... 
brainAreaAtPlay='Amyg';
pt2use=Alzheimer;
pathType=1; %1 is Asyn, 2 is aBeta, 3 is Tau, 4 is TDP  5 is neuron loss 6 is gliosis
secondPath=false;
secondPath2Plot=5;
secondPathCut=2;
snpStatPathCogTestr(pt2use,brainAreaAtPlay, pathTable, cogData,pathType,secondPath,secondPath2Plot,secondPathCut   )


%% analysis 4 looking specifically at the apparent inverse relationship between TDP burden and aSyn in alz
%additionally adding in functionality for neuronal dropout and whatnot for
%all path scores

brainAreaAtPlay='Amyg';
pt2use=over75;
pathCut=1;

figure
pathCogCompTestr(pt2use,brainAreaAtPlay, pathTable, cogData,pathCut   )

%% analysis 5 just showing what brain area has the highest burden in a given disease

pt2use=ALS;
pathType=4; %1 is Asyn, 2 is aBeta, 3 is Tau, 4 is TDP  5 is neuron loss 6 is gliosis

diseaseName='ALS';
brainAreas2Output=15; %num of areas to output max is 22 I think

figure
brainPathBurdenOutputr(pt2use, pathTable, pathType,  brainAreas, diseaseName,brainAreas2Output);

%% analysis 6 let's map TDP spread throughout the brain, this is based on  Brettschneider et al 2013


pt2use=ALS;
disName='ALS';
tdpBrainStagerGPNMB(pt2use, pathTable, disName )

%% analysis 7 outputting total pathological burden for each species


pt2use=ALS;
pathType=4; %1 is Asyn, 2 is aBeta, 3 is Tau, 4 is TDP  5 is neuron loss 6 is gliosis

diseaseName='ALS';

brainRegionCut=0.1; %brain regions must have this percentage of cases with the specified pathBurden to be included
pathCut=1; %minimum pathology to be considered present


brainPathBurdenTotalLoadScore(pt2use, pathTable, pathType,  brainAreas, diseaseName,pathCut,brainRegionCut);


%% analysis 8 now aBeata Spread

pt2use=allDis;
pathType=2;
disName='Alzheimer"s';
plotType='both'; %currently just two options Stacked bar graph ('stackBar')

pathSpreadStandardTestr(pt2use, pathTable,pathType, disName, plotType  )



%% Braak eval 

spreadTable.Braak03;
spreadTable.Braak06;


 BraakEval(pt2use, spreadTable, disName, plotType   )


%% analysis 9 now Tau Spread


%% analysis 10 now Asyn






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




