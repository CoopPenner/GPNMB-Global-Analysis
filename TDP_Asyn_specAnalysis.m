
addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis')
%reading in path scores and other requisite info
pathTable=readtable('/Volumes/PC60/InqueryDatasets/allPatients_updated_pathInfo.xlsx');
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
ofInterest=ALS|corticoBasal;
ParkinsonianDem=corticoBasal|DemLewy;
other= ~neuroPanel;
HC=(contains(globalDx, 'Normal')); 
allCases= ~contains(globalDx,'Fetal'); % just removing fetal cases
allCaseSansAlz= ~contains(globalDx,'Fetal') & ~contains(globalDx,'Alzheimer') ; % just removing fetal cases

overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele



%% analysis 1 looking specifically at the apparent inverse relationship between TDP burden and aSyn in alz
%additionally adding in functionality for neuronal dropout and whatnot for
%all path scores

brainAreaAtPlay='Subcortical';
pt2use=allCases;
pathCut=1;

figure
pathCogCompTestraSynTDPSpec(pt2use,brainAreaAtPlay, pathTable, cogData,pathCut   )


%% hmm ok, seems like this shows up in a lot of different brai areas, let's evaluate every one we have, let's also parcelate on basic characteristics
%additionally adding in functionality for neuronal dropout and whatnot for
%all path scores



pathCOgCompTestraSynTDpAllArea(pt2use,brainAreas, pathTable, cogData,pathCut   )
