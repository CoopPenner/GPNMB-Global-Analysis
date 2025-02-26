addpath('/Users/pennerc/Documents/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis')



%ClinDatTot= readtable('/Users/pennerc/Documents/AllPt_GPNMB_SNP_Status.xlsx');
clinDataTot=readtable('/Volumes/PC60/InqueryDatasets/DetailedBasicInfo_allPatients.xlsx');

%Here I simply walk through and plot all diseases within INDD related to
%GPNMB status 

GPNMBSNP=clinDataTot.rs199347;
ID=clinDataTot.INDDID;
ageAtOnset=(clinDataTot.GlobalAgeOnset);
ageAtDeath=(clinDataTot.AgeatDeath);
disDur=ageAtDeath-ageAtOnset;
globalDx=clinDataTot.GlobalDx;
pdDX=clinDataTot.PDCDx;
CNDRDX=clinDataTot.CNDRDx;
alsDX=clinDataTot.ALSDx;
adDX=clinDataTot.ADCDx;
snpStat=clinDataTot.rs199347;
pathID=clinDataTot.INDDID;
geneCarrier=~cellfun(@isempty,clinDataTot.Mutation_Summary);
SOD1=contains(clinDataTot.Mutation_Summary,'SOD1');



Parkinson= (contains(globalDx, 'Parkinson') |  contains(CNDRDX, 'Parkinson')    );
Alzheimer= (contains(globalDx, 'Alzheimer')  | contains(adDX, 'Alzheimer') | contains(CNDRDX, 'Alzheimer')    );
ALS= (contains(globalDx, 'Amyotrophic'));
DemLewy= (contains(globalDx, 'Dementia with Lewy Bodies')|  contains(CNDRDX, 'Dementia with Lewy Bodies')    );
MCI= (contains(globalDx, 'Mild cognitive impairment')); 
corticoBasal= (contains(globalDx, 'Corticobasal syndrome') | contains(CNDRDX, 'Corticobasal syndrome')    );
bvFTD= (contains(globalDx, 'bvFTD-FTLD')); 
PPA= (contains(globalDx, 'PPA')); 
supraNuc= (contains(globalDx, 'Progressive supranuclear palsy')|  contains(CNDRDX, 'Progressive supranuclear palsy')    );
neuroPanel= Alzheimer | DemLewy  | Parkinson| corticoBasal ;
ParkinsonianPanel= Parkinson |corticoBasal |supraNuc |DemLewy;
ofInterest=ALS|corticoBasal;
ParkinsonianDem=corticoBasal|DemLewy;
other= ~neuroPanel;
HC=(contains(globalDx, 'Normal')); 

normGen= ALS |Alzheimer |DemLewy | MCI | bvFTD | PPA | supraNuc | Parkinson;



overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele



sum(      ~remVals & Parkinson )+ sum(      ~remVals & ALS ) + sum(      ~remVals & Alzheimer )+ sum(      ~remVals & MCI )

sum(~remVals &HC)


remVals= cellfun(@isempty,snpStat);


%% running

% Parkinson's no diff re age of onset or death but the expected diff in the
% allele freq
figure
subplot(1,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, Parkinson, 'Disease Duration', 'Parkinson"s Disease')
subplot(1,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, Parkinson, 'Age At Onset', 'Parkinson"s Disease')
% subplot(2,2,3:4)
% snpPlotterGPNMB2(GPNMBSNP, HC,Parkinson, 'Parkinson"s')




%ALZ nothin
figure
subplot(1,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, Alzheimer, 'Disease Duration', 'Alzheimer"s Disease')
subplot(1,2,2)
 snpPlotterGPNMB(GPNMBSNP,ageAtOnset, Alzheimer, 'Age At Onset', 'Alzheimer"s Disease')
% subplot(2,2,3:4)

figure
% this step runs a permutation 
[pOver, pUnder]=snpPlotterGPNMB2(GPNMBSNP, Parkinson,Alzheimer, 'Alzheimer"s', 'rs199347', 10000,1,'true');


%really clear and interesting trend for age at onset and disease duration

figure
subplot(1,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, ALS, 'Disease Duration', 'ALS')
subplot(1,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, ALS, 'Age At Onset', 'ALS')
% subplot(2,2,3:4)
% [pOver, pUnder]=snpPlotterGPNMB2(GPNMBSNP, HC,ALS, 'ALS', 'rs199347', 10000,1,'true');
% 


%I think there is a clear trend here, we are just very underpowered over
%production has a faster disease progression (ie onset to death)

figure
subplot(1,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, DemLewy, 'Disease Duration', 'Dementia with Lewy Bodies')
subplot(1,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, DemLewy, 'Age At Onset', 'Dementia with Lewy Bodies')

% subplot(2,2,3:4)
% snpPlotterGPNMB2(GPNMBSNP, HC,DemLewy, 'Dementia with Lewy Bodies',1000)
% 

figure
subplot(1,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, MCI, 'Disease Duration', 'MCI')
subplot(1,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, MCI, 'Age At Onset', 'MCI')
% subplot(2,2,3:4)
% snpPlotterGPNMB2(GPNMBSNP, HC,MCI, 'MCI')



%incredibly interesting clear trend in both directions, exact opposite to
%what is seen in ALS
figure
subplot(1,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, corticoBasal, 'Disease Duration', 'corticoBasal Syndrome')
subplot(1,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, corticoBasal, 'Age At Onset', 'corticoBasal Syndrome ')

figure
[pOver, pUnder]=snpPlotterGPNMB2(GPNMBSNP, normGen,corticoBasal, 'CorticoBasal', 'rs199347', 10000,1,'true');

%nothing

figure
subplot(2,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, bvFTD, 'Disease Duration', 'bvFTD')
subplot(2,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, bvFTD, 'Age At Onset', 'bvFTD')
% subplot(2,2,3:4)
% snpPlotterGPNMB2(GPNMBSNP, HC,bvFTD, 'bvFTD')


figure
subplot(2,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, supraNuc, 'Disease Duration', 'Supranuclear Palsy')
subplot(2,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, supraNuc, 'Age At Onset', 'Supranuclear Palsy')
subplot(2,2,3:4)
snpPlotterGPNMB2(GPNMBSNP, HC,supraNuc, 'Supranuclear Palsy')



figure
subplot(2,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, PPA, 'Disease Duration', 'PPA')
subplot(2,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, PPA, 'Age At Onset', 'PPA')
subplot(2,2,3:4)
snpPlotterGPNMB2(GPNMBSNP, HC,PPA, 'PPA')




figure
subplot(1,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, ageOnsetTrend, 'Disease Duration', 'PPA')
subplot(1,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, ageOnsetTrend, 'Age At Onset', 'PPA')
% subplot(2,2,3:4)
% snpPlotterGPNMB2(GPNMBSNP, HC,neuroPanel, 'Collective Neuro Disease')
% 


figure
subplot(2,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, ParkinsonianPanel, 'Disease Duration', 'Parkinsonian Syndromes')
subplot(2,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, ParkinsonianPanel, 'Age At Onset', 'Parkinsonian Syndromes')
subplot(2,2,3:4)
snpPlotterGPNMB2(GPNMBSNP, HC,ParkinsonianPanel, 'Parkinsonian Syndromes')



figure
subplot(2,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, ParkinsonianDem, 'Disease Duration', 'ParkinsonianDem')
subplot(2,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, ParkinsonianDem, 'Age At Onset', 'ParkinsonianDem')
subplot(2,2,3:4)
snpPlotterGPNMB2(GPNMBSNP, HC,ParkinsonianDem, 'ParkinsonianDem')




figure
subplot(2,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, other, 'Disease Duration', 'PPA')
subplot(2,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, other, 'Age At Onset', 'PPA')
subplot(2,2,3:4)
snpPlotterGPNMB2(GPNMBSNP, HC,other, 'Collective Neuro Disease')




%older vs younger hc

figure
snpPlotterGPNMB(GPNMBSNP,ageAtDeath, HC, 'age at Death', 'Healthy Control')








figure
subplot(2,2,1)
snpPlotterGPNMB(GPNMBSNP,disDur, corticoBasal, 'Disease Duration', 'corticoBasal Syndrome')
subplot(2,2,2)
snpPlotterGPNMB(GPNMBSNP,ageAtOnset, corticoBasal, 'Age At Onset', 'corticoBasal Syndrome ')
subplot(2,2,3:4)
snpPlotterGPNMB2(GPNMBSNP, HC,corticoBasal, 'corticoBasal Syndrome')








figure
subplot(2,2,1)
snpPlotterGPNMB(contSNP1,disDur, corticoBasal, 'Disease Duration', 'corticoBasal Syndrome')
subplot(2,2,2)
snpPlotterGPNMB(contSNP1,ageAtOnset, corticoBasal, 'Age At Onset', 'corticoBasal Syndrome ')
subplot(2,2,3:4)
snpPlotterGPNMB2(contSNP1, HC,corticoBasal, 'corticoBasal Syndrome')




figure
subplot(2,2,1)
snpPlotterGPNMB(contSNP3,disDur, corticoBasal, 'Disease Duration', 'corticoBasal Syndrome')
subplot(2,2,2)
snpPlotterGPNMB(contSNP3,ageAtOnset, corticoBasal, 'Age At Onset', 'corticoBasal Syndrome ')
subplot(2,2,3:4)
snpPlotterGPNMB2(contSNP3, HC,corticoBasal, 'corticoBasal Syndrome')





figure
subplot(2,2,1)
snpPlotterGPNMB(contSNP4,disDur, corticoBasal, 'Disease Duration', 'corticoBasal Syndrome')
subplot(2,2,2)
snpPlotterGPNMB(contSNP4,ageAtOnset, corticoBasal, 'Age At Onset', 'corticoBasal Syndrome ')
subplot(2,2,3:4)
snpPlotterGPNMB2(contSNP4, HC,corticoBasal, 'corticoBasal Syndrome')