add2path('')



ClinDatTot= readtable('/Users/pennerc/Documents/AllPt_GPNMB_SNP_Status.xlsx')


%Here I simply walk through and plot all diseases within INDD related to
%GPNMB status 

SNP=clinData.rs199347;
ageAtOnset=(clinData.GlobalAgeOnset);
ageAtDeath=(clinData.AgeatDeath);
disDur=ageAtDeath-ageAtOnset;


Parkinson= (contains(clinData.GlobalDx, 'Parkinson'));
Alzheimer= (contains(clinData.GlobalDx, 'Alzheimer'));
ALS= (contains(clinData.GlobalDx, 'Amyotrophic'));
DemLewy= (contains(clinData.GlobalDx, 'Dementia with Lewy Bodies')); 
MCI= (contains(clinData.GlobalDx, 'Mild cognitive impairment')); 
corticoBasal= (contains(clinData.GlobalDx, 'Corticobasal syndrome')); 
bvFTD= (contains(clinData.GlobalDx, 'bvFTD-FTLD')); 
neuroPanel= Parkinson+ Alzheimer + ALS + DemLewy + MCI + corticoBasal +bvFTD;
other= ~neuroPanel;


snpPlotterGPNMB(SNP,disDur, Parkinson, 'Disease Duration', 'Parkinson"s Disease')



remVals= isnan(ageAtOnset) | isnan(disDur)  


[p,~,stats]=anovan(disDur(~remVals & corticoBasal  ),{SNP(~remVals & corticoBasal  )}, 'display', 'off')

figure
multcompare(stats)







overGPNMB= contains(SNP,'TT'); %the major allele
het=contains(SNP, 'CT');
underGPNMB=contains(SNP,'CC'); %the minor allele

d=figure; % unbelievable how irritating it is to plot this

val2Plot=disDur;

b=bar(1,nanmean(val2Plot(~remVals & overGPNMB )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel('Age at Onset (years')
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Plot(~remVals & het )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Plot(~remVals & underGPNMB )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

hold on
a=scatter(rand(1, sum(overGPNMB))+.5, val2Plot(overGPNMB), 'Marker', 'o' );
a.CData=[.8 .2 .5]; 
b=scatter(rand(1, sum(het))+2.5, val2Plot(het), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, sum(underGPNMB))+4.5, val2Plot(underGPNMB), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

legend({'AA', 'GC','GG'})