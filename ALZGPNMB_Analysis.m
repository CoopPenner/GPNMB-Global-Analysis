%new analysis for GPNMB and ALZ status


clinData= readtable('/Users/pennerc/Documents/AllPt_GPNMB_SNP_Status.xlsx')



SNP=clinData.rs199347;
ageAtOnset=(clinData.GlobalAgeOnset);
ageAtDeath=(clinData.AgeatDeath);
disDur=ageAtDeath-ageAtOnset;
globalDiag=clinData.GlobalDx;
bioSex=clinData.Sex;

ALZ= (contains(clinData.GlobalDx, 'Alzheimer')); %& contains(clinData.GlobalDx, 'Disease Probable') ;
corticoBasal= (contains(clinData.GlobalDx, 'Corticobasal syndrome')); 
HC= (contains(clinData.GlobalDx, 'Normal'));
Parkinson= contains(clinData.GlobalDx, 'Parkinson');
ALS=contains(clinData.GlobalDx, 'Amyotrophic')
Male=contains(bioSex, 'Male');
otherDx= (ALZ+HC+Parkinson)==0;

maj= contains(SNP,'TT');
het=contains(SNP, 'CT');
mini=contains(SNP,'CC');


remVals= isnan(ageAtOnset) | isnan(disDur) 


[p,~,stats]=anovan(disDur(ALZ & ~Male ),{SNP(ALZ & ~Male)})
figure
multcompare(stats)


[n,p]=prop_test( [sum(HC & mini ),sum(ALZ & mini ) ] , [sum(HC ), sum(ALZ )] , false  )


[n,p]=prop_test( [sum(HC & mini ),sum(ALS & mini ) ] , [sum(HC ), sum(ALS )] , false  )










maj= contains(SNP,'TT');
het=contains(SNP, 'CT');
mini=contains(SNP,'CC');

d=figure; % unbelievable how irritating it is to plot this

val2Plot=disDur;

b=bar(1,nanmean(val2Plot(~remVals & maj )))  ;
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


b=bar(5,nanmean(val2Plot(~remVals & maj )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

hold on
a=scatter(rand(1, sum(maj))+.5, val2Plot(maj), 'Marker', 'o' );
a.CData=[.8 .2 .5]; 
b=scatter(rand(1, sum(het))+2.5, val2Plot(het), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, sum(maj))+4.5, val2Plot(maj), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

legend({'AA', 'GC','GG'})