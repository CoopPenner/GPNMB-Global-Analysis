%quick SNP aggregator



clinDataTot=readtable('/Users/pennerc/Documents/GPNMB multi SNP eval.xlsx');


Parkinson= (contains(clinDataTot.GlobalDx, 'Parkinson'));
Alzheimer= (contains(clinDataTot.GlobalDx, 'Alzheimer'));
ALS= (contains(clinDataTot.GlobalDx, 'Amyotrophic'));
DemLewy= (contains(clinDataTot.GlobalDx, 'Dementia with Lewy Bodies')); 
MCI= (contains(clinDataTot.GlobalDx, 'Mild cognitive impairment')); 
corticoBasal= (contains(clinDataTot.GlobalDx, 'Corticobasal syndrome')); 
bvFTD= (contains(clinDataTot.GlobalDx, 'bvFTD-FTLD')); 
PPA= (contains(clinDataTot.GlobalDx, 'PPA')); 
supraNuc= (contains(clinDataTot.GlobalDx, 'Progressive supranuclear palsy')); 
neuroPanel= Alzheimer + ALS + DemLewy  + corticoBasal +bvFTD + PPA +supraNuc;
ParkinsonianPanel= Parkinson +corticoBasal +supraNuc;
ParkinsonianDem=corticoBasal+DemLewy;
other= ~neuroPanel;
HC=(contains(clinDataTot.GlobalDx, 'Normal')); 



contSNP1=clinDataTot.rs17137124; %putative risk for PPA
contSNP2=clinDataTot.rs10139154; %putative risk for ALS
contSNP3=clinDataTot.rs10260404; %putative risk for ALS
contSNP4=clinDataTot.rs10797576; %putative risk for Parkinsons
contSNP5=clinDataTot.rs11136000; %putative risk for ALZ
contSNP6=clinDataTot.rs11158026; %putative risk for Parkinsons
contSNP7=clinDataTot.rs11218343; %putative risk for ALZ
contSNP8=clinDataTot.rs11767557; %putative risk for ALZ
contSNP9=clinDataTot.rs13048019; %putative risk for ALS
GPNMBSNP=clinDataTot.rs199347; %:)


contSNP5=clinDataTot.rs11136000; %putative risk for ALZ
contSNP7=clinDataTot.rs11218343; %putative risk for ALZ
contSNP8=clinDataTot.rs11767557; %putative risk for ALZ
GPNMBSNP=clinDataTot.rs199347; %:)



figure

snpPlotterGPNMB2(contSNP1, HC,Alzheimer, 'Alzheimer',1000)

figure
snpPlotterGPNMB2(contSNP2, HC,Alzheimer, 'Alzheimer',1000)

figure
snpPlotterGPNMB2(contSNP3, HC,Alzheimer, 'Alzheimer',1000)

figure
snpPlotterGPNMB2(contSNP4, HC,Alzheimer, 'Alzheimer',1000)

figure
snpPlotterGPNMB2(contSNP5, HC,Alzheimer, 'Alzheimer',1000)



figure
snpPlotterGPNMB2(contSNP6, HC,Alzheimer, 'Alzheimer',1000)


figure
snpPlotterGPNMB2(contSNP7, HC,Alzheimer, 'Alzheimer',1000)


figure
snpPlotterGPNMB2(contSNP8, HC,Alzheimer, 'Alzheimer',1000)
figure
snpPlotterGPNMB2(contSNP9, HC,Alzheimer, 'Alzheimer',1000)

figure
snpPlotterGPNMB2(contSNP10, HC,Alzheimer, 'Alzheimer',1000)



GPNMBSNP=clinDataTot.rs199347; %:)

figure
subplot(2,2,1)
[pOver, pUnder]=snpPlotterGPNMB2(GPNMBSNP, HC,Alzheimer, 'Alzheimer"s', 'rs199347', 1000,'true');

subplot(2,2,2)
[pOver, pUnder]=snpPlotterGPNMB2(GPNMBSNP, HC,Alzheimer, 'Alzheimer"s', 'rs199347', 1000,'true');


subplot(2,2,3)
[pOver, pUnder]=snpPlotterGPNMB2(GPNMBSNP, HC,Alzheimer, 'Alzheimer"s', 'rs199347', 1000,'true');


subplot(2,2,4)
[pOver, pUnder]=snpPlotterGPNMB2(GPNMBSNP, HC,Alzheimer, 'Alzheimer"s', 'rs199347', 1000,'true');









contSNP5=clinDataTot.rs11136000; %putative risk for ALZ
contSNP7=clinDataTot.rs11218343; %putative risk for ALZ
contSNP8=clinDataTot.rs11767557; %putative risk for ALZ

GPNMBSNP=clinDataTot.rs199347; %:)


figure

subplot(2,2,1)
[pOver, pUnder]=snpPlotterGPNMB2(GPNMBSNP, HC,Alzheimer, 'Alzheimer"s', 'rs199347',10000,4,'false');


subplot(2,2,2)
[pOver, pUnder]=snpPlotterGPNMB2(contSNP5, HC,Alzheimer, 'Alzheimer"s', 'rs11136000',10000,4,'false');


subplot(2,2,3)
[pOver, pUnder]=snpPlotterGPNMB2(contSNP7, HC,Alzheimer, 'Alzheimer"s', 'rs11218343',10000,4,'false');


subplot(2,2,4)
[pOver, pUnder]=snpPlotterGPNMB2(contSNP8, HC,Alzheimer, 'Alzheimer"s', 'rs11767557',10000,4,'false');


