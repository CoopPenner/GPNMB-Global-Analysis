function [outputArg1,outputArg2] = snpPlotterGPNMB2(SNP, HC,DiseaseAtPlay, DiseaseName)
%Generating bar plots and running a simple chi square test of proportions

remVals=cellfun(@isempty, SNP);

overGPNMB= contains(SNP,'TT'); %the major allele
het=contains(SNP, 'CT');
underGPNMB=contains(SNP,'CC'); %the minor allele

b=bar([ (sum(HC & overGPNMB  )  / sum(HC & ~remVals) ) , (sum(DiseaseAtPlay & overGPNMB )/sum(DiseaseAtPlay & ~remVals) );...
 (sum(HC & underGPNMB )/ sum(HC & ~remVals) ) , (sum(DiseaseAtPlay & underGPNMB )/sum(DiseaseAtPlay & ~remVals) ) ])  ;
a=gca;
a.XTickLabel= {'AA (over)', 'GG (under)'};
[~,pOver]=prop_test( [sum(HC & overGPNMB ),sum(DiseaseAtPlay & overGPNMB ) ] , [sum(HC & ~remVals ), sum(DiseaseAtPlay & ~remVals )] , false  );
[~,pUnder]=prop_test( [sum(HC & underGPNMB ),sum(DiseaseAtPlay & underGPNMB ) ] , [sum(HC & ~remVals ), sum(DiseaseAtPlay & ~remVals )] , false  );
title(['rs199347 SNP comparison ', 'AA p=', num2str(pOver), ' GG p=', num2str(pUnder)] )

legend({'HC', DiseaseName})





end