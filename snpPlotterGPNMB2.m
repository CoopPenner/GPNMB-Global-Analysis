function [pOver, pUnder] = snpPlotterGPNMB2(SNP, HC,DiseaseAtPlay, DiseaseName,SNPName,numPerms,numReps,plotOn)
%Generating bar plots and running a simple chi square test of proportions

remVals=cellfun(@isempty, SNP);

SNP(remVals)=[];
HC(remVals)=[];
DiseaseAtPlay(remVals)=[];

overGPNMB= contains(SNP,'TT'); %the major allele
het=contains(SNP, 'CT');
underGPNMB=contains(SNP,'CC'); %the minor allele

[~,pOver]=prop_test( [sum(HC & overGPNMB ),sum(DiseaseAtPlay & overGPNMB ) ] , [sum(HC ), sum(DiseaseAtPlay  )] , false  );
[~,pUnder]=prop_test( [sum(HC & underGPNMB ),sum(DiseaseAtPlay & underGPNMB ) ] , [sum(HC ), sum(DiseaseAtPlay  )] , false  );

pOverPerms=nan(1,numPerms); %initializing
pUnderPerms=nan(1,numPerms);


    for dd=1:numPerms  
       hcShuf= HC(randperm(length(HC))); %HC( randi([1 length(HC)],1,length(HC))  ); %creating random array w/ replacement
       disShuf= DiseaseAtPlay(randperm(length(DiseaseAtPlay))); %DiseaseAtPlay(  randi([1 length(HC)],1,length(HC)) );



       [~,pOverPerms(dd)]=prop_test( [sum(hcShuf & overGPNMB ),sum(disShuf & overGPNMB ) ] , [sum(hcShuf  ), sum(disShuf  )] , false  );
       [~,pUnderPerms(dd)]=prop_test( [sum(hcShuf & underGPNMB ),sum(disShuf & underGPNMB ) ] , [sum(hcShuf ), sum(disShuf  )] , false  );
    end


if plotOn



b=bar([ (sum(HC & overGPNMB  )  / sum(HC ) ) , (sum(DiseaseAtPlay & overGPNMB )/sum(DiseaseAtPlay ) );...
 (sum(HC & underGPNMB )/ sum(HC ) ) , (sum(DiseaseAtPlay & underGPNMB )/sum(DiseaseAtPlay ) ) ])  ;
a=gca;
a.XTickLabel= {'AA', 'GG'};

pOverTitle=GPNMBTitleGen(pOver, pOverPerms, numReps, numPerms );
pUnderTitle=GPNMBTitleGen(pUnder, pUnderPerms, numReps, numPerms );

title([SNPName,' SNP comparison ', 'AA ', pOverTitle, ' GG ', pUnderTitle] )

legend({'HC', DiseaseName})
end

