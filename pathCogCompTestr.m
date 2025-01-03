function [] = pathCogCompTestr(pt2use,brainAreaAtPlay, pathTable, cogTable ,pathCut  )




globalDx=pathTable.GlobalDx;
snpStat=pathTable.rs199347;
pathID=pathTable.INDDID;


overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele

tauName=[brainAreaAtPlay,'Tau'];
aBetaName=[brainAreaAtPlay,'AntibodyPlaques'];
TDPName=[brainAreaAtPlay,'TDP43'];
aSynName=[brainAreaAtPlay,'aSyn'];





cogID=cogTable.ID;
cogSlope=cogTable.cogSlope;
endScorez=cogTable.endScore;

cogScoreAtPlay=nan(1,length(pathID));
endScoreAtPlay=nan(1,length(pathID));

%moving cognitive slope into the same format as path data
  for tt=1:length(pathID)
cogHold=cogSlope( cogID ==pathID(tt)   );
endScoreHold=endScorez(cogID==pathID(tt));
    if ~isempty(endScoreHold)
    cogScoreAtPlay(tt)= cogHold;
    endScoreAtPlay(tt)=endScoreHold;
    end
  end

  [asynCont,tauCont,aBetaCont,TDPCont,neuronLossCont,gliosisCont] = pathScoreGenerate(brainAreaAtPlay,pathTable);




val2Plot=cogScoreAtPlay;
SNP=snpStat;
remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)    ;



path2Plot1=asynCont;
path2Plot2=TDPCont;

cogPathSplitScatBarPlotr(val2Plot, remVals, SNP, pt2use, path2Plot1, path2Plot2, pathCut)

% title(['all Patients cog decline r=', num2str(r(2)), 'p=',num2str(p(2)), ])
legend({'High Asyn Low TDP AA','High TDP Low Asyn AA', '','','High Asyn Low TDP AG','High TDP Low Asyn AG','','',...
    'High Asyn Low TDP GG','High TDP Low Asyn GG','',''})
TDPCont=TDPCont';


%% second figure let's show  inverse relationship between ASYN and TDP 

cutLevel=2;

tdpCut=TDPCont<cutLevel ;

    overGPNMBRatLowTDP=[sum(asynCont(~remVals & overGPNMB & pt2use & tdpCut )==4  )/ sum(((~remVals & overGPNMB & pt2use & tdpCut )))*100,...
sum(asynCont(~remVals & overGPNMB & pt2use & tdpCut )==3  )/ sum(((~remVals & overGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & overGPNMB & pt2use& tdpCut )==2  )/ sum(((~remVals & overGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & overGPNMB & pt2use& tdpCut )==1  )/ sum(((~remVals & overGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & overGPNMB & pt2use& tdpCut )==0  )/ sum(((~remVals & overGPNMB & pt2use & tdpCut)))*100];
    

hetGPNMBRatLowTDP=[sum(asynCont(~remVals & (het|overGPNMB) & pt2use & tdpCut )==4  )/ sum(((~remVals & (het|overGPNMB) & pt2use & tdpCut  )))*100,...
sum(asynCont(~remVals & (het|overGPNMB) & pt2use & tdpCut )==3  )/ sum(((~remVals & (het|overGPNMB) & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & (het|overGPNMB) & pt2use & tdpCut)==2  )/ sum(((~remVals & (het|overGPNMB) & pt2use & tdpCut)))*100,...
sum(asynCont(~remVals & (het|overGPNMB) & pt2use& tdpCut )==1  )/ sum(((~remVals & (het|overGPNMB) & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & (het|overGPNMB) & pt2use& tdpCut )==0  )/ sum(((~remVals & (het|overGPNMB) & pt2use & tdpCut)))*100]; 


underGPNMBRatLowTDP=[sum(asynCont(~remVals & underGPNMB & pt2use & tdpCut )==4  )/ sum(((~remVals & underGPNMB & pt2use & tdpCut   )))*100 ,...
sum(asynCont(~remVals & underGPNMB & pt2use& tdpCut )==3  )/ sum(((~remVals & underGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & underGPNMB & pt2use& tdpCut )==2  )/ sum(((~remVals & underGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & underGPNMB & pt2use & tdpCut)==1  )/ sum(((~remVals & underGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & underGPNMB & pt2use & tdpCut)==0  )/ sum(((~remVals & underGPNMB & pt2use& tdpCut )))*100];



tdpCut=TDPCont>cutLevel ;

    overGPNMBRatHighTDP=[sum(asynCont(~remVals & overGPNMB & pt2use & tdpCut )==4  )/ sum(((~remVals & overGPNMB & pt2use & tdpCut )))*100,...
sum(asynCont(~remVals & overGPNMB & pt2use & tdpCut )==3  )/ sum(((~remVals & overGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & overGPNMB & pt2use& tdpCut )==2  )/ sum(((~remVals & overGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & overGPNMB & pt2use& tdpCut )==1  )/ sum(((~remVals & overGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & overGPNMB & pt2use& tdpCut )==0  )/ sum(((~remVals & overGPNMB & pt2use & tdpCut)))*100];
    

hetGPNMBRatHighTDP=[sum(asynCont(~remVals & (het|overGPNMB) & pt2use & tdpCut )==4  )/ sum(((~remVals & (het|overGPNMB) & pt2use & tdpCut  )))*100,...
sum(asynCont(~remVals & (het|overGPNMB) & pt2use & tdpCut )==3  )/ sum(((~remVals & (het|overGPNMB) & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & (het|overGPNMB) & pt2use & tdpCut)==2  )/ sum(((~remVals & (het|overGPNMB) & pt2use & tdpCut)))*100,...
sum(asynCont(~remVals & (het|overGPNMB) & pt2use& tdpCut )==1  )/ sum(((~remVals & (het|overGPNMB) & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & (het|overGPNMB) & pt2use& tdpCut )==0  )/ sum(((~remVals & (het|overGPNMB) & pt2use & tdpCut)))*100]; 


underGPNMBRatHighTDP=[sum(asynCont(~remVals & underGPNMB & pt2use & tdpCut )==4  )/ sum(((~remVals & underGPNMB & pt2use & tdpCut   )))*100 ,...
sum(asynCont(~remVals & underGPNMB & pt2use& tdpCut )==3  )/ sum(((~remVals & underGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & underGPNMB & pt2use& tdpCut )==2  )/ sum(((~remVals & underGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & underGPNMB & pt2use & tdpCut)==1  )/ sum(((~remVals & underGPNMB & pt2use& tdpCut )))*100,...
sum(asynCont(~remVals & underGPNMB & pt2use & tdpCut)==0  )/ sum(((~remVals & underGPNMB & pt2use& tdpCut )))*100];



figure
% subplot(1,3,1)
% b=bar([1,3], [overGPNMBRatLowTDP;overGPNMBRatHighTDP], 'stacked') ;
subplot(1,2,1)
b=bar([1,3], [hetGPNMBRatLowTDP;hetGPNMBRatHighTDP], 'stacked') ;
subplot(1,2,2)
b=bar([1,3], [underGPNMBRatLowTDP;underGPNMBRatHighTDP], 'stacked') ;











 [r,p]=corrcoef(asynCont(~remVals & pt2use ) ,TDPCont(~remVals & pt2use ),'rows','complete')
% 
 [r,p]=corrcoef(asynCont(~remVals & pt2use & (het |overGPNMB)  ) ,TDPCont(~remVals & pt2use & (het |overGPNMB) ),'rows','complete')
% 
% 
 [r,p]=corrcoef(asynCont(~remVals & pt2use & (overGPNMB)  ) ,TDPCont(~remVals & pt2use & (overGPNMB) ),'rows','complete')
% 
%

 [r,p]=corrcoef(asynCont(~remVals & pt2use & (het)  ) ,TDPCont(~remVals & pt2use & (het) ),'rows','complete')


 [r,p]=corrcoef(asynCont(~remVals & pt2use & (underGPNMB)  ) ,TDPCont(~remVals & pt2use & (underGPNMB) ),'rows','complete')


% 
% 
% 
% [r,p]=corrcoef(asynCont(~remVals & pt2use &(het |underGPNMB) & TDPCont'~=0 ) ,TDPCont(~remVals & pt2use &(het |underGPNMB)& TDPCont'~=0 ),'rows','complete')
% 
% 
% [r,p]=corrcoef(asynCont(~remVals & pt2use &underGPNMB& TDPCont'~=0 ) ,TDPCont(~remVals & pt2use &underGPNMB& TDPCont'~=0 ),'rows','complete')
% 
% 



end