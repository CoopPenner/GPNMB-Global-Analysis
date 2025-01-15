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



%change 2
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

%% second figure, comparing cog scores



%% third figure let's show  inverse relationship between ASYN and TDP 




asynRatCollect=[];

for cut=0:4

    remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)   ;


GPNMBRatHolder=[sum(asynCont(~remVals & TDPCont<=cut    )==4  )/ sum(((~remVals & TDPCont<=cut     )))*100 ,...
sum(asynCont(~remVals & TDPCont<=cut   )==3  )/ sum(((~remVals & TDPCont<=cut  )))*100,...
sum(asynCont(~remVals & TDPCont<=cut   )==2  )/ sum(((~remVals & TDPCont<=cut   )))*100,...
sum(asynCont(~remVals & TDPCont<=cut   )==1  )/ sum(((~remVals & TDPCont<=cut   )))*100,...
sum(asynCont(~remVals  & TDPCont<=cut  )==0  )/ sum(((~remVals & TDPCont<=cut   )))*100];

asynRatCollect=[ asynRatCollect;GPNMBRatHolder]; % going low to high


end

  [r,p]=corrcoef(asynCont(~remVals  ) ,TDPCont(~remVals  ),'rows','complete');


figure 

c=bar(asynRatCollect,'stacked','FaceColor','flat');

c(1).CData=[0.8,0,0];
c(2).CData=[0,.8,.7];
c(3).CData=[0.8,0.3567,0.1];
c(4).CData=[.7,0,.8];
c(5).CData=[0,0.7,0.4];


a=gca; a.XTickLabel={'0 TDP','rare TDP', '1+ TDP', '2+ TDP', '3+ TDP' };

legend({'3+ Asyn',' 2+ Asyn', ' 1+ Asyn', 'rare Asyn', '0 Asyn'}, 'FontSize',14)
title(['Amygdalar Asyn Burden as a function of TDP burden in Alzheimer"s Patients ','r=', num2str(r(2)),' p=',num2str(p(2)) ], 'FontSize',18)
ylabel('Percentage of Samples', 'FontSize',18)

%% test zone
% 
% figure
% 
% remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)   ;
% 
% histogram(asynCont(~remVals)- TDPCont(~remVals) )
% 
% 
% figure
% 
% cutLevel=2;
% 
% 
% remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)|TDPCont<cutLevel     ;
% 
% 
% 
%  overGPNMBRatLowTDP=[sum(asynCont(~remVals & overGPNMB  )==4  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==3  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==2  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==1  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==0  )/ sum(((~remVals & overGPNMB )))*100];
% 
% 
% hetGPNMBRatLowTDP=[sum(asynCont(~remVals & het  )==4  )/ sum(((~remVals & het  )))*100,...
% sum(asynCont(~remVals & het)==3  )/ sum(((~remVals & het  )))*100,...
% sum(asynCont(~remVals & het )==2  )/ sum(((~remVals & het )))*100,...
% sum(asynCont(~remVals & het )==1  )/ sum(((~remVals & het )))*100,...
% sum(asynCont(~remVals & het )==0  )/ sum(((~remVals & het  )))*100]; 
% 
% 
% underGPNMBRatLowTDP=[sum(asynCont(~remVals & underGPNMB  )==4  )/ sum(((~remVals & underGPNMB   )))*100 ,...
% sum(asynCont(~remVals & underGPNMB  )==3  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB  )==2  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB )==1  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB )==0  )/ sum(((~remVals & underGPNMB )))*100];
% 
% 
% 
% 
% remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)|TDPCont>cutLevel     ;
% 
% 
%     overGPNMBRatHighTDP=[sum(asynCont(~remVals & overGPNMB  )==4  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==3  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==2  )/ sum(((~remVals & overGPNMB )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==1  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==0  )/ sum(((~remVals & overGPNMB  )))*100];
% 
% 
% hetGPNMBRatHighTDP=[sum(asynCont(~remVals & (het)   )==4  )/ sum(((~remVals & (het)   )))*100,...
% sum(asynCont(~remVals & (het)  )==3  )/ sum(((~remVals & (het)  )))*100,...
% sum(asynCont(~remVals & (het) )==2  )/ sum(((~remVals & (het) )))*100,...
% sum(asynCont(~remVals & (het)  )==1  )/ sum(((~remVals & (het)  )))*100,...
% sum(asynCont(~remVals & (het)  )==0  )/ sum(((~remVals & (het)  )))*100]; 
% 
% 
% underGPNMBRatHighTDP=[sum(asynCont(~remVals & underGPNMB  )==4  )/ sum(((~remVals & underGPNMB    )))*100 ,...
% sum(asynCont(~remVals & underGPNMB  )==3  )/ sum(((~remVals & underGPNMB )))*100,...
% sum(asynCont(~remVals & underGPNMB  )==2  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB  )==1  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB  )==0  )/ sum(((~remVals & underGPNMB  )))*100];
% 
% 
% 
% figure
%  subplot(1,3,1)
%  b=bar([1,3], [overGPNMBRatLowTDP;overGPNMBRatHighTDP], 'stacked') ;
% subplot(1,3,2)
% b=bar([1,3], [hetGPNMBRatLowTDP;hetGPNMBRatHighTDP], 'stacked') ;
% subplot(1,3,3)
% b=bar([1,3], [underGPNMBRatLowTDP;underGPNMBRatHighTDP], 'stacked') ;
% 
% 
% 
% 
% 
% remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)     ;
% 
% 
% 
% 
% 
% 
% 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (TDPCont+ asynCont'~=0) ) ,TDPCont(~remVals & pt2use & (TDPCont+ asynCont'~=0) ),'rows','complete')
% % 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (het |overGPNMB)  ) ,TDPCont(~remVals & pt2use & (het |overGPNMB) ),'rows','complete')
% % 
% % 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (overGPNMB)  ) ,TDPCont(~remVals & pt2use & (overGPNMB) ),'rows','complete')
% % 
% %
% 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (het)  ) ,TDPCont(~remVals & pt2use & (het) ),'rows','complete')
% 
% 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (underGPNMB)  ) ,TDPCont(~remVals & pt2use & (underGPNMB) ),'rows','complete')
% 
% 
% % 
% % 
% % 
% % [r,p]=corrcoef(asynCont(~remVals & pt2use &(het |underGPNMB) & TDPCont'~=0 ) ,TDPCont(~remVals & pt2use &(het |underGPNMB)& TDPCont'~=0 ),'rows','complete')
% % 
% % 
% % [r,p]=corrcoef(asynCont(~remVals & pt2use &underGPNMB& TDPCont'~=0 ) ,TDPCont(~remVals & pt2use &underGPNMB& TDPCont'~=0 ),'rows','complete')
% % 
% % 
% 


end