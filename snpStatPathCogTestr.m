function [] = snpStatPathCogTestr(pt2use,brainAreaAtPlay, pathTable, cogTable, pathType,secondPath,pathType2,secondPathCut   )
%our third analysis, herein I evaluate how the relationship between
%cognitive decline and pathology is mediated by SNP Status






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




if pathType==1
    path2Plot=asynCont;
    
    elseif pathType==2
    path2Plot=aBetaCont;
    elseif pathType==3
    path2Plot=tauCont;
    elseif pathType==4
    path2Plot=TDPCont;
    elseif pathType==5
    path2Plot=neuronLossCont;
    elseif pathType==6
    path2Plot=gliosisCont;
end



if secondPath && pathType2==1
secondPath2Plot=asynCont;
elseif secondPath && pathType2==2
secondPath2Plot=aBetaCont;

elseif secondPath && pathType2==3
secondPath2Plot=tauCont;

elseif secondPath && pathType2==4
secondPath2Plot=TDPCont;

elseif secondPath && pathType2==5
secondPath2Plot=neuronLossCont;

elseif secondPath && pathType2==6
secondPath2Plot=gliosisCont;

elseif ~secondPath
    secondPath2Plot=nan;
end

val2Plot=cogScoreAtPlay;
SNP=snpStat;
remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)     ;

%% plotting


figure



subplot(1,3,1)
cogPathScatBarPlotr(val2Plot, remVals, overGPNMB,pt2use,path2Plot,secondPath,secondPath2Plot,secondPathCut)
ylabel('cognitive decline slope')
legend({'3+','rare-2+', 'none'})
[r,p]=corrcoef(val2Plot(pt2use& (overGPNMB)& ~remVals ), path2Plot(pt2use& (overGPNMB)& ~remVals), 'rows','complete');
title(['AA (over Production)  SNP & cog decline r=', num2str(r(2)), ' p=',num2str(p(2)), ],  'FontSize', 15)
ylabel('Slope of Cognitive Decline','FontSize',15)


subplot(1,3,2)
cogPathScatBarPlotr(val2Plot, remVals, het,pt2use,path2Plot,secondPath,secondPath2Plot,secondPathCut)
[r,p]=corrcoef(val2Plot(pt2use& (het)& ~remVals ), path2Plot(pt2use& (het)& ~remVals), 'rows','complete');
title(['AG  SNP & cog decline r=', num2str(r(2)), ' p=',num2str(p(2)), ],  'FontSize', 15)
legend({'3+','rare-2+', 'none'},'FontSize',14)


subplot(1,3,3)
cogPathScatBarPlotr(val2Plot, remVals, underGPNMB,pt2use,path2Plot,secondPath,secondPath2Plot,secondPathCut)

[r,p]=corrcoef(val2Plot(pt2use & (underGPNMB) & ~remVals), path2Plot(pt2use& (underGPNMB)& ~remVals), 'rows','complete');

title(['GG (under Production) SNP & cog decline r=', num2str(r(2)), ' p=',num2str(p(2)), ],  'FontSize', 15)
legend({'3+','rare-2+', 'none'},'FontSize',14)



figure

cogPathScatBarPlotr(val2Plot, remVals, true(1,length(underGPNMB))' ,pt2use,path2Plot,secondPath,secondPath2Plot,secondPathCut)

[r,p]=corrcoef(val2Plot(pt2use & ~remVals ), path2Plot(pt2use & ~remVals), 'rows','complete');

title(['all Patients cog decline r=', num2str(r(2)), ' p=',num2str(p(2)), ], 'FontSize', 16)
legend({'3+','rare-2+', 'none'},'FontSize',14)
ylabel('Slope of Cognitive Decline','FontSize',15)








%% experimental section 



% figure
% subplot(1,3,1)
% histogram(asynCont(overGPNMB' & asynCont>0 & TDPCont>0 )-TDPCont(overGPNMB' & asynCont>0 & TDPCont>0 ))
% subplot(1,3,2)
% histogram(asynCont(het' & asynCont>0 & TDPCont>0 )-TDPCont(het' & asynCont>0 & TDPCont>0 ))
% 
% 
% subplot(1,3,3)
% histogram(asynCont(underGPNMB' & asynCont>0 & TDPCont>0 )-TDPCont(underGPNMB' & asynCont>0 & TDPCont>0 ))
% 


end


