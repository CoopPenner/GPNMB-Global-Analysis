function [] = snpStatPathCogTestr(pt2Plot,brainAreaAtPlay, pathTable, cogTable, pathType   )
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

  [asynCont,tauCont,aBetaCont,TDPCont] = pathScoreGenerate(brainAreaAtPlay,pathTable);




if pathType==1
    path2Plot=asynCont;
    elseif pathType==2
    path2Plot=aBetaCont;
    elseif pathType==3
    path2Plot=tauCont;
    elseif pathType==4
    path2Plot=TDPCont;
end





val2Plot=cogScoreAtPlay;
SNP=snpStat;
remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)    ;

%% plotting
figure

subplot(1,3,1)
cogPathScatBarPlotr(val2Plot, remVals, overGPNMB,pt2Plot,path2Plot)
ylabel('cognitive decline slope')
legend({'3+','2+','1+','rare', 'none'})
[r,p]=corrcoef(val2Plot(pt2Plot& (overGPNMB)& ~remVals ), path2Plot(pt2Plot& (overGPNMB)& ~remVals), 'rows','complete');
title(['AA (over Production)  SNP & cog decline r=', num2str(r(2)), ' p=',num2str(p(2)), ])


subplot(1,3,2)
cogPathScatBarPlotr(val2Plot, remVals, het,pt2Plot,path2Plot)
[r,p]=corrcoef(val2Plot(pt2Plot& (het)& ~remVals ), path2Plot(pt2Plot& (het)& ~remVals), 'rows','complete');
title(['AG  SNP & cog decline r=', num2str(r(2)), ' p=',num2str(p(2)), ])
legend({'3+','2+','1+','rare', 'none'})



subplot(1,3,3)
cogPathScatBarPlotr(val2Plot, remVals, underGPNMB,pt2Plot,path2Plot)

[r,p]=corrcoef(val2Plot(pt2Plot & (underGPNMB) & ~remVals), path2Plot(pt2Plot& (underGPNMB)& ~remVals), 'rows','complete');

title(['GG (under Production) SNP & cog decline r=', num2str(r(2)), ' p=',num2str(p(2)), ])
legend({'3+','2+','1+','rare', 'none'})



figure

cogPathScatBarPlotr(val2Plot, remVals, true(1,length(underGPNMB))' ,pt2Plot,path2Plot)

[r,p]=corrcoef(val2Plot(pt2Plot & ~remVals ), path2Plot(pt2Plot & ~remVals), 'rows','complete');

title(['all Patients cog decline r=', num2str(r(2)), 'p=',num2str(p(2)), ])
legend({'3+','2+','1+','rare', 'none'})









end


