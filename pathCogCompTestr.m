function [] = pathCogCompTestr(pt2use,brainAreaAtPlay, pathTable, cogTable   )




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



pathCut=2;
path2Plot1=asynCont;
path2Plot2=TDPCont;

cogPathSplitScatBarPlotr(val2Plot, remVals, SNP, pt2use, path2Plot1, path2Plot2, pathCut)

% title(['all Patients cog decline r=', num2str(r(2)), 'p=',num2str(p(2)), ])
legend({'3+','2+','1+','rare', 'none'})



highTDP=nanmean(val2Plot(overGPNMB & asynCont' <2 & TDPCont'>2  & pt2use ))



highAsyn=nanmean(val2Plot(overGPNMB & asynCont' >2 & TDPCont'<2  & pt2use ))







nanmean(val2Plot(het & asynCont' <1 & TDPCont'>1  & pt2use ))



nanmean(val2Plot(het & asynCont' >1 & TDPCont'<1  & pt2use ))




end