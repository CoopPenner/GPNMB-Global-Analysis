function [] = snpStatPathCogTestr(pt2Plot,brainAreaAtPlay, pathTable, cogTable,pathCut   )
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



  path2Plot=asynCont;

secondaryCut=pathCut;
val2Plot=cogScoreAtPlay;
SNP=snpStat;
remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)   ;



figure


subplot(1,3,1)

b=bar(1,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot==4 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel('end score')
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot==3 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot==2 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];



b=bar(7,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot==1 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.1 0.5 .8];



b=bar(9,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot==0 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.9 0.4 .2];

hold on
scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot & asynCont'==4 ))+.5, val2Plot(~remVals & overGPNMB & pt2Plot & asynCont'==4), 'Marker', 'o', 'MarkerEdgeColor', [.8 .2 .5] );


scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot & asynCont'==3))+2.5, val2Plot(~remVals & overGPNMB & pt2Plot & asynCont'==3), 'Marker', 'o', 'MarkerEdgeColor', [0 0.7 .25] );


scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot & asynCont'==2))+4.5, val2Plot(~remVals & overGPNMB & pt2Plot & asynCont'==2), 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.1 .6]);

scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot & asynCont'==1))+6.5, val2Plot(~remVals & overGPNMB & pt2Plot & asynCont'==1), 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.1 .6]);


scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot & asynCont'==0))+8.5, val2Plot(~remVals & overGPNMB & pt2Plot & asynCont'==0), 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.1 .6]);



legend({'AA (over Production)', 'GC','GG (under Production)'})




% figure
% 
% 
% scatter(asynCont(pt2Plot &overGPNMB ), val2Plot(pt2Plot &overGPNMB), 'MarkerFaceColor','r')
% lsline
% hold on
% scatter(asynCont(pt2Plot &underGPNMB ), val2Plot(pt2Plot &underGPNMB), 'MarkerFaceColor','c')
% lsline
% scatter(asynCont(pt2Plot &het ), val2Plot(pt2Plot &het), 'MarkerFaceColor','g')
% lsline
% 
% [r,p]=corrcoef(asynCont(pt2Plot& (het | underGPNMB) ), val2Plot(pt2Plot& (het | underGPNMB)), 'rows','complete')