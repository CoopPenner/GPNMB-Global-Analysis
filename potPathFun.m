b=1;



[n,p]=corrcoef(TDPCont(~remVals & (het|overGPNMB) & pt2use & (TDPCont'>0  )  ),asynCont(~remVals & pt2use &(het|overGPNMB)   & (TDPCont'>0  ) ),'rows','complete')


[n,p]=corrcoef(asynCont(~remVals  & pt2use & overGPNMB ) ,neuronDropCont(~remVals & pt2use & overGPNMB    ) ,'rows','complete')


[n,p]=corrcoef(TDPCont(~remVals  & pt2use & het  ) ,neuronDropCont(~remVals & pt2use & het ) ,'rows','complete')





[r,p]=corrcoef(val2Plot(~remVals & (het|overGPNMB) & pt2use &TDPCont'>0 ) ,TDPCont(~remVals & (het|overGPNMB) & pt2use &TDPCont'>0   )  )



[r,p]=corrcoef(asynCont(~remVals  & pt2use  ) ,val2Plot(~remVals  & pt2use    ) ,'rows','complete' )


nanmean(val2Plot(overGPNMB & TDPCont'<2 & asynCont'>2)  )


nanmean(val2Plot( overGPNMB &  TDPCont'<2 & asynCont'>2)  )


nanmean(val2Plot(overGPNMB &  TDPCont'>2 & asynCont'<2)  )


nanmean(val2Plot(overGPNMB &  TDPCont'>1 & asynCont'<1)  )



nanmean(val2Plot(overGPNMB &  TDPCont'<1 & asynCont'>1)  )


nanmean(val2Plot(overGPNMB &  TDPCont'>1 & asynCont'<1)  )


figure




remVals= isnan(val2Plot)' | cellfun(@isempty,SNP) |   ;


cogPathScatBarPlotr(val2Plot, remVals, overGPNMB,pt2use,path2Plot)

[r,p]=corrcoef(val2Plot(pt2use & (underGPNMB) & ~remVals), path2Plot(pt2use& (underGPNMB)& ~remVals), 'rows','complete');

title(['GG (under Production) SNP & cog decline r=', num2str(r(2)), ' p=',num2str(p(2)), ])
legend({'3+','2+','1+','rare', 'none'})









b=bar(1,nanmean(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==4 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==3 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==2 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];



b=bar(7,nanmean(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==1 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.1 0.5 .8];



b=bar(9,nanmean(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==0 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.9 0.4 .2];

hold on
scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==4 ))+.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==4), 'Marker', 'o', 'MarkerEdgeColor', [.8 .2 .5] );


scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==3))+2.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==3), 'Marker', 'o', 'MarkerEdgeColor', [0 0.7 .25] );


scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==2))+4.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==2), 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.1 .6]);

scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==1))+6.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==1), 'Marker', 'o', 'MarkerEdgeColor', [0.1 0.5 .8]);


scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==0))+8.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==0), 'Marker', 'o', 'MarkerEdgeColor', [0.9 0.4 .2]);

