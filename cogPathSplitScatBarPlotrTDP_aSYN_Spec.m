function [] = cogPathSplitScatBarPlotrTDP_aSYN_Spec(val2Plot, remVals,pt2Plot,path2Plot1, path2Plot2, pathCut)







b=bar(1,nanmean(val2Plot(~remVals & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.15;
b.CData(1,:) = [.8 .2 .5]; 
a=gca; a.XTickLabel=[];

hold on


b=bar(3.2,nanmean(val2Plot(~remVals & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut ))) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.15;
b.CData(1,:) = [.5 .2 .7]; 
a=gca; a.XTickLabel=[];



b=bar(5.4,nanmean(val2Plot(~remVals &  pt2Plot & path2Plot1'<pathCut & path2Plot2'< pathCut ))) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.15;
b.CData(1,:) = [.2 .8 .5]; 
a=gca; a.XTickLabel=[];




scatter(   rand(1, sum(~remVals &  pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut)  )+.5,...
    val2Plot(~remVals &  pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut) , 'Marker', 'o', 'MarkerEdgeColor', [.8 .2 .5], 'MarkerFaceColor', [.8 .2 .5],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2 );


scatter(   rand(1, sum(~remVals & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut)  )+2.65,...
    val2Plot(~remVals &  pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut) , 'Marker', 'o', 'MarkerEdgeColor', [.5 .2 .7],  'MarkerFaceColor', [.5 .2 .7],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2  );

scatter(   rand(1, sum(~remVals &  pt2Plot & path2Plot1'==0 & path2Plot2'==0)  )+4.8,...
    val2Plot(~remVals &  pt2Plot & path2Plot1'==0 & path2Plot2'==0) , 'Marker', 'o', 'MarkerEdgeColor', [.2 .8 .5],  'MarkerFaceColor', [.2 .8 .5],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2  );



