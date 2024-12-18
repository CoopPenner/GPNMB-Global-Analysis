function [] = cogPathScatBarPlotr(val2Plot, remVals, snpMark,pt2Plot,path2Plot)
%plotting function for scatter bars




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

