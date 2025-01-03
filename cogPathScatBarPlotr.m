function [] = cogPathScatBarPlotr(val2Plot, remVals,snpMark,pt2Plot,path2Plot,secondPath,secondPath2Plot,secondPathCut)
%plotting function for scatter bars



if ~secondPath
b=bar(1,nanmean(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'>=3 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'<3 & path2Plot'>=1   )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'<1 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];



% b=bar(7,nanmeadian(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==1 )))  ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.5;
% b.CData(1,:) = [0.1 0.5 .8];
% 
% 
% 
% b=bar(9,nanmeadian(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==0 )))  ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.5;
% b.CData(1,:) = [0.9 0.4 .2];

hold on
scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'>3 ))+.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'>3), 'Marker', 'o', 'MarkerEdgeColor', [.8 .2 .5], 'MarkerFaceColor', [.8 .2 .5] );


scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'<=3 & path2Plot'>=1))+2.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'<=3 & path2Plot'>=1 ), 'Marker', 'o', 'MarkerEdgeColor', [0 0.7 .25], 'MarkerFaceColor',[0 0.7 .25] );


scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'<1))+4.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'<1), 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.1 .6], 'MarkerFaceColor',[0.3 0.1 .6]);





% scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==1))+6.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==1), 'Marker', 'o', 'MarkerEdgeColor', [0.1 0.5 .8]);
% 
% 
% scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==0))+8.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==0), 'Marker', 'o', 'MarkerEdgeColor', [0.9 0.4 .2]);

elseif secondPath


b=bar(1,nanmedian(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==4 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmedian(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==3 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmedian(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==2 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];



b=bar(7,nanmedian(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==1 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.1 0.5 .8];



b=bar(9,nanmedian(val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==0 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.9 0.4 .2];

hold on
scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==4 & secondPath2Plot'<secondPathCut ))  +.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==4 & secondPath2Plot'<secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [.8 .2 .5],'MarkerFaceColor', [.8 .2 .5] );
scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==4 & secondPath2Plot'>= secondPathCut))+.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==4 & secondPath2Plot'>= secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [.8 .2 .5],'MarkerFaceColor', 'k');


scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==3& secondPath2Plot'<secondPathCut))+2.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==3& secondPath2Plot'<secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [0 0.7 .25],'MarkerFaceColor',[0 0.7 .25] );
scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==3 & secondPath2Plot'>= secondPathCut))+2.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==3 & secondPath2Plot'>= secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [0 0.7 .25],'MarkerFaceColor', 'k');


scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==2& secondPath2Plot'<secondPathCut))+4.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==2& secondPath2Plot'<secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.1 .6],'MarkerFaceColor',[0.3 0.1 .6]);
scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==2 & secondPath2Plot'>= secondPathCut))+4.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==2 & secondPath2Plot'>= secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [0.3 0.1 .6],'MarkerFaceColor','k');

scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==1& secondPath2Plot'<secondPathCut))+6.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==1& secondPath2Plot'<secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [0.1 0.5 .8],'MarkerFaceColor',[0.1 0.5 .8]);
scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==1 & secondPath2Plot'>= secondPathCut))+6.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==1 & secondPath2Plot'>= secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [0.1 0.5 .8],'MarkerFaceColor','k');


scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==0& secondPath2Plot'<secondPathCut))+8.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==0& secondPath2Plot'<secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [0.9 0.4 .2],'MarkerFaceColor',[0.9 0.4 .2]);
scatter(rand(1, sum(~remVals & snpMark & pt2Plot & path2Plot'==0 & secondPath2Plot'>= secondPathCut))+8.5, val2Plot(~remVals & snpMark & pt2Plot & path2Plot'==0 & secondPath2Plot'>= secondPathCut), 'Marker', 'o', 'MarkerEdgeColor', [0.9 0.4 .2],'MarkerFaceColor','k');










end