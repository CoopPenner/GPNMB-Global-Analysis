function [] = cogPathSplitScatBarPlotr(val2Plot, remVals, snpStat,pt2Plot,path2Plot1, path2Plot2, pathCut)




overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele




b=bar(1,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.15;
b.CData(1,:) = [.8 .2 .5]; 
a=gca; a.XTickLabel=[];

hold on


b=bar(3.2,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut ))) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.15;
b.CData(1,:) = [.8 .4 .7]; 
a=gca; a.XTickLabel=[];



% b=bar(5.4,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot1'<pathCut & path2Plot2'< pathCut ))) ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.15;
% b.CData(1,:) = [.8 .6 .9]; 
% a=gca; a.XTickLabel=[];




scatter(   rand(1, sum(~remVals & overGPNMB & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut)  )+.5,...
    val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut) , 'Marker', 'o', 'MarkerEdgeColor', [.8 .2 .5], 'MarkerFaceColor', [.8 .2 .5] );


scatter(   rand(1, sum(~remVals & overGPNMB & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut)  )+2.65,...
    val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut) , 'Marker', 'o', 'MarkerEdgeColor', [.8 .4 .7],  'MarkerFaceColor', [.8 .4 .7]  );

% scatter(   rand(1, sum(~remVals & overGPNMB & pt2Plot & path2Plot1'==0 & path2Plot2'==0)  )+4.8,...
%     val2Plot(~remVals & overGPNMB & pt2Plot & path2Plot1'==0 & path2Plot2'==0) , 'Marker', 'o', 'MarkerEdgeColor', [.8 .6 .9],  'MarkerFaceColor', [.8 .6 .9]  );




%het

b=bar(7.6,nanmean(val2Plot(~remVals & het & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.15;
b.CData(1,:) = [0 0.7 .25]; 
a=gca; a.XTickLabel=[];



b=bar(9.8,nanmean(val2Plot(~remVals & het & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut ))) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.15;
b.CData(1,:) = [0 0.9 .25]; 
a=gca; a.XTickLabel=[];


% b=bar(12,nanmean(val2Plot(~remVals & het & pt2Plot & path2Plot1'==0 & path2Plot2'==0 ))) ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1.15;
% b.CData(1,:) = [0.2 0.9 .5]; 
% a=gca; a.XTickLabel=[];


scatter(   rand(1, sum(~remVals & het & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut)  )+6.95,...
    val2Plot(~remVals & het & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut) , 'Marker', 'o', 'MarkerEdgeColor', [0 0.7 .25], 'MarkerFaceColor', [0 0.7 .25]);


scatter(   rand(1, sum(~remVals & het & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut)  )+9.3,...
    val2Plot(~remVals & het & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut) , 'Marker', 'o', 'MarkerEdgeColor', [0 0.9 .5],  'MarkerFaceColor', [0 0.9 .5] );

% scatter(   rand(1, sum(~remVals & het & pt2Plot & path2Plot1'==0 & path2Plot2'==0)  )+11.45,...
%     val2Plot(~remVals & het & pt2Plot & path2Plot1'==0 & path2Plot2'==0) , 'Marker', 'o', 'MarkerEdgeColor', [0.2 0.9 .5],  'MarkerFaceColor', [0.2 0.9 .5] );
% 

%protective


b=bar(14.2,nanmean(val2Plot(~remVals & underGPNMB & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1;
b.CData(1,:) = [0.2 0.4 .9];
a=gca; a.XTickLabel=[];



b=bar(16.4,nanmean(val2Plot(~remVals & underGPNMB & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut ))) ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1;
b.CData(1,:) = [0.1 0.6 .9];
a=gca; a.XTickLabel=[];



% b=bar(18.6,nanmean(val2Plot(~remVals & underGPNMB & pt2Plot & path2Plot1'==0 & path2Plot2'==0 ))) ;
% b.FaceColor = 'flat';
% b.FaceAlpha=.3;
% b.BarWidth=1;
% b.CData(1,:) = [0 0.8 .9];
% a=gca; a.XTickLabel=[];


scatter(   rand(1, sum(~remVals & underGPNMB & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut)  )+13.65,...
    val2Plot(~remVals & underGPNMB & pt2Plot & path2Plot1'>=pathCut & path2Plot2'<pathCut) , 'Marker', 'o', 'MarkerEdgeColor', [0.2 0.4 .9], 'MarkerFaceColor', [0.2 0.4 .9]);


scatter(   rand(1, sum(~remVals & underGPNMB & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut)  )+15.95,...
    val2Plot(~remVals & underGPNMB & pt2Plot & path2Plot1'<pathCut & path2Plot2'>=pathCut) , 'Marker', 'o', 'MarkerEdgeColor', [0.1 0.6 .9],  'MarkerFaceColor', [0.1 0.6 .9] );

% scatter(   rand(1, sum(~remVals & underGPNMB & pt2Plot & path2Plot1'==0 & path2Plot2'==0)  )+ 18.05,...
%     val2Plot(~remVals & underGPNMB & pt2Plot & path2Plot1'==0 & path2Plot2'==0) , 'Marker', 'o', 'MarkerEdgeColor', [0 0.8 .9],  'MarkerFaceColor', [0 0.8 .9] );
% 
