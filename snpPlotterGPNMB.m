function [] = snpPlotterGPNMB(SNP,val2Plot, pt2Plot, ValName, ptName)
%initializing  figure generation for a host of different neurodegenerative
%conditions in the context of GPNMB SNP status

remVals= isnan(val2Plot) | cellfun(@isempty,SNP) | val2Plot<18  ;

overGPNMB= contains(SNP,'TT'); %the major allele
het=contains(SNP, 'CT');
underGPNMB=contains(SNP,'CC'); %the minor allele


b=bar(1,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel(ValName,'FontSize',15)
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Plot(~remVals & het & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Plot(~remVals & underGPNMB & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

hold on
a=scatter(rand(1, sum(~remVals & overGPNMB & pt2Plot))+.5, val2Plot(~remVals & overGPNMB & pt2Plot), 'Marker', 'o' );
a.CData=[.8 .2 .5]; 
b=scatter(rand(1, sum(~remVals & het & pt2Plot))+2.5, val2Plot(~remVals & het & pt2Plot), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, sum(~remVals & underGPNMB & pt2Plot))+4.5, val2Plot(~remVals & underGPNMB & pt2Plot), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

legend({'AA (over Production)', 'GC','GG (under Production)'},'FontSize',15)


[p,~,~]=anovan(val2Plot(~remVals & pt2Plot  ),{SNP(~remVals & pt2Plot  )}, 'display', 'off');


title( [ptName, 'Patients ', ValName, ' as a Function of  rs199347 Status p=', num2str(p) ])
