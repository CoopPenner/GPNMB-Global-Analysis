function [] = BraakEval(pt2use, pathTable, disName, plotType   )
%Our first basic analysis relating snp Status to levels of path burden,
%currently separated based on asyn pure or asyn copath

%initializing all variables


snpStat=pathTable.rs199347;

Braak03=pathTable.Braak03; Braak06=pathTable.Braak06;
overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele



%converting braak 6 to brakk 3

for dd=1:length(Braak03)
    if isempty(Braak03{dd})
        %in this context 1/2=1 3/4=2 5/6=3
braak6Val=Braak06{dd}; braak6Convert=ceil(braak6Val/2); 
Braak03{dd}=braak6Convert;
    end
end



remVals= isnan(val2Test) | cellfun(@isempty,snpStat)  ;


%% plotting figures
% 
% scatter bar option


if strcmp(plotType,'scatBar')

figure
b=bar(1,nanmean(val2Test(~remVals & overGPNMB & pt2use )   ))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel([brainAreaAtPlay, ' ',pathName])
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Test(~remVals & het & pt2use )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Test(~remVals & underGPNMB & pt2use )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

hold on
a=scatter(rand(1, sum(~remVals & overGPNMB & pt2use))+.5, val2Test(~remVals & overGPNMB & pt2use), 'Marker', 'o' );
a.CData=[.8 .2 .5]; 
b=scatter(rand(1, sum(~remVals & het & pt2use))+2.5, val2Test(~remVals & het & pt2use), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, sum(~remVals & underGPNMB & pt2use))+4.5, val2Test(~remVals & underGPNMB & pt2use), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

legend({'AA (over Production)', 'GC','GG (under Production)'})
title([brainAreaAtPlay, pathName, ' Burden in ', disName, ' Patients' ])


% plotting stack bar option

elseif strcmp(plotType,'stackBar')

figure


   overGPNMBRat=[sum(val2Test(~remVals & overGPNMB & pt2use )==9  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
                        sum(val2Test(~remVals & overGPNMB & pt2use )==8  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
                 sum(val2Test(~remVals & overGPNMB & pt2use )==7  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
                 sum(val2Test(~remVals & overGPNMB & pt2use )==6  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
          sum(val2Test(~remVals & overGPNMB & pt2use )==5  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==4  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==3  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==2  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==1  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==0  )/ sum(((~remVals & overGPNMB & pt2use )))*100];
    

hetGPNMBRat=[sum(val2Test(~remVals & het & pt2use )==9  )/ sum(((~remVals & het & pt2use   )))*100,...
             sum(val2Test(~remVals & het & pt2use )==8  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==7  )/ sum(((~remVals & het & pt2use )))*100,...                  
             sum(val2Test(~remVals & het & pt2use )==6  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==5  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==4  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==3  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==2  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==1  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==0  )/ sum(((~remVals & het & pt2use )))*100]; 


underGPNMBRat=[sum(val2Test(~remVals & underGPNMB & pt2use )==9  )/ sum(((~remVals & underGPNMB & pt2use   )))*100 ,...
                 sum(val2Test(~remVals & underGPNMB & pt2use )==8  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
                 sum(val2Test(~remVals & underGPNMB & pt2use )==7  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==6  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
    sum(val2Test(~remVals & underGPNMB & pt2use )==5  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==4  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==3  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==2  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==1  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==0  )/ sum(((~remVals & underGPNMB & pt2use )))*100];
b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;



legend('9','8','7','6','5','4','3','2','1', 'FontSize', 14)

ylabel('percent of total cases', 'FontSize', 15)


title([ pathName, ' Score in ', disName, ' Patients' ], 'FontSize', 20)

ylim([0,120])

a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};

a.YTick(a.YTick>100)=[];


elseif strcmp(plotType,'both')
figure

subplot(1,2,1)


   overGPNMBRat=[sum(val2Test(~remVals & overGPNMB & pt2use )==9  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
                        sum(val2Test(~remVals & overGPNMB & pt2use )==8  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
                 sum(val2Test(~remVals & overGPNMB & pt2use )==7  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
                 sum(val2Test(~remVals & overGPNMB & pt2use )==6  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
          sum(val2Test(~remVals & overGPNMB & pt2use )==5  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==4  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==3  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==2  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==1  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
        sum(val2Test(~remVals & overGPNMB & pt2use )==0  )/ sum(((~remVals & overGPNMB & pt2use )))*100];
    

hetGPNMBRat=[sum(val2Test(~remVals & het & pt2use )==9  )/ sum(((~remVals & het & pt2use   )))*100,...
             sum(val2Test(~remVals & het & pt2use )==8  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==7  )/ sum(((~remVals & het & pt2use )))*100,...                  
             sum(val2Test(~remVals & het & pt2use )==6  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==5  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==4  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==3  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==2  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==1  )/ sum(((~remVals & het & pt2use )))*100,...
             sum(val2Test(~remVals & het & pt2use )==0  )/ sum(((~remVals & het & pt2use )))*100]; 


underGPNMBRat=[sum(val2Test(~remVals & underGPNMB & pt2use )==9  )/ sum(((~remVals & underGPNMB & pt2use   )))*100 ,...
                 sum(val2Test(~remVals & underGPNMB & pt2use )==8  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
                 sum(val2Test(~remVals & underGPNMB & pt2use )==7  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==6  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
    sum(val2Test(~remVals & underGPNMB & pt2use )==5  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==4  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==3  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==2  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==1  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==0  )/ sum(((~remVals & underGPNMB & pt2use )))*100];
b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;





legend('9','8','7','6','5','4','3','2','1', 'FontSize', 14)

ylabel('percent of total cases', 'FontSize', 15)


title([ pathName, ' Score in ', disName, ' Patients' ], 'FontSize', 20)

ylim([0,120])

a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};

a.YTick(a.YTick>100)=[];



subplot(1,2,2)




b=bar(1,nanmean(val2Test(~remVals & overGPNMB & pt2use )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel([pathName], 'FontSize', 15)
a=gca; a.XTickLabel=[];
hold on

b=bar(3,nanmean(val2Test(~remVals & het & pt2use )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Test(~remVals & underGPNMB & pt2use )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

hold on
a=scatter(rand(1, sum(~remVals & overGPNMB & pt2use))+.5, val2Test(~remVals & overGPNMB & pt2use), 'Marker', 'o' );
a.CData=[.8 .2 .5]; 
b=scatter(rand(1, sum(~remVals & het & pt2use))+2.5, val2Test(~remVals & het & pt2use), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, sum(~remVals & underGPNMB & pt2use) )+4.5, val2Test(~remVals & underGPNMB & pt2use), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

legend({'AA (over Production)', 'GC','GG (under Production)'}, 'FontSize', 15)
title([pathName, ' Score in ', disName, ' Patients' ], 'FontSize', 20)
ylim([0,9])


end

