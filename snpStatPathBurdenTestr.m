function [outputArg1,outputArg2] = snpStatPathBurdenTestr(pt2use,brainAreaAtPlay, pathTable,pathType, disName, plotType, basicData   )
%Our first basic analysis relating snp Status to levels of path burden,
%currently separated based on asyn pure or asyn copath

%initializing all variables
globalDx=pathTable.GlobalDx;
snpStat=pathTable.rs199347;
pathID=pathTable.INDDID;



overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele

% now that we have our path scores as continuous variables, let's output
% categorical for aSyn pure and aSyn copath

[asynCont,tauCont,aBetaCont,TDPCont,neuronDropCont,gliosisCont] = pathScoreGenerate(brainAreaAtPlay,pathTable);

if pathType==1
    val2Test=asynCont;
    pathName='Asyn';
elseif pathType==2
    val2Test=aBetaCont;
        pathName='aBeta';
elseif pathType==3
    val2Test=tauCont;
    pathName='Tau';
elseif  pathType==4
    val2Test=TDPCont;
    pathName='TDP43';
elseif pathType==5
    val2Test=neuronDropCont;
    pathName='NeuronLoss';
elseif pathType==6
    val2Test=gliosisCont;
    pathName='Gliosis';
end

% asynPure=val2Test > 1 & (tauCont <=1 & aBetaCont <=1 & TDPCont <=1); % I consider a positive path score to be > 1 ('rare')
% asynCoPath= val2Test >1 & (tauCont > 1 | aBetaCont >1 | TDPCont >1); % currently not differentiating based on copath subtype, tho this will change depending on power



remVals= isnan(val2Test)' | cellfun(@isempty,snpStat)  ;


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


overGPNMBRat=[sum(val2Test(~remVals & overGPNMB & pt2use )==4  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & overGPNMB & pt2use )==3  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & overGPNMB & pt2use )==2  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & overGPNMB & pt2use )==1  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & overGPNMB & pt2use )==0  )/ sum(((~remVals & overGPNMB & pt2use )))*100];
    

hetGPNMBRat=[sum(val2Test(~remVals & het & pt2use )==4  )/ sum(((~remVals & het & pt2use   )))*100,...
sum(val2Test(~remVals & het & pt2use )==3  )/ sum(((~remVals & het & pt2use )))*100,...
sum(val2Test(~remVals & het & pt2use )==2  )/ sum(((~remVals & het & pt2use )))*100,...
sum(val2Test(~remVals & het & pt2use )==1  )/ sum(((~remVals & het & pt2use )))*100,...
sum(val2Test(~remVals & het & pt2use )==0  )/ sum(((~remVals & het & pt2use )))*100]; 


underGPNMBRat=[sum(val2Test(~remVals & underGPNMB & pt2use )==4  )/ sum(((~remVals & underGPNMB & pt2use   )))*100 ,...
sum(val2Test(~remVals & underGPNMB & pt2use )==3  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==2  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==1  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==0  )/ sum(((~remVals & underGPNMB & pt2use )))*100];
b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;



legend('3+','2+','1+','rare','none')

ylabel('percent of total cases')

title([brainAreaAtPlay,' ', pathName, ' Burden in ', disName, ' Patients' ])

ylim([0,120])

a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};

a.YTick(a.YTick>100)=[];


elseif strcmp(plotType,'both')
figure

subplot(1,2,1)

    overGPNMBRat=[sum(val2Test(~remVals & overGPNMB & pt2use )==4  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & overGPNMB & pt2use )==3  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & overGPNMB & pt2use )==2  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & overGPNMB & pt2use )==1  )/ sum(((~remVals & overGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & overGPNMB & pt2use )==0  )/ sum(((~remVals & overGPNMB & pt2use )))*100];
    

hetGPNMBRat=[sum(val2Test(~remVals & het & pt2use )==4  )/ sum(((~remVals & het & pt2use   )))*100,...
sum(val2Test(~remVals & het & pt2use )==3  )/ sum(((~remVals & het & pt2use )))*100,...
sum(val2Test(~remVals & het & pt2use )==2  )/ sum(((~remVals & het & pt2use )))*100,...
sum(val2Test(~remVals & het & pt2use )==1  )/ sum(((~remVals & het & pt2use )))*100,...
sum(val2Test(~remVals & het & pt2use )==0  )/ sum(((~remVals & het & pt2use )))*100]; 


underGPNMBRat=[sum(val2Test(~remVals & underGPNMB & pt2use )==4  )/ sum(((~remVals & underGPNMB & pt2use   )))*100 ,...
sum(val2Test(~remVals & underGPNMB & pt2use )==3  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==2  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==1  )/ sum(((~remVals & underGPNMB & pt2use )))*100,...
sum(val2Test(~remVals & underGPNMB & pt2use )==0  )/ sum(((~remVals & underGPNMB & pt2use )))*100];
b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;



legend('3+','2+','1+','rare','none')

ylabel('percent of total cases')

title([brainAreaAtPlay,' ', pathName, ' Burden in ', disName, ' Patients' ])

ylim([0,120])

a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};

a.YTick(a.YTick>100)=[];


subplot(1,2,2)




b=bar(1,nanmean(val2Test(~remVals & overGPNMB & pt2use )))  ;
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

c=scatter(rand(1, sum(~remVals & underGPNMB & pt2use) )+4.5, val2Test(~remVals & underGPNMB & pt2use), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

legend({'AA (over Production)', 'GC','GG (under Production)'})
title([brainAreaAtPlay,' ', pathName, ' Burden in ', disName, ' Patients' ])
ylim([0,5])


end

%% statistical testing


%first, ensure that the two tables line up... this is super irritating,
%basically just because INDD is challenging to navigate so I have to output
%things in separate tables...



% Sex=basicData.Sex;
% ageAtDeath=basicData.AgeatDeath;
% ID=basicData.INDDID;
% 
% 
% 
% sexVar=cell(1,length(pathID));
% ageAtDeathVar=nan(1,length(pathID));
%     for dd=1:length(pathID)
%     basicIDX=find( ID==pathID(dd) );
%     sexVar{dd}=Sex{basicIDX}; ageAtDeathVar(dd)=ageAtDeath(basicIDX);
%     end
% 
% 



% ageVar=ageAtDeathVar(~remVals & pt2use)' ;
% IDVar=(pathID(~remVals & pt2use))   ;
% sexVarTable=categorical(sexVar(~remVals & pt2use))' ;
% snpVar= categorical(snpStat(~remVals & pt2use));
% pathVar= val2Test(~remVals & pt2use)';
% 
% varNames=["age","ID", "sex","snp","pathScore"];
% 
% glmeTable = table(ageVar ,IDVar,sexVarTable,snpVar,pathVar,'Variablenames',varNames);
% 
% 
% 
% 
% glme = fitglme(glmeTable,...
% 		'pathScore ~ 1  + sex + age + snp   ',...
% 		'Distribution','Normal','Link','identity','FitMethod','Laplace',...
% 		'DummyVarCoding','reference')
% 
% 






% 
% overGPNMB= contains(snpStat,'TT'); %the major allele
% het=contains(snpStat, 'CT');
% underGPNMB=contains(snpStat,'CC'); %the minor allele
% 
% snpCond=nan(1, length(overGPNMB)   );
% snpCond(overGPNMB & ~remVals &pt2use )=1; snpCond(het & ~remVals &pt2use)=2; snpCond(underGPNMB & ~remVals &pt2use)=3;
% 
% [n,p,groups]=anovan(val2Test(~remVals & pt2use &  TDPCont'~=0 ),{snpCond(~remVals & pt2use &  TDPCont'~=0)})
% 
% figure
% multcompare(groups)
% 
%     [p,n]=prop_test(  [sum(val2Test'==0 & pt2use & ~remVals &  overGPNMB), sum(val2Test'==0 & pt2use & ~remVals & (underGPNMB)) ],...
%     [sum( pt2use &~remVals &overGPNMB), sum( pt2use & ~remVals &(underGPNMB)) ] ,'false')
% 


   % [p,n]=prop_test(  [sum(val2Test'==0 & pt2use & ~remVals &  overGPNMB), sum(val2Test'==0 & pt2use & ~remVals & (underGPNMB)) ],...
   % [sum( pt2use &~remVals &overGPNMB), sum( pt2use & ~remVals &(underGPNMB)) ] ,'false')


% nanmean(TDPCont(~remVals  & pt2use & asynCont'==4))
% nanmean(TDPCont(~remVals  & pt2use & asynCont'==3))
% nanmean(TDPCont(~remVals  & pt2use & asynCont'==2))
% nanmean(TDPCont(~remVals  & pt2use & asynCont'==1))
% nanmean(TDPCont(~remVals  & pt2use & asynCont'==0))
% 
% 
% 
% 
% nanmean(asynCont(~remVals  & pt2use & TDPCont'==4))
% nanmean(asynCont(~remVals  & pt2use & TDPCont'==3))
% nanmean(asynCont(~remVals  & pt2use & TDPCont'==2))
% nanmean(asynCont(~remVals  & pt2use & TDPCont'==1))
% nanmean(asynCont(~remVals  & pt2use & TDPCont'==0))








nanmean(asynCont(~remVals  & pt2use & TDPCont'==4 & neuronDropCont'>0));
nanmean(asynCont(~remVals  & pt2use & TDPCont'==3& neuronDropCont'>0));
nanmean(asynCont(~remVals  & pt2use & TDPCont'==2& neuronDropCont'>0));
nanmean(asynCont(~remVals  & pt2use & TDPCont'==1& neuronDropCont'>0));
nanmean(asynCont(~remVals  & pt2use & TDPCont'==0& neuronDropCont'>0));
 %#ok<NASGU>