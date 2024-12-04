%basic check of rs199347 frequency and cognitive decline

addpath('/Users/pennerc/Documents/MATLAB/GPNMB_Global_Analysis')


%% reading in data
dataRead=true;
%read in metaData


if dataRead

%reading in all of our data

dataTable=readtable('/Volumes/PC60/InqueryDatasets/AllPatients_WithMoca_11_24.xlsx');
end

cogScore=dataTable.MoCATotal;
testDate=dataTable.TestDate;
dx=dataTable.GlobalDx;
ageATTest=dataTable.AgeatTest;
snpAtPlay=dataTable.rs199347;
globalID=dataTable.INDDID;

singleID=unique(globalID);

%initializing collection vectors for real this time
slopeCollect=nan(1, length(singleID));
startScore=nan(1, length(singleID));
year2DecCollect=nan(1, length(singleID));
diagCollect=cell(1, length(singleID));
snpCollect=cell(1, length(singleID));
ageAtTest=nan(1, length(singleID));

for tt= 1:length(singleID) %iterating through each patient in the dataset

relDates=testDate(   globalID==singleID(tt)    ); %all the test dates
relScores=cogScore(globalID==singleID(tt) ); %all the scores
%sigh these are sometimes out of order  
[testYear,testMonth]=ymd(relDates);        
pureTime=datenum(relDates); %convert it to pure time to sort
[sorted,n]=sort(pureTime); %sorting on absolute value of time
testYear=testYear(n); testMonth=testMonth(n); relScores=relScores(n); %then just outputting the reorder

testYear=testYear-testYear(1); %everything relative to start date;
testYear=testYear*12;

testMonth=testMonth-testMonth(1);

allTime=testYear+testMonth;

        allTime((isnan(relScores)))=[];
        relScores(isnan(relScores))=[]; 
        [uniqueA, i,j] = unique(allTime,'first');
        idxToRemove = find(not(ismember(1:numel(allTime),i)))   ;
        allTime(idxToRemove)=[]; relScores(idxToRemove)=[]; %sometimes for some reason there are duplicate times
        allTime=allTime/12; %everything contextualized in years, since that is the standard by which slopes are presented

        if ~isempty(relScores) & length(relScores)~=1  
        scoreVec=diff(relScores);  sigChange=find( (relScores-relScores(1) ) <=-2 ,1); %here we identify when scores sink two points below first measure.      
        P=polyfit(allTime,relScores,1);
        slopeCollect(tt)=P(1);
if P(1)<-30
    b=1;
end


        diagCollect(tt)=unique(dx(globalID==singleID(tt)));
        snpCollect(tt)=unique(snpAtPlay(globalID==singleID(tt)));
        % ageAtTest(tt)=ageAtTest(globalID==globalID(tt)) ;

            if ~isempty(sigChange)      
        year2DecCollect(tt)= allTime(sigChange)  ;
            else
         year2DecCollect(tt)= nan ;
            end
        startScore(tt)=relScores(1); %collect starting point
           
        else
              slopeCollect(tt)=nan;
                diagCollect{tt}='nan';
                snpCollect{tt}='nan';
        
        end

end



%% Plotting some basic clinical metrics

cellfun(@isempty, diagCollect)

Parkinson= (contains(diagCollect, 'Parkinson'));
Alzheimer= (contains(diagCollect, 'Alzheimer'));
ALS= (contains(diagCollect, 'Amyotrophic'));
DemLewy= (contains(diagCollect, 'Dementia with Lewy Bodies')); 
MCI= (contains(diagCollect, 'Mild cognitive impairment')); 
corticoBasal= (contains(diagCollect, 'Corticobasal syndrome')); 
bvFTD= (contains(diagCollect, 'bvFTD-FTLD')); 
PPA= (contains(diagCollect, 'PPA')); 
supraNuc= (contains(diagCollect, 'Progressive supranuclear palsy')); 
neuroPanel= Alzheimer + ALS + DemLewy  +Parkinson+ corticoBasal +bvFTD + PPA +supraNuc;
ParkinsonianPanel= Parkinson +corticoBasal +supraNuc;
ofInterest=Parkinson+ALS+Alzheimer+bvFTD;

ParkinsonianDem=corticoBasal+DemLewy;
other= ~neuroPanel;
HC=(contains(diagCollect, 'Normal')); 


overGPNMB= contains(snpCollect,'TT'); %the major allele
het=contains(snpCollect, 'CT');
underGPNMB=contains(snpCollect,'CC'); %the minor allele


figure
histogram(slopeCollect((overGPNMB) ),'BinWidth',1, 'FaceColor','b')
hold on
histogram(slopeCollect((underGPNMB)),'BinWidth',1, 'FaceColor','r')
histogram(slopeCollect((het)),'BinWidth',1, 'FaceColor','g')
xlabel('MOCA slope')


pt2Plot=[Parkinson];
val2Plot=slopeCollect;
SNP=snpCollect;
remVals= isnan(val2Plot) | cellfun(@isempty,SNP) | startScore<10 ;

figure
b=bar(1,nanmean(val2Plot(~remVals & overGPNMB & pt2Plot )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel('MOCA Slope')
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

legend({'AA (over Production)', 'GC','GG (under Production)'})


pt2Plot=Parkinson



allPt=[ slopeCollect(~remVals & (overGPNMB|het) & pt2Plot),  slopeCollect(~remVals & underGPNMB & pt2Plot)];

allMark=[ones(1, sum(~remVals & (overGPNMB|het)  & pt2Plot)),  2*ones(1, sum(~remVals & underGPNMB & pt2Plot)) ];

[p,~,~]=anovan(allPt,{allMark});



figure

pt2Plot=neuroPanel


remVals= isnan(val2Plot) | cellfun(@isempty,SNP) | startScore<10 ;

compVal=-3;

bar([sum(slopeCollect(~remVals & (overGPNMB) & pt2Plot)<compVal)/ sum(overGPNMB), sum(slopeCollect(~remVals & (het) & pt2Plot)<compVal)/ sum(het),   sum(slopeCollect(~remVals & (underGPNMB) & pt2Plot)<compVal)/ sum(underGPNMB) ])


figure
bar([sum(slopeCollect(~remVals & (overGPNMB) & pt2Plot)>-0)/ sum(overGPNMB), sum(slopeCollect(~remVals & (het) & pt2Plot)>-0)/ sum(het),   sum(slopeCollect(~remVals & (underGPNMB) & pt2Plot)>-.0)/ sum(underGPNMB) ])



%% reorienting data correctly

origData2use=origData(:,12:end);
origNames=origData2use.Properties.VariableNames;