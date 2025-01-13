function [] = brainPathBurdenOutputr(pt2use, pathTable, pathType,  brainAreas,  diseaseName   )



%this will produce a simple overview of the path ratings from each brain
%area for each patient type


%initializing all variables
globalDx=pathTable.GlobalDx;
pathID=pathTable.INDDID;

numSamp=nan(1,length(brainAreas));
avgBurd=nan(1,length(brainAreas));

for dd=1:length(brainAreas)
    brainAreaAtPlay=brainAreas{dd};
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


keepVals= ~isnan(val2Test) & pt2use';

numSamp(dd)= sum( keepVals); % outputting the number of path samples in a given brain area

avgBurd(dd)=mean(val2Test(keepVals)); %outputting mean path burden

end


%now ranking path burden in patients with sufficent samples and outputting
%a simple bar graph 


avgBurd(numSamp< 50)=nan; %just nanning out brain areas with fewer than 20 samples


[highBurd,i]= maxk(avgBurd, 5);




c=bar(highBurd');

a=gca; a.XTickLabel=brainAreas(i);

c.FaceColor='Flat';
c.CData(1,:)= [.8,0,0] ;
c.CData(2,:)= [.6,.2,0];
c.CData(3,:)=[.4,.4,.1];
c.CData(4,:)=[.2,.3,.2];
c.CData(5,:)=[.1,.1,.4];


title([diseaseName, ' Patients ', pathName, ' Burden', ' Across Assesed brain regions' ], 'FontSize',15)

ylabel([pathName, ' pathological ratings'] ,'FontSize', 15)


end