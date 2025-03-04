function [avgPath] = TDPStager(brainAreasAtPlay,pathTable,pathCut,pt2use)

%simple for loop that will be repeated many times


tdpCollect=nan(length(pt2use), length(brainAreasAtPlay)  );



    for dd=1:length(brainAreasAtPlay)
        brainAreaAtPlay=brainAreasAtPlay{dd};
    [asynCont,tauCont,aBetaCont,TDPCont,neuronDropCont,gliosisCont] = pathScoreGenerate(brainAreaAtPlay,pathTable);
    
    
        
remVals= isnan(TDPCont)' | ~pt2use  ;

TDPCont(remVals)=nan;

tdpCollect(:,dd)=(TDPCont);
    end
    
if size(tdpCollect,2)>1
  avgPath=mean(tdpCollect,2);  
  %pathPos=nan(size(avgPath)); pathPos(avgPath>pathCut)=1; pathPos(avgPath<pathCut)=0;
else
    avgPath=tdpCollect;
 % pathPos=nan(size(tdpCollect)); pathPos(tdpCollect>pathCut)=1; pathPos(tdpCollect<pathCut)=0;
end