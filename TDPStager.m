function [outputArg1,outputArg2] = TDPStager(brainAreasAtPlay,pathTable,pathCut,pt2use)

%simple for loop that will be repeated many times


tdpCollect=[];

    for dd=1:length(brainAreasAtPlay)
        brainAreaAtPlay=brainAreas{dd};
    [asynCont,tauCont,aBetaCont,TDPCont,neuronDropCont,gliosisCont] = pathScoreGenerate(brainAreaAtPlay,pathTable);
    
    
        
remVals= isnan(TDPCont)' | ~pt2use  ;

tdpCollect=
o
    
    
    
endc1[\]\