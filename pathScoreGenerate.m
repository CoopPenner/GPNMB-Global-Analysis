function [asynCont,tauCont,aBetaCont,TDPCont,neuronDropCont,gliosisCont] = pathScoreGenerate(brainAreaAtPlay,pathTable)
%Simple function that converts the categorical path data to continuous data

pathID=pathTable.INDDID;
tauName=[brainAreaAtPlay,'Tau'];
aBetaName=[brainAreaAtPlay,'AntibodyPlaques'];
TDPName=[brainAreaAtPlay,'TDP43'];
aSynName=[brainAreaAtPlay,'aSyn'];
neuronDropName=[brainAreaAtPlay, 'NeuronLoss'];
gliosisName=[brainAreaAtPlay,'Gliosis'];



for pathAtPlay=1:6 % I iterate through each path type of interest
holderVar=nan(1,length(pathID));
    if pathAtPlay==1 %first Asyn
        feature2Convert=pathTable.(aSynName);
    elseif pathAtPlay==2 %tau
        feature2Convert=pathTable.(tauName);
    elseif pathAtPlay==3 %aBeta (antibody in INDD)
        feature2Convert=pathTable.(aBetaName);
    elseif pathAtPlay==4 %TDP43
        feature2Convert=pathTable.(TDPName);
    elseif pathAtPlay==5
        try
        feature2Convert=pathTable.(neuronDropName);
        catch
            feature2Convert=nan; %in some brain areas like Ponsthere are no neuronal loss or gliosis markers... not sure why
        end
    elseif pathAtPlay==6 
        try
        feature2Convert=pathTable.(gliosisName);
        catch
            feature2Convert=nan;
        end 
    end

    holderVar=nan(1,length(pathID));
    if length(feature2Convert)>1
    outputVar=InddScoreConvert(holderVar, feature2Convert); %this function converts scores from categorical to continuous
    else
        outputVar=nan;
    end

   if pathAtPlay==1
        asynCont=outputVar;
    elseif pathAtPlay==2
        tauCont=outputVar;
    elseif pathAtPlay==3
        aBetaCont=outputVar;
    elseif pathAtPlay==4 
        TDPCont=outputVar;
   elseif pathAtPlay==5
       neuronDropCont=outputVar;
   elseif pathAtPlay==6 
        gliosisCont=outputVar;
   end
end


