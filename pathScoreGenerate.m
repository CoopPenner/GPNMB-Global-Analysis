function [asynCont,tauCont,aBetaCont,TDPCont] = pathScoreGenerate(brainAreaAtPlay,pathTable)
%Simple function that converts the categorical path data to continuous data

pathID=pathTable.INDDID;
tauName=[brainAreaAtPlay,'Tau'];
aBetaName=[brainAreaAtPlay,'AntibodyPlaques'];
TDPName=[brainAreaAtPlay,'TDP43'];
aSynName=[brainAreaAtPlay,'aSyn'];




for pathAtPlay=1:4 % I iterate through each path type of interest
holderVar=nan(1,length(pathID));
    if pathAtPlay==1 %first Asyn
        feature2Test=pathTable.(aSynName);
    elseif pathAtPlay==2 %tau
        feature2Test=pathTable.(tauName);
    elseif pathAtPlay==3 %aBeta (antibody in INDD)
        feature2Test=pathTable.(aBetaName);
    elseif pathAtPlay==4 %TDP43
        feature2Test=pathTable.(TDPName);
    end

    holderVar=nan(1,length(pathID));
    outputVar=InddScoreConvert(holderVar, feature2Test); %this function converts scores from categorical to continuous

   if pathAtPlay==1
        asynCont=outputVar;
    elseif pathAtPlay==2
        tauCont=outputVar;
    elseif pathAtPlay==3
        aBetaCont=outputVar;
    elseif pathAtPlay==4 
        TDPCont=outputVar;
   end
end
