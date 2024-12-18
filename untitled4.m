function [outputArg1,outputArg2] = snpStatPathBurdenTestr(pt2use,brainAreAtPlay, pathTable,   )
%Our first basic analysis relating snp Status to levels of path burden,
%currently separated based on asyn pure or asyn copath






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

% now that we have our path scores as continuous variables, let's output
% categorical for aSyn pure and aSyn copath


asynPure=asynCont > 1 & (tauCont <=1 & aBetaCont <=1 & TDPCont <=1); % I consider a positive path score to be > 1 ('rare')
asynCoPath= asynCont >1 & (tauCont > 1 | aBetaCont >1 | TDPCont >1); % currently not differentiating based on copath subtype, tho this will change depending on power








end