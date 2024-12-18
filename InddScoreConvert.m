function [outputVar] = InddScoreConvert(outputVar, feature2Test)
%simple conversion of categorical to numerical values of INDD neuropath
%scores

% here I  consider 0 as 0 Rare as 1 1+ as 2, 2+ as 3 3+ as 4



        for dd=1: length(feature2Test) % simply convert the categorical features to numbers
            if ~isempty(feature2Test{dd}) & strcmp(feature2Test{dd}, '0')
                outputVar(dd)=0;
            elseif ~isempty(feature2Test{dd}) & strcmp(feature2Test{dd}, 'Rare')
                outputVar(dd)=1;
            elseif ~isempty(feature2Test{dd}) & strcmp(feature2Test{dd}, '1+')
                outputVar(dd)=2;
            elseif ~isempty(feature2Test{dd}) & strcmp(feature2Test{dd}, '2+')
                outputVar(dd)=3;
            elseif ~isempty(feature2Test{dd}) & strcmp(feature2Test{dd}, '3+')
                outputVar(dd)=4;
            else
                outputVar(dd)=nan; % just for readability 
            end
        end




end