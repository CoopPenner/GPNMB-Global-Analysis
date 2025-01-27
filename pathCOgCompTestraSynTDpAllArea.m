function [] = pathCOgCompTestraSynTDpAllArea(pt2use,allAreas, pathTable, cogTable ,pathCut  )



pValz=nan( length(allAreas),1);
rValz=nan( length(allAreas),1);


  for areaCount=1:length(allAreas)
      brainAreaAtPlay=allAreas{areaCount};
   %   corrValues{areaCount,1}= brainAreaAtPlay;

      [asynCont,tauCont,aBetaCont,TDPCont,neuronLossCont,~] = pathScoreGenerate(brainAreaAtPlay,pathTable);

      remVals=   ~pt2use' | (TDPCont+ asynCont==0)   ;
try
      [r,p]=corrcoef(TDPCont(~remVals  ) ,neuronLossCont(~remVals  ),'rows','complete');

pValz(areaCount)=p(2);
rValz(areaCount)=r(2);



  end






significantIdx = pValz < (.05/length(pValz));

% Extract significant regions and their corresponding values
significantAreas = allAreas(significantIdx);
significantRValues = rValz(significantIdx);
significantPValues = pValz(significantIdx);

% Normalize p-values for color mapping (more significant -> darker color)
normPValues = log10(significantPValues); 
normPValues = abs(normPValues); % Make values positive for color scaling
normPValues = (normPValues - min(normPValues)) / (max(normPValues) - min(normPValues));

% Define a colormap 
numBars = length(significantRValues);
colors = hot(numBars); 
colors = colors(round(normPValues * (size(colors, 1) - 1)) + 1, :);

% Create the bar graph
figure;
barGraph = bar(significantRValues, 'FaceColor', 'flat');

% Apply colors based on p-value significance
for i = 1:numBars
    barGraph.CData(i, :) = colors(i, :);
end

% Customize axes and labels
a = gca;
a.XTick = 1:numBars;  % Ensure there are xticks for all bars
a.XTickLabel = significantAreas; 
xtickangle(45);

xlabel('Brain Region');
ylabel('Correlation (r)');
title('Significant Brain Region Correlations (Bonferroni Corrected)');

% Add colorbar with appropriate label
c = colorbar;
colormap(hot);
caxis([min(normPValues) max(normPValues)]);
c.Ticks = linspace(min(normPValues), max(normPValues), 5);
c.TickLabels = (logspace(log10(min(significantPValues)), log10(max(significantPValues)), 5));
c.Label.String = 'log_{10}(p-values)';

% Adjust figure for readability
set(a, 'FontSize', 12);
grid on;







%% test zone


% pathID=pathTable.INDDID;








% cogID=cogTable.ID;
% cogSlope=cogTable.cogSlope;
% endScorez=cogTable.endScore;
% 
% cogScoreAtPlay=nan(1,length(pathID));
% endScoreAtPlay=nan(1,length(pathID));
% 
% %moving cognitive slope into the same format as path data
%   for tt=1:length(pathID)
% cogHold=cogSlope( cogID ==pathID(tt)   );
% endScoreHold=endScorez(cogID==pathID(tt));
%     if ~isempty(endScoreHold)
%     cogScoreAtPlay(tt)= cogHold;
%     endScoreAtPlay(tt)=endScoreHold;
%     end
%   end







% 
% figure
% 
% remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)   ;
% 
% histogram(asynCont(~remVals)- TDPCont(~remVals) )
% 
% 
% figure
% 
% cutLevel=2;
% 
% 
% remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)|TDPCont<cutLevel     ;
% 
% 
% 
%  overGPNMBRatLowTDP=[sum(asynCont(~remVals & overGPNMB  )==4  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==3  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==2  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==1  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==0  )/ sum(((~remVals & overGPNMB )))*100];
% 
% 
% hetGPNMBRatLowTDP=[sum(asynCont(~remVals & het  )==4  )/ sum(((~remVals & het  )))*100,...
% sum(asynCont(~remVals & het)==3  )/ sum(((~remVals & het  )))*100,...
% sum(asynCont(~remVals & het )==2  )/ sum(((~remVals & het )))*100,...
% sum(asynCont(~remVals & het )==1  )/ sum(((~remVals & het )))*100,...
% sum(asynCont(~remVals & het )==0  )/ sum(((~remVals & het  )))*100]; 
% 
% 
% underGPNMBRatLowTDP=[sum(asynCont(~remVals & underGPNMB  )==4  )/ sum(((~remVals & underGPNMB   )))*100 ,...
% sum(asynCont(~remVals & underGPNMB  )==3  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB  )==2  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB )==1  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB )==0  )/ sum(((~remVals & underGPNMB )))*100];
% 
% 
% 
% 
% remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)|TDPCont>cutLevel     ;
% 
% 
%     overGPNMBRatHighTDP=[sum(asynCont(~remVals & overGPNMB  )==4  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==3  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==2  )/ sum(((~remVals & overGPNMB )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==1  )/ sum(((~remVals & overGPNMB  )))*100,...
% sum(asynCont(~remVals & overGPNMB  )==0  )/ sum(((~remVals & overGPNMB  )))*100];
% 
% 
% hetGPNMBRatHighTDP=[sum(asynCont(~remVals & (het)   )==4  )/ sum(((~remVals & (het)   )))*100,...
% sum(asynCont(~remVals & (het)  )==3  )/ sum(((~remVals & (het)  )))*100,...
% sum(asynCont(~remVals & (het) )==2  )/ sum(((~remVals & (het) )))*100,...
% sum(asynCont(~remVals & (het)  )==1  )/ sum(((~remVals & (het)  )))*100,...
% sum(asynCont(~remVals & (het)  )==0  )/ sum(((~remVals & (het)  )))*100]; 
% 
% 
% underGPNMBRatHighTDP=[sum(asynCont(~remVals & underGPNMB  )==4  )/ sum(((~remVals & underGPNMB    )))*100 ,...
% sum(asynCont(~remVals & underGPNMB  )==3  )/ sum(((~remVals & underGPNMB )))*100,...
% sum(asynCont(~remVals & underGPNMB  )==2  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB  )==1  )/ sum(((~remVals & underGPNMB  )))*100,...
% sum(asynCont(~remVals & underGPNMB  )==0  )/ sum(((~remVals & underGPNMB  )))*100];
% 
% 
% 
% figure
%  subplot(1,3,1)
%  b=bar([1,3], [overGPNMBRatLowTDP;overGPNMBRatHighTDP], 'stacked') ;
% subplot(1,3,2)
% b=bar([1,3], [hetGPNMBRatLowTDP;hetGPNMBRatHighTDP], 'stacked') ;
% subplot(1,3,3)
% b=bar([1,3], [underGPNMBRatLowTDP;underGPNMBRatHighTDP], 'stacked') ;
% 
% 
% 
% 
% 
% remVals= isnan(val2Plot)' | cellfun(@isempty,SNP)| ~pt2use | (TDPCont+ asynCont'==0)     ;
% 
% 
% 
% 
% 
% 
% 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (TDPCont+ asynCont'~=0) ) ,TDPCont(~remVals & pt2use & (TDPCont+ asynCont'~=0) ),'rows','complete')
% % 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (het |overGPNMB)  ) ,TDPCont(~remVals & pt2use & (het |overGPNMB) ),'rows','complete')
% % 
% % 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (overGPNMB)  ) ,TDPCont(~remVals & pt2use & (overGPNMB) ),'rows','complete')
% % 
% %
% 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (het)  ) ,TDPCont(~remVals & pt2use & (het) ),'rows','complete')
% 
% 
%  [r,p]=corrcoef(asynCont(~remVals & pt2use & (underGPNMB)  ) ,TDPCont(~remVals & pt2use & (underGPNMB) ),'rows','complete')
% 
% 
% % 
% % 
% % 
% % [r,p]=corrcoef(asynCont(~remVals & pt2use &(het |underGPNMB) & TDPCont'~=0 ) ,TDPCont(~remVals & pt2use &(het |underGPNMB)& TDPCont'~=0 ),'rows','complete')
% % 
% % 
% % [r,p]=corrcoef(asynCont(~remVals & pt2use &underGPNMB& TDPCont'~=0 ) ,TDPCont(~remVals & pt2use &underGPNMB& TDPCont'~=0 ),'rows','complete')
% % 
% % 
% 


end