%scratch Track Analaysis 


% addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');
% % first upload relevant video
% 
% homePath= '/Volumes/PC60/iMicroVids/KOScratchAssay_6_6_25_allChannels/KO_aB'; %Link the folder containing all videos;
% 
% 
%         %all videos are saved in the following way Genotype_ImageType_Treatment
% % (a space) (number of video).mp4
% % so for example: KO_Scratch_aB (1).mp4   Note the space between aB and (1)
% genoTypePrefix1='KO';
% genoType1Path='/Volumes/PC60/iMicroVids/KOScratchAssay_6_6_25_allChannels/KO_aB';
% 
% genoTypePrefix2='WT';
% genoType2Path='/Volumes/PC60/iMicroVids/KOScratchAssay_6_6_25_allChannels/WT_aB';
% 
% 
% 
% imageTypePrefix='Scratch';
% treatmentPrefix='aB';
% 
% outputSuffix='6_6_25';
% 
% frameRate=5; %one image every 5 minutes
% 
% 
% 


%%  ------------------------- Setup -------------------------
genoTypes = {'KO', 'WT'};
genoPaths = {
    '/Volumes/PC60/iMicroVids/KOScratchAssay_6_6_25_allChannels/KO_aB/ResultsRepo', ...
    '/Volumes/PC60/iMicroVids/WT_ScratchAssay_6_6_25_allChannels_Scale/WT_aB/ResultsRepo'
};

treatmentPrefix = 'aB';
imageTypePrefix = 'Scratch';

structAtPlay = [];  % will become a table

% ------------------------- Loop over genotypes -------------------------
structAtPlay = struct([]);

for dd = 1:length(genoTypes)
    geno = genoTypes{dd};
    folderPath = genoPaths{dd};

    fileList = dir(fullfile(folderPath, '*.mat'));
    fileList = fileList(~startsWith({fileList.name}, '._')); % remove macOS metadata files
    fileList = fileList(contains({fileList.name}, 'scratchMetrics'));

    for f = 1:length(fileList)
        fileName = fileList(f).name;
        fullPath = fullfile(folderPath, fileName);

        loaded = load(fullPath);
        metricsAtPlay = loaded.scratchMet;

        % Extract video number from filename 
        vidNum = regexp(fileName, '\((\d+)\)', 'tokens');
        vidNum = str2double(vidNum{1});

        % Append to structure array
        structAtPlay(end+1).Genotype = geno;
        structAtPlay(end).Treatment = treatmentPrefix;
        structAtPlay(end).VideoNum = vidNum;

        structAtPlay(end).InScratchCount = metricsAtPlay.inScratchCount;
        structAtPlay(end).NearScratchCount = metricsAtPlay.nearScratchCount;
        structAtPlay(end).TotObjects = metricsAtPlay.totObjects;
        structAtPlay(end).AvgDistToScratch = metricsAtPlay.avgDistTot;
        structAtPlay(end).StdDistToScratch = metricsAtPlay.stdDistTot;
        structAtPlay(end).ScratchArea = metricsAtPlay.scratchSize;
    end
end



%% setting up imaging basics


    % Genotype colors
    colorMap = struct('KO', [0.5 0 0.5], ...          % purple
                      'WT', [34,139,34]/255);         % forest green

    figure; hold on;
    set(gcf, 'Color', 'w');
    genotypes = {'KO', 'WT'};
    alphaVal = 0.25;  % transparency for individual lines
    maxLength = 168; %14 hours

    %Pixel to uM transition previously calculated as  % ~1.13 µm/px based
    %on a field of view of 1300 micrometers



%%  ok first let's do average distance to scratch. We will first plot the average distance. We will then plot the cumulative sum of the velocity 
%statistical testing will be done with a glme, output is the interaction
%term between genotype and time passage with average distance as the
%predicted variable. 


figure;

%  containers
traceCatchTot = [];
genoCatchTot = [];
timeCatchTot = [];
vidNumCatchTot = {};
countCatchTot = [];
startValCatchTot=[];
scratchAreaCatchTot = [];
frameRate=5;

maxLength = 14*12; % 14 hours tho some sessions are longer
timeVec =   [0:5:5* (maxLength-1)] ;% Time in minutes

for gg = 1:numel(genotypes)
    geno = genotypes{gg};
    color = colorMap.(geno);

    % Get all entries for this genotype
    genoEntries = structAtPlay(strcmp({structAtPlay.Genotype}, geno));
    
    % Extract vectors
    traces = {genoEntries.AvgDistToScratch};
    totCounts = {genoEntries.TotObjects};
    vidNums = [genoEntries.VideoNum];
    scratchAreas = [genoEntries.ScratchArea];

    subplot(1, 2, gg);
    hold on;

    for jj = 1:numel(traces)

        
        trace = traces{jj}(1:maxLength);
        countz=totCounts{jj}(1:maxLength);

        % Plot individual trace
        plot(timeVec, trace, 'Color', [color alphaVal], 'LineWidth', 1);

        % Store values for GLME
        traceCatchTot = [traceCatchTot, trace'];
        genoCatchTot = [genoCatchTot, repmat({geno}, 1, maxLength)];
        timeCatchTot = [timeCatchTot, timeVec ];
        vidNumCatchTot = [vidNumCatchTot, repmat({sprintf('%s_%d', geno, vidNums(jj))}, 1, maxLength)];
        countCatchTot = [countCatchTot,countz' ];
        scratchAreaCatchTot = [scratchAreaCatchTot, (scratchAreas(jj)* ones(1, maxLength))];
                startValCatchTot=[startValCatchTot, (trace(1)* ones(1, (maxLength)))];

    end

    % Plot mean trace
    allTraces = cellfun(@(x) x(1:maxLength), traces, 'UniformOutput', false);
    traceMat = cell2mat(cellfun(@(x) padarray(x(:)', [0, maxLength - length(x)], NaN, 'post'), allTraces, 'UniformOutput', false)');
    meanTrace = nanmean(traceMat, 1);
    plot(timeVec, meanTrace, '-', 'Color', color, 'LineWidth', 2.5);

    % Format subplot
    xlabel('Time (Minutes)','FontSize',20);
    ylabel('Average Distance to Scratch','FontSize', 20);
    title(sprintf('%s Microglia Treated with %s', geno, treatmentPrefix),'FontSize',30);
    box on;
    xlim([0, maxLength*frameRate]);

 %        subplot(2,2,3:4)
 %        hold on
 %    plot(timeVec, meanTrace, '-', 'Color', color, 'LineWidth', 2.5);
 % xlabel('Time (Minutes)');
 %    ylabel('Average Distance to Scratch');
 %    title(sprintf('%s Microglia Treated with %s', geno, treatmentPrefix));
 %    box on;
 %    xlim([0, maxLength*frameRate]);

end


%% alright time for statistical testing running a glme...


% make model  table
modelData = table( ...
    traceCatchTot(:), ...
    startValCatchTot(:),...
    vidNumCatchTot(:), ...
    countCatchTot(:), ...
    scratchAreaCatchTot(:), ...
    genoCatchTot(:), ...
    timeCatchTot(:), ...
    'VariableNames', {'Dist', 'startVal', 'vidNum', 'TotCells', 'scratchArea', 'genotype', 'time'} ...
);


modelData.vidNum = cell2mat(modelData.vidNum);
modelData.genotype = cell2mat(modelData.genotype);

% Fit GLME  and alt model
glme = fitglme(modelData, ...
    'Dist ~ 1 + TotCells + startVal+ scratchArea + genotype * time + (1|vidNum)', ...
    'Distribution','normal','Link','identity');


pVal = glme.Coefficients.pValue(7);
sgtitle(['Comparative Global Distance From Scratch p= ' num2str(pVal) ],'FontSize',45) 











%% here it is for the ratio of cells in scratch awful this is separate for
%%now but will make a caller function soon 




figure;

%  containers
traceCatchTot = [];
genoCatchTot = [];
timeCatchTot = [];
vidNumCatchTot = {};
countCatchTot = [];
scratchAreaCatchTot = [];
startValCatchTot=[];
frameRate=5;

maxLength = 14*12; % 14 hours tho some sessions are longer
timeVec =   [0:5:5* (maxLength-2)] ;% Time in minutes

for gg = 1:numel(genotypes)
    geno = genotypes{gg};
    color = colorMap.(geno);

    % Get all entries for this genotype
    genoEntries = structAtPlay(strcmp({structAtPlay.Genotype}, geno));
    
    % Extract vectors
    totCounts = {genoEntries.TotObjects};
    traces = {genoEntries.AvgDistToScratch};
    trace2use={genoEntries.NearScratchCount};

    vidNums = [genoEntries.VideoNum];
    scratchAreas = [genoEntries.ScratchArea];
        traceMat=[];
    subplot(1, 2, gg);
    hold on;

    for jj = 1:numel(traces)

        
        trace = trace2use{jj}(2:maxLength);
        countz=totCounts{jj}(2:maxLength);

        traceAtPlay=trace'./countz';

        % Plot individual trace
        plot(timeVec, traceAtPlay, 'Color', [color alphaVal], 'LineWidth', 1);

        % Store values for GLME
        traceCatchTot = [traceCatchTot, traceAtPlay];
        genoCatchTot = [genoCatchTot, repmat({geno}, 1, (maxLength-1))];
        timeCatchTot = [timeCatchTot, timeVec ];
        vidNumCatchTot = [vidNumCatchTot, repmat({sprintf('%s_%d', geno, vidNums(jj))}, 1, (maxLength-1))];
        countCatchTot = [countCatchTot,countz' ];
        scratchAreaCatchTot = [scratchAreaCatchTot, (scratchAreas(jj)* ones(1, (maxLength-1)))  ];
        startValCatchTot=[startValCatchTot, (traceAtPlay(1)* ones(1, (maxLength-1)))];
        traceMat=[traceMat;traceAtPlay];
    end

    % Plot mean trace
    % allTraces = cellfun(@(x) x(1:maxLength), traceAtPlay, 'UniformOutput', false);
    % traceMat = cell2mat(cellfun(@(x) padarray(x(:)', [0, maxLength - length(x)], NaN, 'post'), allTraces, 'UniformOutput', false)');
    meanTrace = nanmean(traceMat, 1);
    plot(timeVec, meanTrace, '-', 'Color', color, 'LineWidth', 2.5);

    % Format subplot
    xlabel('Time (Minutes)','FontSize',20);
    ylabel('Ratio of Detected Cells in Scratch','FontSize', 20);
    title(sprintf('%s Microglia Treated with %s', geno, treatmentPrefix),'FontSize',30);
    box on;
    xlim([0, maxLength*frameRate]);

 %        subplot(2,2,3:4)
 %        hold on
 %    plot(timeVec, meanTrace, '-', 'Color', color, 'LineWidth', 2.5);
 % xlabel('Time (Minutes)');
 %    ylabel('Average Distance to Scratch');
 %    title(sprintf('%s Microglia Treated with %s', geno, treatmentPrefix));
 %    box on;
 %    xlim([0, maxLength*frameRate]);

end


%% alright time for statistical testing running a glme...


% make model  table
modelData = table( ...
    traceCatchTot(:), ...
    startValCatchTot(:),...
    vidNumCatchTot(:), ...
    countCatchTot(:), ...
    scratchAreaCatchTot(:), ...
    genoCatchTot(:), ...
    timeCatchTot(:), ...
    'VariableNames', {'Dist', 'startVal', 'vidNum', 'TotCells', 'scratchArea', 'genotype', 'time'} ...
);


modelData.vidNum = cell2mat(modelData.vidNum);
modelData.genotype = cell2mat(modelData.genotype);

% Fit GLME  and alt model
glme = fitglme(modelData, ...
    'Dist ~ 1 + startVal+ scratchArea + genotype * time + (1|vidNum)', ...
    'Distribution','normal','Link','identity');


pVal = glme.Coefficients.pValue(6);
sgtitle(['Comparative Ratio of Detected Cells in Scratch Over Time p= ' num2str(pVal) ],'FontSize',45) 


%% ok now time to do the same with our individual tracings




genoTypes = {'KO', 'WT'};
genoPaths = {
    '/Volumes/PC60/iMicroVids/KOScratchAssay_6_6_25_allChannels/KO_aB/ResultsRepo', ...
    '/Volumes/PC60/iMicroVids/WT_ScratchAssay_6_6_25_allChannels_Scale/WT_aB/ResultsRepo'
};

treatmentPrefix = 'aB';
imageTypePrefix = 'Scratch';

structAtPlay = [];  % will become a table

% ------------------------- Loop over genotypes -------------------------
structAtPlay = struct([]);

for dd = 1:length(genoTypes)
    geno = genoTypes{dd};
    folderPath = genoPaths{dd};

    fileList = dir(fullfile(folderPath, '*.mat'));
    fileList = fileList(~startsWith({fileList.name}, '._'));  % remove macOS metadata files
    fileList = fileList(contains({fileList.name}, 'trajectoryMetrics'));

    for f = 1:length(fileList)
        fileName = fileList(f).name;
        fullPath = fullfile(folderPath, fileName);

        loaded = load(fullPath);
        metricsAtPlay = loaded.velocitySummary;

        % Loop over each nucleus (i.e. each entry in velocitySummary)
        for nn = 1:length(metricsAtPlay)
            structAtPlay(end+1).Genotype = geno;
            structAtPlay(end).Treatment = treatmentPrefix;
            structAtPlay(end).NucleusID = metricsAtPlay(nn).NucleusID;

            structAtPlay(end).Rate = metricsAtPlay(nn).Rate_um_per_min;
            structAtPlay(end).NetDisplacement = metricsAtPlay(nn).NetDisplacement;
            structAtPlay(end).TotalPath = metricsAtPlay(nn).TotalPath;
            structAtPlay(end).displacementRatio = metricsAtPlay(nn).displacementRatio;
            structAtPlay(end).sharpTurnCount = metricsAtPlay(nn).sharpTurnCount;
        end
    end
end


val2use=[structAtPlay.TotalPath];
netDispMat=[structAtPlay.NetDisplacement];
dispRatMat=[structAtPlay.displacementRatio];




conversionFactor = 0.65 / 1.13;
correctedAverages = nan(numel(structAtPlay), 1);  % Preallocate
numTimePoints=[];
for qq = 1:numel(structAtPlay)
    rate = structAtPlay(qq).Rate;
    
    if isnumeric(rate) && ~isempty(rate)
        corrected = rate * conversionFactor;  % Convert each array
        numTimePoints(qq)=length(corrected);
        correctedAverages(qq) = mean(corrected, 'omitnan');  % Store average
    else
        correctedAverages(qq) = NaN;  % Handle empty/missing cases
    end
end








val2use=[structAtPlay.TotalPath];
netDispMat=[structAtPlay.NetDisplacement];
dispRatMat=[structAtPlay.displacementRatio];



%%

figure

% Define color map
colorMap = struct('WT', [0.13 0.55 0.13], 'KO', [0.6 0.2 0.8]   ); 

% Extract values
val2use = [structAtPlay.displacementRatio];
filter2use=numTimePoints>84 ;
genotypes = {structAtPlay.Genotype};

% Identify groups
isWT = strcmp(genotypes, 'WT');
isKO = strcmp(genotypes, 'KO');

% Make plot
figure; hold on;

% KO Bar
b = bar(1, median(val2use(isKO & filter2use )), 'FaceColor', colorMap.KO, 'FaceAlpha', 0.3, 'BarWidth', 1.2);
scatter(rand(1, sum(isKO& filter2use )) * 0.6 + 0.7, val2use(isKO& filter2use ), ...
    'o', 'MarkerFaceColor', colorMap.KO, 'MarkerEdgeColor', colorMap.KO);

% WT Bar
b = bar(3, median(val2use(isWT& filter2use )), 'FaceColor', colorMap.WT, 'FaceAlpha', 0.3, 'BarWidth', 1.2);
scatter(rand(1, sum(isWT& filter2use )) * 0.6 + 2.7, val2use(isWT& filter2use ), ...
    'o', 'MarkerFaceColor', colorMap.WT, 'MarkerEdgeColor', colorMap.WT);

% Format
xlim([0 4])
set(gca, 'XTick', [1 3], 'XTickLabel', {'KO', 'WT'},'FontSize',15)
ylabel('Displacement Ratio (\mum)','FontSize',20)
title(' Displacement Ratio by Genotype')
box on

anovan(val2use(filter2use), {geno(filter2use)})


val2use = [structAtPlay.displacementRatio];


diffCatch=[];
count=1;
for tt=0:12:160

  filter2use=numTimePoints<=tt+12 & numTimePoints>tt ;

  if sum(isWT& filter2use )<6 | sum(isKO & filter2use)<6

      diffCatch(count)=0;
  else
   diffPoint=median(val2use(isWT& filter2use )) -   median(val2use(isKO & filter2use ));
  end
    diffCatch(count)=diffPoint;
    count=count+1;


end


figure

plot(diffCatch)



%%


genoType = {structAtPlay.Genotype};  
geno = contains((genoType), 'WT');  


anovan(netDispMat(numTimePoints<12), {geno(numTimePoints<12)})




figure
hold on

histogram(val2use(geno==1))
histogram(val2use(geno==0))





anovan(correctedAverages, {geno})
