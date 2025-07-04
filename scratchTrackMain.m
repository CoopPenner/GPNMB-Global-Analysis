%scratchTrack
%This script will upload live cell videos of a nuclear stain and then track
%the dynamics of cells as they move freeely or towards a scratch.... later
%perhaps we will endeavor to capture these morphological shifts

%% initializing relevant variables and loading video
addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');
% first upload relevant video

homePath= '/Volumes/PC60/iMicroVids/KOScratchAssay_6_6_25_allChannels/KO_aB'; %Link the folder containing all videos;


        %all videos are saved in the following way Genotype_ImageType_Treatment
% (a space) (number of video).mp4
% so for example: KO_Scratch_aB (1).mp4   Note the space between aB and (1)
genoTypePrefix='KO';
imageTypePrefix='Scratch';
treatmentPrefix='aB';

outputSuffix='6_6_25';

frameRate=5; %one image every 5 minutes



%Mask Making factors
minSize=250; % this is the minimum size an object can be (some dapi debris always present)
maxSize=1200; % sometimes multiple cells become clumped (as in a swarm) I just remove these
circCut=0.05; % this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot.  
detectSensitivity=.5; % this is the sensitivity of our initial binarization step. I've found .6 is a good middle ground but decrease if there are a lot of artifactual objects
strelRadiusSize=3; %this is the pixel radius of our structuring element. We first apply an erosion of all pixels to remove noise and then proceed with dilation procedure
%a larger strel Radius Size will cause objects that are further apart ( 4
%vs 5 vs 6 pixels) to become aggregated
demoPlot=false; %this will just create a figure with the masks at each of our cleanup stages


% trajectory filtering features
maxJump=30; %this is the upper limit of how much a nucleus can move in one frame
confThreshold=.05; %this is the lower limit of our confidence in tracking a nucleus from one point to another.
areaThreshold=1; %this is how much of a difference in size from frame to frame you'll tolerate this will be protection against fusion or splitting events


%Scratch analysis factors
% Set threshold in pixels ( how far from scratch counts as "near")
proximityThreshold = 5;
%         genoTypePrefix='KO';
% imageTypePrefix='Scratch';


allFilesAtPlay = dir(fullfile(homePath, '*.mp4'));
allFilesAtPlay = allFilesAtPlay(~startsWith({allFilesAtPlay.name}, '._')); % mac keeps these weird ghost files with ._ prefix

allFilesAtPlay= allFilesAtPlay(  contains( {allFilesAtPlay.name}, genoTypePrefix)  & ...
(contains( {allFilesAtPlay.name}, imageTypePrefix))  &...
(contains( {allFilesAtPlay.name}, treatmentPrefix)));


%%start of for Loop through all vids
 for vidNum=1: length(allFilesAtPlay)      
        
        vidName=[homePath,'/', genoTypePrefix,'_',imageTypePrefix,'_',treatmentPrefix,' (',num2str(vidNum),').mp4' ];
 
         %KO_Scratch_aB (1).mp4
         vidObj = VideoReader([homePath,'/', genoTypePrefix,'_',imageTypePrefix,'_',treatmentPrefix,' (',num2str(vidNum),').mp4' ]);
    

    
        %% loading objects into frames
        
        disp('loading video data, this will take a few seconds')
        
        %preallocating objects
        
        numFrames = floor(vidObj.Duration * vidObj.FrameRate);
        frames = cell(numFrames, 1);
        refFrames=cell(numFrames, 1);
        frameIdx = 1;
        while hasFrame(vidObj)
            rgbFrame = readFrame(vidObj);
        
            % Extract the blue channel (DAPI)
            blueChannel = rgbFrame(:,:,3);  % Blue is the 3rd channel in RGB
            frames{frameIdx} = blueChannel;
            refFrames{frameIdx}=rgbFrame;
        
            frameIdx = frameIdx + 1;
        end
        
        
        
        %% ok now that we have our first points identified we will generate masks for every single frame 
        
        
        
        maskCollect=cell(numFrames, 1); %we will collect both masks and info about them
        frameDataCollect=cell(numFrames, 1);
        objNumCollect=nan(1,numFrames);
        
        
        
            if ~demoPlot
                    frames2Assess=length(frames);
                    disp('video loaded, generating masks for each frame of movie')
        
            elseif demoPlot
                frames2Assess=6;
                disp('outputting image of binarized masks')
                figure
        
            end
        
        
                for bb=1: frames2Assess
                    if ~isempty(frames{bb})
                [maskCollect{bb},CC,frameDataCollect{bb}] = scratchMaskMaker(frames{bb},minSize, maxSize, circCut, detectSensitivity,strelRadiusSize,demoPlot);
                objNumCollect(bb)=height(frameDataCollect{bb});
                   
                    end
            
               
            
                end
        
        
            %often the last frame is blank so I just remove it
            frameDataCollect(end)=[];
        
           
            %%  As things are currently implemented we lose all cells that are not present at time 0... this is obviously suboptimal 
        % there are potential solutions but for now I just want to get things
        % started. 
        
        
        baseFrame=frameDataCollect{1};
        
        centroidCatch=cell(1, height(baseFrame));
        confCatch=cell(1, height(baseFrame));
        areaCatch = cell(1, height(baseFrame));  % this will be used to filter out fusing nuclei
        
        
            for nucIdx=1:height(baseFrame)
            
                nucAtPlay=baseFrame(nucIdx,:); %these first identified nuclei will be our seeds
        
                confCatty=nan( height(frameDataCollect),1); 
                centCatty=nan( height(frameDataCollect),2); 
                areaCatty = nan(height(frameDataCollect), 1);  % per-frame area for this nucleus
        
                centCatty(1,1:2)= nucAtPlay.Centroid; %our starting point is the ground truth we take from our first centroid
                updateCentroid=nucAtPlay.Centroid;
                confCatty(1)=1; % our confidence will always be 100% for this first point
                areaCatty(1) = nucAtPlay.Area;
        
                    for frameIdx=1: (height(frameDataCollect)-1) %then we go from frame to frame and trace the seed as it moves, by taking the closest nucleus to it.
        
                        frameComp=frameDataCollect{frameIdx+1}; %this is the next frame in the movie 
                        frameCompCentroid=frameComp.Centroid; %we extract the centroid of the nucleus
        
        
                     
        
                            % Here we get the difference between our frame centroid and all other centroids in the subsequent frame
                            if frameIdx==1
                            centroidDiff=frameCompCentroid-nucAtPlay.Centroid; % if this is the first loop through frames we will use our seed as the start
                            else
                            centroidDiff=frameCompCentroid-updateCentroid; % if we have already identified a secondary Centroid this will be our new reference
                            end
            
        
                                dists = sqrt(sum(centroidDiff.^2, 2));  % Euclidean distances
                                [minDist, mostLikelyLoc] = min(dists); %here we get the location of the closest point (what we presume to be the cell in the next frame)
                                lowThree = mink(dists, 3); % other closest nuclei
           
                            
                            distConfRaw= 1- (minDist/ mean(lowThree(2)));  %now we get the difference of this difference to all other differences, this is how we'll output our confidence
                            %in situations where a microglia is "swarming" it is
                            %virtually impossible to differentiate them
            
                            updateCentroid=frameCompCentroid(mostLikelyLoc,:); % our newly identified point becomes our next comparator
            
                            confCatty(frameIdx+1)=distConfRaw; % populate our collection with confidence and centroid loc (will be used to calc velocity)
                            centCatty(frameIdx+1,1:2)= updateCentroid;
                            areaCatty(frameIdx+1) = frameComp.Area(mostLikelyLoc);  % track area
        
        
        
                    end
               centroidCatch{nucIdx}=centCatty;
               confCatch{nucIdx}=confCatty;
               areaCatch{nucIdx}=areaCatty;
        
            end
        
        %% We will now generate our filters for each trajectory. I am primarily concerned about large jumps in space and shape
        
        
        jumpFilter = cell(size(centroidCatch));
        areaFilter = cell(size(centroidCatch));
        confFilter = cell(size(centroidCatch));
        
        
        for dd=1:length(centroidCatch)
            
            trackatPlay = centroidCatch{dd};  % trajectory to eval
            confAtPlay=confCatch{dd};
            areaAtPlay=areaCatch{dd};
        
            jumpMat = nan(size(trackatPlay, 1), 1);
            areaChange=nan(size(trackatPlay, 1), 1);
            
            for bb = 1:(size(trackatPlay, 1)-1)
                pt1Coor = trackatPlay(bb, :);
                pt2Coor = trackatPlay(bb+1, :);
        
                
                if all(~isnan(pt1Coor)) && all(~isnan(pt2Coor))
                    jumpMat(bb) = norm(pt2Coor - pt1Coor);  % Euclidean distance
                end
        
                if ~isnan(areaAtPlay(bb)) && ~isnan(areaAtPlay(bb+1))
                     areaChange(bb+1) = abs(areaAtPlay(bb+1) - areaAtPlay(bb)) / areaAtPlay(bb); % percent change in area 
                    % areaChange(bb+1) = (areaAtPlay(bb+1) - areaAtPlay(bb)) / areaAtPlay(bb); % percent change in area For now we only worry about increases in size
        
                end
        
            end
        
        % adding in filters for confidence area and location jumps
        jumpPoint = find(jumpMat >= maxJump, 1);
            jumpFilterAtPlay = true(1, length(jumpMat));
            if ~isempty(jumpPoint)
                jumpFilterAtPlay(jumpPoint:end) = false;
            end
            jumpFilter{dd} = jumpFilterAtPlay;
        
            % currently filter for 33% area change
            areaFilterAtPlay = true(1, length(areaAtPlay));
            growthEvent = find( areaChange > areaThreshold , 1);
            areaFilterAtPlay(growthEvent:end) = false;
            areaFilter{dd} = areaFilterAtPlay;
        
            % conf threshold set above
            confFilterAtPlay = true(1, length(confAtPlay));
            lowConfStart = find(confAtPlay < confThreshold, 1);
            if ~isempty(lowConfStart)
                confFilterAtPlay(lowConfStart:end) = false;
            end
            confFilter{dd} = confFilterAtPlay;
        
        
        end
        
        
        
        
        
        
        
        
        %% we have now output all the basic movement data and our confidence in those scores we proceed by defining the scratch boundary
        
if ~exist([homePath,'/maskHold'])
    mkdir([homePath,'/maskHold']); %make a holder folder if it doesn't exist
end

    %once we make the scratch once it will be reloaded for future use
    if ~exist([homePath,'/maskHold/', genoTypePrefix,'_',imageTypePrefix,'_',treatmentPrefix,' (',num2str(vidNum),') MaskOutput.mat' ])
    
         scratchMask = drawScratch(refFrames);

         save([  homePath,'/maskHold','/', genoTypePrefix,'_',imageTypePrefix,'_',treatmentPrefix,' (',num2str(vidNum),') MaskOutput.mat' ],'scratchMask');
    else

        tmp=load([homePath,'/maskHold/', genoTypePrefix,'_',imageTypePrefix,'_',treatmentPrefix,' (',num2str(vidNum),') MaskOutput.mat' ]);
        scratchMask=tmp.scratchMask; %idk why matlab makes loading come out as a structure instead of just a matrix...
    end
        
        close all
            %% alright, finally let's calculate some of our metrics
        
        
            % first velocity, we will take in the direction of our scratch object
            % to be positive and away to be negative. 

        velocitySummary = [];
        speedTotal=cell(1,length(centroidCatch));
        % signedVelocityTotal=cell(1,length(centroidCatch)); % I'm having
        % some trouble with directional velocity, will just output rate for
        % now


        count=1;
        
        for bb = 1:length(centroidCatch)
            coords = centroidCatch{bb};  % Nx2
            jumpy = jumpFilter{bb};      
            conf=confFilter{bb};
            boxedIn=areaFilter{bb};
        
        
            validIdx = jumpy & conf & boxedIn;
            coordsFilt = coords(validIdx, :);
        
            % Skip short (less than 30 min) or invalid tracks
            if size(coordsFilt,1) < 6 || any(isnan(coordsFilt(:)))
                continue
            end
        
        
        
        
        % Euclidean velocity
        diffs = diff(coordsFilt);
        dists = vecnorm(diffs, 2, 2);  % Euclidean step sizes
        velocity_px_per_frame = dists;
        
        
        
          % how I calculate distance
            pxSize=.65;
        % typical camera pixel size is 6.5 uM divide that by ten... for now
        % let's just do 1 since I'm not super clear on this conversion easy
        % enough to fix downstream. 
        velocity_um_per_min = velocity_px_per_frame * pxSize / frameRate;
        % signed_velocity_um_per_min=signed_velocity*pxSize/frameRate;
        
        



% Parameters
framesPerBin = 24;  % 2 hours = 24 frames
numBins = floor(size(coordsFilt, 1) / framesPerBin);  % Max valid bins

% Preallocate output
dispRatioBins = nan(1, 7);  % 7 bins max for 14-hour recording

for binIdx = 1:min(7, numBins)
    % Define range for this bin
    startFrame = (binIdx - 1) * framesPerBin + 1;
    endFrame = binIdx * framesPerBin;

    coordsBin = coordsFilt(startFrame:endFrame, :);

    if any(isnan(coordsBin(:)))
        continue  % skip bin if any NaNs
    end

    % Use first point as origin
    startPos = coordsBin(1, :);
    coordsZeroed = coordsBin - startPos;

    % Net displacement (relative to start of bin)
    netDisp_px = norm(coordsZeroed(end, :));

    % Total path length (sum of step distances)
    diffs_bin = diff(coordsBin, 1, 1);
    totalPath_px = sum(vecnorm(diffs_bin, 2, 2));

    % Displacement ratio
    if totalPath_px > 0
        dispRatioBins(binIdx) = netDisp_px / totalPath_px;
    end
end

% Add to output struct
velocitySummary(count).DisplacementRatioPer2hr = dispRatioBins;






        
            % Net displacement
            netDisp_px = norm(coordsFilt(end,:) - coordsFilt(1,:));
            totalPath_px = sum(dists);
            displacementRatio = netDisp_px / totalPath_px;
        
            %  Directional switches (angle between steps)
            unitVecs = diffs ./ vecnorm(diffs, 2, 2);  % normalize step vectors
            angles = acosd(dot(unitVecs(1:end-1,:), unitVecs(2:end,:), 2));  % angle between steps
        
            % Count sharp turns (>90°)
            sharpTurnCount = sum(angles > 90);
        

            velocitySummary(count).NucleusID=count; %lol
            velocitySummary(count).Rate_um_per_min=velocity_um_per_min; %rate not accounting for directionality
            velocitySummary(count).NetDisplacement=netDisp_px * pxSize;
            velocitySummary(count).TotalPath=totalPath_px*pxSize;
            velocitySummary(count).displacementRatio=displacementRatio;
            velocitySummary(count).sharpTurnCount= sharpTurnCount; 
            count=count+1;
        end
        
          %save output structure

if ~exist([homePath,'/ResultsRepo'])
    mkdir([homePath,'/ResultsRepo']); %make a holder folder if it doesn't exist
end


        save([homePath,'/ResultsRepo/', genoTypePrefix,'_',imageTypePrefix,'_',treatmentPrefix,' (',num2str(vidNum),') trajectoryMetrics.mat'],'velocitySummary' )

%%

        % 
        %         velMat=[];
        % 
        %     for dd=1:length(signedVelocityTotal)
        % 
        %         velAtPlay=signedVelocityTotal{dd};
        % 
        %         if ~isempty(velAtPlay)
        % 
        % velMat=[velMat; [velAtPlay', nan(1, (length(frameDataCollect))-length(velAtPlay)  )]];
        %         end
        % 
        %     end
        % 
        % 
        
        %     if demoPlot
        %  figure 
        %  hold on
        % 
        %  plot(nanmean(velMat))
        % 
        % 
        % 
        % % Compute mean and std
        % velAvg = nanmean(velMat, 1);
        % velStd = nanstd(velMat, 0, 1);
        % 
        % timeVec = (0:length(velAvg)-1) * frameRate;
        % 
        % % Remove any NaNs for fill plotting
        % validIdx = ~isnan(velAvg) & ~isnan(velStd);
        % timeClean = timeVec(validIdx);
        % meanClean = velAvg(validIdx);
        % stdClean = velStd(validIdx);
        % 
        % % Plot
        % figure; hold on;
        % 
        % % Fill shaded region for ±1 std
        % fill([timeClean fliplr(timeClean)], ...
        %      [meanClean + stdClean, fliplr(meanClean - stdClean)], ...
        %      [0.3 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        % 
        % % Plot mean line
        % plot(timeClean, meanClean, 'b-', 'LineWidth', 2);
        % 
        % xlabel('Time (min)');
        % ylabel('Velocity (\mum/min)');
        % title('Average Nucleus Velocity Over Time');
        % legend('±1 SD', 'Mean Velocity', 'Location', 'northeast');
        %     end
        
        %% straightforward Scratch Assay analysis
        
        % number of cells present within 50 um of scratch at T0 and Tend
        
        
        
        % Preallocate structures
        inScratchCount = zeros(length(frameDataCollect), 1);
        nearScratchCount = zeros(length(frameDataCollect), 1);
        totObjects = zeros(length(frameDataCollect), 1);
        avgDistTot= zeros(length(frameDataCollect), 1);
        stdDistTot= zeros(length(frameDataCollect), 1);
        scratchMet=[];
        distFromScratch = bwdist(scratchMask);  % each pixel = distance in pixels from nearest scratch pixel for base mask
        
        for bb = 1:length(frameDataCollect)
            frameData = frameDataCollect{bb};
            objInScratch = 0;
            objNearScratch = 0;
            centroidDistHold=[];
        
        
            for objIdx = 1:height(frameData)
                pixIdxs = frameData.PixelIdxList{objIdx};  % get pixel idx for this object
                pixList = frameData.PixelList{objIdx};  % and straight pixel number
        
                % Get distances of this object’s pixels from scratch
                objDistances = distFromScratch(pixIdxs);
                centroidDistHold(objIdx)=mean(objDistances);
        
        
                %  any pixel inside scratch
                if any(scratchMask(pixIdxs))
                    objInScratch = objInScratch + 1;
                    objNearScratch = objNearScratch + 1;
        
                %  any pixel within the proximity threshold
                elseif any(objDistances <= proximityThreshold)
                    objNearScratch = objNearScratch + 1;
                end
        
            end
        
            inScratchCount(bb) = objInScratch;
            nearScratchCount(bb) = objNearScratch;
            totObjects(bb)=objIdx; %capturing number of identified objects per frame for normalization purposes
            avgDistTot(bb)=nanmean(centroidDistHold); % this is the average distance from the scratch for each identified nucleus
            stdDistTot(bb)=nanstd(centroidDistHold);
        end
        
        scratchMet.inScratchCount=inScratchCount;
        scratchMet.nearScratchCount=nearScratchCount;
        scratchMet.totObjects=totObjects;
        scratchMet.avgDistTot=avgDistTot;
        scratchMet.stdDistTot=stdDistTot;
        scratchMet.scratchSize= (sum(scratchMask(:))) * pxSize^2 ;


        save([homePath,'/ResultsRepo/', genoTypePrefix,'_',imageTypePrefix,'_',treatmentPrefix,' (',num2str(vidNum),') scratchMetrics.mat' ],'scratchMet')
% 
%                             figure
%                         plot([1:length(avgDistTot)],avgDistTot );
%                         title('avg Dist')
% 
% 
%                  figure
%         plot([1:length(nearScratchCount)],nearScratchCount)
% title('nearScratch Dist')
% 

        % 
        % normFactor= mean(  [mean(totObjects(1:3)), mean(totObjects(end-3:end)) ]  );
        % 
        % inScratchCountNorm=inScratchCount./ mean(totObjects(1:3));
        % nearScratchCountNorm=nearScratchCount./ mean(totObjects(1:3));
        % 
        % inScratchCountTotNorm=inScratchCount- mean(totObjects(1:3));
        % 
        % nearScratchCountTotNorm=nearScratchCount- mean(totObjects(1:3));
        % 
         if demoPlot
        
                        figure
                        plot([1:length(inScratchCountNorm)],inScratchCountNorm)
         
                        figure
                        plot([1:length(nearScratchCountNorm)],nearScratchCountNorm)
                      
                        figure
                        plot([1:length(avgDistTot)],avgDistTot );
                        
                        figure
                        plot([1:length(nearScratchCountNorm)],nearScratchCountNorm)
        
                % Plot
                figure; hold on;
                timeVec = (0:length(avgDistTot)-1) * frameRate;
        
                % % Fill shaded region for ±1 std
                % fill([timeVec fliplr(timeVec)], ...
                %      [(avgDistTot + stdDistTot)', fliplr(avgDistTot - stdDistTot)'], ...
                %      [0.3 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                
                % Plot mean line
                plot(timeVec, avgDistTot, 'b-', 'LineWidth', 2);
                
                xlabel('Time (min)');
                ylabel('Avg Dist of all nuclei to scratch');
                title('Avg Dist of all nuclei to scratch over imaging session');
        end
             






        
        %% output image of all trajectories, these should look clean! Reevaluate if not
        if demoPlot
                figure
                imshow(frames{1}); hold on;
                for bb = 1:length(centroidCatch)
                    coords = centroidCatch{bb};
                    conf=confFilter{bb};
                    jumpy=jumpFilter{bb};
                    boxedIn=areaFilter{bb};
                
                % if sum(confFilter)< sum(jumpy)
                %     catchr=1;
                % end
                
                    if sum(jumpy)>=6 && sum(conf)>=6 && sum(boxedIn)>=6 % if at least 10 time points pass the filter
                    plot(coords(jumpy & conf &boxedIn,1), coords(jumpy& conf&boxedIn,2), '-o');  % filtered Tracks
                    end
                end
                
                title('All Trajectories over video duration')
                
           %%     
            % lets make a video of one of our trajectories being generated. 
            
            viableExamp=[];
            for dd=1:length(confFilter)
            
            conf = confFilter{dd};
            jumpy = jumpFilter{dd};
            boxedIn=areaFilter{dd};
            validFrames = find(jumpy & conf & boxedIn);
                if length(validFrames)>12
                viableExamp=[viableExamp,dd];
                end
            
            
            end
            
            
            % ---- Inputs ----
            nucleusID = 135;
            conf = confFilter{nucleusID};
            jumpy = jumpFilter{nucleusID};
            boxedIn=areaFilter{nucleusID};
            
            
            % Combine filters
            validFrames = find(jumpy & conf & boxedIn);
            if length(validFrames) < 12 && sum(jumpy) <12
                error(['Not enough valid timepoints for nucleus ', num2str(nucleusID), ' jump detected early']);
            elseif length(validFrames) < 12 && sum(boxedIn) <12
                error(['Not enough valid timepoints for nucleus ', num2str(nucleusID), ' fusion mischaracterization detected early']);
            elseif length(validFrames) < 12 && sum(conf) <12
                error(['Not enough valid timepoints for nucleus ', num2str(nucleusID), ' ambiguous characterization detected early']);
            end
            
            frameRange = min(validFrames):max(validFrames);
            track = centroidCatch{nucleusID};  % Nx2 [X, Y] coordinates
            
            % ---- Display Setup ----
            figure('Name','Trajectory Playback with Filters');
            tiledlayout(1, 2);
            
            % Left: All nuclei (binarized)
            ax1 = nexttile;
            title(ax1, 'All Nuclei (Binarized)');
            h1 = imshow(maskCollect{1}, 'Parent', ax1);
            
            % Right: Single nucleus + trail
            ax2 = nexttile;
            title(ax2, sprintf('Nucleus %d Filtered Trajectory', nucleusID));
            h2 = imshow(false(size(maskCollect{1})), 'Parent', ax2);  % blank mask
            hold(ax2, 'on');
            trailPlot = plot(ax2, nan, nan, '.', 'Color', 'r');
            
            %animate
            for frameIdx = frameRange
                if ~(jumpy(frameIdx) && conf(frameIdx))
                    continue
                end
            
                frameData = frameDataCollect{frameIdx};
                thisCoord = track(frameIdx, :);
                fullMask = maskCollect{frameIdx};
                
                % Prepare nucleus mask
                if any(isnan(thisCoord))
                    nucleusMask = false(size(fullMask));
                else
                    dists = vecnorm(frameData.Centroid - thisCoord, 2, 2);
                    [~, objIdx] = min(dists);
                    pixelIdx = frameData.PixelIdxList{objIdx};
                    nucleusMask = false(size(fullMask));
                    nucleusMask(pixelIdx) = true;
                end
            
                % ---- LEFT tile: RGB mask with highlighted nucleus ----
                rgbMask = repmat(uint8(fullMask) * 255, 1, 1, 3);  % white binary image
                % Turn the tracked nucleus blue: [R G B]
                rgbMask(:,:,1) = rgbMask(:,:,1) .* uint8(~nucleusMask);  % remove red
                rgbMask(:,:,2) = rgbMask(:,:,2) .* uint8(~nucleusMask);  % remove green
                rgbMask(:,:,3) = max(rgbMask(:,:,3), uint8(nucleusMask) * 255);  % boost blue
            
                set(h1, 'CData', rgbMask);
            
                % ---- RIGHT tile: single nucleus only ----
                set(h2, 'CData', nucleusMask);
            
                % Update trail
                trailCoords = track(validFrames(validFrames <= frameIdx), :);
                validTrail = all(~isnan(trailCoords), 2);
                trailPlot.XData = trailCoords(validTrail, 1);
                trailPlot.YData = trailCoords(validTrail, 2);
            
                drawnow;
                pause(0.1);
            end
        
        
        end
        
  end        
        
        
        




%% code cemetery

%%
% 
% 


% 
% % STEP 2: Distance to scratch edge at each timepoint
% % (Precomputed scratch boundary coordinates)
% 
% [yScratch, xScratch] = find(scratchMask);                  % row, col (i.e., y, x)
% scratchCoords = [xScratch, yScratch];                      % Mx2, in [x, y] format
% 
% % Ensure coordsFilt is valid and rounded
% coordsValid = (coordsFilt);                           % make sure coords are integers
% coordsValid(any(isnan(coordsValid), 2), :) = 1;            % placeholder for NaNs
% 
% % For each nucleus position, find the closest scratch pixel and get distance
% nearestIdx = knnsearch(scratchCoords, coordsValid);        % 1-nearest neighbor index
% nearestScratchPoints = scratchCoords(nearestIdx, :);       % matched scratch points
% 
% % Now compute the distance
% distToScratchVec = vecnorm(coordsValid - nearestScratchPoints, 2, 2);  % Nx1 vector
% 
% 
% % STEP 3: Signed velocity
% deltaDist = diff(distToScratchVec);   % >0 = away, <0 = toward
% signed_velocity = -deltaDist;         % flip sign: toward = positive
% 
% % STEP 4: Mask signed velocity if inside scratch
% inScratch = scratchMask(sub2ind(size(scratchMask), ...
%     round(coordsFilt(:,2)), round(coordsFilt(:,1))));
% inScratch = logical(inScratch);
% signed_velocity(inScratch(1:end-1) | inScratch(2:end)) = NaN;



% 
% figure
% Y = fft(velAvg(1:end-1) - nanmean(velAvg(1:end-1)));
% P2 = abs(Y/length(Y));
% P1 = P2(1:floor(length(Y)/2));
% f = (1/(frameRate*length(velAvg(1:end-1)))) * (0:floor(length(Y)/2)-1);  % in cycles/min
% 
% figure; plot(f, P1); xlabel('Frequency (1/min)'); ylabel('Amplitude');
% title('FFT of Velocity Time Series');
% 




% %%
% nFrames = 10; %how many times ref video will loop
% frameHeight = size(frames{1},1);
% frameWidth = size(frames{1},2);
% 
% figure('Name', 'Annotate Scratch Boundary');
% 
% % Set up two side-by-side axes
% tiledlayout(1,2)
% 
% % --- LEFT: Show dynamic video (loop a few frames to give context) ---
% ax1 = nexttile;
% hImg = imshow(frames{1}, 'Parent', ax1);
% title(ax1, 'Live-cell video (so you can see where scratch is)');
% 
% numLoops = 2;
% pauseTime = 0.1;
% 
% for loopIdx = 1:numLoops
%     for f = 1:nFrames
%         set(hImg, 'CData', frames{f});
%         drawnow;
%         pause(pauseTime);
%     end
% end
% 
% % Hold last frame for drawing
% set(hImg, 'CData', frames{1});
% 
% % now actually drawing the polygon
% ax2 = nexttile;
% imshow(frames{1}, 'Parent', ax2);
% title(ax2, "Draw boundary of scratch (drag corners and don't forget to close shape!");
% 
% scratchROI = drawpolygon(ax2);  % interactive
% scratchMask = createMask(scratchROI);  % binary mask
    
    
    
