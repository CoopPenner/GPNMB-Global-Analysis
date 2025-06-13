%scratchTrack
%This script will upload live cell videos of a nuclear stain and then track
%the dynamics of cells as they move freeely or towards a scratch.... later
%perhaps we will endeavor to capture these morphological shifts

%% initializing relevant variables and loading video
addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');
% first upload relevant video



vidObj = VideoReader('/Volumes/PC60/iMicroVids/WT_ScratchAssay_6_6_25_allChannels_Scale/WT_Scratch_aB (5).mp4');

frameRate=5; %one image every 5 minutes

minSize=50; % this is the minimum size an object can be (some dapi debris always present)
maxSize=300; % sometimes multiple cells become clumped (as in a swarm) I just remove these
circCut=.5; % this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot.  
detectSensitivity=.1; % this is the sensitivity of our initial binarization step. I've found .6 is a good middle ground but decrease if there are a lot of artifactual objects
strelRadiusSize=2; %this is the pixel radius of our structuring element. We first apply an erosion of all pixels to remove noise and then proceed with dilation procedure
%a larger strel Radius Size will cause objects that are further apart ( 4
%vs 5 vs 6 pixels) to become aggregated
demoPlot=false; %this will just create a figure with the masks at each of our cleanup stages
frameRateMin=5; % the sampling rate in minutes


disp('loading video data, this will take a few seconds')




%% loading objects into frames


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
        frames2Assess=20;
        disp('outputting image of binarized masks')

    end



        for bb=1: length(frames)
            if ~isempty(frames{bb})
        [maskCollect{bb},CC] = scratchMaskMaker(frames{bb},minSize, maxSize, circCut, detectSensitivity,strelRadiusSize,demoPlot);
        frameDataCollect{bb} = regionprops("table", CC, ...
        "Area", "Centroid", "EquivDiameter", "Solidity","PixelIdxList","PixelList");
        objNumCollect(bb)=height(frameDataCollect{bb});
           
            end
    
             if demoPlot
                break
            end
    
        end


    %often the last frame is blank so I just remove it
    frameDataCollect(end)=[];

   
    %%  As things are currently implemented we lose all cells that are not present at time 0... this is obviously suboptimal 
% there are potential solutions but for now I just want to get things
% started. 

maxJump=30; %this is the upper limit of how much a nucleus can move in one frame
confThreshold=.1; %this is the lower limit of our confidence in tracking a nucleus from one point to another.
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
    fusionEvent = find(areaChange > areaThreshold, 1);
    areaFilterAtPlay(fusionEvent:end) = false;
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

scratchMask = drawScratch(frames);



%% output image of all trajectories, these should look clean! Reevaluate if not


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
%% lets make a video of one of our trajectories being generated. 



% ---- Inputs ----
nucleusID = 94;
conf = confFilter{nucleusID};
jumpy = jumpFilter{nucleusID};
boxedIn=areaFilter{nucleusID};


% Combine filters
validFrames = find(jumpy & conf & boxedIn);
if length(validFrames) < 4 && sum(jumpy) <4
    error(['Not enough valid timepoints for nucleus ', num2str(nucleusID), ' jump detected early']);
elseif length(validFrames) < 4 && sum(boxedIn) <4
    error(['Not enough valid timepoints for nucleus ', num2str(nucleusID), ' fusion mischaracterization detected early']);
elseif length(validFrames) < 4 && sum(conf) <4
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

% ---- Animate ----
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






    %% alright, finally let's calculate some of our metrics


    % first velocity, we will take in the direction of our scratch object
    % to be positive and away to be negative. 



velocitySummary = table();
velocityTotal=cell(1,length(centroidCatch));

for bb = 1:length(centroidCatch)
    coords = centroidCatch{bb};  % Nx2
    conf = confCatch{bb};
    jumpy = jumpFilter{bb};      
    conf=confFilter{bb};
    boxedIn=areaFilter{bb};


    validIdx = jumpy & conf & boxedIn;
    coordsFilt = coords(validIdx, :);

    % Skip short (less than 40 min) or invalid tracks
    if size(coordsFilt,1) < 4 || any(isnan(coordsFilt(:)))
        continue
    end

    % --- 1. Euclidean velocity
    diffs = diff(coordsFilt);  % (N-1)x2
    dists = sqrt(sum(diffs.^2, 2));  % pixel distances
    velocity_px_per_frame = dists;

  % how I calculate distance

fov_um = 1300;        % width of field of view in microns (approximate)
image_px = height(frames{1});      % width of image in pixels
pxSize = fov_um / image_px;
% ~1.13 µm/px
velocity_um_per_min = velocity_px_per_frame * pxSize / frameRateMin;



    % --- 2. Net displacement
    netDisp_px = norm(coordsFilt(end,:) - coordsFilt(1,:));
    totalPath_px = sum(dists);
    displacementRatio = netDisp_px / totalPath_px;

    % --- 3. Directional switches (angle between steps)
    unitVecs = diffs ./ vecnorm(diffs, 2, 2);  % normalize step vectors
    angles = acosd(dot(unitVecs(1:end-1,:), unitVecs(2:end,:), 2));  % angle between steps

    % Count sharp turns (>90°)
    sharpTurnCount = sum(angles > 90);
    meanAngle = mean(angles);

    % --- 4. Save metrics
    velocitySummary = [velocitySummary; table( ...
        bb, ...
        mean(velocity_um_per_min), ...
        netDisp_px * pxSize, ...
        totalPath_px * pxSize, ...
        displacementRatio, ...
        sharpTurnCount, ...
        meanAngle, ...
        'VariableNames', { ...
            'NucleusID', ...
            'MeanVelocity_um_per_min', ...
            'NetDisplacement_um', ...
            'TotalPath_um', ...
            'DisplacementRatio', ...
            'NumSharpTurns', ...
            'MeanTurnAngle_deg' ...
        })];


        velocityTotal{bb}=velocity_um_per_min;
end

       
      
        velMat=[];

    for dd=1:length(velocityTotal)
        
        velAtPlay=velocityTotal{dd};
        
        if ~isempty(velAtPlay)
        
velMat=[velMat; [velAtPlay', nan(1, (length(frameDataCollect))-length(velAtPlay)  )]];
        end
        
    end


 figure 
 hold on

 plot(nanmean(velMat))



% Compute mean and std
velAvg = nanmean(velMat, 1);
velStd = nanstd(velMat, 0, 1);

timeVec = (0:length(velAvg)-1) * frameRateMin;

% Remove any NaNs for fill plotting
validIdx = ~isnan(velAvg) & ~isnan(velStd);
timeClean = timeVec(validIdx);
meanClean = velAvg(validIdx);
stdClean = velStd(validIdx);

% Plot
figure; hold on;

% Fill shaded region for ±1 std
fill([timeClean fliplr(timeClean)], ...
     [meanClean + stdClean, fliplr(meanClean - stdClean)], ...
     [0.3 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot mean line
plot(timeClean, meanClean, 'b-', 'LineWidth', 2);

xlabel('Time (min)');
ylabel('Velocity (\mum/min)');
title('Average Nucleus Velocity Over Time');
legend('±1 SD', 'Mean Velocity', 'Location', 'northeast');













%% code cemetery

%%
% 
% 
% Y = fft(velAvg - nanmean(velAvg));
% P2 = abs(Y/length(Y));
% P1 = P2(1:floor(length(Y)/2));
% f = (1/(frameRateMin*length(velAvg))) * (0:floor(length(Y)/2)-1);  % in cycles/min
% 
% figure; plot(f, P1); xlabel('Frequency (1/min)'); ylabel('Amplitude');
% title('FFT of Velocity Time Series');





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
    
    
    
