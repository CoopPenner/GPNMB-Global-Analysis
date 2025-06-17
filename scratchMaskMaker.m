function [cleanedMask,CC,nucFeatures] = scratchMaskMaker(frameAtPlay,minSize, maxSize, circCut, detectSensitivity,strelRadiusSize,demoPlot)



%convert to grayscale if necessary
        if size(frameAtPlay,3)==3
        frameAtPlay=rgb2gray(frameAtPlay);
        warning('running algorithm on rgb image, consider using only nuclear channel')
        end

% Binarize the  first frame of the image
nuclearMask = imbinarize(frameAtPlay, 'adaptive', 'Sensitivity', detectSensitivity); % I use a relatively low sensitivity


se = strel('disk',strelRadiusSize);  % now using strel to fill in small gaps between different parts of a given nuclei. 


nuclearMaskOpened = imopen(nuclearMask, se);  % remove specks first
nuclearMaskClosed = imclose(nuclearMaskOpened, se);  % then close gaps



% Apply initial size filtering
cleanedMask = bwareaopen(nuclearMaskClosed, minSize);


% Get connected components and nuclear features
CC = bwconncomp(cleanedMask);
nucFeatures = regionprops("table", CC, ...
    "Area", "Centroid", "EquivDiameter", "Solidity", "Circularity");

% Extract features
areaAtPlay = nucFeatures.Area;
circ = nucFeatures.Circularity;

% Define  nuclei more likely to be artifactual
badIdx = areaAtPlay > maxSize | circ < circCut;
% how I calculate distance


% Set those components to false in the binary mask
labeledMask = labelmatrix(CC);
for ss = find(badIdx)'
    cleanedMask(labeledMask == ss) = 0;
end


% regenerate connected components


CC = bwconncomp(cleanedMask);

nucFeatures = regionprops("table", CC, ...
    "Area", "Centroid", "EquivDiameter", "Solidity", "Circularity","PixelIdxList","PixelList");








if demoPlot

figure

subplot(2,2,1)
imshow(frameAtPlay)

subplot(2,2,2)
imshow(nuclearMask)

subplot(2,2,3)
imshow(nuclearMaskClosed)

subplot(2,2,4)
imshow(cleanedMask)
end



end


%% initial cell selection... this will be the mask that we use to track

% my initial approach will be to prioritize eliminating low intensity
% nuclei and differentiating nearby cells. Ie removing false positives 

% avgFrame = (double(rgb2gray(frames{1})) + double(rgb2gray(frames{2})))  / 2;