function scratchMask = drawScratch(frames)

nFrames = length(frames); 
firstFrame = frames{1};

% Set up figure and layout
hFig = figure('Name', "Draw Scratch Region (don't forget to connect points!", 'Position', [100, 100, 1200, 600]);
tiledlayout(1, 2);

% first a video preview so you can identify the scratch
ax1 = nexttile(1);
hImg = imshow(firstFrame, 'Parent', ax1);
title(ax1, 'Live-cell vid for help identifying scratch');

% Animate a preview loop (non-blocking)
stopPreview = false;
previewTimer = timer( ...
    'ExecutionMode', 'fixedRate', ...
    'Period', 0.1, ...
    'TimerFcn', @(~,~) cycleFrames());

start(previewTimer);

% RIGHT: scratch draw panel
ax2 = nexttile(2);
imshow(firstFrame, 'Parent', ax2);
title(ax2, 'Draw boundary of scratch');
scratchROI = drawpolygon(ax2);
scratchMask = createMask(scratchROI);

% Button: Redraw
uicontrol(hFig, 'Style', 'pushbutton', ...
    'String', '🔁 Redraw Polygon', ...
    'Position', [950, 80, 150, 40], ...
    'Callback', @(~,~) redrawPolygon());

% Button: Confirm
uicontrol(hFig, 'Style', 'pushbutton', ...
    'String', '✅ Confirm Selection', ...
    'Position', [950, 20, 150, 40], ...
    'Callback', @(~,~) confirmSelection());

% % Wait for user to finish
uiwait(hFig);

    % -------- Helper Functions --------
    function cycleFrames()
        persistent idx
        if isempty(idx), idx = 1; end
        if ~stopPreview && ishandle(hImg)
            set(hImg, 'CData', frames{idx});
            idx = mod(idx, nFrames) + 1;
        end
    end

    function redrawPolygon()
        if isvalid(scratchROI)
            delete(scratchROI);
        end
        axes(ax2);  % make sure we’re on the right axis
        scratchROI = drawpolygon(ax2);
        scratchMask = createMask(scratchROI);
    end

    function confirmSelection()
        stopPreview = true;
        stop(previewTimer);
        delete(previewTimer);
        uiresume(hFig);
        close(hFig);
    end

end
