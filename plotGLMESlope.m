function plotGLMESlope(estimate, lowerCI, upperCI,startPt,colorCode)
    % simple function that plots the GLME estimate and its 95% confidence
    % interval, eventually I might do this with residuals in a more
    % thoughtful way...

    % Define x values (time progression)
    x = linspace(0, 120, 100); % Adjust as needed
    
    % Compute the estimated line
    y = (estimate * x) +startPt  ;

    % Compute constant confidence interval range
    ci_range = abs(upperCI - lowerCI) ; 
    y_lower = y - ci_range;
    y_upper = y + ci_range;
    
    % Plot the estimated slope
    plot(x, y, 'b', 'LineWidth', 2,'Color',colorCode);
    
    % Plot confidence interval as shaded region
    fill([x, fliplr(x)], [y_lower, fliplr(y_upper)], colorCode, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
