function [timePassed] = timeOutput(time2use)




    %converting to months
    days2use=hours(time2use)/24; weeks=days2use/7; timePassed=weeks/4;
    if timePassed==0
    timePassed=nan;
    end