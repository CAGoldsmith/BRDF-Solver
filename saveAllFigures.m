function [figureDirectory, success] = saveAllFigures( folderName )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if ischar(folderName)
        curTime = datestr(clock,30);
        figureDirectory = [curTime,folderName];
        mkdir(figureDirectory)
        h = get(0,'children');
        for i=1:length(h)
          saveas(h(i), [figureDirectory,'\figure',num2str(i)], 'fig');
        end
        success = true;
    else
        warning('Input to save_all_figures.m must be a string.')
        success = false;
    end
    
end

