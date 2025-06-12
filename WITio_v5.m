
%-------------------------------------------------------------------------%
% This matlab code requires the plugin described in Holmi, J.T., and
% Lipsanen, H., 2022, WITio: A MATLAB data evaluation toolbox to script
% broader insights into big data from WITec microscopes: SoftwareX, v. 18, p. 101009, 
% doi:10.1016/j.softx.2022.101009.
% 
% code by Christina Cauley June 06 2025
%-------------------------------------------------------------------------%
% Temporarily set user preferences
resetOnCleanup = WITio.tbx.pref.set({'wip_AutoCreateObj', 'wip_AutoCopyObj', 'wip_AutoModifyObj'}, {true, true, true}); % The original values prior to this call are restored when resetOnCleanup-variable is cleared.
%set(0,'DefaultFigureWindowStyle','normal');

WITio.tbx.edit(); % Open this code in Editor
close all; % Close figures

% Define WIP file to open
pathstr = '/Users/christinacauley/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab_Work/Olivine Deformation/Raman/Trip 3'; % Get folder of this script
file = fullfile(pathstr, 'June06_Iki.wip'); % Construct full path of the demo file 

% Open density results file
pathstr2 = '/Users/christinacauley/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab_Work/Olivine Deformation/Raman'; % Get folder of this script
file_densities = fullfile(pathstr2, 'Iki_FI_fitted_2025-06-06.xlsx'); % Construct full path of the demo file 
T = readtable(file_densities);

savePath = fullfile(pathstr2, 'Iki_FI_fitted_2025-06-06_withNotes.xlsx');
%% Use to check for non-ASCII friendly characters you'll
% Get original variable names
varNames = T.Properties.VariableNames;

for i = 1:length(varNames)
    name = varNames{i};
    name = strrep(name, 'σ', 'SE');
    name = strrep(name, '%', 'percent');
    name = strrep(name, char(181), 'u');
    name = strrep(name, '(', '');
    name = strrep(name, ')', '');
    varNames{i} = matlab.lang.makeValidName(name);
end

% Assign cleaned names back to table
T.Properties.VariableNames = varNames;

%% -------------------------------------------------------------------------%
[O_wid, O_wip, O_wit] = WITio.read(file, '-all');
% O_wid.manager  % Open Project Manager - use this command to view interactive Project Manager viewer similar to WITio proprietary 

%% ------------------------------ GET WIN INDEXs -------------------------------------------%
num_obj=length(O_wid);
wid_IndTag = strings(1,num_obj);

% loop though each object in o_Wid to save filenames
for i=1:num_obj
    wid_IndTag(i)=string(O_wid(i).Name);
end

% loop through density table to pull Win index
suffixes = {'', '_05x', '_20x', '_50x'};
colNames = {'point_num', 'x05x_num', 'x20x_num', 'x50x_num'};

for i = 1:height(T)
    % Get base filename string
    if iscell(T.filename)
        baseName = T.filename{i};
    else
        baseName = T.filename(i); 
    end
    baseName = char(baseName);  % Ensure it's a character vector

    % If it contains '_15mw', use it only for point_num
    if contains(baseName, '_15mw')
        % point_num gets full name
        searchName = baseName;
        idx = find(strcmp(wid_IndTag, searchName), 1);
        if isempty(idx)
            T.point_num(i) = NaN;
        else
            T.point_num(i) = idx;
        end

        % Remove '_15mw' for all other suffix matches
        trimmedName = regexprep(baseName, '_15mw', '');

        for j = 2:length(suffixes)  % start at 2 to skip point_num
            searchName = [trimmedName, suffixes{j}];
            idx = find(strcmp(wid_IndTag, searchName), 1);
            if isempty(idx)
                T.(colNames{j})(i) = NaN;
            else
                T.(colNames{j})(i) = idx;
            end
        end
    else
        % Normal case: match all with full baseName
        for j = 1:length(suffixes)
            searchName = [baseName, suffixes{j}];
            idx = find(strcmp(wid_IndTag, searchName), 1);
            if isempty(idx)
                T.(colNames{j})(i) = NaN;
            else
                T.(colNames{j})(i) = idx;
            end
        end
    end
end


if ~ismember('notes_FI', T.Properties.VariableNames)
    T.notes_FI = cell(height(T), 1);  % Correct: preallocate as cell array
    T.notes_FI = repmat({''}, height(T), 1);  % convert to cell array of char
end

%% -------------------------------------------------------------------------%
% Ask user whether to start from beginning or resume
startChoice = menu('Start loop from:', ...
                   'Beginning (row 1)', ...
                   'First row with missing r_um');

if startChoice == 1
    startIdx = 1;
else
    % Find first row where r_um is missing, NaN, or 0
    if iscell(T.r_um)
        isMissing = cellfun(@isempty, T.r_um);
    else
        isMissing = isnan(T.r_um) | T.r_um == 0;
    end
    missingIdx = find(isMissing, 1);

    if isempty(missingIdx)
    disp('No missing, NaN, or zero r_um values found. Script cancelled.');
    return;  % Exit script

    else
    startIdx = missingIdx;
    disp(['Resuming from row ', num2str(startIdx), ' where r_um is missing/NaN/0.']);
    end

end

for i = startIdx:height(T)
    image_num = T.x50x_num(i);
    meta_num = image_num+1;
    point_num = T.point_num(i);

    O_Bitmap = O_wid(image_num); % Get image
    %-------------------------------------------------------------------------%
    % Display image
    figure(2); h1a = O_Bitmap.plot('-scalebar'); % Show positions AND show default bottom-left auto-length scalebar with text
    %figure(2); h1a = O_Bitmap.plot_position(O_wid([point_num])); % plot analytical point
    set(findall(h1a, 'Type', 'line'), 'LineWidth', 3); % Set line color to white and its width to 3
    try
        figure(2); h1a = O_Bitmap.plot_position(O_wid([point_num]));
    catch ME
        warning('Skipping plot_position due to error: %s', ME.message);
        continue;  % Skip to next iteration of the main loop
    end

    %-------------------------------------------------------------------------%
    % Get scaling data 
    O_Text = O_wid(meta_num); % import metadata file
    um_width=str2double(O_Text.Data{13,2});
    um_height=str2double(O_Text.Data{14,2});
    px_width = size(O_Bitmap.Data,1);
    px_height = size(O_Bitmap.Data,2);
    Xscale = um_width/px_width;
    Yscale = um_height/px_height;
    
    %% -------------------------------------------------------------------------%
    %set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    figure(2);  % make full screen
    %set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    drawnow;
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    %% -------------------------------------------------------------------------%
    % Prompt user to choose method for current image
    choice = menu(['Choose radius extraction method for image ', T.filename{i}], ...
                  'Option A: From analytical point outward (50x)', ...
                  'Option B: From analytical point outward (20x)', ...
                  'Option C: Arbitrary edge-to-edge', ...
                  'Stop loop',...
                  'Skip');
    
    if choice == 1
        % ----------------------------- Option A ----------------------------- %
        p0 = [h1a.XData, h1a.YData];
        num_points = 8;
        r_values = zeros(1, num_points);
    
        for j = 1:num_points
            pi = ginput(1);
            lineSec = [p0; pi];
            h = line(lineSec(:,1), lineSec(:,2), 'linewidth', 2, 'DisplayName', ['Segment ' num2str(j)]);
    
            delta_x = (h.XData(end) - h.XData(1)) * Xscale;
            delta_y = (h.YData(end) - h.YData(1)) * Yscale;
            r = sqrt(delta_x^2 + delta_y^2);
            r_values(j) = r;
    
            mid_x = mean(h.XData);
            mid_y = mean(h.YData);
            text(mid_x, mid_y, sprintf('%.2f', r), ...
                'FontSize', 16, 'Color', 'white', 'HorizontalAlignment', 'center');
        end
    
    elseif choice == 2
        % ----------------------------- Option B ----------------------------- %
        close(figure(2));
        image_num = T.x20x_num(i);
        meta_num = image_num+1;
        O_Bitmap = O_wid(image_num); % Get image
        % Display image
        figure(2); h1 = O_Bitmap.plot('-scalebar'); % Show positions AND show default bottom-left auto-length scalebar with text
        figure(2); h1a = O_Bitmap.plot_position(O_wid([point_num])); % Moved figure(4) outside the function call
        set(findall(h1a, 'Type', 'line'), 'LineWidth', 3); % Set line color to white and its width to 3
        % Get scaling data 
        O_Text = O_wid(meta_num); % import metadata file
        um_width=str2double(O_Text.Data{13,2});
        um_height=str2double(O_Text.Data{14,2});
        px_width = size(O_Bitmap.Data,1);
        px_height = size(O_Bitmap.Data,2);
        Xscale = um_width/px_width;
        Yscale = um_height/px_height;
    
        figure(2); % make full screen
        %set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        drawnow;
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

        % ---------------------------------------------------------- %
    
        p0 = [h1a.XData, h1a.YData];
        num_points = 8;
        r_values = zeros(1, num_points);
    
        for j = 1:num_points
            pi = ginput(1);
            lineSec = [p0; pi];
            h = line(lineSec(:,1), lineSec(:,2), 'linewidth', 2, 'DisplayName', ['Segment ' num2str(j)]);
    
            delta_x = (h.XData(end) - h.XData(1)) * Xscale;
            delta_y = (h.YData(end) - h.YData(1)) * Yscale;
            r = sqrt(delta_x^2 + delta_y^2);
            r_values(j) = r;
    
            mid_x = mean(h.XData);
            mid_y = mean(h.YData);
            text(mid_x, mid_y, sprintf('%.2f', r), ...
                'FontSize', 16, 'Color', 'white', 'HorizontalAlignment', 'center');
        end
    
    elseif choice == 3
        % ----------------------------- Option C ----------------------------- %
        num_points = 4;
        d_values = zeros(1, num_points);
    
        for j = 1:num_points
            p0 = ginput(1);
            pi = ginput(1);
            lineSec = [p0; pi];
            h = line(lineSec(:,1), lineSec(:,2), 'linewidth', 2, 'DisplayName', ['Segment ' num2str(j)]);
    
            delta_x = (h.XData(end) - h.XData(1)) * Xscale;
            delta_y = (h.YData(end) - h.YData(1)) * Yscale;
            d = sqrt(delta_x^2 + delta_y^2);
            d_values(j) = d;
    
            mid_x = mean(h.XData);
            mid_y = mean(h.YData);
            text(mid_x, mid_y, sprintf('%.2f', d), ...
                'FontSize', 16, 'Color', 'white', 'HorizontalAlignment', 'left');
        end
    
        r_values = d_values ./ 2;

    elseif choice == 4
        % ----------------------------- stop ----------------------------- %
        %warning('Exiting loop and terminating script.');
        return;
        % ----------------------------- skip ----------------------------- %
    else
        warning('No selection made. Skipping this image.');
        T.r_um(i) = NaN;
        continue;
    end
    
    r_mean = mean(r_values);

    %disp(['Mean distance: ', num2str(r_mean)])
    T.r_um(i) = r_mean;
    %% -------------------------------------------------------------------------%
    pause(3);  % Wait for 3 seconds

    % Ask user for notes on the current image
    notePrompt = {['Enter notes for image ', T.filename{i}]};
    dlgTitle = 'Image Notes';
    numLines = [4 50];  % 4 lines tall, 50 characters wide
    defaultAns = {''};
    userNote = inputdlg('Enter notes for this image:', 'Notes', [4 50]);
    
    if ~isempty(userNote)
        T.notes_FI{i} = userNote{1};  % assign user input
    else
        T.notes_FI{i} = '';           % empty if canceled
end


    % Auto-save the table to Excel after each entry
    writetable(T, savePath);

    close(figure(2));
    close all;
end
