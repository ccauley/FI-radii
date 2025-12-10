
%-------------------------------------------------------------------------%
% This matlab code requires the plugin described in Holmi, J.T., and
% Lipsanen, H., 2022, WITio: A MATLAB data evaluation toolbox to script
% broader insights into big data from WITec microscopes: SoftwareX, v. 18, p. 101009, 
% doi:10.1016/j.softx.2022.101009.
% 
% code by Christina Cauley June 06 2025
% code updated by Christina Cauley Dec 01 2025
%-------------------------------------------------------------------------%
% Temporarily set user preferences
resetOnCleanup = WITio.tbx.pref.set({'wip_AutoCreateObj', 'wip_AutoCopyObj', 'wip_AutoModifyObj'}, {true, true, true}); % The original values prior to this call are restored when resetOnCleanup-variable is cleared.
%set(0,'DefaultFigureWindowStyle','normal');

WITio.tbx.edit(); % Open this code in Editor
close all; % Close figures

% Define WIP file to open
pathstr = '/Users/christinacauley/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab_Work/Olivine Deformation/Raman/WITio_FI_radii/example'; 
file = fullfile(pathstr, 'Input_from_Raman.wip'); 

% Open density results file
pathstr2 = '/Users/christinacauley/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab_Work/Olivine Deformation/Raman/WITio_FI_radii/example'; 
file_densities = fullfile(pathstr2, 'Input_from_Diadfit_fitted_2025-11-19.xlsx'); 
T = readtable(file_densities);
origNames = T.Properties.VariableDescriptions;
%T = readtable(file_densities, 'VariableNamingRule', 'preserve');


newfile_name = 'Output_FI_radii_Diadfit_fitted_2025-11-19.xlsx'
savePath = fullfile(pathstr2, newfile_name);
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
% O_wid.manager;  % Open Project Manager - use this command to view interactive Project Manager viewer similar to WITio proprietary 

%% ------------------------------ GET WIN INDEXs -------------------------------------------%
num_obj=length(O_wid);
wid_IndTag = strings(1,num_obj);

% loop though each object in o_Wid to save filenames
for i=1:num_obj
    wid_IndTag(i)=string(O_wid(i).Name);
end

% loop through density table to pull Win index
suffixes = {'', '_05x', '_20x', '_50x', '_100x'};
colNames = {'point_num', 'x05x_num', 'x20x_num', 'x50x_num', 'x100x_num'};

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
    try
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
    catch ContinueFile
        file_densities_new = fullfile(pathstr2, newfile_name);
        T = readtable(file_densities_new);
        % --- Fix column names if needed ---
        % Case 1: auto-sanitized "r (µm)" → "r__m_"
        if ismember('r__m_', T.Properties.VariableNames)
            T.Properties.VariableNames{'r__m_'} = 'r_um';
            disp('Renamed r__m_ → r_um');
        end
        % Case 2: actual original name
        if ismember('r (µm)', T.Properties.VariableNames)
            T.Properties.VariableNames{'r (µm)'} = 'r_um';
            disp('Renamed r (µm) → r_um');
        end
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
end

for i = startIdx:height(T)

    if ~isnan(T.x50x_num(i))
        image_num = T.x50x_num(i);
    elseif ~isnan(T.x100x_num(i))
        image_num = T.x100x_num(i);
    else
        image_num = T.x20x_num(i);
    end
    
    meta_num  = image_num + 1;
    point_num = T.point_num(i);
    O_Bitmap  = O_wid(image_num);

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
    % -------------------------------------------------------------------------
    % Dynamic menu: hide unavailable options
    menuLabels = {};
    menuMap = [];   % maps displayed index → original option number
    
    % Option A
    if ~isnan(T.x50x_num(i))
        menuLabels{end+1} = 'Option A: From analytical point outward (50x)';
        menuMap(end+1) = 1;
    end
    
    % Option B
    if ~isnan(T.x20x_num(i))
        menuLabels{end+1} = 'Option B: From analytical point outward (20x)';
        menuMap(end+1) = 2;
    end
    
    % Option C
    if ~isnan(T.x100x_num(i))
        menuLabels{end+1} = 'Option C: From analytical point outward (100x)';
        menuMap(end+1) = 3;
    end
    
    % Option D (always shown)
    menuLabels{end+1} = 'Option D: Arbitrary edge-to-edge';
    menuMap(end+1) = 4;
    
    % Stop loop (always shown)
    menuLabels{end+1} = 'Stop loop';
    menuMap(end+1) = 5;
    
    % Skip (always shown)
    menuLabels{end+1} = 'Skip';
    menuMap(end+1) = 6;
    
    % Show dynamic menu
    idx = menu(['Choose radius extraction method for image ', T.filename{i}], ...
               menuLabels);
    
    % Translate menu index → original choice number
    choice = menuMap(idx);
    % -------------------------------------------------------------------------

    if choice == 1
        % ----------------------------- Option A:50x ----------------------------- %
        p0 = [h1a.XData, h1a.YData];
        num_points = 8;
        r_values = zeros(1, num_points);
    
        for j = 1:num_points
            pi = ginput(1);
            lineSec = [p0; pi];
            h = line(lineSec(:,1), lineSec(:,2), 'linewidth', 2, 'DisplayName', ['Segment ' num2str(j)]);
            %h = line(lineSec(:,1), lineSec(:,2), 'linewidth', 1.5);
    
            delta_x = (h.XData(end) - h.XData(1)) * Xscale;
            delta_y = (h.YData(end) - h.YData(1)) * Yscale;
            r = sqrt(delta_x^2 + delta_y^2);
            r_values(j) = r;
    
            mid_x = mean(h.XData);
            mid_y = mean(h.YData);
            %text(mid_x, mid_y, sprintf('%.2f', r), 'FontSize', 16, 'Color', 'white', 'HorizontalAlignment', 'center');
        end
    
    elseif choice == 2
        % ----------------------------- Option B:20x ----------------------------- %
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
        % ----------------------------- Option C:100x ----------------------------- %
        close(figure(2));
        image_num = T.x100x_num(i);
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
            %text(mid_x, mid_y, sprintf('%.2f', r), 'FontSize', 16, 'Color', 'white', 'HorizontalAlignment', 'center');
        end
    % elseif choice == 4
    %     % ----------------------------- Option D ----------------------------- %
    %     num_points = 4;
    %     d_values = zeros(1, num_points);
    % 
    %     for j = 1:num_points
    %         p0 = ginput(1);
    %         pi = ginput(1);
    %         lineSec = [p0; pi];
    %         h = line(lineSec(:,1), lineSec(:,2), 'linewidth', 2, 'DisplayName', ['Segment ' num2str(j)]);
    % 
    %         delta_x = (h.XData(end) - h.XData(1)) * Xscale;
    %         delta_y = (h.YData(end) - h.YData(1)) * Yscale;
    %         d = sqrt(delta_x^2 + delta_y^2);
    %         d_values(j) = d;
    % 
    %         mid_x = mean(h.XData);
    %         mid_y = mean(h.YData);
    %         %text(mid_x, mid_y, sprintf('%.2f', d), 'FontSize', 16, 'Color', 'white', 'HorizontalAlignment', 'left');
    %     end
    % 
    %     r_values = d_values ./ 2

    elseif choice == 4
        % ----------------------------- Option D: Arbitrary edge-to-edge ----------------------------- %
    
        % ------------------------------------------------------------------
        % Ask user which resolution image to use (20x, 50x, 100x)
        % ------------------------------------------------------------------
        resMenuLabels = {};
        resMenuMap = [];
    
        if ~isnan(T.x20x_num(i))
            resMenuLabels{end+1} = 'Use 20x image';
            resMenuMap(end+1) = T.x20x_num(i);
        end
        if ~isnan(T.x50x_num(i))
            resMenuLabels{end+1} = 'Use 50x image';
            resMenuMap(end+1) = T.x50x_num(i);
        end
        if ~isnan(T.x100x_num(i))
            resMenuLabels{end+1} = 'Use 100x image';
            resMenuMap(end+1) = T.x100x_num(i);
        end
    
        % If none exist (should not happen), fallback
        if isempty(resMenuLabels)
            warning('No resolution images available! Using currently loaded image.');
            chosenImage = image_num;
        else
            resIdx = menu('Choose image resolution for Option D:', resMenuLabels);
            chosenImage = resMenuMap(resIdx);
        end

    % ------------------------------------------------------------------
    % Load selected resolution image
    % ------------------------------------------------------------------
    close(figure(2));
    image_num = chosenImage;
    meta_num  = image_num + 1;

    O_Bitmap = O_wid(image_num);
    figure(2); h1 = O_Bitmap.plot('-scalebar');

    % Scaling from metadata
    O_Text = O_wid(meta_num);
    um_width  = str2double(O_Text.Data{13,2});
    um_height = str2double(O_Text.Data{14,2});
    px_width  = size(O_Bitmap.Data,1);
    px_height = size(O_Bitmap.Data,2);
    Xscale = um_width/px_width;
    Yscale = um_height/px_height;

    figure(2);
    drawnow;
    set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]);

    % ------------------------------------------------------------------
    % Perform edge-to-edge radius estimation
    % ------------------------------------------------------------------
    num_points = 4;
    d_values = zeros(1, num_points);

    for j = 1:num_points
        p0 = ginput(1);
        pi = ginput(1);
        lineSec = [p0; pi];
        h = line(lineSec(:,1), lineSec(:,2), 'linewidth', 2, ...
                 'DisplayName',['Segment ' num2str(j)]);

        delta_x = (h.XData(end) - h.XData(1)) * Xscale;
        delta_y = (h.YData(end) - h.YData(1)) * Yscale;
        d = sqrt(delta_x^2 + delta_y^2);
        d_values(j) = d;
    end

    r_values = d_values ./ 2;


    elseif choice == 5
        % ----------------------------- stop ----------------------------- %
        %warning('Exiting loop and terminating script.');
        close all;
        T.Properties.VariableNames{'r_um'} = 'r (µm)';
        for k = 1:numel(origNames)
            if ~isempty(origNames{k})
                T.Properties.VariableNames{k} = origNames{k};
            end
        end
        
        writetable(T, savePath);
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

for k = 1:numel(origNames)
    if ~isempty(origNames{k})
        T.Properties.VariableNames{k} = origNames{k};
    end
end

% After restoring original names, fix radius column name
if ismember('r_um', T.Properties.VariableNames)
    T.Properties.VariableNames{'r_um'} = 'r (µm)';
elseif ismember('r__m_', T.Properties.VariableNames)
    % Handle MATLAB auto-sanitized name
    T.Properties.VariableNames{'r__m_'} = 'r (µm)';
end

writetable(T, savePath);
