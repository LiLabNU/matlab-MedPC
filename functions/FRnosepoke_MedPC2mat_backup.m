function     [data, trialTS, AnimalIDcell] = FRnosepoke_MedPC2mat(MedPCfile, folderDir, dur)

%%% MedPC mapper
% \  D = Completed FR Timestamp
% \  E = Left NosePoke Response Timestamp
% \  F = Tone On Timestamp
% \  G = Shock Frenq
% \  H = Right NosePoke Response Timestamp
% \  I = Sucroce Timestamp
% \  J = Shock Intensity
% \  K = Shock Timestamp
% \  L =
% \  M =
% \  R =
% \  P = Port entry time stamp array
% \  N = Port exit time stamp array


%%% for data extraction
PortEntry = 1;
PortExit = 2;
% Initialize the data_map with direct mappings
data_map = containers.Map({'P:', 'N:', 'E:', 'H:', 'I:'}, ...
    [PortEntry, PortExit, 3, 4, 5]);
% Initialize data_map2 with mappings for the structure fields
data_map2 = containers.Map([PortEntry, 3, 4, 5], ... %
    {'PortEntry', 'ActiveNP', 'InactiveNP', 'Sucrose'});

%%% for timestamps extraction
% Initialize the timestamp_map with direct mappings
timestamp_map = containers.Map({'D:', 'F:'}, ...
    [1,2]);
% Initialize timestamp_map2 with mappings for the structure fields
timestamp_map2 = containers.Map([1, 2], ... %
    {'CompleteNP', 'CueTS'});



if nargin < 3
    % extract the time duration of the session
    for fileNumber = 1:length(MedPCfile)
        cd(folderDir{fileNumber});
        fileName = string(MedPCfile(fileNumber, 1));
        fid = fopen(fileName); %open a medpc file
        while ~feof(fid) %read medpc file
            line = fgets(fid);
            if startsWith(strtrim(line), 'T:')
                parts = strsplit(line, ':'); % Split the line at ':'
                dur(fileNumber) = str2double(strtrim(parts{end})); % Convert the value to double
                break; % Exit the loop once the value is found
            end
        end
        fclose(fid);
    end
    dur = round(max(dur))*100;
end


% Main loop
for fileNumber = 1:length(MedPCfile)
    cd(folderDir{fileNumber});
    fileName = string(MedPCfile(fileNumber, 1));
    fid = fopen(fileName); %open a medpc file
    aline = fgetl(fid);%read line excluding newline character

    data_matrix = [];
    TS_matrix = zeros(1000,size(timestamp_map,1));

    while ~feof(fid)
        if length(aline) == 0
            aline = fgetl(fid);
            continue;
        end

        prefix = aline(1:2);

        if isKey(data_map, prefix) && length(aline) == 2
            col_idx = data_map(prefix);
            aline = fgetl(fid);
            n = 0;

            while length(aline) > 2
                n = n + 1;
                tempdata = regexp(aline, ' ', 'split');
                data_matrix(n, col_idx) = str2num(cell2mat(tempdata(end)));
                aline = fgetl(fid);
            end

        elseif isKey(timestamp_map, prefix) && length(aline) == 2
            col_idx2 = timestamp_map(prefix);
            aline = fgetl(fid);
            n = 0;

            while length(aline) > 2
                n = n + 1;
                tempdata = regexp(aline, ' ', 'split');
                TS_matrix(n, col_idx2) = str2num(cell2mat(tempdata(end)));
                aline = fgetl(fid);
            end

        else
            aline = fgetl(fid);
        end
    end
    fclose(fid);

    % extracting data
    keys = data_map2.keys();
    for k = 1:length(keys)
        col_idx = keys{k};
        length_temp = find(data_matrix(:, col_idx), 1, 'last');

        %         if col_idx == PortEntry && length_temp > 0 && data_matrix(length_temp, PortExit) == 0
        %             data_matrix(length_temp, PortExit) = dur/100;
        %         end



        output_matrix_name = data_map2(col_idx);
        data.(output_matrix_name)(fileNumber,:) = zeros(1,dur);
        % Process entries and exits or single timestamps
        if col_idx == PortEntry
            for numentry = 1:length_temp
                startIdx = round(data_matrix(numentry, col_idx)*100);
                endIdx = round(data_matrix(numentry, PortExit)*100);
                % Ensure the index does not exceed 'dur'
                endIdx = min(endIdx, dur);
                data.(output_matrix_name)(fileNumber, startIdx:endIdx) = 1;
            end
        else
            for numentry = 1:length_temp
                idx = round(data_matrix(numentry, col_idx)*100);
                % Ensure the index does not exceed 'dur'
                if idx <= dur
                    data.(output_matrix_name)(fileNumber, idx) = 1;
                end
            end
        end

        %             % Padding with zeros to ensure length matches 'dur'
        %             currentLength = size(data.(output_matrix_name), 2);
        %             if currentLength < dur
        %                 % Calculate the number of zeros needed to pad
        %                 paddingLength = dur - currentLength;
        %                 % Pad the specific field for the current fileNumber with zeros
        %                 data.(output_matrix_name)(fileNumber, end+1:dur) = zeros(1, paddingLength);
        %             end

    end

    %extracting timestamps
    for key = timestamp_map.keys()
        col_idx2 = timestamp_map(key{:}); % Get column index for the current timestamp type
        ts_type = timestamp_map2(col_idx2); % Get the corresponding field name for this timestamp

        % Check and extract timestamps
        ts_data = TS_matrix(TS_matrix(:, col_idx2) > 0, col_idx2); % Extract non-zero timestamps
        ts_data_rounded = round(ts_data .* 100); % Round timestamps to nearest hundredth

        trialTS.(ts_type){fileNumber} = ts_data_rounded; % Assign rounded timestamps to the correct field
    end

    %extracting names
    fid = fopen(fileName);

    for i = 1:20
        if ~feof(fid)
            line = fgets(fid); % Read the current line
            if contains(line, 'Subject:') % Check if the line contains the string 'Name:'
                % Use regular expression to extract the name
                SubjectID (fileNumber,1) = regexp(line, 'Subject: (.*)', 'tokens');
            elseif contains(line, 'Experiment:') % Check if the line contains the string 'Name:'
                % Use regular expression to extract the name
                ExperimentID (fileNumber,1) = regexp(line, 'Experiment: (.*)', 'tokens');

            elseif contains(line, 'Group:') % Check if the line contains the string 'Name:'
                % Use regular expression to extract the name
                GroupID (fileNumber,1) = regexp(line, 'Group: (.*)', 'tokens');

            elseif contains(line, 'Box:') % Check if the line contains the string 'Name:'
                % Use regular expression to extract the name
                BoxID (fileNumber,1) = regexp(line, 'Box: (.*)', 'tokens');

            end
        end
    end

    AnimalIDcell(fileNumber, 1) = string(SubjectID(fileNumber,1));
    AnimalIDcell(fileNumber, 2) = ExperimentID(fileNumber,1);
    AnimalIDcell(fileNumber, 3) = GroupID(fileNumber,1);
    AnimalIDcell(fileNumber, 4) = BoxID(fileNumber,1);


end
temp = splitlines(AnimalIDcell);
AnimalIDcell = temp(:,:,1);

fclose(fid);
