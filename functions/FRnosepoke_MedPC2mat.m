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
% \  O = First active nosepoke time stamp array


%%% for data extraction
PortEntry = 1;
PortExit = 2;

%%% for timeseries extraction
% Initialize the data_map with direct mappings
data_map = containers.Map({'P:', 'N:', 'E:', 'H:', 'I:'}, ...
    [PortEntry, PortExit, 3, 4, 5]);
% Initialize data_map2 with mappings for the structure fields
data_map2 = containers.Map([PortEntry, 3, 4, 5], ... %
    {'PortEntry', 'ActiveNP', 'InactiveNP', 'Sucrose'});

%%% for timestamps extraction
% Initialize the timestamp_map with direct mappings
timestamp_map = containers.Map({'I:', 'F:', 'E:', 'H:', 'K:','P:','N:','O:'}, ...
    [1,2,3,4,5,6,7,8]);
% Initialize timestamp_map2 with mappings for the structure fields
timestamp_map2 = containers.Map([1,2,3,4,5,6,7,8], ... %
    {'Sucrose', 'Cue', 'ActiveNP', 'InactiveNP', 'Shock', 'PortEntry', 'PortExit', 'FirstActive'});

 if iscell(folderDir)
        folderDir = folderDir{1};
 end

MedPCfile = string(MedPCfile);

if nargin < 3
    % extract the time duration of the session
    for fileNumber = 1:length(MedPCfile)
        cd(folderDir);
        fileName = string(MedPCfile(fileNumber));
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

% extracting binary data
for fileNumber = 1:length(MedPCfile)
    cd(folderDir);
    fileName = string(MedPCfile(fileNumber));
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

        else
            aline = fgetl(fid);
        end
    end
    fclose(fid);

    keys = data_map2.keys();
    for k = 1:length(keys)
        col_idx = keys{k};
        length_temp = find(data_matrix(:, col_idx), 1, 'last');
        
        output_matrix_name = data_map2(col_idx);
        data.(output_matrix_name)(fileNumber,:) = zeros(1,dur);
        if col_idx == PortEntry
            for numentry = 1:length_temp
                startIdx = round(data_matrix(numentry, col_idx)*100);
                endIdx = round(data_matrix(numentry, PortExit)*100);
                endIdx = min(endIdx, dur);
                data.(output_matrix_name)(fileNumber, startIdx:endIdx) = 1;
            end
        else
            for numentry = 1:length_temp
                idx = round(data_matrix(numentry, col_idx)*100);
                if idx <= dur
                    data.(output_matrix_name)(fileNumber, idx) = 1;
                end
            end
        end
    end
end

% extracting timestamps
for fileNumber = 1:length(MedPCfile)
    cd(folderDir);
    fileName = string(MedPCfile(fileNumber));
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

        if isKey(timestamp_map, prefix) && length(aline) == 2
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

    for key = timestamp_map.keys()
        col_idx2 = timestamp_map(key{:});
        ts_type = timestamp_map2(col_idx2);
        
        ts_data = TS_matrix(TS_matrix(:, col_idx2) > 0, col_idx2);
        ts_data_rounded = round(ts_data .* 100);
        
        trialTS.(ts_type){fileNumber} = ts_data_rounded;
    end
end

% extracting names
for fileNumber = 1:length(MedPCfile)
    cd(folderDir);
    fileName = string(MedPCfile(fileNumber));
    fid = fopen(fileName);
    SubjectID = "";
    ExperimentID = "";
    GroupID = "";
    BoxID = "";
    for i = 1:20  % Assuming the relevant information is within the first 20 lines
        if ~feof(fid)
            line = fgets(fid);
            if contains(line, 'Subject:')
                SubjectID(fileNumber,1) = regexp(line, 'Subject: (.*)', 'tokens', 'once');
            elseif contains(line, 'Experiment:')
                ExperimentID(fileNumber,1) = regexp(line, 'Experiment: (.*)', 'tokens', 'once');
            elseif contains(line, 'Group:')
                GroupID(fileNumber,1) = regexp(line, 'Group: (.*)', 'tokens', 'once');
            elseif contains(line, 'Box:')
                BoxID(fileNumber,1) = regexp(line, 'Box: (.*)', 'tokens', 'once');
            end
        end
    end
    fclose(fid);
    
    AnimalIDcell(fileNumber, 1) = string(SubjectID(fileNumber,1));
    AnimalIDcell(fileNumber, 2) = ExperimentID(fileNumber,1);
    AnimalIDcell(fileNumber, 3) = GroupID(fileNumber,1);
    AnimalIDcell(fileNumber, 4) = BoxID(fileNumber,1);



%     if isempty(trialTS.FirstActive{fileNumber})
% 
%         activeTS = trialTS.ActiveNP{fileNumber};     % Vector of timestamps (sorted)
%         sucroseTS = trialTS.Sucrose{fileNumber};     % Vector of timestamps (sorted)       
%        
%         [firstNP_success, firstNP_timeout]=find_first_nosepokes(activeTS,sucroseTS);
%         trialTS.FirstActive{fileNumber} =sort([firstNP_success, firstNP_timeout])';
%     end
end
temp = splitlines(AnimalIDcell);
AnimalIDcell = temp(:,:,1);

