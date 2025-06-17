function     [data, trialTS, AnimalIDcell] = CocaineSA_MedPC2mat(MedPCfile, folderDir, dur)

% MedPCfile and folderDir are cell arrays

%%% MedPC mapper

% \  E = Active nose poke Timestamp
% \  F = Pump on Timestamp
% \  G = pump off Timestamp
% \  H = Inactive NosePoke Response Timestamp
% \  I = Cocaine triggering nose poke Timestamp
% \  K = array for no-cocaine active nosepoke timestamps

%%% for binary data extraction
% PortEntry = 1;
% PortExit = 2;
PumpOn = 1;
PumpOff = 2;
% Initialize the data_map with direct mappings
data_map = containers.Map({'F:', 'G:','E:', 'H:', 'I:', 'K:'}, ...
    [PumpOn, PumpOff, 3, 4, 5, 6]);
% Initialize data_map2 with mappings for the structure fields
data_map2 = containers.Map([PumpOn, 3, 4, 5, 6], ... 
    {'PumpOn', 'ActiveNP',  'InactiveNP', 'CocTrigNP', 'NoCocActiveNP'});

%%% for timestamps extraction
% Initialize the timestamp_map with direct mappings
timestamp_map = containers.Map({'F:', 'G:','E:', 'H:', 'I:', 'K:'}, ...
    [1,2,3,4,5,6]);
% Initialize timestamp_map2 with mappings for the structure fields
timestamp_map2 = containers.Map([1,2,3,4,5,6], ... %
    {'PumpOn', 'PumpOff', 'ActiveNP',  'InactiveNP', 'CocTrigNP', 'NoCocActiveNP'});


if iscell(folderDir)
    folderDir = folderDir{1};
end

MedPCfile = string(MedPCfile);

if nargin < 3
    % extract the time duration of the session
    for fileNumber = 1:length(MedPCfile)
        cd(folderDir);
        fileName = string(MedPCfile(fileNumber, 1));
        fid = fopen(fileName); %open a medpc file
        while ~feof(fid) %read medpc file
            line = fgets(fid);
            if startsWith(strtrim(line), 'B:')
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
    fileName = string(MedPCfile(fileNumber, 1));
    fid = fopen(fileName); %open a medpc file
    aline = fgetl(fid);%read line excluding newline character

    data_matrix = zeros(1,size(data_map,1));

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

    if isempty(data_matrix)
        fprintf('No data found in %s\n', fileName);
        continue;
    end

    keys = data_map2.keys();
    for k = 1:length(keys)
        col_idx = keys{k};
        length_temp = find(data_matrix(:, col_idx), 1, 'last');

        output_matrix_name = data_map2(col_idx);
        data.(output_matrix_name)(fileNumber,:) = zeros(1,dur);
        if col_idx == PumpOn
            for numentry = 1:length_temp
                startIdx = round(data_matrix(numentry, col_idx)*100);
                endIdx = round(data_matrix(numentry, PumpOff)*100);
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
    fileName = string(MedPCfile(fileNumber, 1));
    fid = fopen(fileName); %open a medpc file
    aline = fgetl(fid);%read line excluding newline character

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
    fileName = string(MedPCfile(fileNumber, 1));
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
    %         idx = diff(trialTS.ActiveNP{fileNumber}) >= 6000;
    %         idx = [1;idx];
    %         trialTS.FirstActive{fileNumber} = trialTS.ActiveNP{fileNumber}(logical(idx));
    %     end
end
temp = splitlines(AnimalIDcell);
AnimalIDcell = temp(:,:,1);

