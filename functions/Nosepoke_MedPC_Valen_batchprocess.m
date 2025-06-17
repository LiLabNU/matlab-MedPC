
function [data, trialTS, Totals, Details, Avg, miceIDs, folderName, folderDir] = Nosepoke_MedPC_Valen_batchprocess(keywords, motherDir)

%%
[files, folderDir] = GetFilesFromFolder(0,[],motherDir); %select the mother folder

% put keywords all the subfolders contain
files = string(files);
idx = startsWith(files, keywords, 'IgnoreCase', true) & ~startsWith(files, '_');
files = files(idx);
folderDir = string(folderDir);
folderDir = folderDir(idx);

%% load MedPC files from selected subfolders
for n = 1:length(files)
    folderPath = fullfile(folderDir(n), files(n));
    if isfolder(folderPath)
        [MedPCfile, folder] = GetFilesFromFolder(0,'.txt', folderPath);
        [data{n}, trialTS{n}] = FRnosepoke_MedPC2mat(MedPCfile, folder);
        folderName(n,:) = string(cellfun(@(x) regexp(x, '[^\\]+$', 'match', 'once'), folder{1}, 'UniformOutput', false));
        % get animal IDs
        mice{n} = cellfun(@(x) regexp(x, '(?<= ).*(?=\.txt)', 'match', 'once'), MedPCfile, 'UniformOutput', false);
    else
        error('The specified folderPath (%s) is not a valid folder.', folderPath);
    end
end
[~,I] = max(cell2mat(cellfun(@(x) size(x,1), mice, 'UniformOutput', false)));
refMice = mice{I};

    
%% reorganize orders based on mouse ID
for n = 1:length(files)  
    currMice = mice{n};
    [~, idx] = ismember(refMice, currMice);  
    idxUse = idx(idx > 0);
    newAnimalID = strings(size(refMice));    
    % Assign values based on idx, keeping unmatched entries empty
    newAnimalID(idx > 0) = mice{n}(idxUse);
    mice{n} = newAnimalID;

    fields = fieldnames(data{n});
    for fd = 1:length(fields)
        field = fields{fd};
        newDataField = nan(size(data{I}.(field)));
        id2 = 1;
        for id = 1:length(idx)            
            if idx(id)>0
                newDataField(id,:) = data{n}.(field)(idxUse(id2),:);
                id2 = id2+1;            
            end
        end
        data{n}.(field) = newDataField;
    end

    clear newDataField
    fields = fieldnames(trialTS{n});
    for fd = 1:length(fields)
        field = fields{fd};
        newDataField = cell(size(trialTS{I}.(field)));
        id2 = 1;
        for id = 1:length(idx)
            if idx(id)>0
                newDataField(id) = trialTS{n}.(field)(:,idxUse(id2));
                id2 = id2+1;            
            end
        end
        trialTS{n}.(field)= newDataField;
    end
end

%% calculate
minGap = 20; % 0.2s 

for n = 1:length(files)  
    % calculate totals
    PortEntry = data{n}.PortEntry./100;
    Sucrose = data{n}.Sucrose;
    ActiveNP = data{n}.ActiveNP;
    InactiveNP = data{n}.InactiveNP;
    for i=1:size(PortEntry,1) % for i equal to 1 to the # of medPC files (i.e., animals)
        Totals.Sucrose_n(i,n)=sum(Sucrose(i,:));%set sucrose_n, in the i-th row and first column, to the sum of events in the i-th row of Sucrose
        Totals.ActiveNP_n(i,n)=sum(ActiveNP(i,:));
        Totals.InactiveNP_n(i,n)=sum(InactiveNP(i,:));
        Totals.PortEntry_n(i,n)=sum(PortEntry(i,:));
    end   

    if contains(folderName(n),'FR1_', 'IgnoreCase', true)
        fr = 1;
    elseif contains(folderName(n),'FR3_', 'IgnoreCase', true)
        fr = 3;
    end

    for i=1:size(PortEntry,1)
        SucTS = trialTS{n}.Sucrose{i};   
        ShkTS = trialTS{n}.Shock{i};  
        if ~isempty(SucTS)
            if sum(data{n}.PortEntry(i,:)) > 10 && ~isempty(trialTS{n}.FirstActive{i}) && ~isempty(trialTS{n}.ActiveNP{i})
                [Details.PEdetails{i,n}, SumPE(i,n)] = findBouts(data{n}.PortEntry(i,:), minGap, SucTS, mice{n}(i), 'mean', ShkTS);
                [Details.FRdetails{i,n}, SumFR(i,n)] = findFRduration(trialTS{n}.FirstActive{i}, SucTS, trialTS{n}.ActiveNP{i}, fr);
            end
        end              
        miceIDs(i,n) = mice{n}(i);
    end
end

%% reorganize 
for n = 1:length(files)
    temp = SumPE(:,n);
    fields = fieldnames(temp);
    for fd = 1:length(fields)
        field = fields{fd};
        for j = 1:length(temp)
            if ~isempty(temp(j).(field))
                Avg.SumBouts.(field)(j,n) = temp(j).(field);
            else
                Avg.SumBouts.(field)(j,n) = 0;
            end
        end
    end

    temp = SumFR(:,n);
    fields = fieldnames(temp);
    for fd = 1:length(fields)
        field = fields{fd};
        for j = 1:length(temp)
            if ~isempty(temp(j).(field))
                Avg.SumFR.(field)(j,n) = temp(j).(field);
            else
                Avg.SumFR.(field)(j,n) = 0;
            end
        end
    end
end
end
