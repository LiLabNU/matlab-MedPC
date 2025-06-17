function [Summary_ctrl, Summary_exp] = summarizeSessions(miceInclude, ctrl, exp, Totals, Avg, segmentedSum, sessionNames, motherDir)
% summarizeSessions aggregates session-level metrics for control and experimental groups.
%
%   [Summary_ctrl, Summary_exp] = summarizeSessions(miceInclude, ctrl, exp, Totals, Avg, segmentedSum, sessionNames, motherDir)
%
%   Inputs:
%       miceInclude   - A matrix (or cell array) of mouse IDs (strings) per session.
%       ctrl          - A cell array (or vector) of control mouse IDs.
%       exp           - A cell array (or vector) of experimental mouse IDs.
%       Totals        - Structure containing total counts (fields with '_n' in their names).
%       Avg           - Structure with substructures (e.g., SumBouts, SumFR) of averaged data.
%       segmentedSum  - Structure array with segmented data per session.
%       sessionNames  - Array of session names.
%       motherDir     - String representing the directory; used to extract region and sex info.
%
%   Outputs:
%       Summary_ctrl  - Structure array of control session summaries.
%       Summary_exp   - Structure array of experimental session summaries.

%% Extract Region and Sex from motherDir
% Ensure motherDir is a single string
motherDir = motherDir(1);
pathParts = split(motherDir, filesep);
region = pathParts{end-1};  % Second-to-last part
sex = pathParts{end};       % Last part

%% Remove extraneous field from Avg.SumBouts if it exists
if isfield(Avg.SumBouts, 'mouseID')
    Avg.SumBouts = rmfield(Avg.SumBouts, 'mouseID');
end

%% Initialize Output Structures
numSessions = size(segmentedSum, 1);
Summary_ctrl = struct();
Summary_exp  = struct();

%% Loop through Each Session
for session = 1:numSessions
    % Determine indices for control and experimental mice
    ctrlIdx = ismember(miceInclude(:, session), ctrl);
    expIdx  = ismember(miceInclude(:, session), exp);

    % Assign session name to outputs
    Summary_ctrl(session).sessionNames = sessionNames(session);
    Summary_exp(session).sessionNames = sessionNames(session);

    % Process control group: remove empty strings and count mice
    ctrlMice = miceInclude(ctrlIdx, session);
    ctrlMice = ctrlMice(ctrlMice ~= "");
    Summary_ctrl(session).MiceInclude = ctrlMice;
    Summary_ctrl(session).NumMice = numel(ctrlMice);

    % Process experimental group
    expMice = miceInclude(expIdx, session);
    expMice = expMice(expMice ~= "");
    Summary_exp(session).MiceInclude = expMice;
    Summary_exp(session).NumMice = numel(expMice);

    % Assign region and sex information
    Summary_ctrl(session).Region = string(region);
    Summary_ctrl(session).Sex = string(sex);
    Summary_exp(session).Region = string(region);
    Summary_exp(session).Sex = string(sex);

    %% Add Totals Fields (renamed with '_Total')
    totalFields = fieldnames(Totals);
    for fd = 1:length(totalFields)
        field = totalFields{fd};
        % Create new field name by extracting before '_n' and appending '_Total'
        newFieldName = [extractBefore(field, '_n'), '_Total'];
        Summary_ctrl(session).(newFieldName) = Totals.(field)(ctrlIdx, session);
        Summary_exp(session).(newFieldName) = Totals.(field)(expIdx, session);
    end

    %% Add Avg.SumBouts Fields (append '_SumBouts'), Replace zeros with NaN
    sumBoutsFields = fieldnames(Avg.SumBouts);
    for fd = 1:length(sumBoutsFields)
        field = sumBoutsFields{fd};
        newFieldName = [field, '_SumBouts'];

        val_ctrl = Avg.SumBouts.(field)(ctrlIdx, session);
        val_ctrl(val_ctrl == 0) = NaN;
        Summary_ctrl(session).(newFieldName) = val_ctrl;

        val_exp = Avg.SumBouts.(field)(expIdx, session);
        val_exp(val_exp == 0) = NaN;
        Summary_exp(session).(newFieldName) = val_exp;
    end

    %% Add Avg.SumFR Fields (append '_SumFR'), Replace zeros with NaN
    sumFRFields = fieldnames(Avg.SumFR);
    for fd = 1:length(sumFRFields)
        field = sumFRFields{fd};
        newFieldName = [field, '_SumFR'];

        val_ctrl = Avg.SumFR.(field)(ctrlIdx, session);
        val_ctrl(val_ctrl == 0) = NaN;
        Summary_ctrl(session).(newFieldName) = val_ctrl;

        val_exp = Avg.SumFR.(field)(expIdx, session);
        val_exp(val_exp == 0) = NaN;
        Summary_exp(session).(newFieldName) = val_exp;
    end


    %% Add segmentedSum Fields (append '_Seg'), Replace zero columns with NaN
    segmentedFields = fieldnames(segmentedSum);
    for fd = 1:length(segmentedFields)
        field = segmentedFields{fd};
        newFieldName = [field, '_Seg'];

        % Get control data
        data_ctrl = segmentedSum(session).(field)(ctrlIdx, :);
        % Check columns where all values are 0
        zeroCols_ctrl = all(data_ctrl == 0, 1);
        % Replace entire columns with NaN if all values are zero
        data_ctrl(:, zeroCols_ctrl) = NaN;
        Summary_ctrl(session).(newFieldName) = data_ctrl;

        % Get experimental data
        data_exp = segmentedSum(session).(field)(expIdx, :);
        % Check columns where all values are 0
        zeroCols_exp = all(data_exp == 0, 1);
        % Replace entire columns with NaN if all values are zero
        data_exp(:, zeroCols_exp) = NaN;
        Summary_exp(session).(newFieldName) = data_exp;
    end

end
end
