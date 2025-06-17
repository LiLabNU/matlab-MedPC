function     [allport_P, allport_Q, AnimalIDcell,trialTypes, trialTimeStamps] = PortEntry_MedPC2mat(MedPCfile, folderDir)

column_number1 = 1;
column_number2 = 2;
column_number3 = 3;
column_number4 = 4;
column_number5 = 5;
column_number6 = 6;

% Specify total duration of the session
total_dur = 11045;
total_durus=total_dur*100;
allport_P=zeros(length(MedPCfile),total_durus);
allport_Q=zeros(length(MedPCfile),total_durus);


if iscell(folderDir)
    folderDir = folderDir{1};
end

cd(folderDir)
% Main loop originated by Shan, modified by Hao
for fileNumber = 1:length(MedPCfile)
    fileName = string(MedPCfile(fileNumber, 1));
    fid = fopen(fileName); %open a medpc file
    aline = fgetl(fid);%read line excluding newline character

    data_write = {};
    data_matrix = [];

    while ~feof(fid)%read medpc file
        % to get port entry
        if length(aline) == 0
            aline = fgetl(fid);
        elseif aline(1:2) == 'P:' & length(aline)==2
            n = 1;
            data_write{n,column_number1} = aline(1);
            aline = fgetl(fid); %find P
            while length(aline)>2
                n = n+1;
                tempdata = regexp(aline,' ','split');
                data_write(n,column_number1) = tempdata(end);
                data_matrix(n-1,column_number1) = str2num(cell2mat(tempdata(end)));
                aline = fgetl(fid);
            end

        elseif aline(1:2) == 'N:' & length(aline)==2
            n = 1;
            data_write{n,column_number2} = aline(1);
            aline = fgetl(fid);
            while length(aline)>2
                n = n+1;
                tempdata = regexp(aline,' ','split');
                data_write(n,column_number2) = tempdata(end);
                data_matrix(n-1,column_number2) = str2num(cell2mat(tempdata(end)));
                aline = fgetl(fid);
            end

        elseif aline(1:2) == 'Q:' & length(aline)==2
            n = 1;
            data_write{n,column_number3} = aline(1);
            aline = fgetl(fid);
            while length(aline)>2
                n = n+1;
                tempdata = regexp(aline,' ','split');
                data_write(n,column_number3) = tempdata(end);
                data_matrix(n-1,column_number3) = str2num(cell2mat(tempdata(end)));
                aline = fgetl(fid);
            end

        elseif aline(1:2) == 'R:' & length(aline)==2
            n = 1;
            data_write{n,column_number4} = aline(1);
            aline = fgetl(fid);
            while length(aline)>2
                n = n+1;
                tempdata = regexp(aline,' ','split');
                data_write(n,column_number4) = tempdata(end);
                data_matrix(n-1,column_number4) = str2num(cell2mat(tempdata(end)));
                aline = fgetl(fid);
            end

        elseif aline(1:2) == 'S:' & length(aline)==2
            n = 1;
            data_write{n,column_number5} = aline(1);
            aline = fgetl(fid);
            while length(aline)>2
                n = n+1;
                tempdata = regexp(aline,' ','split');
                data_write(n,column_number5) = tempdata(end);
                data_matrix(n-1,column_number5) = str2num(cell2mat(tempdata(end)));
                aline = fgetl(fid);
            end

        elseif aline(1:2) == 'K:' & length(aline)==2
            n = 1;
            data_write{n,column_number6} = aline(1);
            aline = fgetl(fid);
            while length(aline)>2
                n = n+1;
                tempdata = regexp(aline,' ','split');
                data_write(n,column_number6) = tempdata(end);
                data_matrix(n-1,column_number6) = str2num(cell2mat(tempdata(end)));
                aline = fgetl(fid);
            end

        else
            aline = fgetl(fid);
        end

    end
    fclose(fid);

    %     for P=1:2:size(data_matrix,2)%balance
    %     N=P+1;
    length_P = find(data_matrix(:,column_number1)~=0,1,'last');
    if data_matrix(length_P,column_number2) == 0
        data_matrix(length_P,column_number2) = total_dur;
    end
    %     end

    %     for P=1:2:size(data_matrix,2)%generate file for allport data
    %     N=P+1;
    for numentry=1:find(data_matrix(:,column_number1)~=0,1,'last')
        allport_P(fileNumber,round(data_matrix(numentry,column_number1).*100):round(data_matrix(numentry,column_number2).*100))=1;
    end

    length_Q = find(data_matrix(:,column_number3)~=0,1,'last');
    if data_matrix(length_Q,column_number4) == 0
        data_matrix(length_Q,column_number4) = total_dur;
    end
    %     end

    %     for P=1:2:size(data_matrix,2)%generate file for allport data
    %     N=P+1;
    for numentry=1:find(data_matrix(:,column_number3)~=0,1,'last')
        allport_Q(fileNumber,round(data_matrix(numentry,column_number3).*100):round(data_matrix(numentry,column_number4).*100))=1;
    end

    trialTS = data_matrix(:,column_number5);
    trialtype = data_matrix(:,column_number6);

    if any(trialTS == 0)
        Total = find(trialTS==0,1,'first')-1;
        trialTS = trialTS(1:Total,:);
        trialtype = trialtype(1:Total,:);
    end

    % fid = fopen(fileName);
    % 
    % AnimalInformation = textscan(fid,'%s',6,'HeaderLines',6);
    % %Subject ID
    % SubjectID (fileNumber,1) = string(AnimalInformation{1}(2,1));
    % %Experiment (day1, day2...)
    % ExperimentID (fileNumber,1) = string(AnimalInformation{1}(4,1));
    % % Group (control bs experimental)
    % GroupID (fileNumber,1) = string(AnimalInformation{1}(6,1));
    % 
    % AnimalIDcell(fileNumber, 1) = SubjectID(fileNumber,1);
    % AnimalIDcell(fileNumber, 2) = ExperimentID(fileNumber,1);
    % AnimalIDcell(fileNumber, 3) = GroupID(fileNumber,1);
    % AnimalIDcell(fileNumber, 4) = Box(fileNumber,1);
    % %     totalPortEntry_P(fileNumber,1)= length_P;
    % %     totalPortEntry_Q(fileNumber,1)= length_Q;

    trialTypes{fileNumber} = trialtype;
    trialTimeStamps{fileNumber} = trialTS;

end


% save ('allport.mat','allport');
% save ('AnimalIDcell.mat','AnimalIDcell');
% save ('totalPortEntry.mat','totalPortEntry');
% save ('Trialtype.mat','Trialtype');


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



