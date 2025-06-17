function [firstNP_success, firstNP_timeout] = find_first_nosepokes(nosepokeTS, collectionTS)
% nosepokeTS: sorted array of nosepoke timestamps (ms or 10ms increments)
% collectionTS: sorted array of times the mouse *collects* the sucrose
% If the mouse never collects, we skip the rest of nosepokes after success.

timeWindow = 6000;  % 60 s in ms (use 6000 if each unit is 10ms)
firstNP_success = [];
firstNP_timeout = [];

i = 1;  
nPokes = length(nosepokeTS);

while i <= nPokes
    
    % Candidate is the first NP of a potential chain.
    candidateNP = nosepokeTS(i);
    
    % Check if we have at least i+2 to see if 3rd poke is within 60s
    if (i+2) <= nPokes
        thirdPokeTS = nosepokeTS(i+2);
        
        if (thirdPokeTS - candidateNP) <= timeWindow
            % ============= SUCCESS! =================
            firstNP_success(end+1) = candidateNP; %#ok<AGROW>
            
            % The mouse has triggered sucrose. We ignore further 
            % nosepokes until the mouse collects it (if at all).
            % 1) Mark time of sucrose delivery = time of 3rd poke 
            sucroseDeliveredTS = thirdPokeTS;
            
            % 2) Find next collection that happens AFTER sucroseDeliveredTS
            cIdx = find(collectionTS > sucroseDeliveredTS, 1, 'first');
            if isempty(cIdx)
                % No collection => skip the rest
                break;
            else
                % The mouse collects the sucrose at this time
                timeOfCollection = collectionTS(cIdx);
                
                % 3) Skip all nosepokes that occur before "collection"
                while i <= nPokes && nosepokeTS(i) <= timeOfCollection
                    i = i + 1;
                end
            end
            continue;  % proceed to next iteration
            
        else
            % ============= TIME-OUT ================
            firstNP_timeout(end+1) = candidateNP; %#ok<AGROW>
            
            % Skip all nosepokes that are within 60s of candidateNP
            k = i;
            while (k <= nPokes) && (nosepokeTS(k) - candidateNP <= timeWindow)
                k = k + 1;
            end
            i = k;
            continue;
        end
    else
        % Not enough nosepokes left to form 3 => time-out
        firstNP_timeout(end+1) = candidateNP;
        break;
    end
end
end
