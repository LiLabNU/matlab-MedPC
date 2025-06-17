function [pe] = processPE(cue_onset, win_half, cue_dur, binsize, data, trialtype) 

if nargin < 5
    cue_onsetTS = cue_onset;
else
    cue_onset = cue_onset(cue_onset(:,1)==trialtype,:);
    cue_onsetTS = cue_onset(:,2)*100;
end

totalT = size(cue_onset,1);

Win = [1, 1+win_half/binsize/100];
 
for num = 1:size(data,1)
  
    pe_binary = data(num,:);
    %generate trial*window file
    pe_stru_r = zeros(totalT,win_half*2+1);

    for t = 1:totalT
        if cue_onsetTS(t)+win_half < length(pe_binary)
            pe_stru_r(t,:) = pe_binary(1,floor(cue_onsetTS(t)-win_half):floor(cue_onsetTS(t)+win_half));
        end
    end
    all_pe_stru_r(:,:,num) = pe_stru_r;
    
    %save('Individual Reward CS.mat','all_pe_stru_r');
end
all_pe_stru_r_Day = all_pe_stru_r;

%%%%%%%%%%%%%% generating trial by trial; animal by trial; latency data for
%%%%%%%%%%%%%% each animals
if ~isempty(all_pe_stru_r_Day)
    pe.rawdata.rew=all_pe_stru_r_Day;
    %trial by trial
    pe.txt.rew=mean(all_pe_stru_r_Day,3);
    pe.txt.rew_ms=imresize(pe.txt.rew,[size(pe.txt.rew,1),size(pe.txt.rew,2)/(binsize*100)],'Bilinear');
    
    %animal by trial
    if size(all_pe_stru_r_Day,3) == 1
        pe.mxt.rew = mean(all_pe_stru_r_Day,1);
    else
    pe.mxt.rew=squeeze(mean(all_pe_stru_r_Day,1))';
    end 
    
    pe.mxt.rew_ms=imresize(pe.mxt.rew,[size(pe.mxt.rew,1),size(pe.mxt.rew,2)/(binsize*100)],'Bilinear');
    
    %latency
    pe.latency.rew.txt = zeros(size(all_pe_stru_r_Day,1),size(all_pe_stru_r_Day,3));
    for nm=1:size(all_pe_stru_r_Day,3)
        for n=1:size(all_pe_stru_r_Day,1)
            lat=find(all_pe_stru_r_Day(n,win_half+1:end,nm),1);
            if isempty(lat)
                pe.latency.rew.txt(n,nm)=cue_dur*100;
            else
                pe.latency.rew.txt(n,nm)=lat-1;
            end
        end
    end
    pe.latency.rew.txt = pe.latency.rew.txt/100';
    pe.latency.rew.m = mean(pe.latency.rew.txt);
    pe.DiffScore.mxt.rew_ms = HaoDiffScore(pe.mxt.rew_ms, Win);
    pe.zScore.mxt.rew_ms = HaozScore(pe.mxt.rew_ms, Win);
end
