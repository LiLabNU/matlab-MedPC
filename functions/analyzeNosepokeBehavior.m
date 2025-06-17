function [results, aggregate_results] = analyzeNosepokeBehavior(tone_onset_ts, nosepoke_ts, reward_ts, timestampMultiplier, trial_duration, pretrial_duration, N)
% Analyze mouse behavior in a behavioral experiment
%
% Parameters:
% tone_onset_ts - timestamps of tone onsets
% nosepoke_ts - timestamps of nosepokes
% reward_ts - timestamps of reward deliveries
% timestampMultiplier - factor to adjust the timestamps (default: 100)
% trial_duration - duration of each trial in timestamp units (default: 30*timestampMultiplier)
% pretrial_duration - duration before each trial to consider in timestamp units (default: 30*timestampMultiplier)
% N - interval for aggregating metrics (default: 5 trials)
%
% Returns:
% results - structure containing trial-by-trial analysis
% aggregate_results - structure containing aggregated metrics every N trials

% Default parameter handling
if nargin < 4 || isempty(timestampMultiplier)
    timestampMultiplier = 100;
end
tone_onset_ts = tone_onset_ts/timestampMultiplier;
nosepoke_ts = nosepoke_ts/timestampMultiplier;
reward_ts = reward_ts/timestampMultiplier;

if nargin < 5 || isempty(trial_duration), trial_duration = 30; end
if nargin < 6 || isempty(pretrial_duration), pretrial_duration = 30; end
if nargin < 7 || isempty(N), N = 5; end

n_trials = length(tone_onset_ts);

% Initialize result structures
results = struct('nosepokes_before_tone_count', {}, 'intervals_before_tone', {}, ...
    'nosepokes_during_tone_count', {}, 'intervals_during_tone', {}, ...
    'latencies_to_rewards', {}, 'avg_latency', {}, ...
    'reward_probability_during_cue', {}, 'reward_probability_not_during_cue', {}, ...
    'rewarded_intervals_during_cue', {}, 'rewarded_intervals_not_during_cue', {});

% Process each trial
for trial = 1:n_trials
    trial_start = tone_onset_ts(trial);
    trial_end = trial_start + trial_duration;
    pre_trial_start = max(trial_start - pretrial_duration, 0); % Ensure non-negative start

    % Process nosepokes and rewards
    results(trial) = processTrial(pre_trial_start, trial_start, trial_end, nosepoke_ts, reward_ts);
end

% Aggregate results every N trials, if necessary
if N > 1
    aggregate_results = aggregateMetricsEveryNTrials(results, N);
else
    aggregate_results = []; % No aggregation if N is not greater than 1
end

end

function trialResults = processTrial(pre_trial_start, trial_start, trial_end, nosepoke_ts, reward_ts)
% Process a single trial's nosepoke and reward data
trialResults = struct();

% Segment nosepokes and rewards by timing
pre_trial_nosepokes = nosepoke_ts(nosepoke_ts >= pre_trial_start & nosepoke_ts < trial_start);
trial_nosepokes = nosepoke_ts(nosepoke_ts >= trial_start & nosepoke_ts <= trial_end);
pre_trial_rewards = reward_ts(reward_ts >= pre_trial_start & reward_ts < trial_start);
trial_rewards = reward_ts(reward_ts >= trial_start & reward_ts <= trial_end);

% Count and calculate intervals
trialResults.nosepokes_before_tone_count = numel(pre_trial_nosepokes);
trialResults.nosepokes_during_tone_count = numel(trial_nosepokes);
trialResults.intervals_before_tone = diff(pre_trial_nosepokes);
trialResults.intervals_during_tone = diff(trial_nosepokes);
trialResults.rewarded_intervals_not_during_cue = calculateRewardedIntervals(pre_trial_nosepokes, pre_trial_rewards);
trialResults.rewarded_intervals_during_cue = calculateRewardedIntervals(trial_nosepokes, trial_rewards);

% Calculate reward probabilities
trialResults.reward_probability_not_during_cue = calculateRewardProbability(pre_trial_nosepokes, pre_trial_rewards);
trialResults.reward_probability_during_cue = calculateRewardProbability(trial_nosepokes, trial_rewards);

% Calculate latencies to rewards for nosepokes during the tone
trialResults.latencies_to_rewards = calculateLatencies(trial_nosepokes, trial_rewards);
% Calculate average latency, ignoring NaN values
trialResults.avg_latency = mean(trialResults.latencies_to_rewards, 'omitnan');
end

function rewardedIntervals = calculateRewardedIntervals(nosepoke_ts, reward_ts)
% Calculate intervals between rewarded nosepokes
rewardedIntervals = diff(nosepoke_ts(ismember(nosepoke_ts, reward_ts)));
end

function rewardProbability = calculateRewardProbability(nosepoke_ts, reward_ts)
% Calculate the probability of a reward following a nosepoke
rewardedNosepokes = sum(ismember(nosepoke_ts, reward_ts));
totalNosepokes = numel(nosepoke_ts);
rewardProbability = rewardedNosepokes / max(totalNosepokes, 1); % Avoid division by zero
end

function latencies = calculateLatencies(nosepoke_ts, reward_ts)
% Calculate latencies from each nosepoke to the first subsequent reward
latencies = [];
for i = 1:length(nosepoke_ts)
    subsequent_rewards = reward_ts(reward_ts > nosepoke_ts(i));
    if ~isempty(subsequent_rewards)
        % Find the minimum time difference between nosepoke and subsequent rewards
        latencies(end+1) = min(subsequent_rewards) - nosepoke_ts(i);
    else
        % If there are no subsequent rewards, append NaN or some other placeholder
        latencies(end+1) = NaN; % Use NaN to represent cases with no subsequent reward
    end
end
end

function aggregate_results = aggregateMetricsEveryNTrials(results, N)
% Aggregate metrics over every N trials
n_blocks = ceil(length(results) / N);
aggregate_results = struct('avg_latency', {}, 'sum_nosepokes_during', {}, 'sum_nosepokes_before', {});

for i = 1:n_blocks
    blockStart = (i - 1) * N + 1;
    blockEnd = min(i * N, length(results));
    blockResults = results(blockStart:blockEnd);

    % Aggregate calculations here (implementation depends on desired metrics)
    % Example: aggregate_results(i).avg_latency = mean([blockResults.avg_latency]);
end
end
