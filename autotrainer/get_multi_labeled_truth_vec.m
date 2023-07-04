function [truth_vec] = get_multi_labeled_truth_vec(tab ,total_len_bin,fs ,all_labels,input_type)
%GET_MULTI_LABELED_TRUTH_VEC Summary of this function goes here
%   Detailed explanation goes here

if ~exist('all_labels','var')||isempty(all_labels)
    all_labels = {};
end
if ~iscell(all_labels)
    all_labels = {all_labels};
end

if ~exist('input_type','var')||isempy(input_type)
    input_type  = 'signal';
end

tab_times = tab{:,["TimeStart","TimeEnd"]};
labels_arr = categorical(tab{:,"Label"});

num_labels = numel(all_labels);
truth_vec = false(num_labels,total_len_bin);


%Create labels vector
for ilabel = 1:num_labels
    if num_labels==1
        tab_sub = tab_times;
    else
        label = all_labels{ilabel};
        tab_sub = tab_times(labels_arr==label,:);
    end
    
    [truth_vec_temp, ~] = create_logical_vec_from_table_v2(tab_sub,total_len_bin,fs,input_type);
    %             truth_vec_temp = truth_vec_temp';
    % adjust legnth of truth_vec_temp
    
    truth_vec(ilabel,:) = truth_vec_temp;
end

%inverse bg
bg_ind = find(contains(all_labels,'bg'));
if ~isempty(bg_ind)
    truth_vec(bg_ind,:) = ~any(truth_vec(~contains(all_labels,'bg'),:),1);
end

truth_vec = truth_vec';
end

