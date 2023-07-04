function [truth_vec, error_flag] = create_logical_vec_from_table(tab,total_len,fs,type,total_len_bin)
%create_logical_vec_from_csv(tab,total_len,fs)
%  Create a logical vector from a table
%{
Create a labeld logical vector with the length as described
inputs:
    tab - a table with onset and offset of each syllable 
        - can also be a matrix of onsets and offsets Nx2

    total_len - total length of the output vector in seconds
    fs - sampling rate - defulte is 2.5e5
output:
    truth vec - a logical vector of length total_len*fs
  %} 

error_flag.value = false;
error_flag.msg = [];
error_flag.ind_to_remove = [];

if ~exist('fs','var')||isempty(fs)
   fs = 2.5e5; 
end

if istable(tab)
    tab = table2array(tab);
end
if ~exist('type','var')
    type = 'audio';    
end

%If the table is empty, return an array of zeros
if ~exist('total_len_bin','var')||isempty(total_len_bin)
    total_len_bin = ceil(total_len*fs);
end

if isempty(tab)||(type=="noise")
    truth_vec = false(1,total_len_bin);
    return
end

%Check table for error and remove bad lines
tab(tab>total_len) = total_len;
error_flag = check_input(tab,total_len);
if error_flag.value
    bad_lines = tab(error_flag.ind_to_remove,:);
    tab(error_flag.ind_to_remove,:) = [];
end

syl_onsets = tab(:,1);
syl_offsets = tab(:,2);

temp_mat_bins = floor([syl_onsets,syl_offsets]*fs); %Abs syl onset and offset in bins

num_of_syllables = size(temp_mat_bins,1);
truth_vec = false(1,total_len_bin); %Total number of bins is = lenght of output

for i1 =1:num_of_syllables
    inds = 1+temp_mat_bins(i1,1):temp_mat_bins(i1,2);
    if numel(inds)==0
        inds = temp_mat_bins(i1,1)+1;
    end
      truth_vec(inds)=true  ;
end

function error_flag = check_input(tab,total_len)
error_flag.value = false;
error_flag.msg = [];
ind_to_remove = [];

non_numeric = any(~arrayfun(@isnumeric,tab),2)|any(arrayfun(@isnan,tab),2);
non_positive = any(tab<0,2);
ind = any([non_numeric,non_positive],2);

if sum(ind)   
    ind = find(ind);
    ind_to_remove = [ind_to_remove;ind];
    
    msg = sprintf("All time stamps must be positive numbers at line/s %d",ind);
    error_flag.msg = [error_flag.msg;msg];
end

ind = any((tab(:,2)-tab(:,1))<=0,2);
if sum(ind)
    ind = find(ind);
    ind_to_remove = [ind_to_remove;ind];
    
    msg = sprintf('Invalid time stamps at line/s %d\n',ind);
    error_flag.msg = [error_flag.msg;msg] ;
    
end

%If only the last timestamp is too long, then change it to max, as long
    %as it is still valid
tab(end,end) = min(tab(end,2),total_len);

%Make sure that no time stamps are in range
ind = any(tab> total_len,2);
if sum(ind)
    
    ind = find(ind);
    ind_to_remove = [ind_to_remove;ind];
    msg = sprintf('All values in the table must be less than the total length of the segment- line/s %d\n',ind);
    error_flag.msg = [error_flag.msg;msg] ;
end

%report bad lines to be ingored
error_flag.ind_to_remove = unique(ind_to_remove);
if ~isempty(error_flag.ind_to_remove)
    error_flag.value = true;
end
