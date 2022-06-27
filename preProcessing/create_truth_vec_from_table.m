function [truth_vec] = create_truth_vec_from_table(tab,total_len,fs)
%{
Create a labeld logical vector with the length as described
inputs:
    tab - a table with onset and offset of each syllable and thier absolute
        onset
        - can also be a matrix of onsets and offsets Nx2

    total_len - total length of the output vector in seconds
    fs - sampling rate - defulte is 2.5e5
output:
    truth vec - a logical vector of length total_len*fs
  %} 

if ~exist('fs','var')||isempty(fs)
   fs = 2.5e5; 
end
if istable(tab)
    abs_onsets = tab.('absonset');  %Onset of call (s)
    syl_onsets = tab.('Sonset');    %Onset of syllable within a call (ms)
    syl_offsets = tab.('Soffset');  %Onffset of syllable within a call (ms)
else
    abs_onsets = zeros(size(tab,1),1);
    syl_onsets = tab(:,1);
    syl_offsets = tab(:,2);
end


temp_mat_bins = floor([abs_onsets+syl_onsets/1000,abs_onsets+syl_offsets/1000]*fs); %Abs syl onset and offset in bins

total_len_bin = ceil(total_len*fs);
num_of_syllables = size(temp_mat_bins,1);

truth_vec = false(1,total_len_bin); %Total number of bins is = lenght of output

for i1 =1:num_of_syllables
    inds = 1+temp_mat_bins(i1,1):temp_mat_bins(i1,2);
    if numel(inds)==0
        inds = temp_mat_bins(i1,1)+1;
    end
      truth_vec(inds)=true  ;
end

