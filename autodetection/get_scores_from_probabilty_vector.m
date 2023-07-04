function scores = get_scores_from_probabilty_vector(probability_vectors,tab,fs)
%     scores = get_scores_from_probabilty_vector(probability_vectors,tab,fs)

%Model max - 0.91639

% model_max = 0.91639;
model_max = 1;

onSets = tab.Onsets;
offSets = tab.Offsets;

sz = size(onSets,1);

probability_vectors = reshape(probability_vectors,1,[]);
%option 1
%scoreFun = @(x,t1,t2,fs) mean(x(1+floor((t1*fs)):floor((t2*fs)))); 
%option 2
n_bins = numel(probability_vectors);
scoreFun = @(x,t1,t2,fs) max(x(1+(floor((t1*fs))):(min([ceil((t2*fs)),n_bins])))); 
scores = zeros(sz,1);

%Make sure that time stamps are within the array length
onSets = max(0,onSets);
offSets = min(offSets,n_bins/fs);



for i = 1:sz
%     disp(i)
    scores(i) = scoreFun(probability_vectors,...
        onSets(i),offSets(i),fs);
end
scores = scores./model_max; %Normalize to model max
end