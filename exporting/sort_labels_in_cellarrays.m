function [cellarray_new,LablesMap_new] = sort_labels_in_cellarrays(cellarray_in,LablesMap_old,LablesMap_new)
%Reorder cell array according to a new labels map

%find differences between the two sets
[C,~] = setdiff(LablesMap_old,LablesMap_new);
LablesMap_new = [LablesMap_new, C];

total_num_tags = numel(LablesMap_new);
out_size = [size(cellarray_in,1),total_num_tags];

cellarray_new = cell(out_size);
inds = cellfun((@(x) find(strcmpi(x,LablesMap_old))),LablesMap_new,'UniformOutput',false);
for i1 = 1:numel(inds)
    ind = inds{i1};
    if ~isempty(ind)    
        cellarray_new(:,i1) = cellarray_in(:,ind);
    end
end

