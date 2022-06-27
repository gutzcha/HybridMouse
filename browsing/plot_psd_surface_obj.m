function plot_psd_surface_obj(obj,normalize_method)
%normalize_method : 1- normlize per segment, 2 - normalize all 
%
if ~exist('normalize_method','var')||isempty(normalize_method)
   normalize_method = 1; 
end

%Get all audio segment
allseg = obj.getAllSegments;

%Extract psd from all segment
allseg = cellfun(@get_power_vector,allseg,'UniformOutput',false);

%Reshape to matrix
allseg = cell2mat(allseg);

%Normlize
if normalize_method==1
     allseg = normlize_per_file(allseg);
elseif normalize_method==2
    allseg = normlize_all(allseg);
else
    error('METHOD must be 1 or 2 but it was %d instead',normalize_method)
end


f = figure;
ax = axes(f);
surface(ax,allseg','EdgeColor',"none")
ax.XLabel.String = 'Call Number';
ax.YLabel.String = 'Frequency (kHz)';
ax.FontName = 'Times New Roman';
ax.FontSize = 12;

f = figure;
ax = axes(f);
mean_vec = mean(allseg,1);
std_vec = std(allseg,1);
x = 1:numel(mean_vec);


x2 = [x, fliplr(x)];
inBetween = [mean_vec-std_vec, fliplr(mean_vec+std_vec)];
fill(ax,x2, inBetween, 'g','FaceAlpha',0.1);
hold on
plot(ax,x,mean_vec)

end
function mat = normlize_per_file(mat)
% freq_range = [20:100];
freq_range = [1:125];
max_mat = max(mat(:,freq_range),[],2);
min_mat = min(mat(:,freq_range),[],2);
mat = (mat-min_mat)./(max_mat - min_mat);
end

function mat = normlize_all(mat)
mat = rescale(mat,0,1);
end

