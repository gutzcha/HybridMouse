function hROI = addRoi(ax,source,timeRange,shapeToAdd,timeShift)
hROI = [];
%Delete existing rois
delete(findobj(ax,'Type','images.roi.rectangle'));



%Redraw
tab = source.roiTable;

labels = tab.Label;
if isempty(labels)
    return
end

labelColors = source.allROIs;
if ~exist('timeRange','var')||isempty(timeRange)
    timeRange = source.times;
end

if ~exist('shapeToAdd','var')||isempty(shapeToAdd)
   shapeToAdd = 'rect' ;
end


if ~exist('timeShift','var')||isempty(timeShift)
   timeShift = 0 ;
end

YLIM = ax.YLim;
positions = [tab.TimeStart,repmat(YLIM(1),size(tab.TimeStart)),tab.TimeEnd-tab.TimeStart,repmat(diff(YLIM),size(tab.TimeStart))];
inds = find((timeRange(1)<=positions(:,1))&(timeRange(2)>=(positions(:,1)+positions(:,3))));
positions(:,1) = positions(:,1)+timeShift;
hROI = [];
for i = inds'
    
    thisLabels = labels(i,1);
    col =labelColors{ismember(labelColors.Name,thisLabels),'Color'};
    if size(col,1)>1
        warning('Multiple labels with same name')
        col = col(1,:);
    end
    if isempty(col)
        col = [0 0 0] ;
        warning('Could not find lable "%s", and was added',thisLabels)
        source.addedTypeList = {thisLabels,col,""};
        app.updateROITypeList;
    end
    switch shapeToAdd
        case 'rect'
    hROITemp = drawrectangle('Parent',ax,"Position",positions(i,:),"Label",labels(i,:),"Color",col) ;
    addRoiFunctions(hROITemp,source)
        case 'line'
            hROITemp = xline(ax,positions(i,1)+positions(i,3)/2,"Color",col);            
    end
hROI = [hROI;hROITemp];
end
end


