function roiList = getRoi(ax)
%Get list of ROIs 

while 1
    roiList = [];
    roi =  drawrectangle(ax);
    if  isempty(roi.Position)
        break
    end
    roiList = [roiList;roi]; %#ok
end

