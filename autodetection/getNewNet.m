
function [ret,net,networkpath] = getNewNet(networkpath)
[file,path] = uigetfile('*.mat');
ret = false;
net = [];

if file
    netmatfile = matfile(fullfile(path,file));
else
    disp('No network was chocen')
    sprintf('The current network path is %s\n',networkpath)
    return
end

try
    net = netmatfile.net;
    if ~isa(net,'DAGNetwork')
        net = [];
        disp('No network was chocen')
    else
        networkpath = fullfile(path,file);
        fprintf('The network has been changed to %s\n',networkpath);
        ret = true;
    end
    
catch ME
    warning('Must have a variable names -net- ')
    warning(ME.message)
end

end
