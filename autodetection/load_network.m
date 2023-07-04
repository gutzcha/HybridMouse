function [ret,net,networkpath] = load_network(networkpath,try_to_load)
if nargin == 1
    try_to_load = true;
end
ret = false;
net = [];


try
    net = matfile(networkpath);
    net = net.net;
    ret = true;
    fprintf('Loading network %s\n',networkpath)
catch ME
    if try_to_load
        warning('Unable to find net, please choose net to load:')
        warning(ME.message)
        [ret,net,networkpath] = getNewNet(networkpath);
    end    
        if ~((isa(net,'DAGNetwork'))||isa(net,'nnet.cnn.LayerGraph'))
            warning('Unable to assign network')
            return
        end
    
end
end