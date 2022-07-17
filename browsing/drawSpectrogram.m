function [sh,psd,dataIn,axSpect, axP] = drawSpectrogram(axSpect, axP, dataIn,max_bins_to_plot)
%Draw a spectrogram and a normlized pds
%   Draw plots and update plots accordingly

%{
Inputs:
    axSpect - handel for axis on which to plot the spectrogram
    axP     - handel for axis on which to plot the power spectral density
    dataIn  - an instanse of audioClip class
Outputs:
    sh        - handel for surf plot of spectrogram
    psd       - handel for the line plot psd
    dataOut   - an instanse of audioClip class, updated
       
    
%}

%Initialize outputs
sh = [];
psd = [];
t = 0;
s = 1;
f= 0;
ps = [];

if nargin == 0
    figure1 = figure;
    % Create axes
    axP = axes('Parent',figure1,...
        'Position',[0.132 0.05 0.775 0.19]);
    % Create axes
    axSpect = axes('Parent',figure1,...
        'Position',[0.132 0.32 0.775 0.60]);
    
    return
end



%Disable axis interactivity - no zoom and no pan 
disableDefaultInteractivity(axP)
disableDefaultInteractivity(axSpect)

%Clear any surface plots on these axes
% hSurf = findobj(axSpect.Children,'Type','surface');

%Get data and properties from dataIn
vec = dataIn.vec;

fs = dataIn.fs;
win = dataIn.window;
overlap = dataIn.overlap;
fftWin = dataIn.fft;
th = dataIn.threshold;
clim = dataIn.clim;
times = dataIn.times;
vecLen = numel(vec);

% %Fix win length if the sampling rate was modified
% if fs<5001
%     %Normal bin length is seconds
%     win = 20;
%     fftWin = 20;
%     overlap = 18;
% end

%Fix times
if isempty(times)
    times(1)  = 0;
    times(2) = min(50,vecLen/fs);   %MAX SET AT 50 SEC IF EMPTY
end

if isempty(vec)
    return
end

%Set maximum bins to plot
if ~exist('max_bins_to_plot','var')||isempty(max_bins_to_plot)
   max_bins_to_plot = 2.5e7;
end
if diff(times)*fs>max_bins_to_plot
    warning('Atempted to plot an empty vector or the vector was too long')
    return
end

%Ylims
if isstruct(dataIn)
    if isfield(dataIn,'ylims')
        ylims = dataIn.ylims;
    else
        ylims = [];
    end
elseif isa(dataIn,'audioClip')
    if isprop(dataIn,'ylims')
        ylims = dataIn.ylims;
    else
        ylims = [];
    end
end

%CLim mode
if isstruct(dataIn)
    if isfield(dataIn,'climMode')
        climMode = dataIn.climMode;
    else
        climMode = 'auto';
    end
elseif isa(dataIn,'audioClip')
    if isprop(dataIn,'climMode')
        climMode = dataIn.climMode;
    else
        climMode = 'auto';
    end
end

%make sure vector is vertical

vec = reshape(vec,[],1);


%If it is possible, add window length from the vector from the beggining
%Check if there is a window value
AddedBefore = 0; %Number of padded bins
% AddedAfter = 0;
if ~isempty(win)
    winSec = (win)/fs; %window length in seconds
    %Check if you can go back, if not then zero pad
    startingPointShift = times(1) - winSec;
    if startingPointShift<0
        %If you can't go back, you are too close the the begining of the
        %vector - zero pad it
        AddedBefore = ceil(-startingPointShift*fs);
        vec = [zeros(AddedBefore,1);vec];
    end
    
    
    %Same for the end
    EndingPoint = times(2) + winSec;
    if EndingPoint>vecLen/fs
        %If the vector is too short, add zeros after it
        AddedAfter = ceil((EndingPoint*fs-vecLen));
        vec = [vec;zeros(AddedAfter,1)];
    end
end

times(1) = times(1) - winSec;
times(2) = times(2) + winSec;
segmentTimes = times;

segmentInds = floor(AddedBefore+(1+(segmentTimes(1))*fs):(segmentTimes(2)*fs));

vec = vec(segmentInds);


vec = vec./max(abs(vec)); %Normlize, I dont know why this is important but it is imposible to draw ROI if it is not normlized
vecLen = numel(vec);
%  if fs>50000
%     vec = imgaussfilt((vec),1);
%     vec2 = imgaussfilt((vec),1.5);
%     vec = vec-vec2;
%  end



%Perform fft on audio vector
if vecLen>max(win)
    try
        %Maybe add minithreshhold later
        [s,f,t,ps] = spectrogram(vec,hamming(win),overlap,fftWin,fs,'yaxis');
    catch ME
        warning(ME.message,'Invalid spectrogram parameters')
        [s,f,t,ps] = spectrogram(vec,hamming(win),[],[],fs,'yaxis');
    end
    t = t+times(1); %offset time vector
    %Create plot with set threshold and clim
    
    %Adjust ylims
    [s,f,ylims] = adjustYlims(s,f,ylims);
    try
    axSpect.YLim = ylims;
    catch ME
        
        disp("Unable to set ylims")
        disp(ME.message)
    end
    %     s = imgaussfilt(abs(s),1);
    
    %Convert to power
    s = 10*log10((abs(s)));
    
    %Downsample
    %Allways display a maximum of maxPoint time points
    
    maxPoint = 5000;
    if size(s,2) >maxPoint
        downRatio = ceil(size(s,2)/maxPoint);
        s = (downsample(s',downRatio))';
        vec = downsample(vec,downRatio);
        t = downsample(t,downRatio);
    end
    
    try
        sh = surf(axSpect,t,f,s,'EdgeColor','none');
    catch ME
        disp (ME)
        return
    end
    view(axSpect,2);
    switch climMode
        case 'manual'
            if ~isempty(clim)&&diff(clim)>0
                axSpect.CLim = clim;
            end
        case 'auto'
            axSpect.CLimMode = 'auto';
    end
    
    % colormap(axSpect,'gray')
    colormap(axSpect,'default')
    
    %Create psd plot
    % psSum = smooth(sum(ps),0.01);
    % psSum = psSum./max(psSum);
    %
    % psd = plot(axP,t,psSum);
    vec = vec-mean(vec); %Zero center the vector
    vec = vec./(max(abs(vec)));
    %Remove zeroPadding
    %     vec = vec((AddedBefore+1) :end-AddedAfter);
    tVec = linspace(times(1),times(end),numel(vec));
    psd = plot(axP,tVec,(vec));
    
    
    %return
end




%Make axes tight
if ~isempty(dataIn.times_)
    axSpect.XLim = dataIn.times_;
    axP.XLim = dataIn.times_;
else
    axSpect.XLim = dataIn.times;
    axP.XLim = dataIn.times;
    
    %     axis(axSpect,'tight')
    %     axis(axP,'tight')
end

%Set ylims
if ~isempty(ylims)
    axSpect.YLim = ylims;
end

%Update output values
dataIn.stft = s;
dataIn.cyclicalF = f;
dataIn.timeVec = t;
dataIn.psd = ps;
dataIn.clim = axSpect.CLim;
dataIn.ylims = axSpect.YLim;
dataIn.spectrogram_sample_rate = (numel(t))/(t(end)-t(1));
%Make sure ticks are automatic
axSpect.XTickLabelMode = 'auto';
axP.XTickLabelMode = 'auto';
% zoom(ancestor(axSpect,'Figure'),'reset')
end

function [s,f,ylims] = adjustYlims(s,f,ylims)
%Make sure that ylims are within f limits
if isempty(ylims)
    ylims = [floor(f(1)),ceil(f(end))] ;
    return
end
if ylims(1)<f(1)
    warning('Frequency lower limit if out of bound, setting to min')
    
    ylims(1)  = f(1);
end

if ylims(2)>f(end)
    warning('Frequency higher limit if out of bound, setting to max')
    ylims(2)  = f(end);
end

%find indecies of closes values to ylims
[minfVal,minfInd]= min(abs(f-ylims(1)));
[~,maxfInd] = min(abs(f-ylims(2)));

%recrop s and f
s = s(minfInd:maxfInd,:);
f = f(minfInd:maxfInd,:);
f = f+minfVal;


end

