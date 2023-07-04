function log = fix_missing_electro_files(tab_in)
%FIX_MISSING_ELECTRO_FILES This function creates missing files when the
%electropysiology setup was not working while recording
%   Create dummy and synthetic veriables in case the electrophysiology
%   setup was not working during the recording, and saves the files in the
%   same directory as the origin or as determained by user.
%   Arguments
%       tab_in: a table with a list of files to create
%   Outputs
%       log: list of files that were not created 
sz = size(tab_in,1);
log = cell(sz,1);
for i1 = 1:sz
%     try
        fprintf('Processing file num %d of %d\n', i1, sz)
        all_operations(tab_in(i1,:))
%     catch ME
%         warning('File %d, error msg: %s\n',i1, ME.message)
%        log{i1,1} = ME.message;
%         continue
%     end

%         all_operations(tab_in(i1,:))
end

end

function new_file_name = gen_file_name(file_name, prefix)
    %remove sufix
    [~,file_name] = fileparts(file_name);
    new_file_name = sprintf('%s_%s',file_name,prefix);
end

function save_name = get_save_name(origin_folder,dest_folder,file_name)
    if isempty(dest_folder)
        dest_folder = origin_folder;
    end
    
    if ~exist(dest_folder,'dir')
       mkdir(dest_folder);
    end
    
    save_name = fullfile(dest_folder,[file_name,'.mat']);
end

function new_vec = create_sniffing_vector(input_vec_path, trigger_time)
    input_fs = 2.5e5;
    output_fs = 5000;
    
    %Load audio vector
    input_vec = load_audio_file(input_vec_path);
    
    %Determine a rational approximation to the ratio of the new sample rate
    [P,Q] = rat(output_fs/input_fs); 
    %Resample while applying FIR Antialiasing Lowpass Filter
    new_vec = resample(input_vec,P,Q);
    
    %If trigger time is negative then concatinate zeros to the vector 
    if trigger_time<0
        bins_to_add = floor(output_fs*(trigger_time)*-1);
        new_vec = [zeros(1,bins_to_add), new_vec];
    end
end

function input_vec = load_audio_file(input_vec_path, verbose)
    if ~exist('verbose', 'var') ||isempty(verbose)
       verbose = false; 
    end
    
    %Determain if the audio file is wav or mat
    [~,~,sufix] = fileparts(input_vec_path);
    if verbose
        fprintf('File %s was loaded\n',input_vec_path)
    end
    
    switch lower(sufix)
        
        case '.wav'
            input_vec = load_audio_wav(input_vec_path);
        case '.mat'
            input_vec = load_audio_mat(input_vec_path);
        otherwise
                  error('Unsupported audio format')     
    end
    
    function input_vec = load_audio_mat(input_vec_path)
        %load audio vector
        
        vec_mat = matfile(input_vec_path);
        vec_classes = {whos(vec_mat).class}=="double";
        vec_names = {whos(vec_mat).name};
        
        
        %load first double
        input_vec = eval(sprintf('vec_mat.%s',vec_names{vec_classes(1)}));
        input_vec = reshape(input_vec,1,[]);
    end


 function input_vec = load_audio_wav(input_vec_path)
        input_vec = audioread(input_vec_path);
        input_vec = reshape(input_vec,1,[]);
    end
end



function FrameTimeInSec = create_FramesTimesInSec(total_num_frames, audio_length, video_fps, audio_path)
    
    if isempty(video_fps)||isnan(video_fps)
       %If the FPS is not available, estimate it 
       % Estimate video FPS
       if isempty(audio_length)||isnan(audio_length)
           %If the audio length is also unknow, load the file and get it
           audio_fs = 2.5e5;
           verbose = true;
           input_vec = load_audio_file(audio_path, verbose);
           audio_length = length(input_vec)/audio_fs;
       end
        video_fps = audio_length/total_num_frames;
    end
    

    % Create time vector
    FrameTimeInSec = (1:total_num_frames)./video_fps ;
end

function [ret, msg] = save_sniffing_file(sniffing_vec, save_file_name)
    SniffingDataAt5000Hz = sniffing_vec;
    SniffingSensorAuxChannel = 1;
    
    try
        save(save_file_name,'SniffingDataAt5000Hz','SniffingSensorAuxChannel');
    catch ME
        ret = false;
        msg = ME.message;
    end

end
    


function [ret, msg] = save_opto_file(FramesTimesInSec,...
    StartTriggerTimeInSec, StimulusInsertionFrame, StimulusRemovalFrame, save_file_name)

ret = true;
msg = '';

%list of output variables
    ElectrodeVsRegisteredAreasNum = [];
    FrameNumbersForbehavEventsSynchronization = [0, 0];
    FramesTimesDigitalChannel = 3;
    FramesTimesDigitalChannel_2 = 0;
    FramesTimesDigitalChannel_3 = 0;
    FramesTimesDigitalChannel_4 = 0;
    FramesTimesDigitalChannel_5 = 0;
    FramesTimesInSec = FramesTimesInSec;
    FramesTimesInSecThermalCamera = [];
    FramesTimesInSec_2 = [];
    FramesTimesInSec_3 = [];
    FramesTimesInSec_4 = [];
    FramesTimesInSec_5 = [];
    NonActiveElectroChannels = [];
    OptoStimEndTimesInSec = [];
    OptoStimStartTimesInSec = [];
    OptoStimTimesDigitalChannel = 0;
    StartTriggerTimeInSec = StartTriggerTimeInSec;
    StartingTriggerDigitalChannel = 2;
    StimulusInsertionFrame = StimulusInsertionFrame;
    StimulusRemovalFrame = StimulusRemovalFrame;
    ThermalFramesTimesDigitalChannel = 0;

if StartTriggerTimeInSec<0
    warning('Negative trigger time, setting to 0 and adjusting sniffing vector')
   StartTriggerTimeInSec = 0;  
end

%list of output variables
varialbels_list = {...
    'ElectrodeVsRegisteredAreasNum';...
    'FrameNumbersForbehavEventsSynchronization';...
    'FramesTimesDigitalChannel';...
    'FramesTimesDigitalChannel_2';...
    'FramesTimesDigitalChannel_3';...
    'FramesTimesDigitalChannel_4';...
    'FramesTimesDigitalChannel_5';...
    'FramesTimesInSec';...
    'FramesTimesInSecThermalCamera';...
    'FramesTimesInSec_2';...
    'FramesTimesInSec_3';...
    'FramesTimesInSec_4';...
    'FramesTimesInSec_5';...
    'NonActiveElectroChannels';...
    'OptoStimEndTimesInSec';...
    'OptoStimStartTimesInSec';...
    'OptoStimTimesDigitalChannel';...
    'StartTriggerTimeInSec';...
    'StartingTriggerDigitalChannel';...
    'StimulusInsertionFrame';...
    'StimulusRemovalFrame';...
    'ThermalFramesTimesDigitalChannel'};
try
    save(save_file_name,varialbels_list{:});
catch ME
    ret = false;
    msg = ME.message;
end

end

function all_operations(tab_in)

%extract data from table
src_folder = tab_in.folder_name;

mini_file_name = tab_in.mini_file_name;
mini_file_sufix = tab_in.mini_file_sufix;
mini_file_name = strcat(mini_file_name,".",mini_file_sufix);

file_name = tab_in.file_name;
main_file_sufix = tab_in.main_file_sufix;
file_name = strcat(file_name,".",main_file_sufix);

save_folder_name = tab_in.save_folder_name;
stimulus_insertion_frame = tab_in.stimulus_insertion_frame;
stimulus_removal_frame = tab_in.stimulus_removal_frame;
start_trigger_time_in_sec = tab_in.start_trigger_time_in_sec;
total_number_of_frames = tab_in.total_number_of_frames;
legnth_of_audio_in_sec  = tab_in.legnth_of_audio_in_sec;
video_fps = tab_in.video_fps;

audio_path = fullfile(src_folder, file_name);

%Create optogenetics files
%create FrameTimeInSec variable
FramesTimesInSec = create_FramesTimesInSec(total_number_of_frames, legnth_of_audio_in_sec, video_fps, audio_path);

%get save file names
prefix = 'synthesised_Behavior_and_Optogenetics_TimeStamps';
file_name_opto = gen_file_name(file_name, prefix);
save_name_opto = get_save_name(src_folder,save_folder_name,file_name_opto);

%save files with all variables
save_opto_file(FramesTimesInSec,...
    start_trigger_time_in_sec, stimulus_insertion_frame, stimulus_removal_frame, save_name_opto);

%Create sniffing files
%create sniffing vector
input_vec_path = fullfile(src_folder,mini_file_name);
sniffing_vec = create_sniffing_vector(input_vec_path, start_trigger_time_in_sec);

%get save file names
prefix = 'synthesised_SniffingSensor_Data';
file_name_sniffing = gen_file_name(file_name, prefix);
save_name_sniffing = get_save_name(src_folder,save_folder_name,file_name_sniffing);
% save sniffing file
save_sniffing_file(sniffing_vec, save_name_sniffing);
end

