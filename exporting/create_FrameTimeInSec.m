
function FrameTimeInSec = create_FrameTimeInSec(total_num_frames,insertion_frame,extraction_frame,insertion_time,extraction_time)

% todo the function has to read table and save all variables, add name to
% the end 

% Set Frames of insersion and extraction of stimulus
num_frames_ref = extraction_frame-insertion_frame;

% Set opening and closing times of closet
ref_time = extraction_time-insertion_time;

% Estimate video FPS
video_fps = ref_time/num_frames_ref;

% Create time vector
FrameTimeInSec = (1:total_num_frames).*video_fps ;




