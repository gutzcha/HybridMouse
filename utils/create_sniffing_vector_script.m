
%This script converts audio recordings from the mini-mic and creates a
%sniffing vector

%Step1: load the data to workspace
%Step 2: edit the variables in the script 
%Step 3: run the script


%==================================================================
%========================Edit these line===========================


input_vec = data; %Change data to the variable name in workspace
trigger_time = 0; %Set trigger time
save_file_name = 'Rat3-day2-Free-male.mat' ; %Set file name to save with prefix

% prefix = 'synthesised_SniffingSensor_Data';


%Example: 

%{

trigger_time = 0;
input_vec = data
save_file_name = Rat3-day2-Free-male_synthesised_SniffingSensor_Data.mat

->run
%}

%===========================Stop edit================================
%====================================================================



sniffing_vec = create_sniffing_vector(input_vec,trigger_time);
[ret, msg] = save_sniffing_file(sniffing_vec, save_file_name)


function new_vec = create_sniffing_vector(input_vec, trigger_time)
    input_fs = 2.5e5;
    output_fs = 5000;
    
      
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

function [ret, msg] = save_sniffing_file(sniffing_vec, save_file_name)
    SniffingDataAt5000Hz = sniffing_vec;
    SniffingSensorAuxChannel = 1;
    
    try
        save(save_file_name,'SniffingDataAt5000Hz','SniffingSensorAuxChannel');
        disp('The file was saved')
    catch ME
        ret = false;
        msg = ME.message;
        disp('There was an error')
        disp(msg)
    end

end