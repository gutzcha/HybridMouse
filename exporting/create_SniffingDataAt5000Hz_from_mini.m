function new_vec = create_SniffingDataAt5000Hz_from_mini(input_vec,input_fs,output_fs, trigger_time)
% todo : 

%Determine a rational approximation to the ratio of the new sample rate
[P,Q] = rat(output_fs/input_fs); 
%Resample while applying FIR Antialiasing Lowpass Filter
new_vec = resample(input_vec,P,Q);
