function mat_in_rescales =randomaize_scale(mat_in,scale_range)
    max_bin = max(abs(mat_in),[],1);
    mat_in_rescales = mat_in./max_bin;
    num_files = size(mat_in,2);
    scaling_factors = scale_range*rand(1,num_files);
    mat_in_rescales = mat_in_rescales.*scaling_factors;
end