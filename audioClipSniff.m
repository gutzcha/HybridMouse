classdef audioClipSniff<audioClip
    %AUDIOCLIPSNIFF Summary of this class goes here
    %   Detailed explanation goes here

    properties
        extracted_features
        labels_vector
        t_stamps_vector



        feature_win = hamming(100,'periodic');
        feature_overlap = 98;
        feature_fft = 100;

        reshape_win = 50;
        reshape_overlap = 2;

        class_id_vector


    end

    properties(Dependent)
        feature_extrator
        feature_fs

        resampled_audio_vec

        reshaped_extracted_features
        reshaped_labels_vector
        reshaped_labels_vector_2

        reshaped_t_stamps_vector
        reshaped_resampled_audio_vec

    end


    methods
        function obj = audioClipSniff(audioClipMini)
            if nargin==1
                obj = audioClipSniff.copy_properties(audioClipMini);
                refresh_all_features(obj);
            end
        end

        function refresh_all_features(obj)

            refresh_extracted_features(obj)
            refresh_labels_vector(obj)
            refresh_t_stamps_vector(obj)
        end

        function refresh_extracted_features(obj)
            obj.extracted_features =  extract(obj.feature_extrator,obj.vec);
        end

        function refresh_feaure_fs(obj)
            obj.extracted_features =  extract(obj.feature_extrator,obj.vec);
        end

        function refresh_labels_vector(obj)

            total_len_bin = numel(obj.extracted_features);
            tab = obj.roiTable;
            tab = [tab.TimeStart, tab.TimeEnd];
            input_type = 'audio';

            [truth_vec, ~] = create_logical_vec_from_table_v2(tab,total_len_bin,obj.feature_fs,input_type);
            obj.labels_vector = truth_vec';
        end

        function refresh_t_stamps_vector(obj)

            total_len_bin = numel(obj.extracted_features);
            total_len_sec = obj.audioLen;
            t_stamps_flat = linspace(0,total_len_sec,(total_len_bin+1));
            t_stamps_flat = t_stamps_flat(1:end-1);

            obj.t_stamps_vector = t_stamps_flat;

        end
        % Dependet properties
        function ret = get.feature_extrator(obj)

            ret =  audioFeatureExtractor( ...
                SampleRate=obj.fs, ...
                Window=obj.feature_win, ...
                OverlapLength=obj.feature_overlap, ...
                spectralCentroid=true);
        end
        function ret= get.feature_fs(obj)
            total_len_bin = numel(obj.extracted_features);

            ret = total_len_bin/obj.audioLen;
        end

        function ret = get.resampled_audio_vec(obj)
            new_fs = obj.feature_fs;
            old_fs = obj.fs;
            ret = downsample_audio(obj.vec,old_fs,new_fs);

        end

        function ret = get.reshaped_extracted_features(obj)
            ret = obj.reshape_vector(obj.extracted_features, obj.reshape_win, obj.reshape_overlap)';
        end
        function ret = get.reshaped_labels_vector(obj)
            ret = obj.reshape_vector(obj.labels_vector, obj.reshape_win, obj.reshape_overlap)';
            ret = logical(ret);
            ret = any(ret,2);
        end

%         function ret = get.reshaped_labels_vector_2(obj)
%             tab = obj.roiTable;
%             tab = [tab.TimeStart, tab.TimeEnd];
%             total_len_sec = obj.audioLen;
%             total_len_bin = obj.feature_winreshape_win;
%             reshape_feature_fs = n_tot/total_len_sec;
%             input_type = 'audio';
%             [ret, ~] = create_logical_vec_from_table_v2(tab,total_len_bin,reshape_feature_fs,input_type);
% 
%         end

        function ret = get.reshaped_t_stamps_vector(obj)
            ret = obj.reshape_vector(obj.t_stamps_vector, obj.reshape_win, obj.reshape_overlap)';
        
        end
        
        function ret = get.reshaped_resampled_audio_vec(obj)
            ret = obj.reshape_vector(obj.resampled_audio_vec, obj.reshape_win, obj.reshape_overlap)';
        end

        

    end
    methods(Static)
        function obj_copy = copy_properties(source_obj,p)

            obj_copy = audioClipSniff;
            if ~exist('p','var')||isempty(p)
                p = properties(source_obj);
            end
            if ~iscell(p)
                p = {p};
            end
            num_att = length(p);

            for ip = 1:num_att
                this_att = p{ip};
                mp = findprop(audioClip,this_att);
                % Dont copy constants and dependant properties
                if mp.Dependent || mp.Constant || (mp.SetAccess=="private")
                    continue
                end
                obj_copy.(this_att) = source_obj.(this_att);
            end

        end

        function ret = reshape_vector(vec, win, overlap)
            vec = reshape(vec,[],1);
            v_len = length(vec);
            hop_len = (win-overlap);
            n_hops = floor((v_len-overlap)/hop_len);
            ret = zeros(win,n_hops);
            ind_s=1;
            for is = 1:n_hops

                ind_e = ind_s + win -1;

                v = vec(ind_s:ind_e);
                ret(:,is) = v;
                ind_s = ind_e + 1 - overlap;
            end
        end
    end
end

