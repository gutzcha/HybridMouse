classdef audioClipSniff<audioClip
    %AUDIOCLIPSNIFF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %         parent_audioClip_obj = audioClip
        
        extracted_features
        labels_vector
        all_labels_vector
        t_stamps_vector
        
        
        feature_win_sec = 1;
        feature_overlap_sec = 0.95;
        
        reshape_win_sec = 2;
        reshape_overlap_sec = 1.9;
        
        class_id_vector
        
        n_features = []
        features_len_bin = []
        
        reshaped_features_len_bin = []
        
    end
    
    properties(Dependent)
        feature_win
        feature_overlap
        feature_extrator
        feature_fs
        
        reshape_win
        reshape_overlap
        
        resampled_audio_vec
        
        reshaped_extracted_features
        reshaped_extracted_features_mult
        reshaped_labels_vector_folded
        reshaped_labels_vector_direct
        reshaped_all_labels_vector
        
        reshaped_t_stamps_vector
        reshaped_resampled_audio_vec
        
        labels_before_after_vector
        
    end
    
    
    methods
        
        function obj = audioClipSniff(audioClipMini,config)
            
            if nargin>0
                obj = audioClipSniff.copy_properties(audioClipMini);
            end
            
            if nargin >1
                %                     obj.parent_audioClip_obj = audioClipMini;
                if isfield(config,'feature_win_sec')
                    obj.feature_win_sec = config.feature_win_sec;
                end
                
                if isfield(config,'feature_overlap_sec')
                    obj.feature_overlap_sec = config.feature_overlap_sec;
                end
                
                if isfield(config,'reshape_win_sec')
                    obj.reshape_win_sec = config.reshape_win_sec;
                end
                
                if isfield(config,'reshape_overlap_sec')
                    obj.reshape_overlap_sec = config.reshape_overlap_sec;
                end
            end
            if nargin>0
                refresh_all_features(obj);
            end
        end
        
        function refresh_all_features(obj)
            
            refresh_extracted_features(obj)
            refresh_t_stamps_vector(obj)
            refresh_labels_vector(obj)
            refresh_labels_vector_before_after(obj)
            
        end
        
        function refresh_extracted_features(obj)
            
            
            % Extract features
            
            f_all_flat = extract(obj.feature_extrator,obj.vec);
            
            
            
            total_len_sec = obj.audioLen;
            n_tot = size(f_all_flat,1);
            obj.features_len_bin = n_tot;
            
            
            n_features_ = size(f_all_flat,2);
            
            % Add downsampled version of the audio
            new_vec = obj.resampled_audio_vec;
            f_all_flat = [f_all_flat,new_vec];
            obj.n_features = n_features_ +1 ;
            
            % Normalize the features
            
            % f_all_flat = normalize(f_all_flat,1,'range',[-1,1]);
            f_all_flat = normalize(f_all_flat,1,'zscore');
            obj.extracted_features = f_all_flat;
            %             obj.refresh_feaure_fs;
            
        end
        
        
        %         function ret = refresh_feaure_fs(obj)
        %
        %             total_len_sec = obj.audioLen;
        %             n_tot = size(obj.extracted_features,1);
        %
        %             ret = n_tot/total_len_sec;
        %             obj.feature_fs = ret;
        %         end
        
        function n_tot = refresh_features_len_bin(obj)
            n_tot = size(obj.extracted_features,1);
            obj.features_len_bin = n_tot;
        end
        
        function refresh_labels_vector(obj)
            
            n_tot = size(obj.extracted_features,1);
            tab = obj.roiTable;
            tab = [tab.TimeStart, tab.TimeEnd];
            input_type = 'audio';
            
            [truth_vec, ~] = create_logical_vec_from_table_v2(tab,n_tot,obj.feature_fs,input_type);
            obj.labels_vector = truth_vec';
        end
        
        function refresh_labels_vector_before_after(obj)
            
            n_tot = size(obj.extracted_features,1);
            tab = obj.roiTable;
            
%             subtract_length = true;
            subtract_length = false;
            dt = mean(diff(obj.t_stamps_vector));
            %             k = 5*dt;
            k = 2.5;
            
            if subtract_length
                len_before = @(x,y) [max(x-k*(y-x),0),max(x,0)]; % subtract length
                len_after = @(x,y) [min(y,obj.audioLen), min(y+k*(y-x),obj.audioLen)]; % add length of vocalization
            else
                len_before = @(x,y) [max(x-k,0),max(x,0)]; % subtract constant
                len_after = @(x,y) [min(y,obj.audioLen), min(y+k,obj.audioLen)]; % add const
            end
            
            tab = [tab.TimeStart, tab.TimeEnd];
            tab_before = cell2mat(arrayfun(len_before, tab(:,1), tab(:,2), 'UniformOutput', false));
            tab_after = cell2mat(arrayfun(len_after, tab(:,1), tab(:,2), 'UniformOutput', false));
            
            
            input_type = 'audio';
            
            ret = zeros(n_tot,3);
            [truth_vec, ~] = create_logical_vec_from_table_v2(tab,n_tot,obj.feature_fs,input_type);
            ret(:,2) = truth_vec';
            
            [truth_vec, ~] = create_logical_vec_from_table_v2(tab_before,n_tot,obj.feature_fs,input_type);
            ret(:,1) = truth_vec';
            
            [truth_vec, ~] = create_logical_vec_from_table_v2(tab_after,n_tot,obj.feature_fs,input_type);
            ret(:,3) = truth_vec';
            
            % remover overlaps
            %             ret(logical(ret(:,2)),[1,3]) = 0;
            obj.all_labels_vector = ret;
        end
        
        
        
        function refresh_t_stamps_vector(obj)
            
            total_len_bin = obj.features_len_bin;
            total_len_sec = obj.audioLen;
            t_stamps_flat = linspace(0,total_len_sec,(total_len_bin+1));
            t_stamps_flat = t_stamps_flat(1:end-1);
            
            obj.t_stamps_vector = t_stamps_flat';
            
        end
        % Dependet properties
        function ret = get.feature_extrator(obj)
            
%                         ret =  audioFeatureExtractor( ...
%                             SampleRate=obj.fs, ...
%                             Window=obj.feature_win, ...
%                             OverlapLength=obj.feature_overlap, ...
%                             spectralCentroid=true);
            
            
%             ret =  audioFeatureExtractor( ...
%                 SampleRate=obj.fs, ...
%                 Window=obj.feature_win, ...
%                 OverlapLength=obj.feature_overlap, ...
%                 spectralCentroid=true,...
%                 spectralCrest=true,...
%                 spectralDecrease=true,...
%                 spectralEntropy =true,...
%                 spectralFlatness =true,...
%                 spectralFlux =true,...
%                 spectralKurtosis =true,...
%                 spectralRolloffPoint = true,...
%                 spectralSkewness = true,...
%                 spectralSlope = true,...
%                 spectralSpread = true,...
%                 linearSpectrum = true);
            
               
%             ret =  audioFeatureExtractor( ...
%                 SampleRate=obj.fs, ...
%                 Window=obj.feature_win, ...
%                 OverlapLength=obj.feature_overlap, ...
%                 spectralCentroid=true,...
%                 spectralCrest=true,...
%                 spectralDecrease=true,...
%                 spectralEntropy =true,...
%                 spectralFlatness =true,...
%                 spectralFlux =true,...
%                 spectralKurtosis =true,...
%                 spectralRolloffPoint = true,...
%                 spectralSkewness = true,...
%                 spectralSlope = true,...
%                 spectralSpread = true,...
%                 linearSpectrum = true);
            
                          ret =  audioFeatureExtractor( ...
                            SampleRate=obj.fs, ...
                            Window=obj.feature_win, ...
                            OverlapLength=obj.feature_overlap, ...
                            spectralCentroid=true,...
                            linearSpectrum=true...
                            );
        end
        function ret= get.feature_fs(obj)
            total_len_bin = obj.features_len_bin;
            
            ret = total_len_bin/obj.audioLen;
        end
        
        function ret = get.resampled_audio_vec(obj)
            
            
            n_tot = obj.features_len_bin;
            new_vec = downsample_audio(obj.vec,numel(obj.vec),n_tot);
            ret = reshape0(new_vec,[n_tot,1],'truncate');
            
        end
        function ret = get.feature_win(obj)
            sampleRate = obj.fs;
            ret = hamming(sampleRate*obj.feature_win_sec,'periodic');
        end
        function ret = get.feature_overlap(obj)
            sampleRate = obj.fs;
            ret = floor(sampleRate*obj.feature_overlap_sec);
        end
        
        function ret = get.reshape_win(obj)
            ret = max(1,floor(obj.feature_fs*obj.reshape_win_sec));
        end
        function ret = get.reshape_overlap(obj)
            ret = floor(obj.feature_fs*obj.reshape_overlap_sec);
        end
        
        
        function ret = get.reshaped_extracted_features(obj)
            
            
            % Fold the the data into overlaping windows
            f_all = obj.reshape_vector(obj.extracted_features, obj.reshape_win, obj.reshape_overlap);
            
            % Concatenate all features of each window to 1xdxm vector
            % When
            %   d = number of bins per feature per window
            %   m = number of features
            
            f_all_perm_reshape = [];
            for i1 = 1:obj.n_features
                f_all_perm_reshape = [f_all_perm_reshape;f_all(:,:,i1)]; %#ok
            end
            size(f_all_perm_reshape);
            ret = flip(f_all_perm_reshape,1)';
            
            obj.reshaped_features_len_bin = size(ret,1);
            
        end
        
        function f_all = get.reshaped_extracted_features_mult(obj)
            
            
            % Fold the the data into overlaping windows
            f_all = obj.reshape_vector(obj.extracted_features, obj.reshape_win, obj.reshape_overlap);
            % normalize data
            
            f_all = permute(f_all,[2,3,1]); % final shape - nxtxd
            
            % When
            %   n = number of observations
            %   t = number of bins per feature per window
            %   d = number of features per bin
            
            obj.reshaped_features_len_bin = size(f_all,1);
            
        end
        
        function ret = get.reshaped_labels_vector_folded(obj)
            ret = obj.reshape_vector(obj.labels_vector, obj.reshape_win, obj.reshape_overlap)';
            ret = logical(ret);
            ret = any(ret,2);
        end
        
        function ret_struct = get.reshaped_all_labels_vector(obj)
            ret = squeeze(obj.reshape_vector(obj.labels_before_after_vector, obj.reshape_win, obj.reshape_overlap));
            
            
            % the output in this case is one per time bin
            ret_struct.matrix = ret'; % final dim = nxt
            
            % if you want output one element per segment uncomment this
            % section
            
            f = @(c) mode([repmat(c(c~='o'),[1,3]),'oov']);
            bin_sz = size(ret,2);
            ret_char = repmat('o',[bin_sz,1]);
            for i = 1:bin_sz
                ret_char(i) = f(char(ret(:,i))');
            end
            ret_char = string(ret_char);
            ret_struct.vec = ret_char;
            
            %             ret = squeeze(any(ret,1))';
            %             % handle duplications
            %             ret([1,3],ret(2,:)) = 0;
            %
            
            
        end
        
        function ret = get.reshaped_labels_vector_direct(obj)
            n_after_reshape = obj.reshaped_features_len_bin;
            fs_after_reshape = n_after_reshape/obj.audioLen;
            
            tab = obj.roiTable;
            tab = [tab.TimeStart, tab.TimeEnd];
            input_type = 'audio';
            
            [truth_vec_direct, ~] = create_logical_vec_from_table_v2(tab,n_after_reshape,fs_after_reshape,input_type);
            ret = logical(truth_vec_direct)';
        end
        
        function ret = get.reshaped_t_stamps_vector(obj)
            ret = obj.reshape_vector(obj.t_stamps_vector, obj.reshape_win, obj.reshape_overlap)';
        end
        
        function ret = get.reshaped_resampled_audio_vec(obj)
            ret = obj.reshape_vector(obj.resampled_audio_vec, obj.reshape_win, obj.reshape_overlap)';
        end
        
        function ret = get.labels_before_after_vector(obj)
            labels_mat = obj.all_labels_vector;
            labels_mat = logical(labels_mat);
            keep_overlaping_before_after = labels_mat(:,1)&labels_mat(:,3);
            % what do we do when the tails overlap? maybe mark this as usv?
            % mabye mark this as something else
            sz = size(labels_mat,1);
            ret = repmat("o",[sz,1]);
            ret(labels_mat(:,1)) = "b";
            ret(labels_mat(:,3)) = "a";
            
            ret(keep_overlaping_before_after) = "i";
            
            ret(labels_mat(:,2)) = "v";
            
        end
        % Plots
        function plot_features(obj,truth_vec_pred, method)
            
            
            truth_vec = vertcat(obj.labels_vector);
            if ~exist('truth_vec','var')|| isempty(truth_vec_pred)
                truth_vec_pred = zeros(size(truth_vec));
            end
            
            if ~exist("method",'var')||isempty(method)
                method = 'lines' ;
            end
            max_plot = 100000;
            vec_in = vertcat(obj.t_stamps_vector);
            
            n_tot = sum([obj.features_len_bin]);
            v_range = [1:min(n_tot,max_plot)];
            
            
            vec_in = vec_in(v_range);
            
            
            data = vertcat(obj.extracted_features);
            data = data(v_range,:); %normalize(f_all_flat,1);
            
            figure;
            hold on
            switch method
                case 'lines'
                    k=1;
                    data = data + (1:k:(obj(1).n_features*k));
                    plot(vec_in, data')
                    min_y = min(data,[],'all')*0.8;
                    
                    max_y = max(data,[],'all')*1.2;
                    max_z = 1;
                case 'surface'
                    [X,Y] = meshgrid(vec_in, 1:obj(1).n_features);
                    
                    surface(X,Y, data','EdgeColor','none')
                    min_y = min(Y,[],'all')*0.8;
                    max_y = max(Y,[],'all')*1.2;
                    max_z = max(data,[],'all');
                otherwise
                    error('unidetified method')
            end
            
            
            
            
            
            plotpatch = @(x,w,c) patch( ...
                [x-w, x+w, x+w, x-w],[min_y, min_y, max_y, max_y],[max_z,max_z,max_z,max_z], ...
                c,'FaceAlpha',0.3, ...
                "EdgeColor",'none');
            
            
            truth_vec = truth_vec(v_range);
            x_vert = vec_in(truth_vec);
            dt = mean(diff(vec_in));
            w = 0.5*dt*ones(size(x_vert));
            c = repmat("blue",size(x_vert));
            arrayfun(plotpatch,x_vert,w,c)
            
            if sum(truth_vec_pred)>0
            truth_vec_pred = truth_vec_pred(v_range);
            x_vert = vec_in(truth_vec_pred);
            dt = mean(diff(vec_in));
            w = 0.5*dt*ones(size(x_vert));
            c = repmat("green",size(x_vert));
            arrayfun(plotpatch,x_vert,w,c)
            end
            
            
            hold off
        end
        
        function plot_folded_features(obj,plot_range,include_margin,labels_type)
            if ~exist('plot_range','var')||isempty(plot_range)
                plot_range = obj.times;
            end
            time_vec = obj.reshaped_t_stamps_vector;
            time_vec = time_vec(:,max(1,round(obj.reshape_win/2)));
            [add_vec, ret_ind] =  obj.time2ind(time_vec, plot_range);
            
            if ~exist('labels_type','var')||isempty(labels_type)
                labels_type = 'folded';
            end
            
            if ~exist('include_margin','var')||isempty(include_margin)
                include_margin = false;
            end
            
            data = obj.reshaped_extracted_features(ret_ind,:);
            data = data+add_vec;
            
            switch labels_type
                case 'folded'
                    labels_in = obj.reshaped_labels_vector_folded;
                    
                    
                case 'direct'
                    labels_in = obj.reshaped_labels_vector_direct;
                    
                otherwise
                    error('labels type must be direct or folded')
            end
            
            figure;
            hold on
            x = (1:size(data,2));  % number of features
            
            plot(x', data(:,:)')
            
            
            plotpatch = @(y,w,c) patch( ...
                [x(1), x(end)*1.1, x(end)*1.2, x(1)],[y-w, y-w, y+w, y+w], ...
                c,'FaceAlpha',0.4, ...
                "EdgeColor",'none');
            %
            labels_in = labels_in(ret_ind,:);
            dt = mean(diff(add_vec))/2;
            labels_vec = add_vec(labels_in);
            w = dt*ones(size(labels_vec));
            c = repmat("red",size(labels_vec));
            %          arrayfun(plotpatch,labels_vec,w,c)
            %          labels = labels_in_copy+add_vec;
            if include_margin
                labels_mat = obj.reshaped_all_labels_vector(ret_ind);
                [C,ia,ic] = unique(labels_mat);
                col = 'rbgyk';
                for ik = C'
                    if ik=="o" %|| ik=="v"
                        continue
                    end
                    labels_vec_temp = labels_mat==ik;
                    labels_vec_temp = add_vec(labels_vec_temp);
                    w = dt*ones(size(labels_vec_temp));
                    c = repmat(col(C==ik),size(labels_vec_temp));
                    arrayfun(plotpatch,labels_vec_temp,w,c)
                end
                
                
                
                
                
                
                
                
            end
            hold off
        end
        function s = saveobj(obj)
            s = audioClipSniff.copy_properties(obj,struct);
        end
        function [ds, X_ds, y_ds, X, y_mat, y_vec] = get_ds(obj, remove_bg_flag, keep_percent)
           if ~exist('remove_bg_flag','var')||isempty(remove_bg_flag)
               remove_bg_flag = false;
           end
           
           if ~exist('keep_percent','var')||isempty(keep_percent)
               keep_percent = 0;
           end
           
           [X, y_mat, y_vec] = obj.extract_features(obj, remove_bg_flag, keep_percent) ;
           batch_size = 1;
           
           [ds, X_ds, y_ds] = obj.get_datastore(X, y_mat,batch_size);
            
        end
        function set.reshape_win_sec(obj,val)
        n_files = numel(obj);
        assert(isnumeric(val))
        assert(floor(val)==val)
%         assert(val<obj.features_len_bin)
        
        for i1 = 1:n_files
            obj(i1).reshape_win_sec = val;
        end
        
        end
        
        function set_reshape_overlap_sec(obj,val)
            n_files = numel(obj);
            assert(isnumeric(val))
            assert(floor(val)==val)
            assert(any(val<vertcat(obj.reshape_win_sec)))
            
            for i1 = 1:n_files
                obj(i1).reshape_overlap_sec = val;
            end
            
        end
        
    end
    methods(Static)
        
        
        function obj = loadobj(S)
            obj = audioClipSniff.copy_properties(S,audioClipSniff);
        end
        function obj_copy = copy_properties(source_obj,output,p)
            
            if ~exist('output','var')||isempty(output)
                output = audioClipSniff;
            end
            
            obj_copy = output;
            if ~exist('p','var')||isempty(p)
                if isstruct(source_obj)
                    p = fieldnames(source_obj);
                else
                    p = properties(source_obj);
                end
                
            end
            
            
            if ~iscell(p)
                p = {p};
            end
            num_att = length(p);
            
            if  isstruct(source_obj)
                
                
                for ip = 1:num_att
                    this_att = p{ip};
                    obj_copy.(this_att) = source_obj.(this_att);
                end
            else
                for ip = 1:num_att
                    this_att = p{ip};
                    
                    mp = findprop(source_obj,this_att);
                    % Dont copy constants and dependant properties
                    if mp.Dependent || mp.Constant || (mp.SetAccess=="private")
                        continue
                    end
                    obj_copy.(this_att) = source_obj.(this_att);
                end
            end
        end
        function [ret_time, ret_ind] =  time2ind(time_vec, time_range)
            
            find_closest_val = @(a,v)  min(abs(a-v));
            [~, ind_start] = find_closest_val(time_vec,time_range(1));
            [~, ind_stop] = find_closest_val(time_vec,time_range(2));
            ret_ind = ind_start:ind_stop;
            ret_time = time_vec(ret_ind);
        end
        function ret = reshape_vector(vec, win, overlap)
            % vec = reshape(vec,[],1);
            if ~exist('overlap','var')||isempty(overlap)
                overlap = 0;
            end
            
            [v_len, n_dum] = size(vec);
            if v_len==1 && n_dum>1
                vec = reshape(vec,[],1);
                [v_len, n_dum] = size(vec);
            end
            hop_len = (win-overlap);
            n_hops = floor((v_len-overlap)/hop_len);
            if isnumeric(vec)
                ret = zeros(win,n_hops,n_dum);
            elseif isstring(vec)
                ret = repmat("o",[win,n_hops,n_dum]);
            end
            ind_s=1;
            for is = 1:n_hops
                
                ind_e = ind_s + win -1;
                
                v = vec(ind_s:ind_e,:);
                ret(:,is,:) = v;
                ind_s = ind_e + 1 - overlap;
            end
        end
        function [X, y_mat, y_vec, time_stamps] = extract_features(data, remove_bg_flag, keep_percent)
            
            if ~exist('keep_percent','var')||isempty(keep_percent)
                keep_percent = 0;
            end
            assert(keep_percent<=1 & keep_percent>=0, 'keep_percent must be between 0-1')
            
            % X = vertcat(data.reshaped_extracted_features);
            X = vertcat(data.reshaped_extracted_features_mult);
            % y = vertcat(data.reshaped_all_labels_vector);
            y = vertcat(data.reshaped_all_labels_vector);
            
            y_mat = vertcat(y.matrix); % for training
            y_vec = vertcat(y.vec); % for splitting the data
%             time_stamps = data.
            % remove segments that contain only bg
            if remove_bg_flag
                ind_remove = find(y_vec=="o");
                num2remove = floor(numel(ind_remove)*(1-keep_percent));
                ind_remove = randsample(ind_remove,num2remove);
                
                X(ind_remove,:,:) = [];
                y_mat(ind_remove,:)= [];
                y_vec(ind_remove,:) = [];
            end
            
            
           
            
            y_mat(y_mat=="a") = "o";
            y_mat(y_mat=="b") = "o";
            y_mat(y_mat=="i") = "o";
            
            y_vec(y_vec=="a") = "o";
            y_vec(y_vec=="b") = "o";
            y_vec(y_vec=="i") = "o";
            
            y_vec = categorical(y_vec);
            % y_mat(y_mat=="i") = "v"; % optional : replace time bins labeled "i" with "v"
            % y_vec(y_vec=="i") = "v";
             
            
        end
        function [ds, X_ds, y_ds] = get_datastore(X, y,batch_size)
            % reshape data so that the input will have dxt and not 1xdxt
            X = permute(X,[2,3,1]);
            X_ds = arrayDatastore(X,"ReadSize",batch_size,"IterationDimension",3);
            
            
            if ~isempty(y)
                y = categorical(y);
                y_ds = arrayDatastore(y,"ReadSize",batch_size,"IterationDimension",1);
                ds = combine(X_ds,y_ds);
            else
                y_ds = [];
                ds = [];
            end
            
            
            
            
            
        end
    end
end

