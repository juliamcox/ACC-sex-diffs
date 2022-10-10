function inscopixDffExtract_largeFile(aids, recs, expt, onCluster, neurons_selected)

% Convert, select parameters and run cnmfe for inscopix imaging videos; requires tiff files to be saved in locations {recs} 

% aids is a cell array of animal IDs
% recs is a cell array of recordings 
% expt is same size as aids and indicates cohort (i.e. data location: 'DMS imaging' or 'ACC_DMS_imaging') 

if nargin < 4
    onCluster = 1;
end

autoInit = 1; %use automatic or manual ROI initialization 

%% Add cmmfe files to the path


	
for na = 1:numel(aids)
    for nr = 1:numel(recs{na})
        
        %% Go to folder with data
        
        dataLoc = fullfile(whereAreWe('imaging'), expt{na},  aids{na},  recs{na}{nr});
       
         
        
        %% Has the video been converted?
        if exist(fullfile(dataLoc, 'Y.mat'), 'file') ~= 2
            fprintf('Converting tiff files. . .\n')
            cnmfeDataConvert(dataLoc,recs{na}{nr}); % convert from tif to mat
        end
  
        data = matfile(fullfile(dataLoc,'Y.mat'));
        
        %% Have cnmfe parameters been set? If not, run cnmfe_setUpObj to create sources2D object
        cd(dataLoc)
        variableInfo = who('-file', fullfile(dataLoc,'Y.mat'));
        if ~ismember('min_corr', variableInfo)
            fprintf('Setting up cnmfe. . .\n')
            if onCluster
                return
            else
           
            neuron = cnmfe_setUpObj;
            end
        else
            load(fullfile(dataLoc, 'neuron.mat')); 
        end
        
        
        %% Have neuron locations been selected?
       % neurons_selected=exist(fullfile(dataLoc,'center.mat'));
        
        if onCluster
            %% Initialization
            % -------------------------  Parameters   -------------------------  %
            K = 20;             % maximum number of neurons per patch. when K=[], take as many as possible.
            min_corr = neuron.options.min_corr;     % minimum local correlation for a seeding pixel
            min_pnr = neuron.options.min_pnr;       % minimum peak-to-noise ratio for a seeding pixel
            min_pixel = neuron.options.min_pixel;      % minimum number of nonzero pixels for each neuron
            bd = neuron.options.bd;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
            frame_range = data.zCrop;   % when [], uses all frames
            save_initialization = false;    % save the initialization procedure as a video.
            use_parallel = true;    % use parallel computation for parallel computing
            show_init = false;   % show initialization results
            choose_params = false; % manually choose parameters
            center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)         % set the value as false when the background fluctuation is small (2p)
            debug_on = false;
            save_avi = false;
            patch_par = [4,4]*1; %1;  % divide the optical field into m X n patches and do initialization patch by patch. It can be used when the data is too large
            update_sn=1;
            % -------------------------      MERGING      -------------------------  %
            show_merge = false;  % if true, manually verify the merging step
            merge_thr = 0.65;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
            method_dist = 'max';   % method for computing neuron distances {'mean', 'max'}
            dmin = 5;       % minimum distances between two neurons. it is used together with merge_thr
            dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only.
            merge_thr_spatial = [0.8, 0.4, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)
            
            % -------------------------  Residual   -------------------------  %
            min_corr_res = 0.7;
            min_pnr_res = 6;
            seed_method_res = 'auto';  % method for initializing neurons from the residual
            update_sn = true;
            
            numFrame = data.Ysiz(1,3);
            % -------------------------    COMPUTATION    -------------------------  %
            pars_envs = struct('memory_size_to_use', 16, ...   % GB, memory space you allow to use in MATLAB
                'memory_size_per_patch', 8, ...   % GB, space for loading data within one patch
                'patch_dims', [64, 64],...%GB, patch size
                'batch_frames', size(data,3));  
            
            if ~neurons_selected
                if autoInit
                    neuron.options.seed_method='auto';
                else
                    neuron.options.seed_method = 'manual';
                end
                
                % Check whether data has been cropped, if not, load, crop if necessary and save in data.mat
                
                if ~exist(fullfile(dataLoc, 'data.mat'))
                    Y = data.Y(data.rCrop(1,1):data.rCrop(2,1),data.cCrop(1,1):data.cCrop(2,1),data.zCrop(1,1):data.zCrop(2,1));
                    
                    
                    load info.mat droppedIdx
                    if droppedIdx == 0
                        droppedIdx = [];
                    end
                    
                    if ~isempty(droppedIdx) && numFrame > 5000
                        % Load and crop data
                        Y(:,:,droppedIdx) = [];
                        zCrop = [1 size(Y,3)];
                        numFrame = size(Y,3);
                        fprintf('Removed dropped frames...\n')
                        num2read = numFrame;
                        Ysiz = [numel(data.rCrop(1,1):data.rCrop(2,1)) numel(data.cCrop(1,1):data.cCrop(2,1)) numel(zCrop(1):zCrop(2))];
                    end
                    
                    save('data.mat', 'Y', '-v7.3');
                    
                    clear data.Y;
                else load(fullfile(dataLoc, 'data.mat'))
                    
                end
                
                
                nams = neuron.select_data(fullfile(dataLoc, 'data.mat'));
                
                neuron.getReady(pars_envs);
                
                fprintf('CNMFE initialization. . . \n');
                
                if autoInit
                    show_init = false;
                    [center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
                    neuron.compactSpatial();
                else
                end
                fprintf('CNMFE initialization complete. . .\n');
                
                
                save(fullfile(dataLoc, 'neuron.mat'), 'neuron')
                
                save(fullfile(dataLoc, 'center.mat'), 'center', 'Cn', 'PNR')
            else
                load(fullfile(dataLoc, 'neuron.mat'))
                load(fullfile(dataLoc, 'center.mat'))
            end
            %% estimate the background components
          
            neuron.update_background_parallel(use_parallel);
            neuron_init = neuron.copy();
            
            %%  merge neurons
       
            neuron.merge_neurons_dist_corr(show_merge);
            neuron.merge_high_corr(show_merge, merge_thr_spatial);
           
            %% pick neurons from the residual
            [center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel(10, save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
            %[results, center, Cn, PNR, save_avi] = greedyROI_endoscope(Y, K, options,debug_on, save_avi)
            %% udpate spatial&temporal components, delete false positives and merge neurons
            % update spatial
            if update_sn
                neuron.update_spatial_parallel(use_parallel, true);
                udpate_sn = false;
            else
                neuron.update_spatial_parallel(use_parallel);
            end
            % merge neurons based on correlations
            neuron.merge_high_corr(show_merge, merge_thr_spatial);
            
            for m=1:2
                % update temporal
                neuron.update_temporal_parallel(use_parallel);
                
                % delete bad neurons
                neuron.remove_false_positives();
                
                % merge neurons based on temporal correlation + distances
                neuron.merge_neurons_dist_corr(show_merge);
            end
            %% run more iterations
            neuron.update_background_parallel(use_parallel);
            neuron.update_spatial_parallel(use_parallel);
            neuron.update_temporal_parallel(use_parallel);
            
            K = size(neuron.A,2);
            tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
            neuron.remove_false_positives();
            neuron.merge_neurons_dist_corr(show_merge);
            neuron.merge_high_corr(show_merge, merge_thr_spatial);
            
            if K~=size(neuron.A,2)
                neuron.update_spatial_parallel(use_parallel);
                neuron.update_temporal_parallel(use_parallel);
                neuron.remove_false_positives();
            end
            
            %% save the workspace for future analysis
            neuron.orderROIs('snr');
            save('neuron.mat', 'neuron') 
           % cnmfe_path = neuron.save_workspace();
          %  
           % neuron.save_neurons();
            
         
        end
    end
    
end

end