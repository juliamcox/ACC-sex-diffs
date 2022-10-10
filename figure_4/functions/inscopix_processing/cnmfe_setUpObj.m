function neuron = cnmfe_setUpObj

%Creates sources2D object and allows manual updating of parameters:
%gSig: gaussian kernel for filtering the data (pixels)
%gSiz: neuron diameter (pixels)
%min_corr: min correlation for initialization
%min_pnr: min signal to noise
%bd: removal of border?
%zCrop,rCrop and cCrop if data needs to be trimmed

% load tiff stack as matfile
    data = matfile('Y.mat');
    % extract the size of the imaging data and the sampling rate
    Ysiz = data.Ysiz;
    Fs = data.Fs;
    
    %% Starting parameters
    
    % if Y.mat already contains parameters continue
    variableInfo = who('-file', 'Y.mat');
    if ismember('gSig', variableInfo)
        min_corr = data.min_corr;
        min_pnr  = data.min_pnr;
        bd       = data.bd;
        min_pixel= data.min_pixel;
        rCrop = data.rCrop;
        zCrop = data.zCrop;
        cCrop = data.cCrop;
    else %otherwise load default parameters
        gSig = 3;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
        gSiz = 13;          % pixel, neuron diameter
        min_corr = 0.8;     % minimum local correlation for a seeding pixel
        min_pnr = 15;       % minimum peak-to-noise ratio for a seeding pixel
        min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
        bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
        rCrop = [1; Ysiz(1)];
        cCrop = [1; Ysiz(2)];
        zCrop = [1; Ysiz(3)];
    end
    
    if isempty(rCrop)
        keyboard
    end
    d1       = numel(rCrop(1):rCrop(2));   %height
    d2       = numel(cCrop(1):cCrop(2));   %width
    numFrame = numel(zCrop(1):zCrop(2));   %total number of frames
    
    Ysiz = [numel(rCrop(1):rCrop(2)) numel(cCrop(1):cCrop(2)) numel(zCrop(1):zCrop(2))];
   
    
     %% create Source2D class object for storing results and parameters
    
    neuron = Sources2D();
    nam = 'Y.mat';        
    nam = neuron.select_data(nam); 

    
    %% Other parameters (from demo_batch_1p.m)
    
    % -------------------------    COMPUTATION    -------------------------  %
    %     pars_envs = struct('memory_size_to_use', 8, ...   % GB, memory space you allow to use in MATLAB
    %         'memory_size_per_patch', 0.5, ...   % GB, space for loading data within one patch
    %         'patch_dims', [64, 64],...  %GB, patch size
    %         'batch_frames', 1000);           % number of frames per batch
    % -------------------------      SPATIAL      -------------------------  %
    gSig = 3;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
    gSiz = 13;          % pixel, neuron diameter
    
    ssub = 1;           % spatial downsampling factor
    with_dendrites = true;   % with dendrites or not
    if with_dendrites
        % determine the search locations by dilating the current neuron shapes
        updateA_search_method = 'dilate';
        updateA_bSiz = 5;
        updateA_dist = neuron.options.dist;
    else
        % determine the search locations by selecting a round area
        updateA_search_method = 'ellipse'; %#ok<UNRCH>
        updateA_dist = 5;
        updateA_bSiz = neuron.options.dist;
    end
    spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
    spatial_algorithm = 'hals';
    
    % -------------------------      TEMPORAL     -------------------------  %
    tsub = 1;           % temporal downsampling factor
    deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
        'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
        'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
        'optimize_pars', true, ...  % optimize AR coefficients
        'optimize_b', true, ...% optimize the baseline);
        'max_tau', 100);    % maximum decay time (unit: frame);
    
    nk = 3;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
    % when changed, try some integers smaller than total_frame/(Fs*30)
    detrend_method = 'spline';  % compute the local minimum as an estimation of trend.
    
    % -------------------------     BACKGROUND    -------------------------  %
    bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
    nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
    bg_neuron_factor = 1.4;
    ring_radius = round(bg_neuron_factor * gSiz);  % when the ring model used, it is the radius of the ring used in the background model.
    %otherwise, it's just the width of the overlapping area
    num_neighbors = 50; % number of neighbors for each neuron
    
    % -------------------------      MERGING      -------------------------  %
    show_merge = false;  % if true, manually verify the merging step
    merge_thr = 0.65;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
    method_dist = 'max';   % method for computing neuron distances {'mean', 'max'}
    dmin = 5;       % minimum distances between two neurons. it is used together with merge_thr
    dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only.
    merge_thr_spatial = [0.8, 0.4, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)
    
%     % -------------------------  INITIALIZATION   -------------------------  %
    K = [];             % maximum number of neurons per patch. when K=[], take as many as possible.
%     min_corr = 0.8;     % minimum local correlation for a seeding pixel
%     min_pnr = 8;       % minimum peak-to-noise ratio for a seeding pixel
%     min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
%     bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
    frame_range = zCrop;   % when [], uses all frames
    save_initialization = false;    % save the initialization procedure as a video.
    use_parallel = true;    % use parallel computation for parallel computing
    show_init = true;   % show initialization results
    choose_params = true; % manually choose parameters
    center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
    patch_par = [1 1];
    % set the value as false when the background fluctuation is small (2p)
    
    % -------------------------  Residual   -------------------------  %
    min_corr_res = 0.7;
    min_pnr_res = 6;
    seed_method_res = 'auto';  % method for initializing neurons from the residual
    update_sn = true;
    
    % ----------------------  WITH MANUAL INTERVENTION  --------------------  %
    with_manual_intervention = false;
    
    % -------------------------  FINAL RESULTS   -------------------------  %
    save_demixed = true;    % save the demixed file or not
    kt = 3;                 % frame intervals
    
    % -------------------------    UPDATE ALL    -------------------------  %
    neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
        'gSiz', gSiz, ...
        'ring_radius', ring_radius, ...
        'ssub', ssub, ...
        'search_method', updateA_search_method, ...
        'bSiz', updateA_bSiz, ...
        'dist', updateA_bSiz, ...
        'spatial_constraints', spatial_constraints, ...
        'spatial_algorithm', spatial_algorithm, ...
        'tsub', tsub, ...                       % -------- temporal --------
        'deconv_options', deconv_options, ...
        'nk', nk, ...
        'detrend_method', detrend_method, ...
        'background_model', bg_model, ...       % -------- background --------
        'nb', nb, ...
        'ring_radius', ring_radius, ...
        'num_neighbors', num_neighbors, ...
        'merge_thr', merge_thr, ...             % -------- merging ---------
        'dmin', dmin, ...
        'method_dist', method_dist, ...
        'min_corr', min_corr, ...               % ----- initialization -----
        'min_pnr', min_pnr, ...
        'min_pixel', min_pixel, ...
        'bd', bd, ...
        'center_psf', center_psf,...
        'd1',d1,...
        'd2',d2);
    neuron.Fs = Fs;
    
    %% create Source2D class object for storing results and parameters
    
    Y = data.Y(rCrop(1):rCrop(2),cCrop(1):cCrop(2),1:1000);
    Y = neuron.reshape(Y, 1); 
    T = 1000;
    cnmfe_show_corr_pnr % plot correlation and pnr images 
    f = gca;
    f.Visible = 'on';
    
    x = 1;
    fprintf('min corr: %f \n', min_corr);
    fprintf('min pnr: %f \n', min_pnr);
    fprintf('min pixel %f \n', min_pixel);
    g = input('Change parameters? (y/n)', 's');
    
    if strcmpi(g, 'y')
        min_corr = str2num(cell2mat(inputdlg(sprintf('min corr (current %f)', min_corr))));
        min_pnr = str2num(cell2mat(inputdlg(sprintf('min pnr (current %f', min_pnr))));
        min_pixel = str2num(cell2mat(inputdlg(sprintf('min pixel (current %f)', min_pixel))));
        rCrop = inputdlg({sprintf('row range (current %f:%f) \n Start:', rCrop(1), rCrop(2)); 'End:'});
        cCrop = inputdlg({sprintf('column range (current %f:%f) \n Start:', cCrop(1), cCrop(2)); 'End:'});
        
        rCrop = (cellfun(@(x) str2num(x), rCrop));
        cCrop = (cellfun(@(x) str2num(x), cCrop));
        d1       = numel(rCrop(1):rCrop(2));   %height
        d2       = numel(cCrop(1):cCrop(2));   %width
        
        % update parameters in neuron object
        neuron.updateParams('min_corr',min_corr, 'min_pnr', min_pnr, 'min_pixel', min_pixel,'d1',d1,'d2',d2);      
    else
        x = 0;
    end
    
    g = input('Test initialization? (y/n)', 's');
    if strcmpi(g, 'y')
        Y = data.Y(rCrop(1):rCrop(2),cCrop(1):cCrop(2),1:1000);
        d1 = size(Y,1);
        d2 = size(Y,2);
        
        [center, Cn, pnr] = neuron.initComponents_endoscope(Y, 100, patch_par, false, true);
        x=1;
    else
        x = 0;
    end
    
    while x == 1
        
        g = input('Change parameters? (y/n)', 's');
        
        if strcmpi(g, 'y')
            min_corr = str2num(cell2mat(inputdlg(sprintf('min corr (current %f)', min_corr))));
            min_pnr = str2num(cell2mat(inputdlg(sprintf('min pnr (current %f', min_pnr))));
            min_pixel = str2num(cell2mat(inputdlg(sprintf('min pixel (current %f)', min_pixel))));
            rCrop = inputdlg({sprintf('row range (current %f:%f) \n Start:', rCrop(1), rCrop(2)); 'End:'});
            cCrop = inputdlg({sprintf('column range (current %f:%f) \n Start:', cCrop(1), cCrop(2)); 'End:'});
            
            rCrop = (cellfun(@(x) str2num(x), rCrop));
            cCrop = (cellfun(@(x) str2num(x), cCrop));
            
            % update parameters in neuron object
            neuron.updateParams('min_corr',min_corr, 'min_pnr', min_pnr, 'min_pixel', min_pixel,'d1',d1,'d2',d2);
        else
            x = 0; 
        end
        
        g = input('Test initialization? (y/n)', 's');
        
        if strcmpi(g, 'y')
            Y = data.Y(rCrop(1):rCrop(2),cCrop(1):cCrop(2),1:1000);
            d1 = size(Y,1);
            d2 = size(Y,2);
            
            [center, Cn, pnr] = neuron.initComponents_endoscope(Y, 100, patch_par, false, true);
            
        end
    end
    % save sources2D object
    save('neuron.mat', 'neuron')
    % update Y.mat
    data.Properties.Writable = true; 
    data.gSig  = gSig;
    data.gSize = gSiz;
    data.min_corr = min_corr;
    data.min_pnr = min_pnr;
    data.bd = bd;
    data.min_pixel = min_pixel;
    data.rCrop = rCrop;
    data.cCrop = cCrop;
    data.zCrop = zCrop;
    
    
    