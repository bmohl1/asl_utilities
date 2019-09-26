function [subj_t1_name, task, t1_dir, routine_out] = preproc_asl(subjs)
  %
  % >>MANUAL AC-PC alignment must be completed for T1's and ASL's prior to
  % entering the pipeline.<<
  %-------------------------------------------------------------
  % Purpose: Preprocess ASL data into standard (MNI) space and run basic QC
  % metrics to flag poor quality scans.
  % Author: Brianne Mohl, PhD
  % Version: 0.1 (10.17)
  % -------------------------------------------------------------
  % Optional input variable is a string that specifies a single subject ID
  %
  % This wrapper script utilizes the SPM batch system, matlabbatch. Each
  % step will be saved as a batch file for review and reproducibility.
  % Additionally, a textfile with the parameters will be generated within the
  % participant's ASL directory.
  %
  % Basic workflow is: Coregistering the ASL image to the participant's
  % native T1 image, segmenting the T1 (and normalizing to MNI space),
  % masking the coregistered ASL with the native GM segmentation, applying
  % the non-linear deformation matrix to bring the masked, GM ASL to MNI
  % space. More complete steps are below.

  % ASL preprocessing
  % Convert dicom to NIFTI like normal -
  % dcm_convert_acetate.sh is tailored to create: a) the directory
  % structure from scratch b) convert the dicoms c) read the names of the
  % converted files and copy the T1, ASL (and fMRI) to the newly created
  % folders - All the steps have checks, so that processed files are not
  % re-processed or overwritten Which dicoms get pulled over? - There are
  % two smaller files that are preprocessed files from Siemens for the 2D, 3D
  % qualitativie, and 3D WIP. - In each case, the script looks for the
  % smaller files and then uses the last instance. This may need to change in
  % the future, but there are no other fields in the converted dicom naming
  % scheme to indicate whether the image is perfusion or rCBF.

  % Processing (this script)
  % 1) Since the rCBF map is a single volume, no realignment is completed.
  % This is different from other publications, which must somehow calculate
  % the CBF from multiple maps.
  % 2) The aligned rCBF map (a) is copied purely for redundancy at the
  % time-consuming, AC-PC aligned step.
  % 3) The rCBF map is coregistered to the T1 in native space.

  % -------------------- Side bar about T1 processing-----------------------
  % 4a) The 6 tissue prior maps from SPM (3x3x3 resolution) are used as the
  % basis for the segmentation of GM, WM, CSF, dura, skull, and
  % background/extraskullar elements
  % 4b) This process of using the priors requires non-linear coregistration
  % to MNI space. {The y_file is the forward deformation warp field that takes
  % the subject to standard space.}
  % 4c) The output y_file and GM segmentation are the basis of pulling the
  % voxels in the ASL scan with the highest probability of being white matter
  % forward into the common brain space.

  % Returning to the ASL pipeline
  % 5) The native space GM map is binarized and multiplied to the coregistered (r)
  % ASL. (This step results in background noise from the ASL scan being
  % pulled forward in the analysisâ¦ and can show up at the 2nd level
  % contrast.)
  % 6) The GM (g), coregistered (r) ASL map is transformed into
  % MNI space with the y_file, warp field.
  % 7) The warped (w), GM (g), coregistered (r) ASL map are smoothed with
  % a 6x6x6 kernel to improve SNR.
  % 8) The swgr files can be compared in a 2nd level analysis.

  %% Begin setting processing variables
  tool_dir = fileparts(fileparts(which('preproc_asl')));
  addpath([tool_dir filesep 'general_utilities']);

  [spm_home, template_home] = update_script_paths(tool_dir);
  %Wherever brain masks or other ROIs would be stored
  qc_mask_dir = fullfile(tool_dir, 'templates' , 'asl_qc');

  %Can change to name the directories (e.g., 'asl_3d')
  tasks = {'asl' 'pcaslLow' 'pcaslHigh'};
  %tasks = {'asl'};

  global special_templates subj_t1_dir subj_t1_file ignore_preproc;

  %%
  [special_templates, save_qc, gm_masking,  ignore_preproc,  redo_segment, cancel] = preproc_inputvars_asl_GUI; %allows for non-scripting users to alter the settings easily.
  settings = {};
  settings.gm = gm_masking;
  close(gcf);
  %    if eq(unwarping,1)
  %        prefix = 'u'; % Can set the letters that are expected prior to standard naming scheme on the data (e.g., 'aru' in 'aruPerson1_task1_scanDate.nii')
  %    end
  %%
  if eq(gm_masking,1)
    normedName = 'wg';
  else
    normedName = 'w';
  end

  if eq(cancel,1)
    return
  else
    %% Select the template file(s)
    if eq(special_templates,1)
      global template_file
      disp('Please select a 4D Tissue Probability Map.');
      tempfile = cellstr(spm_select([1,Inf],'image','Select the template to use throughout the analysis','',pwd));
      tempfile = textscan(tempfile{1,1}, '%s', 'Delimiter',',');
      template_file = tempfile{1,1}{1}; %Must have the ",1" removed for accurate handling elsewhere
      settings.template = template_file;
    end
    %% Grab the files
    switch exist ('subjs','var')
      case 1
      [projDir,pth_subjdirs, subjList] = file_selector(subjs);
      otherwise
      [projDir,pth_subjdirs, subjList] = file_selector;
    end

    cd(projDir)

    %% Start setting up the individual's script
    for n = 1:length(pth_subjdirs)
      subj_pth = pth_subjdirs{n};

      %% Redunant??
      if strcmp(fileparts(subj_pth(1,1:end)), filesep);
        [proj_dir subj ~] = fileparts(subj_pth(1,1:end-1));
      else
        [proj_dir subj ~] = fileparts(subj_pth(1,1:end));
      end

      %%
      for t =1:length(tasks)
        task = tasks{t};

        fprintf ('Checking for %s in %s.\n', task, subj)
        cd(projDir);
        input_dir = strcat(proj_dir,filesep,subj);

        if exist ('aslName') && strcmp(task,'asl')
          task = aslName;
        end

        rawDirs = dir(strcat(input_dir,filesep,task,'*'));
        rawDirs = rawDirs([rawDirs.isdir] & ~ismember({rawDirs.name},{'.','..'}));

        for z = 1:length(rawDirs)
          fprintf ('Working with %s in %s\n', subj, rawDirs(z).name);
          raw_dir = strcat(projDir,filesep,subj,filesep,rawDirs(z).name);
          cd(input_dir)

          %% Find the "raw data" to use for searches and processing
          rawFileName = rdir(strcat(raw_dir, filesep, '*',task,'*.*i*')); %Don't use loc_file here, b/c you are literally looking for the raw data file name
          keep = cellfun('isempty',(cellfun(@(x) strfind(x,'.fig'),{rawFileName.name},'UniformOutput', false)));
          rawFileName = rawFileName(keep);
          val=cellfun(@(x) numel(x),{rawFileName.name}); %compare the length of all the nii's
          rawFileName = rawFileName(val==min(val));
          rawFileName = textscan(rawFileName(1).name, '%s','Delimiter',filesep);
          [rawFileName] = char(rawFileName{1,1}{end}); % ensures that the script won't run into dtype issues concatenating later.

          if exist('prefix', 'var')
            loc_file = [prefix,rawFileName];
          else
            loc_file = rawFileName;
          end

          check_if_processed_name = rdir(strcat(raw_dir,filesep,'s',normedName,loc_file));

          %% T1 Coregistration
          [subj_t1_dir subj_t1 t1_ext] = locate_scan_file('t1',subj);
          t1_name = [subj_t1_dir,filesep, subj_t1];
          %% Processing
          if isempty(check_if_processed_name) || eq(ignore_preproc,1)
            clear matlabbatch
            addpath(genpath(spm_home));
            disp('GM segmentation:')
            gmSeg = [subj_t1_dir,filesep,'c1',subj_t1,',1']
            csfSeg = [subj_t1_dir,filesep,'c3',subj_t1,',1'];
            warpFile = [subj_t1_dir,filesep,'y_',subj_t1];
            invWarpFile = [subj_t1_dir,filesep,'iy_',subj_t1];

            %% Segment the T1, if the step has not already been completed.
            if ~exist(warpFile) || eq(redo_segment,1)
              disp('Deconstructing the brain (segmenting). Please wait.')
              %Routine to apply new segmentation to participant's T1
              cd(projDir)
              segmentation_spm12(subj, redo_segment);
            end

            %% Locating ASL raw file
            orig_file_name = rdir(strcat(raw_dir,filesep,loc_file)); %checks to see if the file exists

            %if eq(unwarping,1) && isempty(orig_file_name)
            %     rawFiles = rdir(strcat(raw_dir,filesep,rawFileName));
            %     cd(raw_dir)
            %     unwarp ({rawFiles.name}, subj, 0, '12b')  %The unwarp
            %     assumes multiple volumes, in it's current version.
            %     This would need to be corrected, if we want to
            %     develop an unwarping option.
            %     orig_file_name = dir(strcat(raw_dir,filesep,loc_file));
            %end

            if isempty (orig_file_name)
              try
                orig_file_name = rdir(strcat(raw_dir,filesep,'*.nii'));
              catch
                sprintf('Did not find the file for %s',task);
                break
              end
            end

            if length(orig_file_name) > 1
              val=cellfun(@(x) numel(x),{orig_file_name.name});
              orig_file_name = orig_file_name(val==min(val));
            end

            orig_file_name = char(orig_file_name.name);
            asl_input = [orig_file_name, ',1'];

            [task_dir orig_file_name ext] = fileparts(orig_file_name);

            %% Checking coregistration
            gmSegOut = strcat('g',orig_file_name,ext); % GM masked ASL
            gmSegOut_check = dir(strcat(task_dir,filesep,gmSegOut));
            csfSegOut = strcat('c',orig_file_name,ext); % GM masked ASL
            csfSegOut_check = dir(strcat(task_dir,filesep,csfSegOut));
            coregAsl = strcat(task_dir,filesep,'r',orig_file_name,ext);
            img2norm = coregAsl;

            %% Processing ASL through matlabbatch
            clear matlabbatch
            spm_jobman('initcfg');
            jj = 1; %To accommodate the possibility of leaving out GM masking, introduce counter
            disp ('Full processing... commence.')
            %Coreg ASL to T1
            matlabbatch{jj}.spm.spatial.coreg.estwrite.ref = {t1_name};
            matlabbatch{jj}.spm.spatial.coreg.estwrite.source = {asl_input};
            matlabbatch{jj}.spm.spatial.coreg.estwrite.other = {''};
            matlabbatch{jj}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            matlabbatch{jj}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            matlabbatch{jj}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{jj}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            matlabbatch{jj}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            matlabbatch{jj}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{jj}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            matlabbatch{jj}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
            jj = jj+1;

            %GM mask
            if eq(gm_masking,1)
              matlabbatch{jj}.spm.util.imcalc.input = { coregAsl;  gmSeg };
              matlabbatch{jj}.spm.util.imcalc.output = gmSegOut;
              matlabbatch{jj}.spm.util.imcalc.outdir = {task_dir};
              matlabbatch{jj}.spm.util.imcalc.expression = 'i1.*(i2>.2)';
              matlabbatch{jj}.spm.util.imcalc.var = struct('name', {}, 'value', {});
              matlabbatch{jj}.spm.util.imcalc.options.dmtx = 0;
              matlabbatch{jj}.spm.util.imcalc.options.mask = 0;
              matlabbatch{jj}.spm.util.imcalc.options.interp = 1;
              matlabbatch{jj}.spm.util.imcalc.options.dtype = 4;
              jj=jj+1;
              img2norm = fullfile(task_dir,gmSegOut);

              %Apply CSF mask to coregistered data for comparison
              matlabbatch{jj}.spm.util.imcalc.input = { coregAsl ; csfSeg };
              matlabbatch{jj}.spm.util.imcalc.output = csfSegOut;
              matlabbatch{jj}.spm.util.imcalc.outdir = {task_dir};
              matlabbatch{jj}.spm.util.imcalc.expression = 'i1.*(i2>.2)';
              matlabbatch{jj}.spm.util.imcalc.var = struct('name', {}, 'value', {});
              matlabbatch{jj}.spm.util.imcalc.options.dmtx = 0;
              matlabbatch{jj}.spm.util.imcalc.options.mask = -1;
              matlabbatch{jj}.spm.util.imcalc.options.interp = 1;
              matlabbatch{jj}.spm.util.imcalc.options.dtype = 4;
              jj=jj+1;
            end

            %Apply deformations to coreg ASL GM, so that ASL GM is in
            %MNI-space (checked for good registration, but only pulled
            %GM forward to reduce partial voluming for CBF)
            matlabbatch{jj}.spm.util.defs.comp{1}.def = {warpFile};
            matlabbatch{jj}.spm.util.defs.out{1}.pull.fnames = {img2norm};
            matlabbatch{jj}.spm.util.defs.out{1}.pull.savedir.saveusr = {task_dir};
            matlabbatch{jj}.spm.util.defs.out{1}.pull.interp = 4;
            matlabbatch{jj}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{jj}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            jj=jj+1;

            %Smooth the normalized ASL GM
            matlabbatch{jj}.spm.spatial.smooth.data(1) = cfg_dep('Deformations: Warped Images', substruct('.','val', '{}',{jj-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','warped'));
            matlabbatch{jj}.spm.spatial.smooth.fwhm = [8 8 8];
            matlabbatch{jj}.spm.spatial.smooth.dtype = 0;
            matlabbatch{jj}.spm.spatial.smooth.im = 0;
            matlabbatch{jj}.spm.spatial.smooth.prefix = 's';



            %Retain record
            savefile = strcat(input_dir,filesep,task,'_coreg2MNI_',subj)

            save(savefile, 'matlabbatch');
            spm_jobman('run',matlabbatch)
          end

          %% QC
          %B/c the sw file exists, we can run QC stats on the brain.
          checkFile = rdir(strcat(raw_dir,filesep, 's', normedName,loc_file));
          [task_dir proc_file_name ext] = fileparts(checkFile.name);
          cd (task_dir)
          output = strcat(task_dir, filesep, subj, '_', task, '_qc_summary.txt'); %Note that as of 10.24.17 this is a 95% pASL average.

          if eq (save_qc,1)
            disp('Calculating QC stats')
            base_loc = dir(strcat(task_dir,filesep,'s', normedName,'*'));
            baseimg=spm_vol(strcat(task_dir,filesep,base_loc(1).name));
            basei_data = spm_read_vols(baseimg);

            avgs = load(fullfile(tool_dir,'asl_utilities','qc_avgs.mat'));
            avgs_graph = figure;

            avgs_names = fieldnames(avgs);
            for ag = 1:length(avgs_names)
              errorbar(ag, avgs.(avgs_names{ag}),0,'Marker','*','MarkerSize',10,'MarkerEdgeColor','red'); hold on
              % fieldname = avgs_names{ag}
            end

            fileID = fopen(output,'w+');

            segs = {'gr', 'csf'}; %corresponds to gm and csf
            masks = {'q1','q2','q3','q4','q5','q6','bottom','top','noise'};
            for s = 1:length(segs)
              seg_loc = glob(strcat(task_dir,filesep,segs{s},'*')); %the cube masks sample smoothed, normalized data, but that is not necessary for the GM and CSF, b/c the segments are individualized already
              %CSF is not normalized, so don't get fancy
              if ~isempty (seg_loc);
                segimg = spm_vol(seg_loc{1});
                segi_data=spm_read_vols(segimg);
                fprintf(fileID,'%s %6f %3f %3f \n',segs{s},mean2(segi_data(segi_data~=0)),std2(segi_data(segi_data~=0)), max(segi_data(:)));
                if mean2(segi_data(segi_data~=0)) > 900 %80ml/100g/min
                  fid= strcat(task_dir,filesep,'asl_seg_qcError.txt');
                  fclose(fopen(fid,'w'));
                end
                errorbar(s+length(masks), mean2(segi_data(segi_data~=0)), std2(segi_data(segi_data~=0)),'Marker','s','MarkerEdgeColor',[0 0 .75],'Color',[0 0 .75]); hold on;
              end

            end

            for m = 1:length(masks)
              mask_loc = dir(strcat(qc_mask_dir,filesep,masks{m},'*')); % 6 cubic "chunks" and a spherical noise sample outside the skull
              maskimg = spm_vol(strcat(qc_mask_dir,filesep,mask_loc.name));
              if strfind(masks{m},'noise'); %This will change the definition, so noise is the final mask...
                %% Make the file defintions available, in case the processing has already happened and just the stats are being (re)run
                orig_file = rdir(strcat(raw_dir,filesep,'*.nii')); %provides name, folder, date, isdir...
                if length(orig_file) > 1
                  val=cellfun(@(x) numel(x),{orig_file.name});
                  orig_file = orig_file(val==min(val));
                end
                base_loc = dir(orig_file.name);
                [origDir ofn ext] = fileparts(base_loc(1).name);
                %%
                baseimg=spm_vol(base_loc(1).name); % load the original ASL file
                basei_data = spm_read_vols(baseimg); %load the associated data
              end

              if ~eq(baseimg.dim, maskimg.dim)  % The resizing is essential to make sure the matrix sizes match for the "add_data" line
                    clear matlabbatch;
                                      matlabbatch{1}.spm.spatial.coreg.write.ref = {[base_loc(1).folder,filesep,base_loc(1).name,',1']};
                    matlabbatch{1}.spm.spatial.coreg.write.source = {[mask_loc.folder,filesep,mask_loc.name,',1']};
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
                    spm_jobman('run',matlabbatch);
                maskimg = spm_vol(strcat(qc_mask_dir,filesep,'r',mask_loc.name));
              end
              maski_data = spm_read_vols(maskimg);
              add_data = (basei_data.*maski_data);
              clean_data= add_data((add_data~=0) & ~isnan(add_data));
              fprintf(fileID,'%s %6f %3f %3f\n',masks{m},mean2(clean_data),std2(clean_data),max(clean_data(:)));
              if mean2(clean_data) > 250 %25ml/100g/min for mixed tissue types
                fid= strcat(task_dir,filesep,'asl_chunks_qcError.txt');
                fclose(fopen(fid,'w'));
              end
              errorbar((m), mean2(add_data(add_data~=0)),std2(add_data(add_data~=0)), 'Marker','s','MarkerEdgeColor',[0 0 .75],'Color',[0 0 .75]); hold on; %
            end




            fclose(fileID);

            xlim([0 11]);
            xticks(0:1:11);
            xticklabels(['none',masks,segs,'none']);
            xtickangle(45);
            savefig(fullfile(task_dir,[subj,'_qc_graph.fig'])); % for quick-view in MATLAB
            saveas (avgs_graph,fullfile(task_dir,[subj,'_qc_graph.eps']), 'epsc'); %for overview in preview or offline
            close (gcf);
          end


        end

        %    disp('Recording statistics')
        %    dt = datestr(now,'yyyymmdd');
        %    log_file = [projDir filesep dt '_aslCortSummary.txt']
        %    asl_inputa = {[raw_dir filesep aslOuta ] };
        %    roi_input = {[subj_t1_dir filesep 't1_' roi '_' subj '.img']};
        %    subjb = [subj_prefix, 'b'];
        %    log_roi_batch( asl_inputa, roi_input, log_file, 'subject', subj);

        %cmd=['chmod -R 777' ' ' subj_prefix '*']; %needed to make sure that server default setting don't lock out other users from further analysis
        %system(cmd);

      end

    end
    cd(projDir)
    disp ('Process completed')
  end
end
