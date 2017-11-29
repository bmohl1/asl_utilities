function [subj_asl_dir, subj_asl_file, asl_ext] = locate_asl(subj,aslDirName)
dir_check = strfind(pwd, subj(1:4));

if isempty(dir_check);
    cwd = [pwd, filesep, subj];
else
    cwd = pwd;
end
check = textscan(cwd,'%s','Delimiter','/');
ix = strfind(check{1,1},subj); %most specific criterion to find the subject
if isempty(arrayfun(@(x) isempty(x),ix));
    ix = strfind(check{1,1},subj(1:3)); %Catch to find the study prefix in the filepath
end
ix = ~cellfun('isempty',ix); %Find all instances, where the study or subj ID were located
if sum(ix) > 0;
    ix = find(ix==1); %boolean mask of binary ix
    ix = (max(ix)); %finds the last instance in the filepath
    proj_dir = fullfile(filesep,check{1,1}{1:ix-1});
else
    proj_dir = fullfile(filesep,check{1,1}{1:end});
end

if ~exist ('aslDirname','var');
    aslDirName = 'asl'; %default value used in other scripts.
end

x = strcat(proj_dir,filesep,subj,'*',filesep,'*',aslDirName,'*',filesep,'*.*i*');
[asl_file] = rdir(x);
if length(asl_file) < 1 ;
    x = strcat(proj_dir,filesep,subj,'*',filesep,'*',upper(aslDirName),'*',filesep,'*','.nii');
    [asl_file] = rdir(x);
    if length(asl_file) < 1 ;
        disp('Cannot auto-detect ASL, please select the file.')
            [asl_file] = cellstr(spm_select([1,Inf],'file','Select the first persons ASL file','',pwd));
    end
end

tmp = arrayfun(@(x) strfind(x.name,'brain'), asl_file, 'UniformOutput',false);
asl_file = asl_file(find(cellfun(@isempty,tmp)));

val=cellfun(@(x) numel(x),{asl_file.name}); %compare the length of all the nii's

asl_file =  asl_file(val==min(val)); %Find the shortest length, because it will not have been pre-processed
[asl_file] = asl_file.name; %partial name, can't use as output filepath
asl_name = (strcat(asl_file,',1'));
[subj_asl_dir subj_asl_file asl_ext] = fileparts(asl_file);
subj_asl_file = strcat(subj_asl_file,asl_ext);
end