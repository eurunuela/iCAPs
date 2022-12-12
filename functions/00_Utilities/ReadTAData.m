%% This function reads the input data for total activation; it supports
% reading from a 4D NIFTI, 3D NIFTI, or IMG/HDR files
%
% Inputs:
% - Path is the path towards where the data lies
% - param is a structure containing relevant TA parameters; here, we will
% require the fields 'Folder_functional' (name of the functional folder; []
% if directly lying in Path), 'TA_func_prefix' (string with the prefix for
% functional data to read), 'Folder_GM' (name of the folder with the
% probabilistic map; [] if directly lying in Path), and 'TA_GM_prefix'
% (string with the prefix of the probabilistic map to read)
%
% Outputs:
% - fData is the functional data in 4D format (X x Y x Z x T)
% - pData is the probabilistic GM data as a 3D volume
% - fHeader is a header with .dim and .mat entries for functional data
% - pHeader is the same type of header for the probabilistic GM data

% 17.01.2017 (DZ):  - removing dependence on subject index (soi input)
%               while including modifications for multiple subjects
%               - adding option of c1 masks that are not in functional
%               resolution: map to functional space and save the resulting
%               map (then it won?t be needed in the CheckTAInputs anymore)


function [fData,pData,fHeader,pHeader] = ReadTAData(Path,sidx,param,fid)

    if ~isfield(param,'Folder_functional')
        param.Folder_functional=[];
    end
    if ~isfield(param,'Folder_GM')
        param.Folder_GM=[];
    end

    % compatibility with input as cell of strings for every subject or for a
    % simple string with the same subfoldeer for all subjects
    if ~isempty(param.Folder_functional) && iscell(param.Folder_functional)
        if length(param.Folder_functional)~=param.n_subjects
            error('param.Folder_functional: wrong number of functional subfolders');
        end
        param.Folder_functional=param.Folder_functional{sidx};
    end
    if ~isempty(param.Folder_GM) && iscell(param.Folder_GM)
        if length(param.Folder_GM)~=param.n_subjects
            error('param.Folder_GM: wrong number of gray matter subfolders');
        end
        param.Folder_GM=param.Folder_GM{sidx};
    end
    if ~isempty(param.TA_func_prefix) && iscell(param.TA_func_prefix)
        if length(param.TA_func_prefix)~=param.n_subjects
            error('param.TA_func_prefix: wrong number of functional prefix names');
        end
        param.TA_func_prefix=param.TA_func_prefix{sidx};
    end
    if ~isempty(param.TA_gm_prefix) && iscell(param.TA_gm_prefix)
        if length(param.TA_gm_prefix)~=param.n_subjects
            error('param.TA_gm_prefix: wrong number of functional subfolders');
        end
        param.TA_gm_prefix=param.TA_gm_prefix{sidx};
    end



    % Checks that the input folder and the folders containing the data
    % exist
    if ~exist(fullfile(Path,param.Folder_functional))
        error('There is no functional folder...');
    end

    if ~exist(fullfile(Path,param.Folder_GM))
        error('There is no GM map folder...');
    end


    %%%%%%%%%%%%%%%%% Reads the functional data %%%%%%%%%%%%%%%%%

    % Path towards the functional data (already checked for existence
    % before)
    FPath = fullfile(Path,param.Folder_functional);

    % Looks for files with a NIFTI end
    tmp_data = cellstr(spm_select('List',FPath,['^' param.TA_func_prefix '.*\.' 'nii' '$']));

    if ~strcmp(tmp_data,'')
        tmp_data = fullfile(FPath,tmp_data);
        % Number of NIFTI files found
        n_files = length(tmp_data);
    else
        n_files = 0;
    end

    switch n_files

        % If there is only one NIFTI file...
        case 1
            WriteInformation(fid,['Reading data for subject ',Path,': ',num2str(n_files),' NIFTI file']);

            % ... if it only contains one volume, then there is a problem
            % (only one provided frame). Else, we read all the data
            Vh = spm_vol(tmp_data);
            n_vol = length(Vh);

            if n_vol == 1
                error('Only one NIFTI volume provided...');
            else
                fHeader = Vh(1);
                for i = 1:n_vol
                    fData(:,:,:,i) = spm_read_vols(Vh(i));
                end
            end

            % If there is no NIFTI file...
        case 0

            WriteInformation(fid,['Reading data for subject ',Path,': no NIFTI file']);

            % ... we look for img/hdr files instead. If there are such
            % files (and an equal amount of img/hdr), we read them
            tmp_hdr0 = cellstr(spm_select('List',FPath,['^' param.TA_func_prefix '.*\.' 'hdr' '$']));
            tmp_img0 = cellstr(spm_select('List',FPath,['^' param.TA_func_prefix '.*\.' 'img' '$']));
            tmp_hdr = fullfile(FPath,tmp_hdr0);
            tmp_img = fullfile(FPath,tmp_img0);

            if ~strcmp(tmp_hdr0,'') & ~strcmp(tmp_img0,'')
                if (length(tmp_hdr) == length(tmp_img) && is_same(tmp_hdr,tmp_img))

                    n_vol = length(tmp_hdr);

                    fHeader = spm_vol(tmp_hdr{1});

                    for i = 1:n_vol
                        fData(:,:,:,i) = spm_read_vols(spm_vol(tmp_hdr{i}));
                    end
                else
                    WriteInformation(fid,['Reading data for subject ',Path,': not the same amount of IMG/HDR files']);
                    error('Wrong balance between hdr and img files or mismatches names...');
                end
            else
                WriteInformation(fid,['Reading data for subject ',Path,': no HDR/IMG files']);
                error('No functional files provided...');
            end

        otherwise
            WriteInformation(fid,['Reading data for subject ',Path,': ',num2str(n_files),' NIFTI files']);

            fHeader = spm_vol(tmp_data{1});

            for i = 1:n_files
                fData(:,:,:,i) = spm_read_vols(spm_vol(tmp_data{i}));
            end
    end

    %%%%%%%%%%%%%%%%% Reads the probabilistic map data %%%%%%%%%%%%%%%%%

    % Path towards the functional data (already checked for existence
    % before)
    PPath = fullfile(Path,param.Folder_GM);

    % Looks for files with a NIFTI end
    tmp_data = cellstr(spm_select('List',PPath,['^' param.TA_gm_prefix '.*\.' 'nii' '$']));

    %%% if the functional GM file starts with f and does not exist, check if
    %%% there is a structural GM file without that prefix
    if isempty(tmp_data) && param.TA_gm_prefix(1)=='f'
        param.TA_gm_prefix=param.TA_gm_prefix(2:end);
        tmp_data = cellstr(spm_select('List',PPath,['^' param.TA_gm_prefix '.*\.' 'nii' '$']));
    end

    if ~strcmp(tmp_data,'')
        tmp_data = fullfile(PPath,tmp_data);
        % Number of NIFTI files found
        n_files = length(tmp_data);
    else
        n_files = 0;
    end

    switch n_files

        % If there is only one NIFTI file...
        case 1
            WriteInformation(fid,['Reading probabilistic GM data for subject ',Path,': ',num2str(n_files),' NIFTI file']);

            pHeader = spm_vol(tmp_data{1});
            pData = spm_read_vols(pHeader);

            %%% Dani: checking if functional c1 file exists, if not, map structural
            %%% GM file to functional space and save the resulting map
            if ~isequal(pHeader.dim,fHeader.dim)
                WriteInformation(fid,'Different data dimensions for GM map and functional data: converting GM map...');
                pData=mapVTV(pData,pHeader,fHeader);
                pHeader.dim=fHeader.dim;
                pHeader.fname=strrep(tmp_data{1},param.TA_gm_prefix,['f' param.TA_gm_prefix]); % saving GM mask in functional resolution with prefix 'f'
                spm_write_vol(pHeader,pData);
            end

            % If there is no NIFTI file...
        case 0
            WriteInformation(fid,['Reading data for subject ',Path,': no NIFTI file']);

            % ... we look for img/hdr files instead. If there are such
            % files (and an equal amount of img/hdr), we read them
            tmp_hdr0 = cellstr(spm_select('List',PPath,['^' param.TA_gm_prefix '.*\.' 'hdr' '$']));
            tmp_img0 = cellstr(spm_select('List',PPath,['^' param.TA_gm_prefix '.*\.' 'img' '$']));
            tmp_hdr = fullfile(PPath,tmp_hdr0);
            tmp_img = fullfile(PPath,tmp_img0);

            if (~strcmp(tmp_hdr0,'') & ~strcmp(tmp_img0,''))
                if (length(tmp_hdr) == length(tmp_img) && length(tmp_hdr) == 1 && is_same(tmp_hdr,tmp_img))
                    pHeader = spm_vol(tmp_hdr{1});
                    pData = spm_read_vols(pHeader);


                else
                    WriteInformation(fid,['Reading data for subject ',Path,': not the same amount of IMG/HDR files']);
                    error('Wrong balance between hdr and img files or mismatches names...');
                end
            else
                WriteInformation(fid,['Reading data for subject ',Path,': no HDR/IMG files']);
                error('No functional files provided...');
            end

        otherwise
            error('More than one file with the chosen GM map prefix: please check your inputs...');
    end
end