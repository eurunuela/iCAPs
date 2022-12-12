%% This function creates a binary mask containing only within-brain voxels
% that are retained for total activation. This is determined on the basis
% of the combination of the three segmentation image outputs: we keep
% voxels within GM, WM and CSF (but not outside of the brain)
%
% Inputs:
% - pMap is a 3D volume with the probabilistic gray matter information
% - param contains TA relavant parameters; here, we access the field 'T_gm' 
% (the threshold (between 0 and 1) past which we consider that we are
% in gray matter)
%
% Outputs:
% - mask is a 1D logical vector with the elements to retain
function [mask,mask_3D] = CreateTAMask(param,fid)

%     try
        if isfield(param,'is_morpho') && param.is_morpho
            
            % We take the 3D data and threshold it
            mask_3D = param.GM_map;
            mask_3D(mask_3D <= param.T_gm) = 0;
            mask_3D(mask_3D > param.T_gm) = 1;
            
            % We then perform closure followed by opening, to 'fill the
            % holes'
            mask_3D = imopen(imclose(logical(mask_3D),ones(param.n_morpho_voxels,...
                param.n_morpho_voxels,param.n_morpho_voxels)),...
                ones(param.n_morpho_voxels,param.n_morpho_voxels,...
                param.n_morpho_voxels));
            
            % mask is our 1Ded data
            mask_3D=logical(mask_3D);
            mask = mask_3D(:);

            WriteInformation(fid,['Created mask with a threshold of ',num2str(param.T_gm),' and opening/closure...']);
        else
            mask_3D = param.GM_map;
            mask_3D(mask_3D <= param.T_gm) = 0;
            mask_3D(mask_3D > param.T_gm) = 1;
            mask_3D=logical(mask_3D);
            mask = mask_3D(:);

            WriteInformation(fid,['Created mask with a threshold of ',num2str(param.T_gm),'...']);
        end
%     catch
%         mask = param.GM_map(:);
%         mask(mask <= param.T_gm) = 0;
%         mask(mask > param.T_gm) = 1;
%         mask = logical(mask);
% 
%         WriteInformation(fid,['Created mask with a threshold of ',num2str(param.T_gm),'...']);
%     end

    
end