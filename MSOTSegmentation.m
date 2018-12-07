% MSOTSegmentation
% M.Naser
% Surface segmentation
% Read PA image file provided and save the generated binary surf_mask to the provided results folder 
%%
% Example: [completed] = MSOTSegmentation('results_folder','Image file name'); 
%%
% INPUTS
% results_folder - string of FULL path for folder to save the mask (i.e., FluenceCorrection folder)
% Image file name - Reconstruced PA Image file (e.g., ReconMB1.mat)
%%
% OUTPUTS
% completed - integer indicating if action was successfully completed (1) or not (0)
%% 

function [completed] = MSOTSegmentation(results_folder,image_file)
completed = 0;
%%Hardcoded****************************************************************
seg_prog = 'C:\Program Files\ITK-SNAP 3.6\bin\ITK-SNAP'; % Use ITK-SNAP program for manual segmentation
nifti_image = [results_folder,filesep,'Nifti_Image.nii']; % Converted matlab matrices of the PA images into Nifti format to be read by ITK-SNAP
nifti_mask = [results_folder,filesep,'Untitled.nii.gz']; % saved generated binary mask by ITK-SNAP 
[mask_file, num]= latestfile(results_folder,'SurfMask',1); % Save the generated binary mask into matlab file SurfMask
Lfile = 1; filesused{Lfile} = image_file;
%%%%%**********************************************************************
load(image_file); 
wave_index = find(abs(Wavelengths - 800) == min(abs(Wavelengths - 800)),1); % Choose 800 nm wavelength for the segmentation

for L = 1:size(Image,4)
    T = squeeze(Image(:,:,1,L,1,wave_index));
    IMM(:,:,L) = T;
    T = T - min(T(:));
    T = uint8(255*T/max(T(:))); 
    IM(:,:,L) = T;
end
%%%%***********************************************************************
if(num > 1)
    [mask_file, num]= latestfile(results_folder,'SurfMask',0);
    load(mask_file,'surf_mask');
    
    while 1
        MSOTSegmCheck(IMM,surf_mask,3,3);
        ButtonName = questdlg('Segmentation?','Segmentation Question','Ok','Manual','Ok');
        switch ButtonName
            case 'Ok'
                ButtonName2 = questdlg('Saving a New file or skip?','Saving Question','Save','Do not save','Do not save');
                switch ButtonName2
                    case 'Save'
                        [mask_file, num]= latestfile(results_folder,'SurfMask',1);
                        savestamp = clock;
                        save(mask_file,'surf_mask','x','y','ZPositions','filesused','savestamp');
                    case 'Do not save'
                end
                completed = 1;
                close(gcf);
                break;
                
            case 'Manual'
                index = input('Input frames index between brackets like [5 17 18 19], as an example, or press [Enter] for all frames:');
                surf_mask = MSOTSurfSegmManual(IMM,surf_mask,index);
                close(gcf);
        end
    end
    %%%%*******************************************************************
else
    %%% There is no current segmentation -> generating new segmentation
    %%% using ITK-SNAP*****************************************************
    niftiwrite(IM, nifti_image);
    command = ['"',seg_prog,'"',' -g ',nifti_image];
    system(command);
    
    surf_mask = logical(niftiread(nifti_mask));
    while 1
        MSOTSegmCheck(IMM,surf_mask,3,3);
        ButtonName = questdlg('Segmentation?','Segmentation Question','Ok','Manual','Ok');
        switch ButtonName
            case 'Ok'
                savestamp = clock;
                save(mask_file,'surf_mask','x','y','ZPositions','filesused','savestamp');
                completed = 1;
                close(gcf);
                break;
                
            case 'Manual'
                index = input('Input frames index between brackets like [5 17 18 19], as an example, or press [Enter] for all frames:');
                surf_mask = MSOTSurfSegmManual(IMM,surf_mask,index);
                close(gcf);
        end
    end
end
%%%%////////////////////////////////////////////////////////////////////////
