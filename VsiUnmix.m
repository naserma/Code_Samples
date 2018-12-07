% Vevo_Unmix
% M.Naser
% Linear Unmixing, Local Fluence Correction, Surface segmentation, SNR
% Mask, FEM mesh generation
%%
% Example: [completed] = VsiUnmix('Results folder',param,'Image Raw folder','Noise Raw folder','DICOM folder'); % param = []; defauls param.fileds are give in Generate_Unmix_Param function; 
%%
% INPUTS
% results_folder - string of FULL path for folder containing image data, surface mask, SNR mask, Tumor mask, etc. and to save unmixing results 
% varargin{1} - (optional) Structre with different fields to determine the choices used for the whole unmixing process
% varargin{2} - (optional) string of FULL path for folder containing Image raw data - used to extract ImageData file used for unmixing
% varargin{3} - (optional) string of FULL path for folder containing raw noise data to be extracted
% varargin{4} - (optional) string of FULL path for folder containing DICOM file used for Tumor mask generation
%%
% OUTPUTS
% completed - integer indicating if action was successfully completed (1) or not (0)
% varargout{1} - (optional) Structre with different fields to determine the choices used for the whole unmixing process
%%

function [completed,varargout] = VsiUnmix(results_folder,varargin)
completed = 1;
%% Hard coding
file_spectrum = 'Spectrum';
file_flux = 'SurfFlux1';
%%
%%%Set default parameters
param = Generate_Unmix_Param;
%%
%%%Overwrite default fields using optional user input:(e.g., exp param.new = 1);
if (nargin > 1)&&(~isempty(varargin{1}))
    if(isfield(varargin{1},'image'))
        param.image = varargin{1}.image;
    end
    if(isfield(varargin{1},'new'))
        param.new = varargin{1}.new;
    end
    if(isfield(varargin{1},'unmix'))
        param.unmix = varargin{1}.unmix;
    end
    if(isfield(varargin{1},'LFC'))
        param.LFC = varargin{1}.LFC;
    end
    if(isfield(varargin{1},'chrom'))
        param.chrom = varargin{1}.chrom;
    end
    if(isfield(varargin{1},'surf_mask'))
        param.surf_mask = varargin{1}.surf_mask;
    end
    if(isfield(varargin{1},'SNR_mask'))
        param.SNR_mask = varargin{1}.SNR_mask;
    end
    if(isfield(varargin{1},'ROI'))
        param.ROI = varargin{1}.ROI;
    end
    if(isfield(varargin{1},'waves'))
        param.waves = varargin{1}.waves;
    end
    if(isfield(varargin{1},'frames'))
        param.frames = varargin{1}.frames;
    end
    if(isfield(varargin{1},'segm_check'))
        param.segm_check = varargin{1}.segm_check;
    end
    if(isfield(varargin{1},'FEM'))
        param.FEM = varargin{1}.FEM;
    end 
    if(isfield(varargin{1},'No_waves'))
        param.No_waves = varargin{1}.No_waves;
    end
    if(isfield(varargin{1},'SNR_th'))
        param.SNR_th = varargin{1}.SNR_th;
    end
    if(isfield(varargin{1},'Tumor_mask'))
        param.Tumor_mask = varargin{1}.Tumor_mask;
    end
    if(isfield(varargin{1},'a'))
        param.a = varargin{1}.a;
    end
    if(isfield(varargin{1},'b'))
        param.b = varargin{1}.b;
    end
    if(isfield(varargin{1},'beta'))
        param.beta = varargin{1}.beta;
    end
    if(isfield(varargin{1},'maxIT'))
        param.maxIT = varargin{1}.maxIT;
    end
end
varargout{1} = param;
%%%%%Displaying used unmixing parameters***********************************
disp('Unmix_param used:');
disp('==================');
disp(['Unmix_param.image = ',num2str(param.image)]);
disp(['Unmix_param.new = ',num2str(param.new)]);
disp(['Unmix_param.unmix = ',num2str(param.unmix)]);
disp('Unmix_param.chrom = '), disp(param.chrom);
disp(['Unmix_param.waves = ',num2str(param.waves)]);
disp(['Unmix_param.frames = ',num2str(param.frames)]);
disp(['Unmix_param.surf_mask = ',num2str(param.surf_mask)]);
disp(['Unmix_param.segm_check = ',num2str(param.segm_check)]);
disp(['Unmix_param.ROI = ',num2str(param.ROI)]);
disp(['Unmix_param.SNR_mask = ',num2str(param.SNR_mask)]);
if(param.SNR_mask)
    disp(['Unmix_param.No_waves = ',num2str(param.No_waves)]);
    disp(['Unmix_param.SNR_th = ',num2str(param.SNR_th)]);
end
disp(['Unmix_param.Tumor_mask = ',num2str(param.Tumor_mask)]);
disp(['Unmix_param.FEM = ',num2str(param.FEM)]);
if(param.unmix)&&(param.LFC)
    disp(['Unmix_param.a = ',num2str(param.a)]);
    disp(['Unmix_param.b = ',num2str(param.b)]);
    disp(['Unmix_param.beta = ',num2str(param.beta)]);
    disp(['Unmix_param.maxIT = ',num2str(param.maxIT)]);
    disp(['Unmix_param.d = ',num2str(param.d)]);
end
disp('==================');
%%%%% End Displaying used unmixing parameters***********************************
%%
if(isempty(results_folder))
    completed = 0;
    disp('results folder is empty:');
    return;
end
%% ImageData file creation or loading************************************
Lfile=1;
if(param.image)
    if(nargin > 2)&&(~isempty(varargin{2}))
        [comp] = VsiSaveRaw(varargin{2},results_folder);
        if(comp==0)
            completed = 0; disp('ImageData was not created sucessfully:'); return;
        else
            [data_file, num]= latestfile(results_folder,'ImageData',0);
            filesused{Lfile} = data_file; Lfile = Lfile+1;
            load(data_file);
            disp('ImageData was created sucessfully:');
        end
    else
        completed = 0; disp('ImageData was not created sucessfully - Image Raw is empty or not provided:'); return;
    end
else
    [data_file, num]= latestfile(results_folder,'ImageData',0);
    if(num>0)
        filesused{Lfile} = data_file; load(data_file); Lfile = Lfile+1;
        disp('ImageData was loaded sucessfully:');
    elseif(nargin > 2)&&(~isempty(varargin{2}))
        disp('ImageData was not found - A new file will be created:');
        [comp] = VsiSaveRaw(varargin{2},results_folder);
        if(comp==0)
            completed = 0; disp('ImageData was not found and could not be created:'); return;
        else
            [data_file, num]= latestfile(results_folder,'ImageData',0);
            filesused{Lfile} = data_file; Lfile = Lfile+1;
            load(data_file);
            disp('ImageData was created sucessfully:');
        end
    else
        disp('ImageData was not found and could not be created:');
        completed=0; return;
    end   
end
%%%%End ImageData creation or loading**************************************
%% Tumor Mask creation****************************************************
[data_file, num]= latestfile(results_folder,'TumorSegment',0);
if((param.Tumor_mask)&&(param.new))||((param.Tumor_mask)&&(num==0))
    if(nargin > 4)&&(~isempty(varargin{4}))
        [comp] =  VsiDICOMSegment(varargin{4},results_folder,results_folder);
        if(comp)
            disp('Tumor Mask file was created successfully:');
        else
            disp('Tumor Mask file was failed:');
        end
    else
        disp('Tumor Mask file was not created - DICOM folder not provided or empty:');
    end
elseif((param.Tumor_mask)&&(num>0))
    disp('Existing Tumor Mask will be used');
end
%%%%Tumor Mask Mask creation**************************************
%% SurfaceMask
[data_file, num]= latestfile(results_folder,'SurfMask',0);
if((param.surf_mask)&&(param.new))||(((param.surf_mask)||(param.SNR_mask)||(param.FEM)||((param.LFC)&&(param.unmix)))&&(num==0))
    disp('Creating new SurfMask file:');
    [comp] = VsiSurfaceSegment(results_folder,results_folder);
    if(comp)
        Segm_Check(results_folder,4,5);
        while 1
            choice = input('Check Segmentations[1], Manual Segmentation[2], Or Segmentation is Ok and proceed:[enter]?');
            if(isempty(choice))
                break;
            elseif(choice==1)
                Segm_Check(results_folder,4,5);
            elseif(choice==2)
                index = input('Input frames index between brackets like [5 17 18 19] as an example: ');
                Surf_Segm_Manual(results_folder,index,0);
            end
        end
        [data_file, num]= latestfile(results_folder,'SurfMask',0);
        filesused{Lfile} = data_file; Lfile = Lfile+1;
        load(data_file,'PA_m');
        disp('SurfMask file was created sucessfully:');
    else
        disp('SurfMask file was failed - no mask will be applied:');
        PA_m = true(size(RawdataP,1),size(RawdataP,2),size(RawdataP,3));
        param.Surf_mask=0;
    end    
elseif(((param.surf_mask)||(param.SNR_mask)||(param.FEM)||((param.LFC)&&(param.unmix)))&&(num>0))
    filesused{Lfile} = data_file; Lfile = Lfile+1;
    load(data_file,'PA_m');
    disp('SurfMask file was loaded successfully:');
    if(param.segm_check)
        Segm_Check(results_folder,4,5);
    end
else
    PA_m = true(size(RawdataP,1),size(RawdataP,2),size(RawdataP,3));
    param.surf_mask=0;
end
%%%%End Surface Mask creation or loading**************************************
%% SNRMask ***********************************************************
[data_file, num]=latestfile(results_folder,'ROIMask',0);
if((param.SNR_mask)&&(param.new))||(((param.SNR_mask)||((param.LFC)&&(param.unmix)))&&(num==0))
    [data_file, num2]=latestfile(results_folder,'SNR',0);
    if(num2>0)
        disp('SNR file was loaded successfully:');
        [comp] = VsiSNRMask(results_folder,results_folder);
        if(comp>0)
             [data_file, num3]=latestfile(results_folder,'ROIMask',0);
             filesused{Lfile} = data_file; Lfile = Lfile+1;
             load(data_file,'SNR_mask','TGC','param_SNR');
             param.ROI = param_SNR.ROI;
             param.SNR_th = param_SNR.th;
             param.No_waves = param_SNR.No;
             disp('SNR Mask was created successfully:');
        else
            disp('SNR Mask was failed - No SNR Mask will be used');
            param.SNR_mask=0;
            SNR_mask = true(size(RawdataP));
        end
    else
        disp('SNR file was not found - creating a new file:');
        if(nargin > 3)&&(~isempty(varargin{3}))
            [comp] = VsiSNR(results_folder,varargin{3},results_folder);
            if(comp>0)
                [comp] = VsiSNRMask(results_folder,results_folder);
                if(comp>0)
                    [data_file, num4]=latestfile(results_folder,'ROIMask',0);
                    filesused{Lfile} = data_file; Lfile = Lfile+1;
                    load(data_file,'SNR_mask','TGC','param_SNR');
                    param.ROI = param_SNR.ROI;
                    param.SNR_th = param_SNR.th;
                    param.No_waves = param_SNR.No;
                    disp('SNR Mask was created successfully:');
                else
                    disp('SNR Mask was failed - No SNR Mask will be used');
                    param.SNR_mask=0;
                    SNR_mask = true(size(RawdataP));
                end 
            else
                disp('Could not create SNR images - No SNR Mask will be used:');
                param.SNR_mask=0;
                SNR_mask = true(size(RawdataP));
            end
        else
            disp('No ImageData and/or NoiseData folders provided - No SNR Mask will be used:');
            param.SNR_mask=0;
            SNR_mask = true(size(RawdataP));
            TGC = ones(size(RawdataP,1),1);
        end
    end
elseif(((param.SNR_mask)||((param.LFC)&&(param.unmix)))&&(num>0))
    filesused{Lfile} = data_file; Lfile = Lfile+1;
    load(data_file,'SNR_mask','TGC','param_SNR');
    param.ROI = param_SNR.ROI;
    param.SNR_th = param_SNR.th;
    param.No_waves = param_SNR.No;
    param.SNR_mask=1;
    disp('SNR Mask was loaded successfully:');
else
    SNR_mask = true(size(RawdataP));
    param.SNR_mask=0;
    TGC = ones(size(RawdataP,1),1);
end

%%%%End SNR Mask creation or loading**************************************
%% Frames selected for unmixing********************************************
frames = sort(unique(param.frames));
Frames_total = 1:size(RawdataP,3);
index_frames=[];
if(isempty(frames))
    index_frames=1:size(RawdataP,3);
else
    Lw=0;
    for L=1:length(frames)
        ind=find(Frames_total==frames(L));
        if(~isempty(ind))
            Lw=Lw+1; index_frames(Lw) = ind;
        end
    end
    
    if(~isempty(index_frames))
        index_frames=sort(index_frames);
        RawdataP = RawdataP(:,:,index_frames,:);
        PA_m = PA_m(:,:,index_frames);
        SNR_mask = SNR_mask(:,:,index_frames,:);
    else
        disp(['No frames selected in the range: 1 to ',num2str(size(RawdataP,3))]); completed=0; return;
    end
end
%%% End Frames selected for unmixing***************************************
%% FEM Mesh Generation****************************************************
[data_file, num]=latestfile(results_folder,'FEMMesh',0);
if((param.FEM)&&(param.new))||(((param.FEM)||((param.LFC)&&(param.unmix)))&&(num==0))
    disp('Creating new FEM Mesh file:');
    %%%%Tissue-surface irradiation choice**********************************
    if(exist('ParamP','var'))
        if(isfield(ParamP,'TranName'))
            disp(['Transducer name: ', ParamP.TranName]);
            if(strcmp(ParamP.TranName,'LZ-201'))
                choice = 1;
            else
                choice = [];
            end
        end
    else
        disp('Transducer type could not be found:');
        choice = input('Use look-up-table for LZ-201 with narrow-field (~1 cm) fiber-bundle for tissue-surface irradiation[1], Uniform tissue-surface irradiation[enter]?');
    end

    if(isempty(choice))
        disp('Uniform tissue-surface irradiation will be used...');
        file_name_flux = [];
    else
        disp('look-up-table for LZ-201 with narrow-field (~1 cm) fiber-bundle for tissue-surface irradiation will be used...');
        file_name_flux = file_flux;
    end
    %%%%%******************************************************************
    index_frames_local=1:length(index_frames);
    [comp] = FEM_mesh_gen(results_folder, index_frames, [],file_name_flux);
    if(comp)
        disp('FEM Mesh file was created successfully:');
        if(param.LFC==1)
            [data_file, num]=latestfile(results_folder,'FEMMesh',0);
            filesused{Lfile} = data_file; Lfile = Lfile+1;
        end
    else
        disp('FEM Mesh file was failed - Local fluence correction cannot be used: Linear Unmixing will be used');
        param.LFC=0;
    end
elseif(((param.FEM)||((param.LFC)&&(param.unmix)))&&(num>0))
    disp('FEM Mesh file was found:');
    load(data_file,'meshes','mesh_frames','mesh_ROI','index_ROI_pa','A_pa_mesh','A_mesh_pa','flux_file');
    meshes1 = meshes; 
    mesh_frames1 = mesh_frames;
    index_ROI_pa1 = index_ROI_pa;
    A_pa_mesh1 = A_pa_mesh;
    A_mesh_pa1 = A_mesh_pa;
    clear meshes mesh_frames index_ROI_pa A_pa_mesh A_mesh_pa;
    
    ind=sort(unique([mesh_frames1 index_frames]))';
    indx=zeros(length(ind),3);
    indx(:,1)=ind;
    for L=1:size(indx,1)
        i=find(mesh_frames1==indx(L,1));
        j=find(index_frames==indx(L,1));
        if(~isempty(i))
            indx(L,2)=i;
        elseif(~isempty(j))
            indx(L,3)=j;
        end
    end
    indxx=find(indx(:,3)~=0);
    mesh_frames=indx(:,1)';
    if(~isempty(indxx))
        disp('FEM meshes will be created for missing frames:');
        [comp,meshes2,mesh_frames2,mesh_ROI,index_ROI_pa2,A_pa_mesh2,A_mesh_pa2,surf_file,flux_file] = FEM_mesh_gen(results_folder, indx(indxx,1)', [],flux_file,0);
        meshes=cell(1,size(indx,1)); index_ROI_pa=cell(1,size(indx,1)); A_pa_mesh=cell(1,size(indx,1)); A_mesh_pa=cell(1,size(indx,1));
        meshes(indx(:,2)~=0) = meshes1; meshes(indx(:,3)~=0) = meshes2;
        index_ROI_pa(indx(:,2)~=0) = index_ROI_pa1; index_ROI_pa(indx(:,3)~=0) = index_ROI_pa2;
        A_pa_mesh(indx(:,2)~=0) = A_pa_mesh1; A_pa_mesh(indx(:,3)~=0) = A_pa_mesh2;
        A_mesh_pa(indx(:,2)~=0) = A_mesh_pa1; A_mesh_pa(indx(:,3)~=0) = A_mesh_pa2;
        disp('A new file for FEM mesh was created:');
        [data_file, num]=latestfile(results_folder,'FEMMesh',1);
        save(data_file,'meshes','mesh_frames','mesh_ROI','index_ROI_pa','A_pa_mesh','A_mesh_pa','surf_file','flux_file');
    end
    index_frames_local=zeros(size(index_frames));
    for L=1:length(index_frames_local)
        index_frames_local(L)=find(mesh_frames==index_frames(L));
    end
    disp(index_frames_local);
    filesused{Lfile} = data_file; Lfile = Lfile+1;
end
%%%%End FEM Mesh Generation ***********************************************
%% Chromophores Spectra Generation ****************************************
chrom_data = chromophore_spectrum(file_spectrum,param.chrom);
filesused{Lfile} = file_spectrum; Lfile = Lfile+1;
for L=1:length(chrom_data)
    s_wave(:,L) = interp1(chrom_data{L}.wavelengths',chrom_data{L}.mua',Wavelengths);
end
%%%%End  Chromophores Spectra Generation **********************************
%% Wavelengths selected for unmixing***************************************
waves=unique(param.waves);
clear index_waves;
Lw=0;
if(~isempty(waves))
    for L=1:length(waves)
        ind=find(Wavelengths==waves(L));
        if(~isempty(ind))
            Lw=Lw+1; index_waves(Lw) = ind;
        end
    end
    if(Lw==0)
        disp(['Provided Wavelengths [',num2str(waves),'] are not matching with the PA acquisitions Wavelengths [', num2str(Wavelengths),'] - Unmix can not be performed:']);
        completed = 0; return;
    elseif(length(index_waves)<length(chrom_data))
        disp(['Cannot unmix with wavelengths [', num2str(Wavelengths(index_waves)),'] < ',num2str(length(chrom_data)),' :']);
        completed = 0; return;
    elseif(length(index_waves) < length(Wavelengths))
        RawdataP = RawdataP(:,:,:,index_waves);
        Wavelengths = Wavelengths(index_waves);
        SNR_mask = SNR_mask(:,:,:,index_waves);
        s_wave = s_wave(index_waves,:);
    end
else
    index_waves=1:length(Wavelengths);
end
%%% End Wavelengths selected for unmixing**********************************
%% Linear Unmixing*********************************************************
if(param.LFC==0)&&(param.unmix==1)
    disp(['Unmixing using wavelengths: ',num2str(Wavelengths)]);
    if(isempty(param.frames))
        disp('Unmixing all imaging frames');
    else
        disp(['Unmixing frames: ',num2str(index_frames)]);
    end
    
    op_param.SNR_mask = param.SNR_mask;
    op_param.ROI = param.ROI; 
    op_param.ROI_so2 = param.ROI;
    op_param.all_waves = 1;
    No_waves = param.No_waves;
    [chrom,mask_ROI,mask_SNR]=Unmix_Linear(RawdataP,WidthP,DepthP,PA_m,SNR_mask,op_param,s_wave,No_waves,op_param.all_waves);
    savestamp = clock;
    Unmix_param = param;
    [data_file, num]= latestfile(results_folder,'UnmixedUncorr',1);
    save(data_file,'chrom','chrom_data','Wavelengths','WidthP','DepthP','ElevP','mask_ROI','mask_SNR','op_param','Unmix_param','s_wave','index_waves','index_frames','filesused','savestamp');
end
%%% End Linear Unmixing****************************************************
%% Local Fluence Correction************************************************
if(param.LFC==1)&&(param.unmix==1)
    disp(['Unmixing using wavelengths: ',num2str(Wavelengths)]);
    if(isempty(param.frames))
        disp('Unmixing all imaging frames');
    else
        disp(['Unmixing frames: ',num2str(index_frames)]);
    end
    
    if(isempty(param.ROI))
        param.ROI = [-10 10 10];
    end
    op_param.SNR_mask = param.SNR_mask;
    op_param.ROI = param.ROI; 
    op_param.ROI_so2 = param.ROI;
    op_param.ROI_obj = [param.ROI(1)+param.d, param.ROI(2)-param.d, param.ROI(3)];
    op_param.all_waves = 1;
    No_waves = param.No_waves;
    Unmix_param = param;
    [data_file, num]=latestfile(results_folder,'FEMMesh',0);
    if(num==0)
        disp('No FEM file found:'); completed = 0; return;
    else
        load(data_file);
        meshes=meshes(index_frames_local); index_ROI_pa=index_ROI_pa(index_frames_local);
        A_pa_mesh=A_pa_mesh(index_frames_local); A_mesh_pa=A_mesh_pa(index_frames_local);
        mesh_frames=mesh_frames(index_frames_local);
    end
    [chrom,mask_ROI,mask_SNR,phi_final,a_best,b_best] = Unmix_LFC(RawdataP,WidthP,DepthP,meshes,index_ROI_pa,A_pa_mesh,A_mesh_pa,SNR_mask,SNR_mask,...
                                                      TGC,op_param,Wavelengths,s_wave,param.a,param.b,mesh_frames,No_waves,op_param.all_waves,param.beta,param.maxIT);
    savestamp = clock;
    [data_file, num]= latestfile(results_folder,'UnmixedCorr',1);
    save(data_file,'chrom','chrom_data','Wavelengths','WidthP','DepthP','ElevP','mask_ROI','mask_SNR','op_param','Unmix_param','s_wave','index_waves','index_frames','phi_final','a_best','b_best','filesused','savestamp');
end
%%% End Local Fluence Correction*******************************************
%% End of the VsiUnmix***************************************************
%%%%///////////////////////////////////////////////////////////////////////