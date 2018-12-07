% VsiSurfaceSegment
% M.Naser
% Segments surface of VisualSonic B-mode rat image (tissue surfcae extraction)

function [completed,varargout] = VsiSurfaceSegment(image_folder,varargin)
%% Example: [completed] = VsiSurfaceSegment('C:\MNaser\Data1\ExtractedData','C:\MNaser\Data1\ExtractedData')
%%
%% INPUTS
%% image_folder - string of FULL path for folder containing extracted image data
%% varargin - (optional) string of FULL path for destintation folder; WILL ALSO save result and supress optional output
%%
%% OUTPUTS
%% completed - integer indicating if action was successfully completed (1) or not (0)
%% varargout{1} - (optional) PA_m 3D PA binary mask 1 for voxels inside tissue volume and zeros elsewhere
%% varargout{2} - (optional) WidthP 1D array of the PA Width 
%% varargout{3} - (optional) DepthP 1D array of the PA Depth 
%% varargout{4} - (optional) ElevP 1D array of the PA Elevation 
%% varargout{5} - (optional) param the segementation parametes used by for the filtering and surface extraction 

completed = 1;

%%%Hardcoded Parameters****************************************************
param.dp=10; 
param.wp=10;
param.ep=1;
param.th=1; 
param.yc_min=-inf; 
param.yc_max=inf; 
Nr=4; 
Nc=4; 
line_th=20;
%%%%Load ImageData*********************************************************

[image_file,savenum] = latestfile(image_folder,'ImageData');
if savenum == 0
    completed = 0;
    return;
else
    load(image_file)
end
%%%%***********************************************************************
%% Generate binary mask for the PA images using matched-size B-Mode images
index_f=1:size(RawdataB,3);
[BM_s,PA_s,PA_m,param] = segm_Bmode_PA_mask(RawdataB,WidthB,DepthB,ElevB,WidthP,DepthP,param,[],Nr,Nc,index_f,line_th);

filesused = image_file;
savestamp = clock;

if (nargin == 2)&&(~isempty(varargin{1}))
    [savefile] = latestfile(varargin{1},'SurfMask',1);
    fnc = ['save(''',savefile,''',''filesused'',''savestamp'',''PA_m'',''WidthP'',''DepthP'',''ElevP'',''param'')'];
    eval(fnc);
else
    varargout{1} = PA_m;
    varargout{2} = WidthP;
    varargout{3} = DepthP;
    varargout{4} = ElevP;
    varargout{5} = param;
end
%%%%/////////////////////////////////////////////////////////////////////////
