% VsiSNRMask
% M.Naser
% Provides SNR mask based on imaging SNR

function [completed,varargout] = VsiSNRMask(image_folder,varargin)
%% Example: [completed] = VsiSNRMask('C:\MNaser\Data1\ExtractedData','C:\MNaser\Data1\ExtractedData')
%%
%% INPUTS
%% image_folder - string of FULL path for folder containing extracted image data
%% varargin{1} - (optional) string of FULL path for destintation folder; WILL ALSO save result and supress optional output
%%
%% OUTPUTS
%% completed - integer indicating if action was successfully completed (1) or not (0)
%% varargout{1} - (optional) SNR_mask 4D binary mask 1 for voxels having SNR>thershold at (x,y,z, lambda) and 0 elsewhere
%% varargout{2} - (optional) WidthP 1D array of the PA Width 
%% varargout{3} - (optional) DepthP 1D array of the PA Depth 
%% varargout{4} - (optional) ElevP 1D array of the PA Elevation 
%% varargout{5} - (optional) TGC 1D array showing the time gain control used estimated from the noise data  
%% varargout{6} - (optional) param_SNR structre used for estimating the SNR mask
%% varargout{7} - (optional) SNR_mask_out 4D binary mask 1 for voxels having SNR>thershold at (x,y,z, lambda) and 0 elsewhere outside the ROI
%%
completed=1;

[data_file, num] = latestfile(image_folder,'SNR',0);
if(num==0)
    completed=0;
    return;
else
    load(data_file);
    fileused{1} = data_file;
    [data_file, num]= latestfile(image_folder,'ImageData',0);
    if(num==0)
        completed=0;
        return;
    end
    load(data_file);
    fileused{2} = data_file;
    BM=mask_BM_to_PA(RawdataB,WidthB,DepthB,WidthP,DepthP);
    clear RawdataB;
    [data_file, num]= latestfile(image_folder,'SurfMask',0);
    if(num==0)
        completed=0;
        return;
    end
    load(data_file);
    fileused{3} = data_file;
    [data_file,num]=latestfile(image_folder,'TumorSegment',0); 
    if(num==0)
        MaskT=[]; SegPlanes=[];
    else
        load(data_file,'MaskT','SegPlanes');
    end
    
    SNR(repmat(PA_m,[1,1,1,size(SNR,4)])==0)=0;
    
    screen_pos = set_screen([1 1],50);
    h=figure('position',screen_pos{1,1});
    th=9;
    No=2;
    ROI=[-10+0.5*(min(WidthP)+max(WidthP)) 10+0.5*(min(WidthP)+max(WidthP)) 10];
    tx=6;ty=6; tr=6;
    
    while 1
        SNR_mask= (SNR>th);
        SNR_mask_out=SNR_mask;
        for L=1:5:size(SNR_mask,3)
            box=PA_m(:,:,L); box(:,(WidthP<ROI(1))|(WidthP>ROI(2)))=0; box(DepthP>(min(DepthP(sum(box,2)>0))+ROI(3)),:)=0;
            SNR_mask(:,:,L,:)=SNR_mask(:,:,L,:).*repmat(box,[1 1 1 size(SNR_mask,4)]);
            SNR_mask_out(:,:,L,:)=SNR_mask_out(:,:,L,:).*(repmat(box,[1 1 1 size(SNR_mask,4)])==0);
            so2_mask=(squeeze(sum(SNR_mask(:,:,L,:),4))>=2);
            so2_mask_out=(squeeze(sum(SNR_mask_out(:,:,L,:),4))>=2);
            se=strel('disk',tr); temp=imerode(box,se); S=logical(box-temp);
            Lf=find(SegPlanes==L);
            if(~isempty(Lf))
                se=strel('disk',tr); temp=imerode(MaskT(:,:,Lf,1),se); TS=logical(MaskT(:,:,Lf,1)-temp);
            else
                TS=[];
            end
            figure(h);
            hs=subplot(3,5,1); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L),colormap(gray),[],[],S,so2_mask,[],TS,so2_mask_out);colormap(hs,color_s);colorbar off;title([num2str(L),', Press any key']);
            hs=subplot(3,5,11); imagesc(WidthP,DepthP,RawdataP(:,:,L,1)); title(num2str(Wavelengths(1))); colormap(hs,hot); caxis([0 min(1000,max(max(RawdataP(:,:,L,1))))]);
            S=get_border(PA_m(:,:,L),tx,ty);
            hs=subplot(3,5,6); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L),colormap(gray),[],[],S,[],[],TS);colormap(hs,color_s);colorbar off;
            if(L+1<=size(SNR_mask,3))
                box=PA_m(:,:,L+1); box(:,(WidthP<ROI(1))|(WidthP>ROI(2)))=0; box(DepthP>(min(DepthP(sum(box,2)>0))+ROI(3)),:)=0;
                SNR_mask(:,:,L+1,:)=SNR_mask(:,:,L+1,:).*repmat(box,[1 1 1 size(SNR_mask,4)]);
                SNR_mask_out(:,:,L+1,:)=SNR_mask_out(:,:,L+1,:).*(repmat(box,[1 1 1 size(SNR_mask,4)])==0);
                so2_mask=(squeeze(sum(SNR_mask(:,:,L+1,:),4))>=2);
                so2_mask_out=(squeeze(sum(SNR_mask_out(:,:,L+1,:),4))>=2);
                se=strel('disk',tr); temp=imerode(box,se); S=logical(box-temp);
                Lf=find(SegPlanes==L+1);
                if(~isempty(Lf))
                    se=strel('disk',tr); temp=imerode(MaskT(:,:,Lf,1),se); TS=logical(MaskT(:,:,Lf,1)-temp);
                else
                    TS=[];
                end
                hs=subplot(3,5,2); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L+1),colormap(gray),[],[],S,so2_mask,[],TS,so2_mask_out);colormap(hs,color_s);colorbar off;title(num2str(L+1));
                hs=subplot(3,5,12); imagesc(WidthP,DepthP,RawdataP(:,:,L+1,1)); title(num2str(Wavelengths(1))); colormap(hs,hot); caxis([0 min(1000,max(max(RawdataP(:,:,L,1))))]);
                S=get_border(PA_m(:,:,L+1),tx,ty);
                hs=subplot(3,5,7); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L+1),colormap(gray),[],[],S,[],[],TS);colormap(hs,color_s);colorbar off;
            end
            if(L+2<=size(SNR_mask,3))
                box=PA_m(:,:,L+2); box(:,(WidthP<ROI(1))|(WidthP>ROI(2)))=0; box(DepthP>(min(DepthP(sum(box,2)>0))+ROI(3)),:)=0;
                SNR_mask(:,:,L+2,:)=SNR_mask(:,:,L+2,:).*repmat(box,[1 1 1 size(SNR_mask,4)]);
                SNR_mask_out(:,:,L+2,:)=SNR_mask_out(:,:,L+2,:).*(repmat(box,[1 1 1 size(SNR_mask,4)])==0);
                so2_mask=(squeeze(sum(SNR_mask(:,:,L+2,:),4))>=2);
                so2_mask_out=(squeeze(sum(SNR_mask_out(:,:,L+2,:),4))>=2);
                se=strel('disk',tr); temp=imerode(box,se); S=logical(box-temp);
                Lf=find(SegPlanes==L+2);
                if(~isempty(Lf))
                    se=strel('disk',tr); temp=imerode(MaskT(:,:,Lf,1),se); TS=logical(MaskT(:,:,Lf,1)-temp);
                else
                    TS=[];
                end
                hs=subplot(3,5,3); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L+2),colormap(gray),[],[],S,so2_mask,[],TS,so2_mask_out);colormap(hs,color_s);colorbar off;title(num2str(L+2));
                hs=subplot(3,5,13); imagesc(WidthP,DepthP,RawdataP(:,:,L+2,1)); title(num2str(Wavelengths(1))); colormap(hs,hot); caxis([0 min(1000,max(max(RawdataP(:,:,L,1))))]);
                S=get_border(PA_m(:,:,L+2),tx,ty);
                hs=subplot(3,5,8); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L+2),colormap(gray),[],[],S,[],[],TS);colormap(hs,color_s);colorbar off;
            end
            if(L+3<=size(SNR_mask,3))
                box=PA_m(:,:,L+3); box(:,(WidthP<ROI(1))|(WidthP>ROI(2)))=0; box(DepthP>(min(DepthP(sum(box,2)>0))+ROI(3)),:)=0;
                SNR_mask(:,:,L+3,:)=SNR_mask(:,:,L+3,:).*repmat(box,[1 1 1 size(SNR_mask,4)]);
                SNR_mask_out(:,:,L+3,:)=SNR_mask_out(:,:,L+3,:).*(repmat(box,[1 1 1 size(SNR_mask,4)])==0);
                so2_mask=(squeeze(sum(SNR_mask(:,:,L+3,:),4))>=2);
                so2_mask_out=(squeeze(sum(SNR_mask_out(:,:,L+3,:),4))>=2);
                se=strel('disk',tr); temp=imerode(box,se); S=logical(box-temp);
                Lf=find(SegPlanes==L+3);
                if(~isempty(Lf))
                    se=strel('disk',tr); temp=imerode(MaskT(:,:,Lf,1),se); TS=logical(MaskT(:,:,Lf,1)-temp);
                else
                    TS=[];
                end
                hs=subplot(3,5,4); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L+3),colormap(gray),[],[],S,so2_mask,[],TS,so2_mask_out);colormap(hs,color_s);colorbar off;title(num2str(L+3));
                hs=subplot(3,5,14); imagesc(WidthP,DepthP,RawdataP(:,:,L+3,1)); title(num2str(Wavelengths(1))); colormap(hs,hot); caxis([0 min(1000,max(max(RawdataP(:,:,L,1))))]);
                S=get_border(PA_m(:,:,L+3),tx,ty);
                hs=subplot(3,5,9); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L+3),colormap(gray),[],[],S,[],[],TS);colormap(hs,color_s);colorbar off;
            end
            if(L+4<=size(SNR_mask,3))
                box=PA_m(:,:,L+4); box(:,(WidthP<ROI(1))|(WidthP>ROI(2)))=0; box(DepthP>(min(DepthP(sum(box,2)>0))+ROI(3)),:)=0;
                SNR_mask(:,:,L+4,:)=SNR_mask(:,:,L+4,:).*repmat(box,[1 1 1 size(SNR_mask,4)]);
                SNR_mask_out(:,:,L+4,:)=SNR_mask_out(:,:,L+4,:).*(repmat(box,[1 1 1 size(SNR_mask,4)])==0);
                so2_mask=(squeeze(sum(SNR_mask(:,:,L+4,:),4))>=2);
                so2_mask_out=(squeeze(sum(SNR_mask_out(:,:,L+4,:),4))>=2);
                se=strel('disk',tr); temp=imerode(box,se); S=logical(box-temp);
                Lf=find(SegPlanes==L+4);
                if(~isempty(Lf))
                    se=strel('disk',tr); temp=imerode(MaskT(:,:,Lf,1),se); TS=logical(MaskT(:,:,Lf,1)-temp);
                else
                    TS=[];
                end
                hs=subplot(3,5,5); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L+4),colormap(gray),[],[],S,so2_mask,[],TS,so2_mask_out);colormap(hs,color_s);colorbar off;title(num2str(L+4));
                hs=subplot(3,5,15); imagesc(WidthP,DepthP,RawdataP(:,:,L+4,1)); title(num2str(Wavelengths(1))); colormap(hs,hot); caxis([0 min(1000,max(max(RawdataP(:,:,L,1))))]);
                S=get_border(PA_m(:,:,L+4),tx,ty);
                hs=subplot(3,5,10); [overlay,color_s] = Overlay_MN(WidthP,DepthP,BM(:,:,L+4),colormap(gray),[],[],S,[],[],TS);colormap(hs,color_s);colorbar off;
            end
            pause;
        end
        ButtonName = questdlg('Does the mask look ok?','SNR mask Question','Yes','No','Yes');
        switch ButtonName
            case 'Yes'
                param_SNR.ROI=ROI;
                param_SNR.th=th;
                param_SNR.No=2;
                
                SNR_mask=logical(SNR_mask);
                SNR_mask_out=logical(SNR_mask_out);
                close(h);
                savestamp = clock;
                if (nargin > 1)&&(~isempty(varargin{1}))
                    SNR_file = latestfile(varargin{1},'ROIMask',1);
                    save(SNR_file,'fileused','savestamp','SNR_mask','WidthP','DepthP','ElevP','TGC','param_SNR','SNR_mask_out');
                else
                    varargout{1}=SNR_mask;
                    varargout{2}=WidthP;
                    varargout{3}=DepthP;
                    varargout{4}=ElevP;
                    varargout{5}=TGC;
                    varargout{6}=param_SNR;
                    varargout{7}=SNR_mask_out;
            end
                break;
            case 'No' 
                if(isempty(ROI))
                    ROI=[-10+0.5*(min(WidthP)+max(WidthP)) 10+0.5*(min(WidthP)+max(WidthP)) 10];
                end
                [ROI(1),ROI(2),ROI(3),th,No]=mask_data(ROI(1),ROI(2),ROI(3),th, No);
        end
    end
end   
%%%%%%/////////////////////////////////////////////////////////////////////////////
function M = get_border(maskV,tx,ty)
M = imfill(maskV,'holes');
mask = M;
temp = zeros(size(M));
for LL = 1:ceil(ty)
    temp(LL+1:end,:,:) = temp(LL+1:end,:,:)+M(1:end-LL,:,:);
    temp(1:end-LL,:,:) = temp(1:end-LL,:,:)+M(LL+1:end,:,:);
end
for LL = 1:ceil(tx)
    temp(:,LL+1:end,: )= temp(:,LL+1:end,:)+M(:,1:end-LL,:);
    temp(:,1:end-LL,:) = temp(:,1:end-LL,:)+M(:,LL+1:end,:);
end
M = M+temp;
M(M>0) = 1;
M = logical(M-mask);
%%%%%///////////////////////////////////////////////////////////////////////
function [wmin,wmax,dp,th,No] = mask_data(wmin,wmax,dp,th,No)
prompt = {'Lateral min:','Lateral max:','Depth:','SNR threshold','Waves No'};
name = 'Inputs for mask';
numlines = 1;
defaultanswer = {num2str(wmin),num2str(wmax),num2str(dp),num2str(th),num2str(No)};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
pointsource_data = inputdlg(prompt,name,numlines,defaultanswer,options);
wmin = str2double(pointsource_data(1));
wmax = str2double(pointsource_data(2));
dp = str2double(pointsource_data(3));
th = str2double(pointsource_data(4));
No = str2double(pointsource_data(5));
% %%%%/////////////////////////////////////////////////////////////////////////////