% MSOTSegmCheck
% M.Naser
% Surface segmentation Check
%%
% Example: MSOTSegmCheck(image,mask,Nr,Nc); 
%%
% INPUTS
% image - PA image (1:N, 1:M, 1:L), 3D image (i.e.,PA(Depth,Width,Elevation))
% mask - 3D binary matrix - the same dimension as PA
% Nr and Nc # of subplots in the y and x direction respectively
%%
% OUTPUTS
%% 
function MSOTSegmCheck(image,mask,Nr,Nc)
%%%%%**********************************************************************
screen_pos = set_screen([1 1],50);
se = strel('disk',4);
Lf = 0; %% Frame counter

for LLL = 1:ceil(size(mask,3)/(Nr*Nc))
    h = figure('position',screen_pos{1,1});
    for L = 1:Nr
        for LL = 1:Nc
            Lf = Lf+1;
            if(Lf <= size(mask,3))
                me = imerode(mask(:,:,Lf), se)>0;
                S = (mask(:,:,Lf) - me) > 0;
                hs = subplot(Nr,Nc,(L-1)*Nc+LL);
                temp = image(:,:,Lf);
                thmin = nanmean(temp(:)) - 2*nanstd(temp(:));
                thmax = nanmean(temp(:)) + 2*nanstd(temp(:));
                BM = min(max(image(:,:,Lf),thmin),thmax);
                [overlay,color_s] = Overlay_MN(1:size(image,2),1:size(image,1),BM,colormap(gray),[],[],[],S,[]);colormap(hs,color_s);colorbar off;
                title(num2str(Lf));
                drawnow;
            else
                break;
            end
       end
    end
   suptitle('Press Any Key to Proceed:');
   pause;
   close(h);
 end 
%%%%///////////////////////////////////////////////////////////////////////