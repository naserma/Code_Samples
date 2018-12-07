% VsiSNRMask
% M.Naser
% Provides SNR map based on provided PA images and matched noise image

function [completed,varargout]= VsiSNR(image_folder,noise_folder,varargin)
%% Example: [completed] = VsiSNR('C:\MNaser\Data1\ExtractedData','C:\MNaser\Data1\RawNoise','C:\MNaser\Data1\ExtractedData')
%%
%% INPUTS
%% image_folder - string of FULL path for folder containing extracted image data
%% noise_folder - string of FULL path for folder containing raw noise data to be extracted
%% varargin{1} - (optional) string of FULL path for destintation folder; WILL ALSO save result and supress optional output
%%
%% OUTPUTS
%% completed - integer indicating if action was successfully completed (1) or not (0)
%% varargout{1} - (optional) SNR_mask 4D matrix for the SNR in dB (x,y,z, lambda)
%% varargout{2} - (optional) WidthP 1D array of the PA Width 
%% varargout{3} - (optional) DepthP 1D array of the PA Depth 
%% varargout{4} - (optional) ElevP 1D array of the PA Elevation 
%% varargout{5} - (optional) TGC 1D array showing the time gain control used estimated from the noise data  
%% varargout{6} - (optional) Noise 2D array showing the Noise mean values
%% varargout{7} - (optional) sig 2D array showing the Noise std values  
%%
completed=1;

if(isempty(noise_folder))
    completed=0;
    return;
else
    [data_file, num]= latestfile(image_folder,'ImageData',0);
    if(num==0)
        completed=0;
        return;
    else
        load(data_file,'RawdataP','WidthP','DepthP');
        [comp , ~, ~, ~, ~, ~, ~, ~, N , ~, WidthN, DepthN , ~, ElevP , ~, Evalues] = VsiSaveRaw(noise_folder);
        if(comp)
        E=squeeze(mean(Evalues,1));
        N_size=size(N);
        if(length(N_size)>3)
            for L=1:size(N,3)
                for LL=1:size(N,4)
                    if(E(L,LL)>0)
                        N(:,:,L,LL)=N(:,:,L,LL)*E(L,LL);
                    end
                end
            end
        else
            for L=1:size(N,3)
                if(E(L)>0)
                    N(:,:,L)=N(:,:,L)*E(L);
                end
            end
        end
        N=reshape(N,[size(N,1) size(N,2) size(N,3)*size(N,4)]);
        N=mask_BM_to_PA(N,WidthN,DepthN,WidthP,DepthP);
        Nlog=log(N);
        mu_log=mean(Nlog,3);
        sig_log=std(Nlog,1,3);
        mu=exp(mu_log+((sig_log.^2)/2));
        sig=((exp(sig_log.^2)-1).*exp(mu_log*2+(sig_log.^2))).^0.5;
        Noise=medfilt2(mu,[6 3]); sig=medfilt2(sig,[6 3]);
        [SNR,TGC]=SNR_Mask_Gen(Noise,sig,RawdataP,ones(size(RawdataP,1),size(RawdataP,2),size(RawdataP,3)));
        
        fileused{1} = data_file;
        noise_file = findfile(noise_folder,'raw','pamode');
        fileused{2} = fullfile(noise_folder,noise_file);
        savestamp = clock;
        
        if (nargin > 2)&&(~isempty(varargin{1}))
            SNR_file = latestfile(varargin{1},'SNR',1);
            save(SNR_file,'fileused','savestamp','SNR','WidthP','DepthP','ElevP','TGC','Noise','sig');
        else
            varargout{1}=SNR;
            varargout{2}=WidthP;
            varargout{3}=DepthP;
            varargout{4}=ElevP;
            varargout{5}=TGC;
            varargout{6}=Noise;
            varargout{7}=sig;
        end
        else
            completed=0;
            return;
        end
    end
end
%%%%/////////////////////////////////////////////////////////////////////////////