% MSOTUpdateOP
% M.Naser
% To get the tissue optical properties updated into the FEM mesh
%%
% Example: mesh = MSOTUpdateOP(mesh, a, b, chrom, s_waves, waves)
%%
% INPUTS
% FEM mesh
% a, b - the reduced scattering parametes (i.e., mus = a * (wave (nm)/500)^-b).
% chrom -  chrom concentrations (N X M) - N is the number of nodes in the FEM mesh and M is the number of chrom (i.e., HbO2, Hb, ICG, etc.).
% s_waves - Exteinction coeff of chrom (NW X M) and NW is the number of wavlengths (i.e. length(waves)).
% waves - 1D array containing the wavelengths in nm (NW X 1)
%%
% OUTPUTS
% mesh with updated field mesh.mus, mesh.mua, and mesh.kappa with sizes (N X NW)
%% 
function mesh = MSOTUpdateOP(mesh, a, b, chrom, s_waves, waves)
%%%***********************************************
Nw = length(waves);
mus = zeros(size(chrom, 1), Nw);
mua = zeros(size(chrom, 1), Nw);
waves = waves/500;

for L = 1:Nw
    mus(:,L) = a.*(waves(L).^(-b));
    mua(:,L) = (s_waves(L,:)*(chrom'))';
end
mesh.mus = mus;
mesh.mua = mua;
mesh.kappa = 1./(3*(mesh.mus+mesh.mua));
%%%%///////////////////////////////////////////////////////////////////////