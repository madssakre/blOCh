clc
clear 
close all
% This script loads the data from Dr. Guerin (Harvard Medical School and 
% the A. A. Martinos Center for Biomedical Imaging ) made publicly available
% at 
%
% https://ptx.martinos.org/index.php/Small_flip-angle_spokes_with_SAR_constraints
%
% (Last checked 2017, April 26)
%
% ... and arranges it properly for the blOCh conventions.
%
% Please download these data yourself.
%
% Dr. Guerin contributed to the inclusion of strict SAR constraints in
% blOCh, and his work deserves proper citation too.
%
% This script as an example script is based on Dr. Guerins own "loading" scripts 
% of his data sets and comes without warranty for correct use.
% 
% 



b0file = ['b0map.mat'];
b1file = ['b1maps.mat'];
lsarfile = ['vop_1gram_2percent.txt'];
gsarfile = ['glob_sar_mat.txt'];

nchannels = 8;


%%
% global SAR matrix
gsar_temp = load(gsarfile,'-ascii');
Qmat = complex(gsar_temp(:,1:2:end),gsar_temp(:,2:2:end));
UnitQmat = 'W/kg/V^2';

save('global_sar.mat','-mat','Qmat','UnitQmat');

%%

% local SAR matrices
file=fopen(lsarfile,'r');
tline=fgets(file); nvop=tline_to_param(tline,18,size(tline,2)-1);
tline=fgets(file); umax=tline_to_param(tline,51,size(tline,2)-1);
Qmat=zeros(nvop+1,nchannels,nchannels);
for n=1:nvop
    if n==1
        tline=fgets(file); tline=fgets(file);
    else
        tline=fgets(file); tline=fgets(file); tline=fgets(file);
    end
    for i=1:nchannels
        tmp=fscanf(file,'%e',2*nchannels);
        Qmat(n,i,:)=( tmp(1:2:2*nchannels-1)+1j*tmp(2:2:2*nchannels) );
    end
   
end
fclose(file);

Qmat = permute(Qmat,[2,3,1]);

UnitQmat = 'W/kg/V^2';

save('local_sar_vops.mat','-mat','Qmat','UnitQmat');

%%
% B1 maps

b1 = load(b1file);

nx=size(b1.b1maps,1);
ny=size(b1.b1maps,2);
nz=size(b1.b1maps,3);

b1maps=conj( permute( b1.b1maps,[1,2,3,5,4] ) );

B1_nom_amp_unit = 'T/V';
B1_nom_amp_val = max(abs(b1maps(:)));
B1_type = 'absolute_valued';

B1map = b1maps./max(abs((b1maps(:))));
B1map = B1map(end:-1:1,:,:,:,:);
R = [nx,ny,nz];Rv = 1;

D = [-100,-100,-2.5;100,100,2.5]*1e-3 - [34.5060,0,0;34.5060,0,0]*1e-3;
Dv = [0;0];
%% 
B0map = load(b0file);

v0map = B0map.b0map;
v0map = v0map(end:-1:1,:);
T1map = ones(size(v0map))*1e16;
T2map = ones(size(v0map))*1e16-1;

save('external.mps','-mat','D','Dv','R','Rv','B1map','B1_nom_amp_unit','B1_nom_amp_val','B1_type','v0map','T1map','T2map')
