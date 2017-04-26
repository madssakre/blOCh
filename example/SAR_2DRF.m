clear
close all
clc
% This is an example script for optimizing a 2DRF pulse with strict local
% and global SAR constraints, as well as peak and average RF power, as
% described in:
% 
%   "Local SAR, global SAR, and power-constrained large-flip-
%   angle pulses with optimal control and virtual observation points."
%   Vinding, M.S., et al.,
%   Magn Reson Med. 2017, 77, 1, 374-384.
%   doi: 10.1002/mrm.26086
%
%  (Please cite this paper)
%
% The example also exploits B1maps, a B0 maps and SAR matrices made publicly available by 
% Dr. Guerin (Harvard Medical School and the A. A. Martinos Center for Biomedical Imaging) at 
%
% https://ptx.martinos.org/index.php/Small_flip-angle_spokes_with_SAR_constraints
% 
% (be sure to cite Dr. Guerin)
%
% To use these data you need to download the data and run the adjacent
% script "Load_external_data.m" which will arrange the data for the blOCh
% conventions.
%
% For further information about this example or anything related to blOCh
% please contact me

%%
% make the spatial grid for optimization
%
M0 = ''; % Default equilibrium magnetization
Md = ''; % Default deisred magnetization described by the following:
TH = [0.02,0.01]; % desired magnetization square box size [m]
c_TH = [-34.5060e-3,0]; % center of box
R = 10; % spatial grid size [#]
Mdflip = [45,0]; % desired flip angle and flip phase [degrees]

L = [0.2,0.2]; % FOX [m]
r0 = [-34.5060e-3,0,0]; % centering of FOX

MPS = './external.mps'; % mat file with B1, B0, T1, and T2 maps (blOCh convention)

spc = blOCh__spc([],'2DAx',M0,Md,MPS,'pTx',8,'Show',1,'L',L,'R',[R,R],'Mdflip',Mdflip,'TH',TH,'c_TH',c_TH,'r0',r0);
% make the spatial grid for simulation
spc2 = blOCh__spc([],'2DAx',M0,Md,MPS,'pTx',8,'Show',1,'L',L,'R',[R,R]*3,'Mdflip',Mdflip,'TH',TH,'c_TH',c_TH,'r0',r0);

%%
% make the temporal dependent trajectory
% 
% in this case 2d inward spiral based on the Glover'1999 approach
%
khr= blOCh__khr([],spc,'2Dspii',[],'Show',1,'Sys','7T','dt',10e-6);
%%

% chose the external optimization function
%
% in this case the quasi-Newton L-BFGS method, using matlabs fmincon as
% optimizer
Method = {'blOCh__QN_LBFGS_SAR'};
% besure to include the script in folder that matlab looks in.

Par.Constr.Dutycycle = 0.1; % intended duty cycle

Par.Constr.SARlocal.Type = 'nlc';  % non linear constraint
Par.Constr.SARlocal.Lim = 8; % limit value
Par.Constr.SARlocal.Nvop  = -1; % use all available VOPs
Par.Constr.SARlocal.Unit = 'W/kg'; % unit of limit
Par.Constr.SARlocal.FreeIterations = 10; % number of iterations free from haulting

Par.Constr.SARlocal.Filename = 'local_sar_vops.mat'; % a file with VOPs (blOCh convention)

Par.Constr.SARglobal.Type = 'nlc';
Par.Constr.SARglobal.Lim = 3;
Par.Constr.SARglobal.Unit = 'W/kg';
Par.Constr.SARglobal.FreeIterations = 10;
Par.Constr.SARglobal.Filename = 'global_sar.mat';

Par.Constr.RFpeak.Type = 'bnd'; % bound as bound with MATLABs fmincon
Par.Constr.RFpeak.Lim = 5000;
Par.Constr.RFpeak.Unit = 'W';
Par.Constr.RFpeak.FreeIterations = 10;

Par.Constr.RFave.Type = 'nlc';
Par.Constr.RFave.Lim = 100;
Par.Constr.RFave.Unit = 'W';
Par.Constr.RFave.FreeIterations = 10;

% what to save during and after optimization. See describtion of subfunction 
% 'Save_Job" in blOCh__opt.m

Save = '111111';


MaxIter = 3; % number of iterations

% How to cope with constraints:
%
% Ignore:  completely ignore it
% Monitor: just monitor the value
% Stop:    monitor and stop if the value is above the limit, but only after 
%          FreeIterations
% Constr:  monitor it and constrain it if need be.



Par.Constr.SARlocal.Cope = 'Constr'; 
Par.Constr.SARglobal.Cope = 'Constr';
Par.Constr.RFpeak.Cope = 'Constr';
Par.Constr.RFave.Cope = 'Monitor';





% now run the optimization
opt = blOCh__opt([],spc,khr,Method,Par,'dt',10e-6,'Init','Random','Handover',0,'TolFun',1e-6,...
    'TolX',1e-6,'MaxIter',MaxIter,'Save',Save,'Mask',0,'deriv_check','off','par_Ncores',1,'par_Type',1,'Show',1);
%%
% finally simulate on a finer grid (spc2)
sim = blOCh__sim([],spc2,khr,opt,'Show',1);
