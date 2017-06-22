clear
close all
clc
% This is an example script for optimizing a 2DRF pulse with GRAPE as
% described in:
%
%  "Fast numerical design of spatial-selective rf pulses in MRI using Krotov and quasi-Newton based optimal control methods"
%  Vinding, M.S., et al.,
%  The Journal of Chemical Physics, 2012, 137, 054203
%  doi: 10.1063/1.4739755
%
% For further information about this example or anything related to blOCh
% please contact me

%%
% make the spatial grid for optimization
%
% init magnetization
M0 = ''; % Default equilibrium magnetization

% target magnetization
% as a box
Md = ''; 
TH = [0.05,0.05]; % desired magnetization square box size [m]
c_TH = [0,0]; % center of box
% or your png file of a gray-scale ROI
%
% Md = 'some_image.png';

% and described by the following:

R = 20; % spatial grid size [#]
Mdflip = [45,0]; % desired flip angle and flip phase [degrees]
L = [0.2,0.2]; % FOX [m]
r0 = [0,0,0]; % centering of FOX

% defaut maps
MPS = ''; % mat file with B1, B0, T1, and T2 maps (blOCh convention)
% or see SAR example how to load pTx maps etc.


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
% in this case the GRAPE Khaneja version
Method = {'GRAPE_Khaneja'};

% what to save during and after optimization. See describtion of subfunction 
% 'Save_Job" in blOCh__opt.m

Save = '000000';

MaxIter = 30; % number of iterations

% parameters according to paper.
Par.lambda = 1e-7;

Par.epsilon = 5e-3;

% Init = './somepulse.txt';
% Columns in this txt-file should be
% [u(s=1),...,u(s=pTx),v(s=1),...,v(s=pTx)] with u and v according to
% Vinding_MRM2017, and each row is a time point.

Init = 'Random';

% now run the optimization
opt = blOCh__opt([],spc,khr,Method,Par,'dt',10e-6,'Init',Init,'Handover',0,'TolFun',1e-7,...
    'TolX',1e-6,'MaxIter',MaxIter,'Save',Save,'Mask',0,'deriv_check','off','par_Ncores',1,'par_Type',1,'Show',0);
%%
% finally simulate on a finer grid (spc2)
sim = blOCh__sim([],spc2,khr,opt,'Show',1);
