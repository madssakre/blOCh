
function opt = blOCh__QN_LBFGS_SAR(spc,khr,opt)
% function opt = blOCh__QN_LBFGS_SAR(spc,khr,opt)
%
%     This script performs a Low-memory BFGS optimization with exact
%     gradients. It has the ability to include strict constraints on local
%     SAR, global SAR, peak and average RF power.
% 
%     Please cite: 
%
%     Local SAR, global SAR, and power-constrained large-flip-angle pulses with optimal control and virtual observation points.
%     Vinding, M.S., et al.,
%     Magn Reson Med. 2017, 77, 1, 374-384.
%     doi: 10.1002/mrm.26086
%
%     Copyright (C) 2017  Mads Sloth Vinding
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

global history opt2
tic_Opt_Wall_Clock_Time = tic;


Par.Grad= 'Schirmer';


Par.Constr.Dutycycle = 0.1;

Par.Constr.SARlocal.Type = 'nlc';
Par.Constr.SARlocal.Lim = 8;
Par.Constr.SARlocal.Nvop  = -1;
Par.Constr.SARlocal.Unit = 'W/kg';
Par.Constr.SARlocal.FreeIterations = 10;
Par.Constr.SARlocal.Cope = 'Ignore';
Par.Constr.SARlocal.Filename = 'default';

Par.Constr.SARglobal.Type = 'nlc';
Par.Constr.SARglobal.Lim = 3;
Par.Constr.SARglobal.Unit = 'W/kg';
Par.Constr.SARglobal.FreeIterations = 10;
Par.Constr.SARglobal.Cope = 'Ignore';
Par.Constr.SARglobal.Filename = 'default';

Par.Constr.RFpeak.Type = 'bnd';
Par.Constr.RFpeak.Lim = 1000;
Par.Constr.RFpeak.Unit = 'Hz';
Par.Constr.RFpeak.FreeIterations = 10;
Par.Constr.RFpeak.Cope = 'Ignore';

Par.Constr.RFave.Type = 'nlc';
Par.Constr.RFave.Lim = 1;
Par.Constr.RFave.Unit = '(rad/s)^2';
Par.Constr.RFave.FreeIterations = 10;
Par.Constr.RFave.Cope = 'Ignore';


Par.RFrobu_num = 1;
Par.RFrobu_min = 0.9;
Par.RFrobu_max = 1.1;
Par.FWHM = 0.1;
Par.RFrobu_distribution = 'none';
Par.RFrobu_show = 0;





opt.Par = blOCh__khr('Get_NewPar',[],[],[],opt.Par,Par);


[opt.RFrobu_val,opt.RFrobu_weight] = RFrobu_4_QN_LBFGS_SAR(opt.Par.RFrobu_num,opt.Par.RFrobu_min,opt.Par.RFrobu_max,opt.Par.RFrobu_distribution,opt.Par.FWHM,opt.Par.RFrobu_show);


opt.k = 0;

opt = Arrange_Constraints_4_QN_LBFGS_SAR(spc,khr,opt);
opt = blOCh__opt('Allocate_Variables',[],[],[],[],spc,khr,opt);

opt.gx = opt.g(1,:);
opt.gy = opt.g(2,:);
opt.gz = opt.g(3,:);





x0 = real(blOCh__opt('Rearrange_controls',[],[],[],[],opt.u,opt.v,opt.gx,opt.gy,opt.gz));

x0(isinf(x0)) = 0;
x0(isnan(x0)) = 0;

opt.maxRF = opt.w1m;

if strcmp(opt.Par.Constr.RFpeak.Cope,'Ignore')
    xlimRF = ones(1,opt.N*spc.pTx*2)*Inf;
else
    if strcmp(opt.Par.Constr.RFpeak.Type,'bnd')
        xlimRF = ones(1,opt.N*spc.pTx*2);
    elseif strcmp(opt.Par.Constr.RFpeak.Type,'nlc')
        xlimRF = ones(1,opt.N*spc.pTx*2)*Inf;
    end
    
end











options=optimset('TolX',opt.TolX,'TolFun',opt.TolFun,'Display','off','MaxIter',opt.MaxIter,'SubproblemAlgorithm','cg',...
    'MaxFunEvals',opt.MaxFunEvals,'Algorithm','interior-point','Hessian',{'lbfgs',10},'ScaleProblem','none',...
    'GradObj','on','DerivativeCheck',opt.deriv_check,'FinDiffType','central','LargeScale','off','GradConstr','on','OutputFcn',@Output_QN_LBFGS_SAR);



opt.idx_u = [1,spc.pTx*opt.N];
opt.idx_v = [spc.pTx*opt.N+1,spc.pTx*opt.N*2];

xlim = xlimRF;

temp = zeros(1,opt.N);
temp(opt.mon) = 1;

temp = repmat(temp,[1,spc.pTx*2]);
xlim = xlim.*temp;

xlim(isnan(xlim)) = eps;

x0_RF = x0(opt.idx_u(1):opt.idx_v(2))./opt.maxRF;
x0 = x0_RF;

opt.Go = true;
opt2 = opt;

[x,fval,opt.exitflag,opt.output] =fmincon(@Iteration_QN_LBFGS_SAR,x0,[],[],[],[],-xlim,xlim,@Constr2_QN_LBFGS_SAR2,options,spc,khr,opt);



blOCh__opt('Display_Message',[],[],[],[],opt.output.message,1);
opt.fig = opt2.fig;
opt.Go = false;
opt.k = 1;
opt.ksafe = history.safeiter;
opt.Par.Constr.SARlocal.Viol = history.Constr.SARlocal.Viol;
opt.Par.Constr.SARglobal.Viol = history.Constr.SARglobal.Viol;
opt.Par.Constr.RFpeak.Viol = history.Constr.RFpeak.Viol;
opt.Par.Constr.RFave.Viol = history.Constr.RFave.Viol;

ksafe = ones(6,1).*opt.ksafe;
if strcmp(opt.Par.Constr.SARlocal,'Stop')
    if opt.Par.Constr.SARlocal.Viol > 0
        ksafe(1) = ksafe(1)-1; if ksafe(1) < 1;ksafe(1) = 1;
        end
    end
end
if strcmp(opt.Par.Constr.SARglobal,'Stop');if opt.Par.Constr.SARglobal.Viol > 0
        ksafe(2) = ksafe(2)-1; if ksafe(2) < 1;ksafe(2) = 1;end;end;end
if strcmp(opt.Par.Constr.RFpeak,'Stop');if opt.Par.Constr.RFpeak.Viol > 0
        ksafe(3) = ksafe(3)-1; if ksafe(3) < 1;ksafe(3) = 1;end;end;end
if strcmp(opt.Par.Constr.RFave,'Stop');if opt.Par.Constr.RFave.Viol > 0
        ksafe(4) = ksafe(4)-1; if ksafe(4) < 1;ksafe(4) = 1;end;end;end


opt.ksafe = min(ksafe);

opt.Fun = history.Fun(1:opt.ksafe);
opt.con.SARlocal   = zeros(1,opt.ksafe);
opt.con.SARglobal   = zeros(1,opt.ksafe);
opt.con.RFpeak   = zeros(1,opt.ksafe);
opt.con.RFave   = zeros(1,opt.ksafe);

if isfield(history,'Durations')
    opt.Durations = history.Durations;
else
    opt.Durations = -1;
end





opt.uo = zeros(spc.pTx,opt.N,opt.ksafe);
opt.vo = zeros(spc.pTx,opt.N,opt.ksafe);


for p = 1:opt.ksafe
    
    opt.k = p;
    x_prog = Rescalex(history.x(p,:),opt);
    
    
    
    [opt.u,opt.v] = blOCh__opt('Rearrange_controls',[],[],[],[],x_prog,opt.N,spc.pTx);
    opt.uo(:,:,p) = opt.u;
    opt.vo(:,:,p) = opt.v;
    
    
    
    
    if ~strcmp(opt.Par.Constr.SARlocal.Cope,'Ignore')
        
        rf = complex(opt.u,opt.v);
        rf_long = rf(:);
        
        [c1,dc1,opt.con.SARlocal(p),opt.Par.Constr.SARlocal.Viol]   = Estim_SARlocal(opt.Par.Constr.SARlocal.Cope,  rf_long,opt.maxRF,opt.Par.Constr.SARlocal.QmatCor,opt.Par.Constr.SARlocal.Lim,opt.Par.Constr.SARlocal.LimCor,opt.Par.Constr.SARlocal.NQ);
    else
        opt.con.SARlocal(p) = history.con(p).SARlocal;
    end
    opt.con.SARglobal(p) = history.con(p).SARglobal;
    opt.con.RFpeak(p)   = history.con(p).RFpeak;
    opt.con.RFave(p)   = history.con(p).RFave;
    
    
    
end


opt.Durations(opt.ksafe+1:end) = [];
clearvars -global history
opt.Opt_Wall_Clock_Time = toc(tic_Opt_Wall_Clock_Time);
opt.Go = false;
end


function opt = Arrange_Constraints_4_QN_LBFGS_SAR(spc,khr,opt)



if ~strcmp(opt.Par.Constr.SARlocal.Cope,'Ignore')
    
    if strcmp(spc.B1_nom_amp_unit,'T/V')
        
        [opt,existQmat,existUnitQmat] = Load_SARlocal(spc,khr,opt);
        
        if existQmat && existUnitQmat
            val_units = {'W/kg/V^2'};
            
            if ~isempty(validatestring(opt.Par.Constr.SARlocal.Unit,val_units))
                opt.Par.Constr.SARlocal.LimCor = opt.Par.Constr.SARlocal.Lim./opt.Par.Constr.Dutycycle;
                opt.Par.Constr.SARlocal.UnitCor = 'W/kg';
            else
                Msg = sprintf('Arrange_Constraints_4_QN_LBFGS_SAR: LocSARLim''s unit (%s) is incompatible with B1 map unit (%s)',opt.Par.Constr.SARlocal.Lim,spc.B1_nom_amp_unit);
                blOCh__opt('Display_Message',[],[],[],[],Msg);
            end
        else
            
            opt.Par.Constr.SARlocal.LimCor = NaN;
            opt.Par.Constr.SARlocal.UnitCor = 'W/kg';
            opt.Par.Constr.SARlocal.NQ = 0;
            opt.Par.Constr.SARlocal.Nvop = 0;
            
            Msg = sprintf('Arrange_Constraints_4_QN_LBFGS_SAR: Something wrong in local SAR compilation');blOCh__opt('Display_Message',[],[],[],[],Msg);
        end
        
        
    else
        Msg = sprintf('Arrange_Constraints_4_QN_LBFGS_SAR: The B1 map unit (%s) is incompatible with the local SAR control',spc.B1_nom_amp_unit);blOCh__opt('Display_Message',[],[],[],[],Msg);
        opt.Par.Constr.SARlocal.LimCor = NaN;
        opt.Par.Constr.SARlocal.UnitCor = 'W/kg';
        opt.Par.Constr.SARlocal.NQ = 0;
        opt.Par.Constr.SARlocal.Nvop = 0;
    end
    
else
    opt.Par.Constr.SARlocal.LimCor = NaN;
    opt.Par.Constr.SARlocal.UnitCor = 'W/kg';
    opt.Par.Constr.SARlocal.NQ = 0;
    opt.Par.Constr.SARlocal.Nvop = 0;
end

if ~strcmp(opt.Par.Constr.SARglobal.Cope,'Ignore')
    
    if strcmp(spc.B1_nom_amp_unit,'T/V')
        
        [opt,existQmat,existUnitQmat] = Load_SARglobal(spc,khr,opt);
        
        if existQmat && existUnitQmat
            val_units = {'W/kg'};
            
            if ~isempty(validatestring(opt.Par.Constr.SARglobal.Unit,val_units))
                opt.Par.Constr.SARglobal.LimCor = opt.Par.Constr.SARglobal.Lim./opt.Par.Constr.Dutycycle;
                opt.Par.Constr.SARglobal.UnitCor = 'W/kg';
            else
                Msg = sprintf('Arrange_Constraints_4_QN_LBFGS_SAR: LocSARLim''s unit (%s) is incompatible with B1 map unit (%s)',opt.Par.Constr.SARglobal.Lim,spc.B1_nom_amp_unit);
                blOCh__opt('Display_Message',[],[],[],[],Msg);
            end
        else
            
            opt.Par.Constr.SARglobal.LimCor = NaN;
            opt.Par.Constr.SARglobal.UnitCor = 'W/kg';
            
            Msg=sprintf('Arrange_Constraints_4_QN_LBFGS_SAR: Something wrong in global SAR compilation');blOCh__opt('Display_Message',[],[],[],[],Msg);
        end
        
        
    else
        Msg = sprintf('Arrange_Constraints_4_QN_LBFGS_SAR: The B1 map unit (%s) is incompatible with the global SAR control',spc.B1_nom_amp_unit);blOCh__opt('Display_Message',[],[],[],[],Msg);
        opt.Par.Constr.SARglobal.LimCor = NaN;
        opt.Par.Constr.SARglobal.UnitCor = 'W/kg';
    end
    
else
    opt.Par.Constr.SARglobal.LimCor = NaN;
    opt.Par.Constr.SARglobal.UnitCor = 'W/kg';
end


%%
% RF peak
if strcmp(opt.Par.Constr.RFpeak.Type,'nlc') || strcmp(opt.Par.Constr.RFpeak.Type,'bnd')
    
    switch opt.Par.Constr.RFpeak.Unit
        
        case 'W'
            
            opt.Par.Constr.RFpeak.Power = 2;
            
            if strcmp(spc.B1_nom_amp_unit,'T/V')
                
                opt.Par.Constr.RFpeak.Unit_cor_factor = 8*50*khr.gamma^2/spc.f_B1_val^2; % (rad/s)^2/W
                
                opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Unit_cor_factor;
                opt.Par.Constr.RFpeak.UnitCor = '(rad/s)^2';
                
                
                opt.Par.Constr.RFpeak.Unit_cor_factor_VV_2_W = khr.gamma^2*spc.B1_nom_amp_val^2./opt.Par.Constr.RFpeak.Unit_cor_factor;
            else
                blOCh__opt('Display_Message',[],[],[],[],sprintf('blOCh__opt: Par.Constr.RFpeak.Unit''s value (%s) is incompatible with B1 map unit (%s)',opt.Par.Constr.RFpeak.Unit,spc.B1_nom_amp_unit),1);
                
            end
            
        case 'V'
            
            opt.Par.Constr.RFpeak.Power = 1;
            if strcmp(spc.B1_nom_amp_unit,'T/V')
                
                opt.Par.Constr.RFpeak.Unit_cor_factor = khr.gamma/spc.f_B1_val;
                
                opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Unit_cor_factor;
                opt.Par.Constr.RFpeak.UnitCor = 'rad/s';
            else
                blOCh__opt('Display_Message',[],[],[],[],sprintf('blOCh__opt: Par.Constr.RFpeak.Unit''s value (%s) is incompatible with B1 map unit (%s)',opt.Par.Constr.RFpeak.Unit,spc.B1_nom_amp_unit),1);
                
            end
        case 'rad/s'
            
            opt.Par.Constr.RFpeak.Power = 1;
            
            opt.Par.Constr.RFpeak.Unit_cor_factor = 1;
            
            opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Unit_cor_factor;
            opt.Par.Constr.RFpeak.UnitCor = 'rad/s';
        case '(rad/s)^2'
            
            opt.Par.Constr.RFpeak.Power = 2;
            
            opt.Par.Constr.RFpeak.Unit_cor_factor = 1;
            
            opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Unit_cor_factor;
            opt.Par.Constr.RFpeak.UnitCor = '(rad/s)^2';
        case 'T'
            
            opt.Par.Constr.RFpeak.Power = 1;
            
            opt.Par.Constr.RFpeak.Unit_cor_factor = khr.gamma;
            opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Unit_cor_factor;
            opt.Par.Constr.RFpeak.UnitCor = 'rad/s';
        case 'Hz'
            
            opt.Par.Constr.RFpeak.Power = 1;
            opt.Par.Constr.RFpeak.Unit_cor_factor = 2*pi;
            opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Unit_cor_factor;
            opt.Par.Constr.RFpeak.UnitCor = 'rad/s';
    end
    
    if opt.Par.Constr.RFpeak.Power == 1
        opt.w1m = opt.Par.Constr.RFpeak.LimCor./sqrt(2);
    else
        opt.w1m = sqrt(opt.Par.Constr.RFpeak.LimCor)./sqrt(2);
    end
else
    opt.Par.Constr.RFpeak.LimCor = NaN;
    opt.Par.Constr.RFpeak.UnitCor = 'rad/s';% MSV 150811
    %     opt.Par.Constr.RFpeak.LimCor = {1,'rad/s'};
    opt.w1m = 2*pi*10e3;
end
% for fmincon's, fmin's and xlim's

if ~isfield(opt,'w1m')
    opt.w1m = 2*pi*10e3;
    
end
if strcmp(opt.Par.Constr.RFave.Type,'nlc')
    
    switch opt.Par.Constr.RFave.Unit
        
        case 'W'
            
            opt.Par.Constr.RFave.Power = 2;
            
            if strcmp(spc.B1_nom_amp_unit,'T/V')
                
                
                opt.Par.Constr.RFave.Unit_cor_factor = 8*50*khr.gamma^2/spc.f_B1_val^2;
                
                opt.Par.Constr.RFave.LimCor = opt.Par.Constr.RFave.Lim*opt.Par.Constr.RFave.Unit_cor_factor./opt.Par.Constr.Dutycycle;
                opt.Par.Constr.RFave.UnitCor = '(rad/s)^2';
                
            else
                blOCh__opt('Display_Message',[],[],[],[],sprintf('blOCh__3_9__opt: AvLim''s unit (%s) is incompatible with B1 map unit (%s)',opt.Par.Constr.RFave.Unit,spc.B1_nom_amp_unit),1);
                
            end
            
        case '(rad/s)^2'
            
            opt.Par.Constr.RFave.Power = 2;
            
            opt.Par.Constr.RFave.Unit_cor_factor = 1;
            
            opt.Par.Constr.RFave.LimCor = opt.Par.Constr.RFave.Lim*opt.Par.Constr.RFave.Unit_cor_factor./opt.Par.Constr.Dutycycle;
            opt.Par.Constr.RFave.UnitCor = '(rad/s)^2';
        case 'T^2'
            
            opt.Par.Constr.RFave.Power = 2;
            opt.Par.Constr.RFave.Unit_cor_factor = (khr.gamma)^2;
            opt.Par.Constr.RFave.LimCor = opt.Par.Constr.RFave.Lim*opt.Par.Constr.RFave.Unit_cor_factor./opt.Par.Constr.Dutycycle;
            opt.Par.Constr.RFave.UnitCor = '(rad/s)^2';
            
        case 'Hz^2'
            
            opt.Par.Constr.RFave.Power = 2;
            opt.Par.Constr.RFave.Unit_cor_factor = (2*pi)^2;
            opt.Par.Constr.RFave.LimCor = opt.Par.Constr.RFave.Lim*opt.Par.Constr.RFave.Unit_cor_factor./opt.Par.Constr.Dutycycle;
            opt.Par.Constr.RFave.UnitCor = '(rad/s)^2';
        otherwise
            blOCh__opt('Display_Message',[],[],[],[],sprintf(' The RFave unit %s is not supported currently\n',opt.Par.Constr.RFave.Unit),1);
            return
    end
else
    opt.Par.Constr.RFave.LimCor = NaN;
    opt.Par.Constr.RFave.UnitCor = '(rad/s)^2';
    opt.Par.Constr.RFave.Unit_cor_factor = 1;
    opt.Par.Constr.RFave.Power = 2;
end

%%


blOCh__opt('Display_Message',[],[],[],[],sprintf('Arrange_Constraints_4_QN_LBFGS_SAR:'));
blOCh__opt('Display_Message',[],[],[],[],sprintf('    SAR local Limit w. dutycycle (%f): %1.2e/%f %s  ~ %1.2e %s',...
    opt.Par.Constr.Dutycycle,...
    opt.Par.Constr.SARlocal.Lim,opt.Par.Constr.Dutycycle,opt.Par.Constr.SARlocal.Unit,...
    opt.Par.Constr.SARlocal.LimCor,opt.Par.Constr.SARlocal.UnitCor),1);

blOCh__opt('Display_Message',[],[],[],[],sprintf('    SAR global Limit w. dutycycle (%f): %1.2e/%f %s  ~ %1.2e %s',...
    opt.Par.Constr.Dutycycle,...
    opt.Par.Constr.SARglobal.Lim,opt.Par.Constr.Dutycycle,opt.Par.Constr.SARglobal.Unit,...
    opt.Par.Constr.SARglobal.LimCor,opt.Par.Constr.SARglobal.UnitCor),1);

blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF Peak Limit: %1.2e %s ~ %1.2e %s',...
    opt.Par.Constr.RFpeak.Lim,opt.Par.Constr.RFpeak.Unit,opt.Par.Constr.RFpeak.LimCor,opt.Par.Constr.RFpeak.UnitCor),1);

blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF Ave Limit w. dutycycle (%f): %1.2e/%f %s ~ %1.2e %s ~ %1.2e %s',...
    opt.Par.Constr.Dutycycle,...
    opt.Par.Constr.RFave.Lim,opt.Par.Constr.Dutycycle,opt.Par.Constr.RFave.Unit,...
    opt.Par.Constr.RFave.Lim/opt.Par.Constr.Dutycycle,opt.Par.Constr.RFave.Unit,...
    opt.Par.Constr.RFave.LimCor,opt.Par.Constr.RFave.UnitCor),1);



end

function x_prog = Rescalex(x,opt)

x_prog = x.*opt.maxRF;
end

function [J,Grad] = Iteration_QN_LBFGS_SAR(x,spc,khr,opt)

global M_T Duration
tic_One_Iteration = tic;
x_prog = Rescalex(x,opt);

[u,v] = blOCh__opt('Rearrange_controls',[],[],[],[],x_prog,opt.N,spc.pTx);

Gu = zeros(spc.pTx,opt.N);
Gv = zeros(spc.pTx,opt.N);

for rfr = 1:opt.Par.RFrobu_num
    
    opt.u = u.*opt.RFrobu_val(rfr);
    opt.v = v.*opt.RFrobu_val(rfr);
    
    opt = SubIteration_QN_LBFGS_SAR(spc,opt,0,0);
    
    [Gu_,Gv_] = blOCh__opt('RF_Update_Terms',[],[],[],[],spc,opt,[]);
    
    
    
    
    Gu = Gu+Gu_.*opt.RFrobu_weight(rfr);
    Gv = Gv+Gv_.*opt.RFrobu_weight(rfr);
end

Gu = -opt.maxRF.*Gu./(spc.P);
Gv = -opt.maxRF.*Gv./(spc.P);



M_T = opt.M_t(:,end);

opt.Fun = blOCh__opt('Efficiency_QN_LBFGS',[],[],[],[],spc,opt,M_T);



J = -opt.Fun;

if nargout == 2
    
    
    
    Grad = blOCh__opt('Rearrange_controls',[],[],[],[],Gu,Gv);
    
    
    
end
Duration = toc(tic_One_Iteration);
end

function [Gu,Gv] = RF_Update_Terms(spc,opt,n)


if ~isempty(n)
    if ismember(n,opt.mon)
        
        switch opt.Par.Grad
            case 'Schirmer'
                [Gu,Gv] =    blOCh__opt('Get_pTx_u_v_Schirmer_Gradients',[],[],[],[],opt.u,opt.v,opt.M_t,opt.L_t,opt.g,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.w0,opt.dt,n);
                
                
            otherwise
                Gu = zeros(spc.pTx,1);
                Gv = zeros(spc.pTx,1);
        end
        
        
    else
        Gu = zeros(spc.pTx,1);
        Gv = zeros(spc.pTx,1);
    end
else
    Gu = zeros(spc.pTx,opt.N);
    Gv = zeros(spc.pTx,opt.N);
    
    
    switch opt.Par.Grad
        
        case 'Schirmer'
            
            if opt.par_Ncores > 1 && opt.par_Type == 1
                % Type 2 parallelizes Gu,Gv per n step
                
                parfor (n = 1:opt.N,opt.par_Ncores)
                    if ismember(n,opt.mon)
                        [Gu(:,n),Gv(:,n)] =    blOCh__opt('Get_pTx_u_v_Schirmer_Gradients',[],[],[],[],opt.u,opt.v,opt.M_t,opt.L_t,opt.g,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.w0,opt.dt,n,opt.par_Ncores,opt.par_Type);
                    end
                end
                
                
            elseif opt.par_Ncores > 1 && opt.par_Type == 2
                % Type 2 parallelizes Gu,Gv per pTx channel
                for n = 1:opt.N
                    if ismember(n,opt.mon)
                        [Gu(:,n),Gv(:,n)] =    blOCh__opt('Get_pTx_u_v_Schirmer_Gradients',[],[],[],[],opt.u,opt.v,opt.M_t,opt.L_t,opt.g,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.w0,opt.dt,n,opt.par_Ncores,opt.par_Type);
                    end
                end
                
                
            else
                for n = 1:opt.N
                    if ismember(n,opt.mon)
                        [Gu(:,n),Gv(:,n)] =    blOCh__opt('Get_pTx_u_v_Schirmer_Gradients',[],[],[],[],opt.u,opt.v,opt.M_t,opt.L_t,opt.g,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.w0,opt.dt,n,opt.par_Ncores,opt.par_Type);
                        
                    end
                end
                
                
            end
            
            
    end
    
end
Gu(isnan(Gu)) = 0;
Gv(isnan(Gv)) = 0;



end

function [con,opt,Grad] = Constr1_QN_LBFGS_SAR(spc,khr,opt)


x = blOCh__opt('Rearrange_controls',[],[],[],[],opt.u,opt.v);



[c,ceq,dc,dceq,con,opt]=Constr2_QN_LBFGS_SAR2(x,spc,khr,opt);


if nargout == 3
    Grad = dc;
end

end



function [c,ceq,dc,dceq,con,opt]=Constr2_QN_LBFGS_SAR2(x,spc,khr,opt)

ceq = [];
dceq   = [];
x = Rescalex(x,opt);


[rfx_,rfy_] = blOCh__opt('Rearrange_controls',[],[],[],[],x,opt.N,spc.pTx);
gx_ = opt.g(1,:);
gy_ = opt.g(2,:);

gz_ = opt.g(3,:);




rf = complex(rfx_,rfy_);
rf_long = rf(:);

%% SAR
% local
if ~strcmp(opt.Par.Constr.SARlocal.Cope,'Ignore') % this assures that SARlocal is not estimated during optimiztion unless there is a need for it
    
    [c1,dc1,con.SARlocal,opt.Par.Constr.SARlocal.Viol]   = Estim_SARlocal(opt.Par.Constr.SARlocal.Cope,  rf_long,opt.maxRF,opt.Par.Constr.SARlocal.QmatCor,opt.Par.Constr.SARlocal.Lim,opt.Par.Constr.SARlocal.LimCor,opt.Par.Constr.SARlocal.NQ);
else
    c1 = [];
    dc1 = [];
    con.SARlocal = 0;
    opt.Par.Constr.SARlocal.Viol = 0;
end
% global
if ~strcmp(opt.Par.Constr.SARglobal.Cope,'Ignore')
    [c2,dc2,con.SARglobal,opt.Par.Constr.SARglobal.Viol] = Estim_SARglobal(opt.Par.Constr.SARglobal.Cope,rf_long,opt.maxRF,opt.Par.Constr.SARglobal.QmatCor,opt.Par.Constr.SARglobal.Lim,opt.Par.Constr.SARglobal.LimCor,opt.N,spc.pTx);
else
    c2 = [];
    dc2 = [];
    con.SARglobal = 0;
    opt.Par.Constr.SARglobal.Viol = 0;
end
%% RF
% peak
[c3,dc3,con.RFpeak,opt.Par.Constr.RFpeak.Viol] = Estim_RFpeak(opt.Par.Constr.RFpeak.Cope,opt.Par.Constr.RFpeak.Type,opt.Par.Constr.RFpeak.Power,opt.Par.Constr.RFpeak.Lim,opt.Par.Constr.RFpeak.LimCor,opt.Par.Constr.RFpeak.Unit_cor_factor,rfx_,rfy_,opt.maxRF);
% ave
[c4,dc4,con.RFave,opt.Par.Constr.RFave.Viol] = Estim_RFave(opt.Par.Constr.RFave.Cope,opt.Par.Constr.RFave.Power,opt.Par.Constr.RFave.Lim,opt.Par.Constr.RFave.LimCor,opt.Par.Constr.RFave.Unit_cor_factor,rf,opt.maxRF,opt.N,spc.pTx);

%% Compile

if ~strcmp(opt.Par.Constr.SARlocal.Cope,'Constr')
    
    c1 = [];
    dc1 = [];
end

if ~strcmp(opt.Par.Constr.SARglobal.Cope,'Constr')
    
    c2 = [];
    dc2 = [];
end

if (strcmp(opt.Par.Constr.RFpeak.Cope,'Constr') && strcmp(opt.Par.Constr.RFpeak.Type,'bnd'))
    
    c3 = [];
    dc3 = [];
elseif (strcmp(opt.Par.Constr.RFpeak.Cope,'Constr') && strcmp(opt.Par.Constr.RFpeak.Type,'nlc'))
    
else
    
    c3 = [];
    dc3 = [];
end


if ~strcmp(opt.Par.Constr.RFave.Cope,'Constr')
    
    c4 = [];
    dc4 = [];
end

c= [c1;c2;c3;c4];

[A1,B1] = size(dc1);
[A2,B2] = size(dc2);
[A3,B3] = size(dc3);
[A4,B4] = size(dc4);

Arf = opt.N*spc.pTx*2;

Atot = Arf;



if A1 == 0
    dc1_ = [];
else
    dc1_ = [dc1;zeros(Atot-Arf,B1)];
end
if A2 == 0
    dc2_ = [];
else
    dc2_ = [dc2;zeros(Atot-Arf,B2)];
end
if A3 == 0
    dc3_ = [];
else
    dc3_ = [dc3;zeros(Atot-Arf,B3)];
end
if A4 == 0
    dc4_ = [];
else
    dc4_ = [dc4;zeros(Atot-Arf,B4)];
end


dc = sparse([dc1_,dc2_,dc3_,dc4_]);


if (strcmp(opt.Par.Constr.RFpeak.Cope,'Constr') && strcmp(opt.Par.Constr.RFpeak.Type,'nlc')) || ...
        strcmp(opt.Par.Constr.RFave.Cope,'Constr') || ...
        strcmp(opt.Par.Constr.SARlocal.Cope,'Constr') || ...
        strcmp(opt.Par.Constr.SARglobal.Cope,'Constr')
    opt.con = con;
    Print_QN_LBFGS_SAR(opt,'con',[],con);
end


end




function [c,dc,SARlocal,Viol] = Estim_SARlocal(Cope,rf,maxRF,Qmat,Lim,LimCor,NQ)

if strcmp(Cope,'Ignore')
    
    c = [];
    dc = [];
    SARlocal = 0;
    Viol = 0;
else
    
    c = zeros(NQ,1);
    if strcmp(Cope,'Constr')
        dc = zeros(length(rf)*2,NQ);
    end
    for s = 1:NQ
        
        temp = rf'*Qmat{s};
        
        if strcmp(Cope,'Constr')
            dc(1:length(rf),s)     =2.*real(temp);
            dc(length(rf)+1:end,s) = -2.*imag(temp);
        end
        
        c(s) = real(temp*rf)   - LimCor;
        
    end
    if strcmp(Cope,'Constr')
        %         dc_ = zeros(size(dc));
        %
        %         for s = 1:NQ
        %             dc_(:,s) = blOCh__opt('Rearrange_controls',[],[],[],[],dc(:,s).',opt.N,spc.pTx);
        %         end
        %         dc = dc_.*maxRF;
        dc = dc.*maxRF;
    else
        dc = [];
    end
    
    [SARlocal,idx] = max(c+LimCor);
    
    if LimCor-SARlocal < 0
        Viol = 1;
    else
        Viol = 0;
    end
    
    if strcmp(Cope,'Stop') || strcmp(Cope,'Monitor')
        c = [];
    end
    
    
end
end

function [c,dc,SARglobal,Viol] = Estim_SARglobal(Cope,rf,maxRF,Qmat,Lim,LimCor,N,pTx)

if strcmp(Cope,'Ignore')
    SARglobal = 0;
    Viol = 0;
    c = [];
    dc = [];
else
    dc = zeros(length(rf)*2,1);
    
    
    temp = rf'*Qmat;
    
    c = real(temp*rf)   - LimCor; %
    
    if strcmp(Cope,'Constr')
        dc(1:length(rf),1)     =2.*real(temp);
        dc(length(rf)+1:end,1) = -2.*imag(temp);
        %         dc_ = blOCh__opt('Rearrange_controls',[],[],[],[],dc.',N,pTx);
        dc = dc(:).*maxRF;
    else
        dc = [];
    end
    
    SARglobal = c+LimCor;
    
    
    if LimCor-SARglobal < 0
        Viol = 1;
    else
        Viol = 0;
    end
    if strcmp(Cope,'Stop') || strcmp(Cope,'Monitor')
        c = [];
    end
    
end
end

function [c,dc,RFpeak,Viol] = Estim_RFpeak(Cope,Type,Power,Lim,LimCor,Unit_cor_factor,rfx,rfy,maxRF)

if strcmp(Cope,'Ignore')
    
    c = [];
    dc = [];
    Viol = 0;
    RFpeak = 0;
else
    
    if Power == 2
        %[x(idx_u(1):idx_u(2)),x(idx_v(1):idx_v(2))]
        rfx_ = permute(rfx,[2,1]); rfx_ = rfx_(:).';
        rfy_ = permute(rfy,[2,1]); rfy_ = rfy_(:).';
        
        c = abs(complex(rfx_,rfy_)).^2 - LimCor;
        if strcmp(Cope,'Constr') && strcmp(Type,'nlc')
            dc = [2*diag(rfx_);2*diag(rfy_)].*maxRF;
        else
            dc = [];
        end
    else
        rfx_ = permute(rfx,[2,1]); rfx_ = rfx_(:).';
        rfy_ = permute(rfy,[2,1]); rfy_ = rfy_(:).';
        c = abs(complex(rfx_,rfy_)).^2 - LimCor.^2;
        if strcmp(Cope,'Constr') && strcmp(Type,'nlc')
            dc = [2*diag(rfx_);2*diag(rfy_)].*maxRF;
        else
            dc = [];
        end
    end
    c = c(:);
    
    
    if Power == 2
        RFpeak = max(c)+LimCor;
        RFpeak = RFpeak./Unit_cor_factor;
    else
        RFpeak = max(c)+LimCor^2;
        RFpeak = RFpeak./Unit_cor_factor;
    end
    
    if Lim-RFpeak < 0
        Viol = 1;
    else
        Viol = 0;
    end
    if strcmp(Cope,'Stop') || strcmp(Cope,'Monitor')
        c = [];
    end
    
end
end

function [c,dc,RFave,Viol] = Estim_RFave(Cope,Power,Lim,LimCor,Unit_cor_factor,rf,maxRF,N,pTx)

if strcmp(Cope,'Ignore')
    
    Viol = 0;
    RFave = 0;
    c = [];
    dc = [];
else
    
    if Power == 2
        c = sum(abs(rf).^2,2)./N - LimCor;
        
        if strcmp(Cope,'Constr')
            tmp = zeros(N*pTx,pTx);
            k = 1;
            for n = 1:N
                tmp(k:pTx+k-1,:) = diag(rf(:,n));
                k = k+pTx;
            end
            
            
            dc_ = 2.*[real(tmp);imag(tmp)];
            dc = zeros(size(dc_));
            for n = 1:pTx
                dc(:,n) = blOCh__opt('Rearrange_controls',[],[],[],[],dc_(:,n).',N,pTx);
            end
            dc = dc.*maxRF./N;
            
        else
            dc =[];
        end
        RFave = max(c)+LimCor;
        
        if LimCor-RFave < 0
            
            Viol = 1;
        else
            Viol = 0;
        end
        
        
        RFave = RFave./Unit_cor_factor; % This conversion must happen after Viol is made to take dutycycle into account
        
        if strcmp(Cope,'Stop') || strcmp(Cope,'Monitor')
            c = [];
        end
    end
    
end
end



function [stop] = Output_QN_LBFGS_SAR(x,optimValues,state,spc,khr,opt)
% function [stop] = Output_QN_LBFGS_SAR(x,optimValues,state,spc,khr,opt)
%
%   This script is an output function for the QN_LBFGS_SAR functions


global M_T history Duration opt2

stop = false;

switch state
    case 'init'
        %               disp('init')
        %              iter =  optimValues.iteration
        history.x = zeros(opt2.MaxIter,length(x)); % This repairs the bug with history.x being stored in matlabs weird background. So now you don't have to "clear all" if you change number of controls between two optimizations
        
        history.x(1,:) = x;
        
        history.fval(1)= optimValues.fval;
        
        history.safeiter  = optimValues.iteration+1;
        
        
        [opt2.u,opt2.v] = blOCh__opt('Rearrange_controls',[],[],[],[],x,opt2.N,spc.pTx);
        opt2.gx = opt.g(1,:);
        opt2.gy = opt.g(2,:);
        
        opt2.gz = opt.g(3,:);
        
        
        
        
        
        opt2.g = [opt2.gx;opt2.gy;opt2.gz];
        
        
        opt2.M_T = M_T;
        history.M_T{1} = M_T;
        
        
        Fun = blOCh__opt('Efficiency_QN_LBFGS',[],[],[],[],spc,opt,M_T);
        
        
        [con,opt2] = Constr1_QN_LBFGS_SAR(spc,khr,opt2);
        
        
        
        
        
        history.con(1) = con;
        
        %       history.Durations(1) = 0;
        opt2.con = history.con;
        history.Fun(1) = Fun;
        opt2 = Plot_Progress_QN_LBFGS_SAR(spc,khr,opt2);
        opt2.dFun = [];
        opt2 = Progress_QN_LBFGS_SAR(spc,khr,opt2);
        %         disp('here')
        Print_QN_LBFGS_SAR(opt2,'all',Fun,con);
        opt2.dFun = 0;
        history.Constr.SARlocal.Viol = opt2.Par.Constr.SARlocal.Viol;
        history.Constr.SARglobal.Viol = opt2.Par.Constr.SARglobal.Viol;
        history.Constr.RFpeak.Viol = opt2.Par.Constr.RFpeak.Viol;
        history.Constr.RFave.Viol = opt2.Par.Constr.RFave.Viol;
        
    case 'iter'
        %       disp('iter')
        if optimValues.iteration > 0
            iter = optimValues.iteration+1;
            %             disp('now here')
            %             iter
            %             opt2.k = iter;
            
            history.x(iter,:) = x;
            
            history.fval(iter)= optimValues.fval;
            history.M_T{iter} = M_T;
            
            history.safeiter  = iter;
            
            
            
            
            
            
            [opt2.u,opt2.v] = blOCh__opt('Rearrange_controls',[],[],[],[],x,opt2.N,spc.pTx);
            opt2.gx = opt.g(1,:);
            opt2.gy = opt.g(2,:);
            
            opt2.gz = opt.g(3,:);
            
            opt2.g = [opt2.gx;opt2.gy;opt2.gz];
            opt2.M_T = M_T;
            Fun = blOCh__opt('Efficiency_QN_LBFGS',[],[],[],[],spc,opt,M_T);
            
            
            [con,opt2] = Constr1_QN_LBFGS_SAR(spc,khr,opt2);
            
            
            
            history.Fun(iter) = Fun;
            
            history.con(iter) = con;
            opt2.con = history.con;
            opt2.Fun = history.Fun;
            opt2 = Plot_Progress_QN_LBFGS_SAR(spc,khr,opt2);
            opt2.dFun = history.Fun(iter)-history.Fun(iter-1);%1; %abs(history.fval(iter))-abs(history.fval(iter-1));
            opt2.xintermed = x;
            opt2 = Progress_QN_LBFGS_SAR(spc,khr,opt2);
            Print_QN_LBFGS_SAR(opt2,'all',Fun,con);
            history.dFun(iter) = opt2.dFun;
            history.Constr.SARlocal.Viol = opt2.Par.Constr.SARlocal.Viol;
            history.Constr.SARglobal.Viol = opt2.Par.Constr.SARglobal.Viol;
            history.Constr.RFpeak.Viol = opt2.Par.Constr.RFpeak.Viol;
            history.Constr.RFave.Viol = opt2.Par.Constr.RFave.Viol;
            history.Durations(iter) = Duration;
            if ~opt2.Go
                stop = true;
            end
            
        end
    case 'done'
        %         disp(1)
        %         hold off
    otherwise
end
end

function opt = SubIteration_QN_LBFGS_SAR(spc,opt,Option,Penalty)

if nargin == 2
    Option = 1; % This propagates state vectors and updates controls
end
R11 = cell(1,opt.N);R12 = cell(1,opt.N);R13 = cell(1,opt.N);
R21 = cell(1,opt.N);R22 = cell(1,opt.N);R23 = cell(1,opt.N);
R31 = cell(1,opt.N);R32 = cell(1,opt.N);R33 = cell(1,opt.N);
opt.g = [opt.gx;opt.gy;opt.gz];
for n = 1:opt.N
    [R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n}] = blOCh__opt('Get_Rotator',[],[],[],[],opt.u,opt.v,opt.g,opt.w0,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.dt,n);
end



for n = 1:opt.N
    opt.M_t(:,n+1) = blOCh__opt('Rotate',[],[],[],[],opt.M_t(:,n),R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n},'Forward');
end
for n = opt.N:-1:1
    opt.L_t(:,n) = blOCh__opt('Rotate',[],[],[],[],opt.L_t(:,n+1),R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n},'Backward');
end


if Option
    [Gu,Gv] = RF_Update_Terms(spc,opt,[]);
    
    for n = 1:opt.N
        opt = Update_GRAPE_Khaneja_new(spc,opt,Gu,Gv,n);
    end
    
end
end


function Print_QN_LBFGS_SAR(opt,Type,Fun,con)

if nargin == 2
    Type = 'all';
end


switch Type
    case 'all'
        
        %         opt.k
        if opt.k == 1
            
            Msg = sprintf('\nIter    Fun            Constraint (Value & Violation Boolean)');
            blOCh__opt('Display_Message',[],[],[],[],Msg,1);
            
            Msg = sprintf(  '                       SARlocal    SARglobal   RFpeak      RFave     ');
            blOCh__opt('Display_Message',[],[],[],[],Msg,1);
            
            
            
            
            Msg = sprintf('% 5i   %1.3e      %1.3e %i %1.3e %i %1.3e %i %1.3e %i  ',opt.k-1,Fun,con.SARlocal,opt.Par.Constr.SARlocal.Viol,con.SARglobal,opt.Par.Constr.SARglobal.Viol,con.RFpeak,opt.Par.Constr.RFpeak.Viol,con.RFave,opt.Par.Constr.RFave.Viol);
            blOCh__opt('Display_Message',[],[],[],[],Msg,1);
            
            
            
        else
            Msg = sprintf('% 5i   %1.3e      %1.3e %i %1.3e %i %1.3e %i %1.3e %i  ',opt.k-1,Fun,con.SARlocal,opt.Par.Constr.SARlocal.Viol,con.SARglobal,opt.Par.Constr.SARglobal.Viol,con.RFpeak,opt.Par.Constr.RFpeak.Viol,con.RFave,opt.Par.Constr.RFave.Viol);
            blOCh__opt('Display_Message',[],[],[],[],Msg,1);
            
        end
        
        
    case 'con'
        Msg = sprintf('                       %1.3e %i %1.3e %i %1.3e %i %1.3e %i  ',con.SARlocal,opt.Par.Constr.SARlocal.Viol,con.SARglobal,opt.Par.Constr.SARglobal.Viol,con.RFpeak,opt.Par.Constr.RFpeak.Viol,con.RFave,opt.Par.Constr.RFave.Viol);
        
        blOCh__opt('Display_Message',[],[],[],[],Msg,1);
end

end

function opt = Progress_QN_LBFGS_SAR(spc,khr,opt)


Stop = zeros(4,1);
if strcmp(opt.Par.Constr.SARlocal.Cope,'Stop') && opt.Par.Constr.SARlocal.Viol && opt.k > opt.Par.Constr.SARlocal.FreeIterations
    Stop(1) = 1;
end
if strcmp(opt.Par.Constr.SARglobal.Cope,'Stop') && opt.Par.Constr.SARlocal.Viol && opt.k > opt.Par.Constr.SARglobal.FreeIterations
    Stop(2) = 1;
end
if strcmp(opt.Par.Constr.RFpeak.Cope,'Stop') && opt.Par.Constr.RFpeak.Viol && opt.k > opt.Par.Constr.RFpeak.FreeIterations
    Stop(3) = 1;
end
if strcmp(opt.Par.Constr.RFave.Cope,'Stop') && opt.Par.Constr.RFave.Viol && opt.k > opt.Par.Constr.RFave.FreeIterations
    Stop(4) = 1;
end


% save intermediate
if opt.k-1 > 0
    if mod(opt.k-1,opt.ksaveintermediate) == 0
        
        %     isfield(opt,'xintermed')
        tempsave = opt.Save; % take backup
        tempk = opt.k;
        opt.k = opt.k-1;
        %
        
        opt.Save = '000000'; % reset
        opt.Save(2) = tempsave(2); % get initial intention back for data saving
        opt.Save(3) = tempsave(3); % get initial intention back for controls saving
        % thus, if intermediate saving of data and controls were intended it
        % will be saved. But other (bundle, pics, scripts etc.) will not, at
        % least not before the end if intended.
        try blOCh__opt('Save_Job',[],[],[],[],spc,khr,opt);
        catch me;
            blOCh__opt('Display_Message',[],[],[],[],['Progress_QN_LBFGS_SAR: Save_Job: ',me.message],2);
        end
        
        opt.Save = tempsave;
        opt.k = tempk;
        
    end
end

if sum(Stop)>0
    opt.Go = false;
    opt.ksafe = opt.k;
    return
else
    opt.Go = true;
    opt.ksafe = opt.k;
    opt.k = opt.k+1;
end




end

function opt = Update_GRAPE_Khaneja_new(spc,opt,Gu,Gv,n)


switch opt.Par.Grad
    case {'1st','type1'}
        opt.u(:,n) = opt.u(:,n)+opt.Par.epsilon.*Gu(:,n)./(2*opt.Par.lambda*spc.P);
        opt.v(:,n) = opt.v(:,n)+opt.Par.epsilon.*Gv(:,n)./(2*opt.Par.lambda*spc.P);
    case {'Schirmer'}
        opt.u(:,n) = opt.u(:,n)+opt.Par.epsilon.*Gu(:,n)./(2*opt.Par.lambda);
        opt.v(:,n) = opt.v(:,n)+opt.Par.epsilon.*Gv(:,n)./(2*opt.Par.lambda);
        
end
opt.g(:,n) = [opt.gx(n);opt.gy(n);opt.gz(n)];
end
%%
function [val,wei] = RFrobu_4_QN_LBFGS_SAR(Num,Min,Max,Dist,FWHM,Show)


if abs(abs(Max)-abs(Min)) <= eps
    blOCh__opt('Display_Message',[],[],[],[],'RFrobu_4_QN_LBFGS_SAR: RF robustness limits equal, applying eps',1);
    Max = Max + eps;
    Min = Min - eps;
end
if Num == 1
    val = 1;
    wei = 1;
    %     return
else
    val = linspace(Min,Max,Num);
end

switch Dist
    case 'lorentzian'
        % 1./(pi./FWHM/2).*(FWHM/2)^2.
        wei = 1./(1+((val-1)/(FWHM/2)).^2);
        wei = wei./sum(wei);
        
        
        %         wei = 1./(pi./FWHM/2).*(FWHM/2)^2.*1./((val-1).^2+(FWHM/2)^2);
        % %         ?
        %         wei = wei./sum(wei);
        
        blOCh__opt('Display_Message',[],[],[],[],sprintf('RFrobu_4_QN_LBFGS_SAR:'));
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF inhomogeneity distribution:     %s',Dist),1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    Number of RF inhomogeneity points: %i',Num),1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF inhomogeneity FWHM [%%%%]:         %f',FWHM*100),1);
        
        
    case 'evenly'
        
        wei = ones(1,Num)./Num;
        
        blOCh__opt('Display_Message',[],[],[],[],sprintf('RFrobu_4_QN_LBFGS_SAR:'));
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF inhomogeneity distribution:     %s',Dist),1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    Number of RF inhomogeneity points: %i',Num),1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF inhomogeneity minimum [%%%%]:      %f',Min*100),1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF inhomogeneity maximum [%%%%]:      %f',Max*100),1);
        
    case 'none'
        wei = 1;
        val = 1;
    otherwise
        wei = ones(1,Num)./Num;
        blOCh__opt('Display_Message',[],[],[],[],'RFrobu_4_QN_LBFGS_SAR: RF robustness weighting distribution unknown, assigning a distribution with even weights',1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('RFrobu_4_QN_LBFGS_SAR:'));
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF inhomogeneity distribution:     %s',Dist),1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    Number of RF inhomogeneity points: %i',Num),1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF inhomogeneity minimum [%%%%]:      %f',Min*100),1);
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    RF inhomogeneity maximum [%%%%]:      %f',Max*100),1);
end

switch Dist
    case 'none'
        
        
        
    otherwise
        blOCh__opt('Display_Message',[],[],[],[],sprintf('    Value [%%%%]\t\tWeight'),1);
        for n=1:Num
            blOCh__opt('Display_Message',[],[],[],[],sprintf('    %1.4f\t\t%1.4f',val(n)*100,wei(n)),1);
        end
        if Show
            
            
            figure
            plot(val,wei,'r')
            xlabel('RF value')
            ylabel('RF inhomogeneity weight')
            title('RF robustness')
            drawnow
        end
end







end


%%

function opt = Plot_Progress_QN_LBFGS_SAR(spc,khr,opt)

if opt.Show
    if ~isfield(opt,'h')
        
        mp = get(0, 'MonitorPositions');
        if size(mp,1) > 1
            mp = mp(1,:);
        end
        ratio = mp(4)/mp(3);
        
        opt.fig = figure('color','white');
        set(opt.fig,'Units','pixels')
        set(opt.fig,'Position',[mp(1),mp(2),round(mp(3)*0.9),round(mp(4)*0.9)])
        set(opt.fig,'Units','normalized')
        
        opt.h.ax(1) = axes;
        opt.h.ax(2) = axes;
        opt.h.ax(3) = axes;
        opt.h.ax(4) = axes;
        opt.h.ax(5) = axes;
        opt.h.ax(6) = axes;
        opt.h.ax(7) = axes;
        opt.h.ax(8) = axes;
        opt.h.ax(9) = axes;
        % 		disp('here')
        h.Hax56789 = 0.1;
        h.Wax56789 = 0.95;
        
        h.Yax9 = 0.05;
        h.Yax8 = h.Yax9+h.Hax56789;
        h.Yax7 = h.Yax8+h.Hax56789;
        h.Yax6 = h.Yax7+h.Hax56789;
        h.Yax5 = h.Yax6+h.Hax56789;
        
        h.Yax1234 = h.Yax5+h.Hax56789+0.05;
        h.Hax1234 = 1-h.Yax1234;
        h.Wax1234 = (1-0.05*5)/4;
        
        h.Xax1 = 0.05;
        h.Xax2 = h.Xax1+h.Wax1234+0.05;
        h.Xax3 = h.Xax2+h.Wax1234+0.05;
        h.Xax4 = h.Xax3+h.Wax1234+0.05;
        
        set(opt.h.ax(1),'Units','normalized')
        set(opt.h.ax(1),'Position',[h.Xax1 h.Yax1234 h.Wax1234 h.Hax1234])
        
        set(opt.h.ax(2),'Units','normalized')
        set(opt.h.ax(2),'Position',[h.Xax2 h.Yax1234 h.Wax1234 h.Hax1234])
        
        set(opt.h.ax(3),'Units','normalized')
        set(opt.h.ax(3),'Position',[h.Xax3 h.Yax1234 h.Wax1234 h.Hax1234])
        
        set(opt.h.ax(4),'Units','normalized')
        set(opt.h.ax(4),'Position',[h.Xax4 h.Yax1234 h.Wax1234 h.Hax1234])
        
        set(opt.h.ax(5),'Units','normalized')
        set(opt.h.ax(5),'Position',[0.05*ratio h.Yax5 h.Wax56789 h.Hax56789])
        
        set(opt.h.ax(6),'Units','normalized')
        set(opt.h.ax(6),'Position',[0.05*ratio h.Yax6 h.Wax56789 h.Hax56789])
        
        set(opt.h.ax(7),'Units','normalized')
        set(opt.h.ax(7),'Position',[0.05*ratio h.Yax7 h.Wax56789 h.Hax56789])
        
        set(opt.h.ax(8),'Units','normalized')
        set(opt.h.ax(8),'Position',[0.05*ratio h.Yax8 h.Wax56789 h.Hax56789])
        %
        set(opt.h.ax(9),'Units','normalized')
        set(opt.h.ax(9),'Position',[0.05*ratio h.Yax9 h.Wax56789 h.Hax56789])
        
        
        
    end
    % 	opt.h.ax
    
    if isfield(opt,'M_T')
        %         if strcmp(opt.OptMagnType,'x,y')
        %             Mtemp = [real(opt.M_T).';imag(opt.M_T).';zeros(size(opt.M_T)).'];
        
        %             h.M_T = blOCh__list2grid(Mtemp(:),[spc.R,spc.Rv],spc.idxtot,3);
        %         else
        h.M_T = blOCh__spc('list2grid',[],[],[],[],opt.M_T,[spc.R,spc.Rv],spc.idxtot,3);
        %         end
    else
        if isfield(opt,'M_t')
            h.M_T = blOCh__spc('list2grid',[],[],[],[],opt.M_t(:,end),[spc.R,spc.Rv],spc.idxtot,3);
        else
            disp('No Magn to plot')
            h.M_T = blOCh__spc('list2grid',[],[],[],[],spc.Md,[spc.R,spc.Rv],spc.idxtot,3);
        end
    end
    
    
    h.u = opt.u;
    h.v = opt.v;
    h.g = opt.g;
    % 	h
    for n = 1:9
        
        h.Type = n;
        h=Plotting(spc,khr,opt,h);
    end
    
    drawnow
end
end

function h=Plotting(spc,khr,opt,h)

h = what2plot(spc,khr,opt,h);

if h.Type <= 4
    switch spc.Dim(1:2)
        
        case {'1D','0+'}
            
            
            axes(opt.h.ax(h.Type))
            plot(h.ordinate)
            set(opt.h.ax(h.Type),'YLim',[-1,1])
            set(opt.h.ax(h.Type),'XTick',h.Dim1axis)
            set(opt.h.ax(h.Type),'XTickLabel',h.Dim1ticklabel)
            xlabel(h.Dim2label)
        case {'1+','2D','2+','3D','3+'}
            % 			opt
            % isfield(opt,'h')
            % opt.fig;
            % gcf
            % opt.h.ax
            % gca
            
            axes(opt.h.ax(h.Type))
            imagesc(h.ordinate,[h.ordinate_mn,h.ordinate_mx])
            set(opt.h.ax(h.Type),'XTick',h.Dim2axis,'YTick',h.Dim1axis)
            set(opt.h.ax(h.Type),'XTickLabel',h.Dim2ticklabel,'YTickLabel',h.Dim1ticklabel)
            set(opt.h.ax(h.Type),'Ydir','normal')
            xlabel(h.Dim2label)
            ylabel(h.Dim1label)
            
    end
else
    axes(opt.h.ax(h.Type))
    if ~isfield(opt,'gSRok')
        gSRok = true;
    else
        gSRok = opt.gSRok;
    end
    if gSRok
        plot(h.ordinate.','g'), axis tight
    else
        plot(h.ordinate.','r'), axis tight
        
    end
    if h.Type < 9
        set(opt.h.ax(h.Type),'XTick',[])
    end
end


end

function h = what2plot(spc,khr,opt,h)


switch h.Type
    
    case 1
        
        
        h.Colormapno = 1;
        %       'Mx [M0]';
        ordinate = h.M_T(:,:,:,:,1);
        h.ordinate_mn = -1;
        h.ordinate_mx = 1;
    case 2
        h.Colormapno = 1;
        %       'My [M0]';
        ordinate = h.M_T(:,:,:,:,2);
        h.ordinate_mn = -1;
        h.ordinate_mx = 1;
    case 3
        %       'Mz [M0]';
        ordinate = h.M_T(:,:,:,:,3);
        h.ordinate_mn = -1;
        h.ordinate_mx = 1;
    case 4
        
        h.Colormapno = 1;
        ordinate = abs(complex(h.M_T(:,:,:,:,1),h.M_T(:,:,:,:,2)));
        h.ordinate_mn = 0;
        h.ordinate_mx = 1;
        %       '|Mxy| [M0]';
        
    case 5
        ordinate = abs(complex(h.u,h.v));
    case 6
        ordinate = angle(complex(h.u,h.v));
    case 7
        ordinate = h.g(1,:);
    case 8
        ordinate = h.g(2,:);
    case 9
        ordinate = h.g(3,:);
        
end
if h.Type <= 4
    switch spc.Dim
        case '0+1D'
            h.Dim1ticklabel = linspace(spc.Dv(1),spc.Dv(2),5);
            h.Dim1axis = linspace(1,spc.Rv,5);
            h.Dim2label = ' v [Hz]';
            h.ordinate = squeeze(ordinate);
        case '1DSI'
            h.Dim1ticklabel = linspace(spc.D(1,3),spc.D(2,3),5);
            h.Dim1axis = linspace(1,spc.R(3),5);
            h.Dim2label = ' z [m]';
            h.ordinate = squeeze(ordinate);
        case '1DAP'
            h.Dim1ticklabel = linspace(spc.D(1,2),spc.D(2,2),5);
            h.Dim1axis = linspace(1,spc.R(2),5);
            h.Dim2label = ' y [m]';
            h.ordinate = squeeze(ordinate);
        case '1DRL'
            h.Dim1ticklabel = linspace(spc.D(1,1),spc.D(2,1),5);
            h.Dim1axis = linspace(1,spc.R(1),5);
            h.Dim2label = ' x [m]';
            h.ordinate = squeeze(ordinate);
        case '1+1DSI'
            h.Dim1ticklabel = linspace(spc.D(1,3),spc.D(2,3),5);
            h.Dim2ticklabel = linspace(spc.Dv(1),spc.Dv(2),5);
            h.Dim1axis = linspace(1,spc.R(3),5);
            h.Dim2axis = linspace(1,spc.Rv,5);
            h.Dim1label = ' z [m]';
            h.Dim2label = ' v [Hz]';
            h.ordinate = squeeze(ordinate);
        case '1+1DAP'
            h.Dim1ticklabel = linspace(spc.D(1,3),spc.D(2,3),5);
            h.Dim2ticklabel = linspace(spc.Dv(1),spc.Dv(2),5);
            h.Dim1axis = linspace(1,spc.R(2),5);
            h.Dim2axis = linspace(1,spc.Rv,5);
            h.Dim1label = ' y [m]';
            h.Dim2label = ' v [Hz]';
            h.ordinate = squeeze(ordinate);
        case '1+1DRL'
            h.Dim1ticklabel = linspace(spc.D(1,1),spc.D(2,1),5);
            h.Dim2ticklabel = linspace(spc.Dv(1),spc.Dv(2),5);
            h.Dim1axis = linspace(1,spc.R(1),5);
            h.Dim2axis = linspace(1,spc.Rv,5);
            h.Dim1label = ' x [m]';
            h.Dim2label = ' v [Hz]';
            h.ordinate = squeeze(ordinate);
        case '2DAx'
            h.Dim2ticklabel = linspace(spc.D(1,1),spc.D(2,1),5);
            h.Dim1ticklabel = linspace(spc.D(1,2),spc.D(2,2),5);
            h.Dim2axis = linspace(1,spc.R(1),5);
            h.Dim1axis = linspace(1,spc.R(2),5);
            h.Dim1label = ' y [m]';
            h.Dim2label = ' x [m]';
            h.ordinate = squeeze(ordinate);
        case '2DCo'
            h.Dim2ticklabel = linspace(spc.D(1,1),spc.D(2,1),5);
            h.Dim1ticklabel = linspace(spc.D(1,3),spc.D(2,3),5);
            h.Dim2axis = linspace(1,spc.R(1),5);
            h.Dim1axis = linspace(1,spc.R(3),5);
            
            h.Dim2label = ' x [m]';
            h.Dim1label = ' z [m]';
            h.ordinate = squeeze(ordinate);
        case '2DSa'
            h.Dim1ticklabel = linspace(spc.D(1,2),spc.D(2,2),5);
            h.Dim2ticklabel = linspace(spc.D(1,3),spc.D(2,3),5);
            h.Dim1axis = linspace(1,spc.R(2),5);
            h.Dim2axis = linspace(1,spc.R(3),5);
            h.Dim2label = ' y [m]';
            h.Dim1label = ' z [m]';
            h.ordinate = squeeze(ordinate);
        case '2+1DAx'
            h.freqNo = round(spc.Rv/2);
            h.Dim1ticklabel = linspace(spc.D(1,1),spc.D(2,1),5);
            h.Dim2ticklabel = linspace(spc.D(1,2),spc.D(2,2),5);
            h.Dim1axis = linspace(1,spc.R(1),5);
            h.Dim2axis = linspace(1,spc.R(2),5);
            h.Dim2label = ' y [m]';
            h.Dim1label = ' x [m]';
            h.ordinate = squeeze(ordinate(:,:,:,h.freqNo));
        case '2+1DCo'
            h.freqNo = round(spc.Rv/2);
            h.Dim2ticklabel = linspace(spc.D(1,1),spc.D(2,1),5);
            h.Dim1ticklabel = linspace(spc.D(1,3),spc.D(2,3),5);
            h.Dim2axis = linspace(1,spc.R(1),5);
            h.Dim1axis = linspace(1,spc.R(3),5);
            
            h.Dim2label = ' x [m]';
            h.Dim1label = ' z [m]';
            h.ordinate = squeeze(ordinate(:,:,:,h.freqNo));
        case '2+1DSa'
            h.freqNo = round(spc.Rv/2);
            h.Dim1ticklabel = linspace(spc.D(1,2),spc.D(2,2),5);
            h.Dim2ticklabel = linspace(spc.D(1,3),spc.D(2,3),5);
            h.Dim1axis = linspace(1,spc.R(2),5);
            h.Dim2axis = linspace(1,spc.R(3),5);
            h.Dim2label = ' y [m]';
            h.Dim1label = ' z [m]';
            h.ordinate = squeeze(ordinate(:,:,:,h.freqNo));
        case '3D'
            h.SlNo = round(spc.R(3)/2);
            h.Dim1ticklabel = linspace(spc.D(1,1),spc.D(2,1),5);
            h.Dim2ticklabel = linspace(spc.D(1,2),spc.D(2,2),5);
            h.Dim1axis = linspace(1,spc.R(1),5);
            h.Dim2axis = linspace(1,spc.R(2),5);
            h.Dim2label = ' x [m]';
            h.Dim1label = ' y [m]';
            h.ordinate = squeeze(ordinate(:,:,h.SlNo,1));
        case '3+1D'
            h.SlNo = round(spc.R(3)/2);
            h.freqNo = round(spc.Rv/2);
            h.Dim1ticklabel = linspace(spc.D(1,1),spc.D(2,1),5);
            h.Dim2ticklabel = linspace(spc.D(1,2),spc.D(2,2),5);
            h.Dim1axis = linspace(1,spc.R(1),5);
            h.Dim2axis = linspace(1,spc.R(2),5);
            h.Dim2label = ' x [m]';
            h.Dim1label = ' y [m]';
            h.ordinate = squeeze(ordinate(:,:,h.SlNo,h.freqNo));
    end
else
    h.ordinate = ordinate;
    
    
    
end

end


%%

function [opt,existQmat,existUnitQmat] = Load_SARlocal(spc,khr,opt)
Msg = [];

if exist(opt.Par.Constr.SARlocal.Filename,'file')
    
    SAR = load(opt.Par.Constr.SARlocal.Filename,'-mat');
    
else
    SAR.Qmat = Load_default_SARlocal_Qmat;
    SAR.UnitQmat = 'W/kg/V^2';
end
existQmat = 0;
existUnitQmat = 0;
if isfield(SAR,'Qmat');
    [x1,x2,x3] = size(SAR.Qmat);
    if x1 > x2 && x2 == x3
        if x2 == spc.pTx
            SAR.Qmat = permute(SAR.Qmat,[2,3,1]);
            existQmat = 1;
        else
            Msg = 'Load_SARlocal: Q of wrong size';blOCh__opt('Display_Message',[],[],[],[],Msg);
        end
        
    elseif x2 > x1 && x1 == x3
        if x1 == spc.pTx
            SAR.Qmat = permute(SAR.Qmat,[1,3,2]);
            existQmat = 1;
        else
            Msg = 'Load_SARlocal: Q of wrong size';blOCh__opt('Display_Message',[],[],[],[],Msg);
        end
    elseif x3 > x1 && x1 == x2
        if x1 == spc.pTx
            existQmat = 1;
        else
            Msg = 'Load_SARlocal: Q of wrong size';blOCh__opt('Display_Message',[],[],[],[],Msg);
        end
    else
        Msg = 'Load_SARlocal: Q of wrong size';blOCh__opt('Display_Message',[],[],[],[],Msg);
    end
    
else
    existQmat = 0;
    Msg = 'Load_SARlocal: The file must include a field named Q';blOCh__opt('Display_Message',[],[],[],[],Msg);
end



if isfield(SAR,'UnitQmat');
    try validateattributes(SAR.UnitQmat,{'char'},{'nonempty'});
        
        val_units = {'W/kg/V^2'};
        if ~isempty(validatestring(SAR.UnitQmat,val_units))
            UnitQmat = SAR.UnitQmat;
            existUnitQmat = 1;
            
            UnitQmat_cor_factor = spc.f_B1_val^2./khr.gamma^2;
            
        else
            Msg = 'Load_SARlocal: The SARmat did not include a UnitQmat or a proper Q_unit (''W/kg/V^2''). Cannot continue.';blOCh__opt('Display_Message',[],[],[],[],Msg);
            existUnitQmat = 0;
        end
    catch me;
        Msg = ['Load_SARlocal: ',me.message,' The SARmat did not include a UnitQmat or a proper UnitQmat (''W/kg/V^2''). Cannot continue.'];blOCh__opt('Display_Message',[],[],[],[],Msg);
        existUnitQmat = 0;
    end
else
    
end

if existQmat && existUnitQmat
    if opt.Par.Constr.SARlocal.Nvop <= 0
        
        
        opt.Par.Constr.SARlocal.NQ = size(SAR.Qmat,3);
        Msg = sprintf('Load_SARlocal: Importing all %i local SAR Q matrices.',opt.Par.Constr.SARlocal.NQ);blOCh__opt('Display_Message',[],[],[],[],Msg);
        
    else
        if opt.Par.Constr.SARlocal.Nvop < size(SAR.Qmat,3)
            opt.Par.Constr.SARlocal.NQ = opt.Par.Constr.SARlocal.Nvop;
            Msg = sprintf('Load_SARlocal: Importing %i of %i local SAR Q matrices.',opt.Par.Constr.SARlocal.Nvop,size(SAR.Qmat,3));blOCh__opt('Display_Message',[],[],[],[],Msg);
            
        else
            opt.Par.Constr.SARlocal.NQ = size(SAR.Qmat,3);
            Msg = sprintf('Load_SARlocal: Importing all %i local SAR Q matrices.',opt.Par.Constr.SARlocal.NQ);blOCh__opt('Display_Message',[],[],[],[],Msg);
        end
    end
    
    I = speye(opt.N,opt.N);
    LocU = cell(opt.Par.Constr.SARlocal.NQ,1);
    for s = 1:opt.Par.Constr.SARlocal.NQ
        LocU{s} = kron(I,SAR.Qmat(:,:,s)./opt.N*UnitQmat_cor_factor);
        
    end
    
    opt.Par.Constr.SARlocal.QmatCor = LocU;
    opt.Par.Constr.SARlocal.UnitQmatCor = 'W/kg/(rad/s)^2';
    
    
    opt.Par.Constr.SARlocal.Q_unit_cor_factor = UnitQmat_cor_factor; % [V^2 s^2 / rad^2]
    opt.Par.Constr.SARlocal.Q_unit_cor_invfactor = 1/UnitQmat_cor_factor;
    
    
else
    opt.Par.Constr.SARlocal.QmatCor = [];
    opt.Par.Constr.SARlocal.UnitQmatCor = 'W/kg/(rad/s)^2';
    
end

end

function [opt,existQmat,existUnitQmat] = Load_SARglobal(spc,khr,opt)
Msg = [];
if exist(opt.Par.Constr.SARglobal.Filename,'file')
    
    SAR = load(opt.Par.Constr.SARglobal.Filename,'-mat');
    
else
    SAR.Qmat = Load_default_SARglobal_Qmat;
    SAR.UnitQmat = 'W/kg/V^2';
end
existQmat = 0;
existUnitQmat = 0;
if isfield(SAR,'Qmat');
    [x1,x2,x3] = size(SAR.Qmat);
    if x1 == x2 && x2 == spc.pTx && x3 == 1
        existQmat = 1;
    else
        Msg = 'Load_SARglobal: Q of wrong size';
    end
    
else
    existQmat = 0;
    Msg = 'Load_SARglobal: The file must include a field named Q';blOCh__opt('Display_Message',[],[],[],[],Msg);
end



if isfield(SAR,'UnitQmat')
    try validateattributes(SAR.UnitQmat,{'char'},{'nonempty'});
        
        val_units = {'W/kg/V^2'};
        if ~isempty(validatestring(SAR.UnitQmat,val_units))
            UnitQmat = SAR.UnitQmat;
            existUnitQmat = 1;
            
            UnitQmat_cor_factor = spc.f_B1_val^2./khr.gamma^2;
            
        else
            Msg = 'Load_SARglobal: The SARmat did not include a Q_unit or a proper Q_unit (''W/kg/V^2''). Cannot continue.';blOCh__opt('Display_Message',[],[],[],[],Msg);
            existUnitQmat = 0;
        end
    catch me;
        Msg = ['Load_SARglobal: ',me.message,' The SARmat did not include a Q_unit or a proper Q_unit (''W/kg/V^2''). Cannot continue.'];blOCh__opt('Display_Message',[],[],[],[],Msg);
        existUnitQmat = 0;
    end
else
    
end

if existQmat && existUnitQmat
    
    
    I = speye(opt.N,opt.N);
    
    
    opt.Par.Constr.SARglobal.UnitQmatCor = 'W/kg/(rad/s)^2';
    
    
    opt.Par.Constr.SARglobal.QmatCor = kron(I,SAR.Qmat./opt.N*UnitQmat_cor_factor);
    
    
    opt.Par.Constr.SARglobal.UnitQmat_cor_factor = UnitQmat_cor_factor; % [V^2 s^2 / rad^2]
    opt.Par.Constr.SARglobal.UnitQmat_cor_invfactor = 1/UnitQmat_cor_factor;
    
    
else
    opt.Par.Constr.SARglobal.QmatCor = [];
    opt.Par.Constr.SARglobal.UnitQmatCor = 'W/kg/(rad/s)^2';
    
end

end


function Qmat = Load_default_SARlocal_Qmat
% This is a Qmat generator making a representative Qmat array.
locQmatr = [7.92e-04 -4.64e-05 1.49e-05 -1.49e-06 2.74e-05 1.29e-04 5.08e-04 -1.92e-06 -4.64e-05 9.20e-04 -4.69e-05 3.67e-05 7.15e-05 -2.65e-04 -9.66e-05 -7.79e-05 1.49e-05 -4.69e-05 7.60e-04 -3.79e-05 1.11e-05 6.53e-05 4.84e-05 2.88e-05 -1.49e-06 3.67e-05 -3.79e-05 7.29e-04 2.05e-05 -1.50e-04 4.41e-04 3.63e-05 2.74e-05 7.15e-05 1.11e-05 2.05e-05 7.09e-04 -5.84e-05 2.07e-04 1.25e-06 1.29e-04 -2.65e-04 6.53e-05 -1.50e-04 -5.84e-05 1.33e-03 -5.62e-04 2.35e-05 5.08e-04 -9.66e-05 4.84e-05 4.41e-04 2.07e-04 -5.62e-04 6.31e-03 2.28e-04 -1.92e-06 -7.79e-05 2.88e-05 3.63e-05 1.25e-06 2.35e-05 2.28e-04 8.75e-04
    1.11e-03 7.45e-06 1.31e-04 -3.94e-05 2.12e-04 -2.65e-04 -1.19e-04 1.76e-04 7.45e-06 6.24e-04 -5.98e-05 -7.63e-05 2.53e-05 1.41e-04 -5.64e-05 6.50e-05 1.31e-04 -5.98e-05 7.14e-04 2.78e-05 -7.06e-04 -4.36e-04 -5.85e-05 1.49e-04 -3.94e-05 -7.63e-05 2.78e-05 6.46e-04 1.05e-05 -7.20e-05 -2.82e-05 -9.06e-05 2.12e-04 2.53e-05 -7.06e-04 1.05e-05 3.88e-03 1.30e-03 2.00e-04 -3.78e-04 -2.65e-04 1.41e-04 -4.36e-04 -7.20e-05 1.30e-03 1.72e-03 1.42e-04 6.89e-05 -1.19e-04 -5.64e-05 -5.85e-05 -2.82e-05 2.00e-04 1.42e-04 5.82e-04 -3.69e-05 1.76e-04 6.50e-05 1.49e-04 -9.06e-05 -3.78e-04 6.89e-05 -3.69e-05 1.15e-03
    1.26e-03 -1.93e-05 2.60e-04 -7.62e-06 5.94e-06 -4.42e-04 -1.89e-04 5.24e-05 -1.93e-05 3.39e-04 -9.54e-05 -9.81e-05 -3.85e-05 1.25e-04 -8.66e-05 1.66e-05 2.60e-04 -9.54e-05 4.38e-04 4.46e-05 -6.91e-04 -3.90e-04 -4.82e-05 1.32e-04 -7.62e-06 -9.81e-05 4.46e-05 3.45e-04 -1.16e-04 -8.27e-05 7.39e-05 -6.53e-05 5.94e-06 -3.85e-05 -6.91e-04 -1.16e-04 3.44e-03 1.42e-03 2.17e-04 -2.51e-04 -4.42e-04 1.25e-04 -3.90e-04 -8.27e-05 1.42e-03 1.68e-03 1.43e-04 1.67e-04 -1.89e-04 -8.66e-05 -4.82e-05 7.39e-05 2.17e-04 1.43e-04 3.43e-04 -4.40e-05 5.24e-05 1.66e-05 1.32e-04 -6.53e-05 -2.51e-04 1.67e-04 -4.40e-05 7.76e-04
    2.00e-03 -6.49e-05 2.63e-04 8.77e-05 -1.25e-04 -2.74e-04 -1.33e-04 8.25e-05 -6.49e-05 6.31e-04 -2.09e-05 -8.77e-05 -6.08e-05 4.79e-05 -8.24e-05 6.64e-05 2.63e-04 -2.09e-05 6.64e-04 9.60e-05 -6.06e-04 -2.91e-04 1.19e-05 2.94e-05 8.77e-05 -8.77e-05 9.60e-05 6.36e-04 -7.24e-05 -1.62e-05 2.69e-05 -3.65e-05 -1.25e-04 -6.08e-05 -6.06e-04 -7.24e-05 3.50e-03 1.52e-03 2.35e-04 -1.71e-04 -2.74e-04 4.79e-05 -2.91e-04 -1.62e-05 1.52e-03 1.66e-03 1.15e-04 1.21e-04 -1.33e-04 -8.24e-05 1.19e-05 2.69e-05 2.35e-04 1.15e-04 6.13e-04 -1.97e-05 8.25e-05 6.64e-05 2.94e-05 -3.65e-05 -1.71e-04 1.21e-04 -1.97e-05 8.75e-04
    4.58e-04 -1.13e-04 1.33e-04 -8.14e-05 -7.86e-05 -8.85e-05 7.40e-05 -1.03e-04 -1.13e-04 8.45e-04 1.29e-04 -2.18e-05 4.27e-05 -2.04e-05 7.29e-06 2.19e-04 1.33e-04 1.29e-04 2.68e-03 -4.31e-04 4.02e-04 -1.03e-04 1.92e-03 3.88e-04 -8.14e-05 -2.18e-05 -4.31e-04 6.47e-04 -8.19e-05 6.32e-05 -4.13e-04 -7.49e-05 -7.86e-05 4.27e-05 4.02e-04 -8.19e-05 3.93e-04 5.72e-05 4.26e-04 8.86e-05 -8.85e-05 -2.04e-05 -1.03e-04 6.32e-05 5.72e-05 4.90e-04 -1.48e-04 -8.29e-05 7.40e-05 7.29e-06 1.92e-03 -4.13e-04 4.26e-04 -1.48e-04 2.68e-03 2.48e-04 -1.03e-04 2.19e-04 3.88e-04 -7.49e-05 8.86e-05 -8.29e-05 2.48e-04 7.79e-04
    6.69e-04 -3.63e-04 3.48e-05 -9.18e-05 -6.55e-05 -7.08e-05 -2.58e-04 3.36e-05 -3.63e-04 2.66e-03 -1.65e-05 -3.00e-05 3.80e-04 2.03e-04 1.93e-03 1.08e-04 3.48e-05 -1.65e-05 8.78e-04 -9.54e-05 4.57e-05 1.33e-04 1.19e-04 5.10e-05 -9.18e-05 -3.00e-05 -9.54e-05 5.82e-04 -1.40e-04 -7.34e-05 -1.59e-04 -9.14e-05 -6.55e-05 3.80e-04 4.57e-05 -1.40e-04 3.92e-04 7.21e-05 4.96e-04 7.64e-05 -7.08e-05 2.03e-04 1.33e-04 -7.34e-05 7.21e-05 7.04e-04 1.54e-04 -1.20e-04 -2.58e-04 1.93e-03 1.19e-04 -1.59e-04 4.96e-04 1.54e-04 2.71e-03 1.15e-04 3.36e-05 1.08e-04 5.10e-05 -9.14e-05 7.64e-05 -1.20e-04 1.15e-04 5.46e-04
    6.55e-04 -3.20e-04 -2.39e-05 -1.04e-04 -3.37e-05 -1.28e-05 -2.71e-04 4.60e-06 -3.20e-04 2.91e-03 1.60e-05 -4.10e-06 5.26e-04 3.13e-04 1.86e-03 1.35e-04 -2.39e-05 1.60e-05 5.76e-04 -4.12e-05 -2.14e-05 1.78e-04 -7.39e-05 2.97e-05 -1.04e-04 -4.10e-06 -4.12e-05 5.07e-04 -1.92e-04 -1.17e-04 -3.95e-05 -1.13e-04 -3.37e-05 5.26e-04 -2.14e-05 -1.92e-04 4.03e-04 7.50e-05 4.82e-04 8.28e-05 -1.28e-05 3.13e-04 1.78e-04 -1.17e-04 7.50e-05 1.05e-03 -1.70e-04 -2.71e-05 -2.71e-04 1.86e-03 -7.39e-05 -3.95e-05 4.82e-04 -1.70e-04 2.61e-03 6.82e-05 4.60e-06 1.35e-04 2.97e-05 -1.13e-04 8.28e-05 -2.71e-05 6.82e-05 4.10e-04
    7.80e-04 -1.04e-04 4.78e-06 1.13e-04 3.63e-04 -2.80e-04 -1.04e-04 -1.48e-04 -1.04e-04 4.55e-04 -7.41e-05 -2.08e-05 -2.75e-04 8.65e-05 3.52e-05 -1.07e-05 4.78e-06 -7.41e-05 4.32e-04 5.67e-06 -1.66e-04 -1.10e-04 -1.31e-04 -8.44e-06 1.13e-04 -2.08e-05 5.67e-06 7.22e-04 6.71e-04 1.64e-05 -3.10e-05 -1.00e-04 3.63e-04 -2.75e-04 -1.66e-04 6.71e-04 3.83e-03 1.44e-04 2.21e-04 -3.66e-04 -2.80e-04 8.65e-05 -1.10e-04 1.64e-05 1.44e-04 8.64e-04 1.73e-04 1.62e-04 -1.04e-04 3.52e-05 -1.31e-04 -3.10e-05 2.21e-04 1.73e-04 5.38e-04 2.04e-04 -1.48e-04 -1.07e-05 -8.44e-06 -1.00e-04 -3.66e-04 1.62e-04 2.04e-04 1.09e-03
    6.78e-04 -2.55e-04 -6.38e-05 -4.84e-05 -8.49e-06 -2.12e-04 2.73e-05 3.29e-05 -2.55e-04 3.15e-03 -2.53e-04 -2.45e-05 4.95e-04 1.79e-03 4.51e-04 -1.41e-04 -6.38e-05 -2.53e-04 6.31e-04 1.07e-04 -1.06e-04 -1.17e-04 -7.74e-05 -9.19e-06 -4.84e-05 -2.45e-05 1.07e-04 7.25e-04 -1.74e-04 -1.36e-05 -8.07e-05 -5.40e-05 -8.49e-06 4.95e-04 -1.06e-04 -1.74e-04 5.20e-04 3.65e-04 1.57e-04 1.78e-05 -2.12e-04 1.79e-03 -1.17e-04 -1.36e-05 3.65e-04 2.20e-03 -1.55e-04 -8.69e-05 2.73e-05 4.51e-04 -7.74e-05 -8.07e-05 1.57e-04 -1.55e-04 1.51e-03 -1.51e-04 3.29e-05 -1.41e-04 -9.19e-06 -5.40e-05 1.78e-05 -8.69e-05 -1.51e-04 5.45e-04
    1.39e-03 -3.36e-05 2.10e-04 -1.16e-05 -8.22e-05 -6.05e-04 -1.85e-04 4.27e-06 -3.36e-05 4.81e-04 5.68e-05 -9.17e-05 -3.80e-05 7.29e-05 -9.44e-05 5.78e-08 2.10e-04 5.68e-05 5.14e-04 6.60e-05 -5.50e-04 -4.49e-04 -4.31e-05 9.55e-05 -1.16e-05 -9.17e-05 6.60e-05 6.11e-04 -1.76e-04 -7.02e-05 -7.96e-06 1.17e-06 -8.22e-05 -3.80e-05 -5.50e-04 -1.76e-04 3.21e-03 1.26e-03 2.21e-04 -1.70e-04 -6.05e-04 7.29e-05 -4.49e-04 -7.02e-05 1.26e-03 1.67e-03 1.05e-04 2.22e-04 -1.85e-04 -9.44e-05 -4.31e-05 -7.96e-06 2.21e-04 1.05e-04 5.29e-04 -4.64e-05 4.27e-06 5.78e-08 9.55e-05 1.17e-06 -1.70e-04 2.22e-04 -4.64e-05 8.64e-04
    8.75e-04 -1.79e-04 -7.39e-05 -2.26e-04 -1.06e-04 -3.21e-04 1.03e-04 -6.66e-05 -1.79e-04 6.83e-04 -5.70e-05 2.59e-05 -3.25e-04 1.03e-04 -3.08e-05 -3.51e-04 -7.39e-05 -5.70e-05 6.20e-04 7.79e-05 1.90e-04 -2.08e-05 -1.00e-04 1.19e-04 -2.26e-04 2.59e-05 7.79e-05 9.99e-04 8.44e-04 2.24e-04 -1.21e-04 1.32e-04 -1.06e-04 -3.25e-04 1.90e-04 8.44e-04 3.16e-03 2.79e-04 -4.91e-05 8.96e-04 -3.21e-04 1.03e-04 -2.08e-05 2.24e-04 2.79e-04 8.29e-04 -4.81e-05 7.55e-05 1.03e-04 -3.08e-05 -1.00e-04 -1.21e-04 -4.91e-05 -4.81e-05 6.33e-04 4.95e-05 -6.66e-05 -3.51e-04 1.19e-04 1.32e-04 8.96e-04 7.55e-05 4.95e-05 1.73e-03
    2.18e-03 -3.95e-05 3.31e-04 1.11e-04 -2.51e-04 9.72e-05 -9.19e-05 1.98e-04 -3.95e-05 5.52e-04 -8.27e-06 -6.95e-05 -2.58e-05 5.67e-05 -1.78e-05 7.50e-05 3.31e-04 -8.27e-06 6.16e-04 6.30e-05 -5.93e-04 -2.22e-04 9.15e-05 2.03e-05 1.11e-04 -6.95e-05 6.30e-05 5.60e-04 -6.71e-05 2.55e-05 4.07e-05 -7.09e-05 -2.51e-04 -2.58e-05 -5.93e-04 -6.71e-05 3.12e-03 1.43e-03 1.77e-04 -1.35e-04 9.72e-05 5.67e-05 -2.22e-04 2.55e-05 1.43e-03 1.53e-03 1.08e-04 5.68e-05 -9.19e-05 -1.78e-05 9.15e-05 4.07e-05 1.77e-04 1.08e-04 4.61e-04 -6.99e-05 1.98e-04 7.50e-05 2.03e-05 -7.09e-05 -1.35e-04 5.68e-05 -6.99e-05 6.43e-04
    7.32e-04 -2.82e-04 -2.18e-05 -5.04e-05 -3.89e-05 -1.28e-04 1.67e-04 -5.71e-05 -2.82e-04 2.63e-03 -2.99e-04 1.63e-05 4.82e-04 1.81e-03 2.83e-05 -2.29e-04 -2.18e-05 -2.99e-04 6.49e-04 1.16e-04 -1.41e-04 -2.14e-04 -1.05e-04 2.03e-05 -5.04e-05 1.63e-05 1.16e-04 6.34e-04 -1.70e-04 -2.29e-05 -1.74e-04 -7.06e-05 -3.89e-05 4.82e-04 -1.41e-04 -1.70e-04 5.22e-04 3.91e-04 1.52e-04 2.49e-05 -1.28e-04 1.81e-03 -2.14e-04 -2.29e-05 3.91e-04 2.23e-03 8.95e-05 -2.22e-04 1.67e-04 2.83e-05 -1.05e-04 -1.74e-04 1.52e-04 8.95e-05 8.53e-04 -7.17e-05 -5.71e-05 -2.29e-04 2.03e-05 -7.06e-05 2.49e-05 -2.22e-04 -7.17e-05 4.40e-04
    6.85e-04 -1.32e-04 -7.41e-05 -2.04e-04 -1.79e-04 -2.74e-04 1.37e-04 -1.79e-04 -1.32e-04 4.92e-04 -8.99e-05 1.14e-04 -3.43e-04 1.36e-04 -1.08e-04 -3.48e-04 -7.41e-05 -8.99e-05 4.98e-04 1.27e-04 1.07e-04 -4.02e-05 -1.63e-04 1.29e-04 -2.04e-04 1.14e-04 1.27e-04 9.84e-04 6.76e-04 1.56e-04 -1.77e-04 2.74e-05 -1.79e-04 -3.43e-04 1.07e-04 6.76e-04 2.76e-03 1.96e-04 -6.81e-05 1.02e-03 -2.74e-04 1.36e-04 -4.02e-05 1.56e-04 1.96e-04 6.51e-04 -7.51e-05 1.25e-04 1.37e-04 -1.08e-04 -1.63e-04 -1.77e-04 -6.81e-05 -7.51e-05 6.05e-04 1.17e-04 -1.79e-04 -3.48e-04 1.29e-04 2.74e-05 1.02e-03 1.25e-04 1.17e-04 1.95e-03
    7.47e-04 -1.89e-06 6.34e-05 2.75e-05 4.28e-04 -1.48e-04 -7.12e-05 -7.36e-05 -1.89e-06 3.76e-04 -3.97e-05 3.95e-06 -1.21e-04 1.24e-04 -9.57e-05 -4.62e-05 6.34e-05 -3.97e-05 3.25e-04 3.86e-05 -5.25e-04 -2.91e-04 -1.02e-04 1.50e-04 2.75e-05 3.95e-06 3.86e-05 3.26e-04 5.21e-05 -8.92e-05 -2.01e-05 -1.14e-04 4.28e-04 -1.21e-04 -5.25e-04 5.21e-05 3.67e-03 4.67e-04 1.81e-04 -4.54e-04 -1.48e-04 1.24e-04 -2.91e-04 -8.92e-05 4.67e-04 1.26e-03 1.02e-04 1.54e-04 -7.12e-05 -9.57e-05 -1.02e-04 -2.01e-05 1.81e-04 1.02e-04 4.68e-04 6.81e-05 -7.36e-05 -4.62e-05 1.50e-04 -1.14e-04 -4.54e-04 1.54e-04 6.81e-05 1.04e-03
    5.85e-04 3.81e-05 4.12e-05 -6.62e-05 2.21e-04 -2.16e-04 -8.83e-05 -6.99e-05 3.81e-05 3.43e-04 -4.19e-06 -4.15e-05 -7.48e-05 6.47e-05 -1.19e-04 -7.89e-05 4.12e-05 -4.19e-06 2.98e-04 4.79e-05 -5.75e-04 -2.83e-04 -8.99e-05 1.11e-04 -6.62e-05 -4.15e-05 4.79e-05 1.96e-04 -1.40e-04 -1.05e-04 1.89e-05 -5.28e-05 2.21e-04 -7.48e-05 -5.75e-04 -1.40e-04 3.31e-03 7.42e-04 2.04e-04 -3.32e-04 -2.16e-04 6.47e-05 -2.83e-04 -1.05e-04 7.42e-04 1.38e-03 5.39e-05 2.06e-04 -8.83e-05 -1.19e-04 -8.99e-05 1.89e-05 2.04e-04 5.39e-05 4.54e-04 2.56e-05 -6.99e-05 -7.89e-05 1.11e-04 -5.28e-05 -3.32e-04 2.06e-04 2.56e-05 8.28e-04
    8.94e-04 -1.85e-04 -2.02e-04 -1.06e-04 1.98e-04 4.23e-04 -1.78e-04 3.57e-06 -1.85e-04 1.87e-03 1.73e-05 -1.50e-06 4.50e-04 -3.71e-04 -4.46e-04 -1.11e-04 -2.02e-04 1.73e-05 6.37e-04 2.95e-05 -6.77e-05 -3.82e-04 3.22e-04 -7.78e-05 -1.06e-04 -1.50e-06 2.95e-05 5.26e-04 3.20e-05 -1.52e-04 1.02e-04 8.17e-05 1.98e-04 4.50e-04 -6.77e-05 3.20e-05 7.62e-04 -5.41e-05 -1.89e-04 -6.19e-05 4.23e-04 -3.71e-04 -3.82e-04 -1.52e-04 -5.41e-05 2.43e-03 -1.53e-03 1.71e-04 -1.78e-04 -4.46e-04 3.22e-04 1.02e-04 -1.89e-04 -1.53e-03 2.59e-03 -1.65e-04 3.57e-06 -1.11e-04 -7.78e-05 8.17e-05 -6.19e-05 1.71e-04 -1.65e-04 5.64e-04
    5.32e-04 -3.00e-04 -3.94e-05 -6.83e-05 5.43e-05 -1.44e-04 -1.06e-04 -1.04e-05 -3.00e-04 3.10e-03 -1.39e-04 3.07e-05 4.98e-04 1.46e-03 6.29e-04 -8.39e-05 -3.94e-05 -1.39e-04 4.74e-04 1.74e-05 -1.20e-04 -8.95e-05 -8.52e-09 2.21e-05 -6.83e-05 3.07e-05 1.74e-05 4.58e-04 -2.03e-04 -7.40e-05 -5.59e-05 -5.63e-05 5.43e-05 4.98e-04 -1.20e-04 -2.03e-04 4.25e-04 2.72e-04 1.52e-04 8.34e-05 -1.44e-04 1.46e-03 -8.95e-05 -7.40e-05 2.72e-04 1.99e-03 -4.31e-04 5.40e-05 -1.06e-04 6.29e-04 -8.52e-09 -5.59e-05 1.52e-04 -4.31e-04 1.78e-03 -1.32e-04 -1.04e-05 -8.39e-05 2.21e-05 -5.63e-05 8.34e-05 5.40e-05 -1.32e-04 3.72e-04
    4.20e-04 -3.45e-05 -2.81e-05 -8.61e-06 2.08e-04 -2.17e-04 -2.92e-05 -9.09e-05 -3.45e-05 2.91e-04 -5.85e-05 2.27e-05 -3.00e-04 1.99e-05 5.77e-05 3.09e-05 -2.81e-05 -5.85e-05 2.53e-04 -8.65e-07 -2.88e-04 -4.04e-05 -1.23e-04 -4.84e-06 -8.61e-06 2.27e-05 -8.65e-07 5.76e-04 1.81e-04 3.46e-05 -5.08e-05 2.46e-05 2.08e-04 -3.00e-04 -2.88e-04 1.81e-04 3.49e-03 -7.41e-05 2.39e-04 -3.98e-04 -2.17e-04 1.99e-05 -4.04e-05 3.46e-05 -7.41e-05 5.88e-04 6.11e-05 1.68e-04 -2.92e-05 5.77e-05 -1.23e-04 -5.08e-05 2.39e-04 6.11e-05 2.54e-04 1.58e-04 -9.09e-05 3.09e-05 -4.84e-06 2.46e-05 -3.98e-04 1.68e-04 1.58e-04 1.11e-03
    5.16e-04 4.23e-05 1.57e-05 -7.74e-05 4.09e-04 -1.97e-04 -1.13e-04 -2.27e-04 4.23e-05 3.53e-04 -3.24e-05 -5.45e-06 -2.19e-04 9.63e-05 -9.76e-05 -8.11e-05 1.57e-05 -3.24e-05 1.91e-04 1.85e-05 -4.62e-04 -1.22e-04 -7.77e-05 1.12e-04 -7.74e-05 -5.45e-06 1.85e-05 1.34e-04 -2.42e-05 -6.79e-05 7.07e-06 -1.84e-05 4.09e-04 -2.19e-04 -4.62e-04 -2.42e-05 3.31e-03 8.24e-05 1.79e-04 -4.40e-04 -1.97e-04 9.63e-05 -1.22e-04 -6.79e-05 8.24e-05 1.04e-03 9.55e-05 2.20e-04 -1.13e-04 -9.76e-05 -7.77e-05 7.07e-06 1.79e-04 9.55e-05 4.04e-04 1.69e-04 -2.27e-04 -8.11e-05 1.12e-04 -1.84e-05 -4.40e-04 2.20e-04 1.69e-04 8.30e-04
    5.81e-04 -1.48e-06 -4.21e-05 -4.22e-05 -2.25e-05 3.73e-05 -1.98e-05 -9.51e-05 -1.48e-06 5.80e-04 -2.32e-05 -2.95e-05 -7.64e-06 -5.35e-05 2.84e-05 -1.44e-04 -4.21e-05 -2.32e-05 1.23e-03 2.66e-04 3.96e-04 -3.39e-04 -3.66e-04 -9.42e-05 -4.22e-05 -2.95e-05 2.66e-04 9.68e-04 3.00e-04 -3.05e-04 -3.75e-05 5.43e-04 -2.25e-05 -7.64e-06 3.96e-04 3.00e-04 9.17e-04 -2.44e-04 -1.39e-04 1.06e-04 3.73e-05 -5.35e-05 -3.39e-04 -3.05e-04 -2.44e-04 8.13e-04 9.63e-05 -3.95e-04 -1.98e-05 2.84e-05 -3.66e-04 -3.75e-05 -1.39e-04 9.63e-05 9.85e-04 2.78e-04 -9.51e-05 -1.44e-04 -9.42e-05 5.43e-04 1.06e-04 -3.95e-04 2.78e-04 2.49e-03
    4.45e-04 -1.02e-04 1.45e-04 -5.67e-05 -6.97e-05 -6.08e-05 1.07e-04 -1.33e-04 -1.02e-04 6.63e-04 6.34e-05 -2.84e-05 -1.07e-05 -4.22e-05 -4.69e-05 1.47e-04 1.45e-04 6.34e-05 2.67e-03 -4.18e-04 4.06e-04 -8.77e-05 1.43e-03 4.51e-04 -5.67e-05 -2.84e-05 -4.18e-04 6.24e-04 -3.95e-05 6.15e-05 -3.44e-04 -4.47e-05 -6.97e-05 -1.07e-05 4.06e-04 -3.95e-05 4.70e-04 7.75e-05 2.65e-04 7.68e-05 -6.08e-05 -4.22e-05 -8.77e-05 6.15e-05 7.75e-05 4.67e-04 -1.34e-04 -1.71e-05 1.07e-04 -4.69e-05 1.43e-03 -3.44e-04 2.65e-04 -1.34e-04 2.10e-03 -4.01e-06 -1.33e-04 1.47e-04 4.51e-04 -4.47e-05 7.68e-05 -1.71e-05 -4.01e-06 9.15e-04
    6.10e-04 -4.56e-04 1.44e-04 -1.41e-04 -1.22e-04 -4.93e-05 -2.29e-04 -5.73e-05 -4.56e-04 2.08e-03 5.50e-05 5.64e-05 2.52e-04 3.64e-05 1.10e-03 2.46e-04 1.44e-04 5.50e-05 2.05e-03 -4.88e-04 2.45e-04 2.69e-06 1.13e-03 1.22e-04 -1.41e-04 5.64e-05 -4.88e-04 7.72e-04 -1.25e-04 9.27e-06 -4.02e-04 -5.56e-05 -1.22e-04 2.52e-04 2.45e-04 -1.25e-04 5.13e-04 1.50e-05 5.89e-04 6.25e-05 -4.93e-05 3.64e-05 2.69e-06 9.27e-06 1.50e-05 6.56e-04 6.39e-05 -2.20e-04 -2.29e-04 1.10e-03 1.13e-03 -4.02e-04 5.89e-04 6.39e-05 2.67e-03 1.75e-04 -5.73e-05 2.46e-04 1.22e-04 -5.56e-05 6.25e-05 -2.20e-04 1.75e-04 7.32e-04
    9.47e-04 -3.61e-04 3.01e-05 -7.16e-05 -1.00e-04 -1.84e-05 8.31e-05 -4.82e-05 -3.61e-04 2.40e-03 -2.74e-04 5.64e-05 4.83e-04 1.53e-03 -1.11e-05 -1.81e-04 3.01e-05 -2.74e-04 6.97e-04 1.41e-04 -1.73e-04 -2.67e-04 -1.53e-04 1.99e-05 -7.16e-05 5.64e-05 1.41e-04 6.93e-04 -1.32e-04 -9.10e-05 -1.87e-04 -3.93e-05 -1.00e-04 4.83e-04 -1.73e-04 -1.32e-04 5.75e-04 3.96e-04 1.98e-04 -9.95e-06 -1.84e-05 1.53e-03 -2.67e-04 -9.10e-05 3.96e-04 2.17e-03 4.82e-05 -2.89e-04 8.31e-05 -1.11e-05 -1.53e-04 -1.87e-04 1.98e-04 4.82e-05 7.77e-04 -2.51e-05 -4.82e-05 -1.81e-04 1.99e-05 -3.93e-05 -9.95e-06 -2.89e-04 -2.51e-05 5.03e-04
    8.06e-04 -2.66e-05 -2.44e-05 5.85e-05 7.45e-04 -7.59e-05 4.06e-05 -6.43e-05 -2.66e-05 4.25e-04 -2.80e-05 7.18e-06 -1.67e-04 -3.17e-05 -3.69e-05 -9.65e-05 -2.44e-05 -2.80e-05 1.82e-04 -3.64e-07 -3.46e-04 -5.56e-05 5.00e-06 1.48e-04 5.85e-05 7.18e-06 -3.64e-07 4.92e-04 2.05e-04 5.68e-05 -2.51e-05 -1.45e-04 7.45e-04 -1.67e-04 -3.46e-04 2.05e-04 3.11e-03 1.28e-04 1.44e-04 -5.63e-04 -7.59e-05 -3.17e-05 -5.56e-05 5.68e-05 1.28e-04 5.07e-04 4.87e-05 1.59e-04 4.06e-05 -3.69e-05 5.00e-06 -2.51e-05 1.44e-04 4.87e-05 3.02e-04 1.83e-04 -6.43e-05 -9.65e-05 1.48e-04 -1.45e-04 -5.63e-04 1.59e-04 1.83e-04 9.44e-04
    5.08e-04 -4.90e-04 1.89e-04 -1.53e-04 -1.24e-04 -3.19e-05 -3.13e-04 4.16e-05 -4.90e-04 2.20e-03 -1.46e-04 1.47e-04 2.73e-04 7.28e-05 1.31e-03 1.65e-04 1.89e-04 -1.46e-04 1.66e-03 -4.34e-04 1.72e-04 8.78e-05 7.61e-04 1.68e-04 -1.53e-04 1.47e-04 -4.34e-04 6.16e-04 -2.25e-04 -6.56e-05 -3.96e-04 -1.43e-04 -1.24e-04 2.73e-04 1.72e-04 -2.25e-04 3.79e-04 2.92e-05 5.76e-04 4.03e-05 -3.19e-05 7.28e-05 8.78e-05 -6.56e-05 2.92e-05 5.44e-04 1.37e-04 -2.00e-04 -3.13e-04 1.31e-03 7.61e-04 -3.96e-04 5.76e-04 1.37e-04 2.57e-03 1.08e-04 4.16e-05 1.65e-04 1.68e-04 -1.43e-04 4.03e-05 -2.00e-04 1.08e-04 5.14e-04
    5.65e-04 -4.74e-04 5.69e-05 -1.61e-04 -3.56e-05 2.13e-05 -3.90e-04 1.04e-04 -4.74e-04 2.52e-03 1.15e-04 4.91e-05 3.72e-04 1.48e-04 1.42e-03 1.72e-04 5.69e-05 1.15e-04 9.77e-04 -2.72e-04 1.10e-04 3.14e-04 3.65e-05 5.36e-05 -1.61e-04 4.91e-05 -2.72e-04 7.54e-04 -1.44e-04 -1.86e-04 -9.04e-05 -1.05e-04 -3.56e-05 3.72e-04 1.10e-04 -1.44e-04 4.40e-04 7.11e-05 4.47e-04 1.16e-04 2.13e-05 1.48e-04 3.14e-04 -1.86e-04 7.11e-05 8.33e-04 -1.63e-05 -4.32e-05 -3.90e-04 1.42e-03 3.65e-05 -9.04e-05 4.47e-04 -1.63e-05 2.12e-03 -2.89e-05 1.04e-04 1.72e-04 5.36e-05 -1.05e-04 1.16e-04 -4.32e-05 -2.89e-05 5.19e-04
    5.36e-04 -5.23e-04 1.13e-04 -1.66e-04 -6.86e-05 -4.93e-05 -3.71e-04 1.11e-04 -5.23e-04 2.19e-03 -8.38e-05 1.34e-04 3.10e-04 -3.33e-05 1.40e-03 1.15e-04 1.13e-04 -8.38e-05 1.24e-03 -3.56e-04 1.59e-04 2.45e-04 2.33e-04 1.22e-04 -1.66e-04 1.34e-04 -3.56e-04 7.02e-04 -2.21e-04 -1.59e-04 -2.71e-04 -1.49e-04 -6.86e-05 3.10e-04 1.59e-04 -2.21e-04 3.62e-04 4.05e-05 5.34e-04 9.11e-05 -4.93e-05 -3.33e-05 2.45e-04 -1.59e-04 4.05e-05 6.31e-04 1.16e-04 -1.26e-04 -3.71e-04 1.40e-03 2.33e-04 -2.71e-04 5.34e-04 1.16e-04 2.24e-03 4.12e-05 1.11e-04 1.15e-04 1.22e-04 -1.49e-04 9.11e-05 -1.26e-04 4.12e-05 5.55e-04
    4.50e-04 -3.52e-04 -8.16e-05 -3.84e-05 2.85e-05 -2.89e-05 -2.27e-04 8.30e-05 -3.52e-04 3.04e-03 -1.07e-04 4.79e-06 5.18e-04 1.21e-03 8.33e-04 1.61e-05 -8.16e-05 -1.07e-04 3.65e-04 -2.85e-05 -6.62e-05 -7.80e-06 -7.75e-06 -1.74e-05 -3.84e-05 4.79e-06 -2.85e-05 4.85e-04 -1.95e-04 -1.17e-04 -1.24e-05 -9.93e-05 2.85e-05 5.18e-04 -6.62e-05 -1.95e-04 3.41e-04 2.12e-04 1.85e-04 1.09e-04 -2.89e-05 1.21e-03 -7.80e-06 -1.17e-04 2.12e-04 1.85e-03 -5.55e-04 3.43e-05 -2.27e-04 8.33e-04 -7.75e-06 -1.24e-05 1.85e-04 -5.55e-04 2.05e-03 -1.24e-04 8.30e-05 1.61e-05 -1.74e-05 -9.93e-05 1.09e-04 3.43e-05 -1.24e-04 3.87e-04
    5.83e-04 -3.26e-05 1.78e-04 -6.38e-05 3.49e-05 3.14e-05 3.07e-04 -1.81e-04 -3.26e-05 7.31e-04 1.64e-05 -1.99e-04 -8.43e-05 -2.43e-04 3.62e-04 -3.51e-04 1.78e-04 1.64e-05 2.06e-03 -2.55e-04 3.90e-04 -1.44e-04 -4.01e-04 -1.60e-04 -6.38e-05 -1.99e-04 -2.55e-04 6.90e-04 8.72e-05 9.01e-05 -9.57e-05 2.76e-04 3.49e-05 -8.43e-05 3.90e-04 8.72e-05 6.98e-04 -3.29e-05 -8.73e-06 1.25e-05 3.14e-05 -2.43e-04 -1.44e-04 9.01e-05 -3.29e-05 6.67e-04 -2.58e-04 2.64e-04 3.07e-04 3.62e-04 -4.01e-04 -9.57e-05 -8.73e-06 -2.58e-04 3.09e-03 -1.02e-03 -1.81e-04 -3.51e-04 -1.60e-04 2.76e-04 1.25e-05 2.64e-04 -1.02e-03 1.32e-03
    5.99e-04 -3.50e-04 1.55e-04 -1.20e-04 -1.56e-04 -1.03e-04 -9.95e-05 -1.77e-04 -3.50e-04 1.35e-03 5.48e-06 7.45e-05 1.36e-04 6.20e-05 3.27e-04 3.37e-04 1.55e-04 5.48e-06 2.20e-03 -4.60e-04 2.82e-04 -6.37e-05 1.47e-03 9.69e-05 -1.20e-04 7.45e-05 -4.60e-04 6.30e-04 -1.12e-04 1.18e-04 -4.51e-04 -9.02e-05 -1.56e-04 1.36e-04 2.82e-04 -1.12e-04 4.32e-04 9.47e-05 5.53e-04 8.13e-05 -1.03e-04 6.20e-05 -6.37e-05 1.18e-04 9.47e-05 5.75e-04 -4.44e-05 -1.43e-04 -9.95e-05 3.27e-04 1.47e-03 -4.51e-04 5.53e-04 -4.44e-05 2.49e-03 2.38e-04 -1.77e-04 3.37e-04 9.69e-05 -9.02e-05 8.13e-05 -1.43e-04 2.38e-04 6.79e-04
    8.57e-04 -3.04e-04 -2.56e-04 -1.59e-04 1.72e-04 2.86e-04 -2.10e-04 6.64e-05 -3.04e-04 2.10e-03 5.41e-05 -3.53e-06 4.41e-04 -1.14e-04 -3.21e-04 -3.50e-06 -2.56e-04 5.41e-05 6.29e-04 1.52e-05 -7.61e-05 -2.68e-04 3.10e-04 -1.24e-04 -1.59e-04 -3.53e-06 1.52e-05 5.29e-04 -4.63e-05 -1.05e-04 1.56e-04 4.58e-05 1.72e-04 4.41e-04 -7.61e-05 -4.63e-05 6.86e-04 2.16e-05 -1.75e-04 -7.65e-06 2.86e-04 -1.14e-04 -2.68e-04 -1.05e-04 2.16e-05 1.98e-03 -1.51e-03 2.23e-04 -2.10e-04 -3.21e-04 3.10e-04 1.56e-04 -1.75e-04 -1.51e-03 2.40e-03 -2.67e-04 6.64e-05 -3.50e-06 -1.24e-04 4.58e-05 -7.65e-06 2.23e-04 -2.67e-04 5.51e-04
    5.59e-04 -2.59e-04 1.16e-04 -1.16e-04 -1.45e-04 -1.08e-04 1.06e-05 -2.26e-04 -2.59e-04 8.45e-04 1.09e-04 3.05e-06 5.22e-05 -2.95e-05 -6.55e-05 3.25e-04 1.16e-04 1.09e-04 2.44e-03 -4.57e-04 3.89e-04 -5.28e-05 1.47e-03 3.02e-04 -1.16e-04 3.05e-06 -4.57e-04 5.52e-04 -5.88e-05 1.28e-04 -4.76e-04 -6.39e-05 -1.45e-04 5.22e-05 3.89e-04 -5.88e-05 3.96e-04 8.20e-05 3.73e-04 1.01e-04 -1.08e-04 -2.95e-05 -5.28e-05 1.28e-04 8.20e-05 5.21e-04 -1.91e-04 -6.82e-05 1.06e-05 -6.55e-05 1.47e-03 -4.76e-04 3.73e-04 -1.91e-04 1.98e-03 1.49e-04 -2.26e-04 3.25e-04 3.02e-04 -6.39e-05 1.01e-04 -6.82e-05 1.49e-04 7.52e-04
    4.93e-04 -4.12e-04 1.83e-04 -1.14e-04 -1.74e-04 -1.13e-04 -2.08e-04 -1.44e-04 -4.12e-04 1.51e-03 -1.92e-04 1.62e-04 1.89e-04 1.36e-04 5.84e-04 2.67e-04 1.83e-04 -1.92e-04 1.90e-03 -4.53e-04 1.97e-04 -7.78e-05 1.27e-03 9.19e-05 -1.14e-04 1.62e-04 -4.53e-04 4.62e-04 -1.17e-04 1.18e-04 -5.50e-04 -1.16e-04 -1.74e-04 1.89e-04 1.97e-04 -1.17e-04 2.89e-04 7.67e-05 6.26e-04 8.33e-05 -1.13e-04 1.36e-04 -7.78e-05 1.18e-04 7.67e-05 3.98e-04 2.00e-06 -1.80e-04 -2.08e-04 5.84e-04 1.27e-03 -5.50e-04 6.26e-04 2.00e-06 2.44e-03 1.77e-04 -1.44e-04 2.67e-04 9.19e-05 -1.16e-04 8.33e-05 -1.80e-04 1.77e-04 4.68e-04
    8.50e-04 -4.70e-05 -1.35e-04 -8.85e-05 -1.00e-04 -1.25e-04 1.82e-04 -1.06e-04 -4.70e-05 6.31e-04 -3.87e-05 1.66e-04 -2.84e-04 6.05e-05 -3.57e-05 -2.72e-04 -1.35e-04 -3.87e-05 6.28e-04 5.03e-05 7.02e-05 -1.42e-05 -1.05e-04 5.15e-05 -8.85e-05 1.66e-04 5.03e-05 1.42e-03 3.84e-04 1.77e-04 -2.33e-04 -9.87e-05 -1.00e-04 -2.84e-04 7.02e-05 3.84e-04 2.29e-03 1.53e-04 3.16e-05 1.22e-03 -1.25e-04 6.05e-05 -1.42e-05 1.77e-04 1.53e-04 5.64e-04 -8.66e-05 8.65e-05 1.82e-04 -3.57e-05 -1.05e-04 -2.33e-04 3.16e-05 -8.66e-05 6.86e-04 1.00e-04 -1.06e-04 -2.72e-04 5.15e-05 -9.87e-05 1.22e-03 8.65e-05 1.00e-04 1.89e-03
    7.91e-04 -3.62e-04 5.18e-04 2.13e-05 -2.68e-04 -1.32e-04 -3.33e-04 1.37e-04 -3.62e-04 1.95e-03 -1.19e-03 4.67e-04 -2.79e-05 2.48e-04 -9.59e-05 -1.45e-04 5.18e-04 -1.19e-03 2.05e-03 -3.42e-04 -1.20e-04 -1.85e-04 -1.31e-04 1.62e-04 2.13e-05 4.67e-04 -3.42e-04 9.42e-04 -1.52e-04 6.40e-05 -3.91e-04 -2.27e-04 -2.68e-04 -2.79e-05 -1.20e-04 -1.52e-04 5.21e-04 1.05e-04 4.51e-04 6.35e-05 -1.32e-04 2.48e-04 -1.85e-04 6.40e-05 1.05e-04 7.68e-04 2.02e-04 -2.15e-04 -3.33e-04 -9.59e-05 -1.31e-04 -3.91e-04 4.51e-04 2.02e-04 1.64e-03 1.19e-04 1.37e-04 -1.45e-04 1.62e-04 -2.27e-04 6.35e-05 -2.15e-04 1.19e-04 8.93e-04
    5.84e-04 -1.72e-04 -1.55e-05 -1.29e-04 9.02e-05 -2.16e-04 -5.39e-05 -1.20e-04 -1.72e-04 2.87e-04 1.94e-05 -1.67e-05 -2.95e-04 9.48e-05 4.64e-05 -7.65e-05 -1.55e-05 1.94e-05 2.65e-04 6.26e-05 -4.99e-05 -1.66e-04 -9.52e-05 5.33e-05 -1.29e-04 -1.67e-05 6.26e-05 5.38e-04 4.68e-04 3.26e-05 -1.31e-04 6.97e-06 9.02e-05 -2.95e-04 -4.99e-05 4.68e-04 3.24e-03 1.28e-04 1.24e-04 -2.67e-04 -2.16e-04 9.48e-05 -1.66e-04 3.26e-05 1.28e-04 6.95e-04 1.38e-04 1.84e-04 -5.39e-05 4.64e-05 -9.52e-05 -1.31e-04 1.24e-04 1.38e-04 3.86e-04 9.19e-05 -1.20e-04 -7.65e-05 5.33e-05 6.97e-06 -2.67e-04 1.84e-04 9.19e-05 1.33e-03
    4.79e-04 -5.05e-04 2.12e-04 -1.30e-04 -8.48e-05 -8.50e-05 -3.76e-04 1.59e-04 -5.05e-04 1.97e-03 -4.27e-04 2.01e-04 2.40e-04 -3.28e-05 1.19e-03 3.44e-05 2.12e-04 -4.27e-04 1.41e-03 -4.18e-04 1.18e-04 1.66e-04 3.04e-04 1.47e-04 -1.30e-04 2.01e-04 -4.18e-04 6.88e-04 -2.35e-04 -1.08e-04 -3.88e-04 -2.10e-04 -8.48e-05 2.40e-04 1.18e-04 -2.35e-04 3.05e-04 1.06e-04 6.09e-04 7.99e-05 -8.50e-05 -3.28e-05 1.66e-04 -1.08e-04 1.06e-04 5.82e-04 1.12e-04 -1.67e-04 -3.76e-04 1.19e-03 3.04e-04 -3.88e-04 6.09e-04 1.12e-04 2.35e-03 5.79e-05 1.59e-04 3.44e-05 1.47e-04 -2.10e-04 7.99e-05 -1.67e-04 5.79e-05 5.50e-04
    5.35e-04 -5.01e-04 2.33e-04 -1.43e-04 -1.60e-04 -1.11e-04 -3.28e-04 9.52e-07 -5.01e-04 1.92e-03 -4.79e-04 2.08e-04 2.37e-04 8.89e-05 9.66e-04 1.09e-04 2.33e-04 -4.79e-04 1.83e-03 -4.64e-04 1.62e-04 -4.54e-05 8.69e-04 1.40e-04 -1.43e-04 2.08e-04 -4.64e-04 6.05e-04 -1.68e-04 3.66e-05 -5.12e-04 -1.69e-04 -1.60e-04 2.37e-04 1.62e-04 -1.68e-04 3.16e-04 7.97e-05 6.77e-04 8.17e-05 -1.11e-04 8.89e-05 -4.54e-05 3.66e-05 7.97e-05 4.98e-04 6.09e-05 -2.22e-04 -3.28e-04 9.66e-04 8.69e-04 -5.12e-04 6.77e-04 6.09e-05 2.63e-03 1.19e-04 9.52e-07 1.09e-04 1.40e-04 -1.69e-04 8.17e-05 -2.22e-04 1.19e-04 5.58e-04
    3.06e-04 -3.00e-04 -5.47e-05 -3.77e-05 -2.49e-05 7.82e-05 -2.75e-04 -7.55e-06 -3.00e-04 2.54e-03 6.10e-05 3.08e-07 4.63e-04 4.96e-04 1.26e-03 1.16e-04 -5.47e-05 6.10e-05 2.44e-04 -9.27e-05 2.68e-05 1.06e-04 -7.89e-05 4.16e-05 -3.77e-05 3.08e-07 -9.27e-05 2.42e-04 -1.24e-04 -1.38e-04 4.50e-05 -2.29e-05 -2.49e-05 4.63e-04 2.68e-05 -1.24e-04 2.14e-04 5.26e-05 2.76e-04 6.23e-05 7.82e-05 4.96e-04 1.06e-04 -1.38e-04 5.26e-05 1.15e-03 -4.98e-04 -1.17e-05 -2.75e-04 1.26e-03 -7.89e-05 4.50e-05 2.76e-04 -4.98e-04 2.10e-03 3.04e-05 -7.55e-06 1.16e-04 4.16e-05 -2.29e-05 6.23e-05 -1.17e-05 3.04e-05 1.68e-04
    6.25e-04 -4.86e-04 4.39e-04 -6.71e-05 -2.04e-04 -1.82e-04 -3.74e-04 5.29e-05 -4.86e-04 1.71e-03 -1.34e-03 4.24e-04 8.34e-05 2.06e-04 9.28e-05 -1.10e-04 4.39e-04 -1.34e-03 1.77e-03 -4.31e-04 -3.75e-05 -1.24e-04 3.69e-05 1.49e-04 -6.71e-05 4.24e-04 -4.31e-04 7.30e-04 -2.09e-04 2.66e-05 -4.56e-04 -2.59e-04 -2.04e-04 8.34e-05 -3.75e-05 -2.09e-04 3.32e-04 1.50e-04 5.30e-04 9.02e-05 -1.82e-04 2.06e-04 -1.24e-04 2.66e-05 1.50e-04 6.55e-04 1.16e-04 -2.43e-04 -3.74e-04 9.28e-05 3.69e-05 -4.56e-04 5.30e-04 1.16e-04 1.84e-03 2.27e-05 5.29e-05 -1.10e-04 1.49e-04 -2.59e-04 9.02e-05 -2.43e-04 2.27e-05 7.23e-04
    3.12e-04 -5.71e-05 -4.15e-05 -4.12e-05 1.29e-04 2.31e-04 1.42e-04 -6.35e-05 -5.71e-05 5.64e-04 -3.37e-05 2.38e-05 1.26e-04 -8.57e-05 -3.15e-04 -7.52e-05 -4.15e-05 -3.37e-05 1.96e-04 4.67e-06 1.03e-05 -6.69e-05 7.36e-05 -3.41e-05 -4.12e-05 2.38e-05 4.67e-06 2.17e-04 6.35e-05 -1.92e-04 2.33e-04 8.78e-05 1.29e-04 1.26e-04 1.03e-05 6.35e-05 3.43e-04 -9.79e-05 9.75e-07 -1.71e-05 2.31e-04 -8.57e-05 -6.69e-05 -1.92e-04 -9.79e-05 1.45e-03 -8.33e-04 -1.62e-06 1.42e-04 -3.15e-04 7.36e-05 2.33e-04 9.75e-07 -8.33e-04 3.35e-03 1.17e-04 -6.35e-05 -7.52e-05 -3.41e-05 8.78e-05 -1.71e-05 -1.62e-06 1.17e-04 3.75e-04
    6.93e-04 -8.10e-05 -1.73e-04 -1.17e-04 2.28e-04 5.20e-04 -1.11e-04 -7.26e-05 -8.10e-05 1.39e-03 3.77e-05 2.48e-05 4.11e-04 -3.97e-04 -4.64e-04 -1.79e-04 -1.73e-04 3.77e-05 5.06e-04 2.17e-05 -1.62e-05 -3.69e-04 2.15e-04 -9.49e-05 -1.17e-04 2.48e-05 2.17e-05 3.72e-04 5.52e-05 -1.98e-04 6.74e-05 8.82e-05 2.28e-04 4.11e-04 -1.62e-05 5.52e-05 7.08e-04 -8.80e-05 -1.96e-04 -9.33e-05 5.20e-04 -3.97e-04 -3.69e-04 -1.98e-04 -8.80e-05 2.37e-03 -1.18e-03 1.30e-04 -1.11e-04 -4.64e-04 2.15e-04 6.74e-05 -1.96e-04 -1.18e-03 2.63e-03 -1.44e-04 -7.26e-05 -1.79e-04 -9.49e-05 8.82e-05 -9.33e-05 1.30e-04 -1.44e-04 4.87e-04
    5.43e-04 -5.34e-05 1.56e-04 -1.64e-04 1.46e-04 5.10e-05 3.23e-04 -1.79e-04 -5.34e-05 6.86e-04 7.77e-05 -1.82e-04 -1.04e-05 -2.82e-04 2.12e-04 -3.15e-04 1.56e-04 7.77e-05 1.50e-03 -2.35e-04 3.79e-04 -3.73e-05 -4.17e-04 -2.84e-04 -1.64e-04 -1.82e-04 -2.35e-04 5.82e-04 1.58e-04 3.04e-05 -1.77e-05 3.51e-04 1.46e-04 -1.04e-05 3.79e-04 1.58e-04 6.68e-04 -2.20e-05 -1.47e-05 -7.51e-05 5.10e-05 -2.82e-04 -3.73e-05 3.04e-05 -2.20e-05 6.51e-04 -1.19e-04 2.02e-04 3.23e-04 2.12e-04 -4.17e-04 -1.77e-05 -1.47e-05 -1.19e-04 3.20e-03 -7.56e-04 -1.79e-04 -3.15e-04 -2.84e-04 3.51e-04 -7.51e-05 2.02e-04 -7.56e-04 1.18e-03
    9.25e-04 2.19e-04 9.70e-06 -1.35e-04 5.31e-04 -1.67e-05 -2.02e-04 -3.89e-04 2.19e-04 6.57e-04 -3.17e-05 -2.89e-05 1.25e-04 1.84e-04 -2.32e-04 -3.32e-04 9.70e-06 -3.17e-05 2.57e-04 5.11e-05 -5.63e-05 -9.41e-05 1.36e-05 -2.42e-05 -1.35e-04 -2.89e-05 5.11e-05 2.35e-04 -1.11e-04 -1.22e-04 4.19e-05 6.48e-05 5.31e-04 1.25e-04 -5.63e-05 -1.11e-04 1.14e-03 -3.20e-04 -1.09e-04 -2.65e-04 -1.67e-05 1.84e-04 -9.41e-05 -1.22e-04 -3.20e-04 2.38e-03 1.99e-05 -2.24e-04 -2.02e-04 -2.32e-04 1.36e-05 4.19e-05 -1.09e-04 1.99e-05 6.83e-04 2.27e-04 -3.89e-04 -3.32e-04 -2.42e-05 6.48e-05 -2.65e-04 -2.24e-04 2.27e-04 6.56e-04
    5.86e-04 -6.11e-05 -5.99e-05 -1.85e-04 -2.36e-04 -1.75e-04 1.69e-04 -2.44e-04 -6.11e-05 3.24e-04 -3.51e-05 1.61e-04 -2.42e-04 9.10e-05 -7.72e-05 -3.31e-04 -5.99e-05 -3.51e-05 4.20e-04 8.32e-05 2.75e-05 -9.04e-05 -2.41e-04 5.28e-05 -1.85e-04 1.61e-04 8.32e-05 9.76e-04 4.27e-04 9.86e-05 -2.10e-04 -1.51e-04 -2.36e-04 -2.42e-04 2.75e-05 4.27e-04 2.15e-03 1.20e-04 -5.16e-05 1.00e-03 -1.75e-04 9.10e-05 -9.04e-05 9.86e-05 1.20e-04 5.05e-04 -2.75e-05 1.20e-04 1.69e-04 -7.72e-05 -2.41e-04 -2.10e-04 -5.16e-05 -2.75e-05 5.33e-04 9.73e-05 -2.44e-04 -3.31e-04 5.28e-05 -1.51e-04 1.00e-03 1.20e-04 9.73e-05 2.01e-03
    4.16e-04 -1.45e-04 -8.95e-06 -1.19e-04 -5.56e-05 -3.25e-04 3.24e-05 -1.68e-04 -1.45e-04 2.27e-04 -2.69e-05 2.27e-05 -3.17e-04 1.40e-04 -3.56e-05 -2.14e-04 -8.95e-06 -2.69e-05 1.68e-04 1.05e-04 3.19e-05 -9.52e-05 -1.44e-04 1.13e-04 -1.19e-04 2.27e-05 1.05e-04 4.52e-04 6.79e-04 9.71e-05 -1.03e-04 6.50e-05 -5.56e-05 -3.17e-04 3.19e-05 6.79e-04 2.99e-03 1.79e-04 5.07e-07 2.82e-04 -3.25e-04 1.40e-04 -9.52e-05 9.71e-05 1.79e-04 5.90e-04 8.87e-05 1.23e-04 3.24e-05 -3.56e-05 -1.44e-04 -1.03e-04 5.07e-07 8.87e-05 2.85e-04 5.69e-05 -1.68e-04 -2.14e-04 1.13e-04 6.50e-05 2.82e-04 1.23e-04 5.69e-05 1.45e-03
    4.96e-04 -3.79e-04 -7.14e-05 -5.77e-05 4.05e-06 4.06e-05 -3.59e-04 1.10e-04 -3.79e-04 2.77e-03 1.44e-05 2.34e-05 4.90e-04 7.69e-04 1.03e-03 8.83e-05 -7.14e-05 1.44e-05 3.96e-04 -6.82e-05 -2.48e-05 1.33e-04 -6.23e-05 -1.00e-05 -5.77e-05 2.34e-05 -6.82e-05 5.30e-04 -2.12e-04 -1.70e-04 9.26e-06 -1.22e-04 4.05e-06 4.90e-04 -2.48e-05 -2.12e-04 3.29e-04 1.23e-04 2.41e-04 1.21e-04 4.06e-05 7.69e-04 1.33e-04 -1.70e-04 1.23e-04 1.47e-03 -5.24e-04 5.59e-05 -3.59e-04 1.03e-03 -6.23e-05 9.26e-06 2.41e-04 -5.24e-04 2.14e-03 -1.24e-04 1.10e-04 8.83e-05 -1.00e-05 -1.22e-04 1.21e-04 5.59e-05 -1.24e-04 3.82e-04
    2.22e-03 -1.67e-04 1.28e-04 -2.85e-04 2.95e-04 1.38e-03 -1.64e-04 -5.64e-05 -1.67e-04 9.32e-04 5.16e-05 7.00e-06 -1.65e-04 9.59e-05 1.06e-04 -1.36e-04 1.28e-04 5.16e-05 7.69e-04 8.39e-05 9.75e-07 2.61e-05 -2.55e-05 -5.10e-05 -2.85e-04 7.00e-06 8.39e-05 6.86e-04 -5.67e-05 -2.77e-04 -1.50e-05 -4.93e-05 2.95e-04 -1.65e-04 9.75e-07 -5.67e-05 8.51e-04 2.97e-04 6.06e-05 4.13e-05 1.38e-03 9.59e-05 2.61e-05 -2.77e-04 2.97e-04 2.02e-03 -5.62e-05 -1.92e-04 -1.64e-04 1.06e-04 -2.55e-05 -1.50e-05 6.06e-05 -5.62e-05 5.32e-04 -1.06e-04 -5.64e-05 -1.36e-04 -5.10e-05 -4.93e-05 4.13e-05 -1.92e-04 -1.06e-04 5.45e-04
    1.12e-03 -5.09e-04 1.11e-04 -1.08e-04 -1.63e-04 1.27e-04 -1.64e-06 -9.34e-05 -5.09e-04 2.21e-03 -2.52e-04 7.95e-05 4.48e-04 1.29e-03 3.00e-05 -1.93e-04 1.11e-04 -2.52e-04 6.34e-04 1.33e-04 -2.09e-04 -2.96e-04 -1.98e-04 -2.55e-05 -1.08e-04 7.95e-05 1.33e-04 7.30e-04 -1.15e-04 -1.77e-04 -1.61e-04 -4.18e-05 -1.63e-04 4.48e-04 -2.09e-04 -1.15e-04 5.94e-04 4.16e-04 2.42e-04 8.36e-08 1.27e-04 1.29e-03 -2.96e-04 -1.77e-04 4.16e-04 2.13e-03 1.91e-05 -3.09e-04 -1.64e-06 3.00e-05 -1.98e-04 -1.61e-04 2.42e-04 1.91e-05 7.00e-04 3.25e-06 -9.34e-05 -1.93e-04 -2.55e-05 -4.18e-05 8.36e-08 -3.09e-04 3.25e-06 4.43e-04
    6.48e-04 -2.21e-05 -1.10e-04 1.63e-04 -2.27e-04 -1.13e-04 -8.69e-05 -7.71e-05 -2.21e-05 2.88e-03 2.58e-04 7.88e-04 -3.66e-04 1.19e-04 3.35e-04 -1.02e-04 -1.10e-04 2.58e-04 6.61e-04 -3.89e-05 1.66e-05 1.27e-04 8.57e-05 1.14e-04 1.63e-04 7.88e-04 -3.89e-05 8.13e-04 -2.49e-04 -1.17e-04 -2.83e-06 -1.79e-04 -2.27e-04 -3.66e-04 1.66e-05 -2.49e-04 3.64e-04 1.02e-04 9.41e-05 1.07e-04 -1.13e-04 1.19e-04 1.27e-04 -1.17e-04 1.02e-04 7.27e-04 2.18e-04 4.92e-06 -8.69e-05 3.35e-04 8.57e-05 -2.83e-06 9.41e-05 2.18e-04 6.75e-04 3.19e-05 -7.71e-05 -1.02e-04 1.14e-04 -1.79e-04 1.07e-04 4.92e-06 3.19e-05 4.23e-04
    5.90e-04 -1.79e-04 -5.93e-05 -1.75e-04 3.99e-05 -1.63e-04 -3.80e-05 -1.61e-04 -1.79e-04 2.98e-04 5.63e-05 -2.28e-05 -2.71e-04 7.02e-05 4.16e-05 -1.13e-04 -5.93e-05 5.63e-05 2.82e-04 7.87e-05 1.10e-05 -1.65e-04 -1.14e-04 2.43e-05 -1.75e-04 -2.28e-05 7.87e-05 5.22e-04 5.20e-04 1.11e-05 -1.58e-04 2.57e-05 3.99e-05 -2.71e-04 1.10e-05 5.20e-04 2.97e-03 1.29e-04 7.33e-05 -1.52e-05 -1.63e-04 7.02e-05 -1.65e-04 1.11e-05 1.29e-04 6.60e-04 9.63e-05 1.78e-04 -3.80e-05 4.16e-05 -1.14e-04 -1.58e-04 7.33e-05 9.63e-05 4.00e-04 8.77e-05 -1.61e-04 -1.13e-04 2.43e-05 2.57e-05 -1.52e-05 1.78e-04 8.77e-05 1.44e-03
    6.87e-04 -5.62e-04 3.45e-04 -1.16e-04 -2.00e-04 -1.46e-04 -3.68e-04 4.40e-05 -5.62e-04 1.86e-03 -1.20e-03 3.08e-04 1.42e-04 1.39e-04 3.12e-04 -8.13e-05 3.45e-04 -1.20e-03 1.87e-03 -5.34e-04 2.71e-05 -1.46e-04 2.16e-04 1.13e-04 -1.16e-04 3.08e-04 -5.34e-04 8.13e-04 -2.00e-04 6.66e-06 -4.49e-04 -1.99e-04 -2.00e-04 1.42e-04 2.71e-05 -2.00e-04 3.43e-04 1.32e-04 5.89e-04 9.28e-05 -1.46e-04 1.39e-04 -1.46e-04 6.66e-06 1.32e-04 6.90e-04 7.74e-05 -2.84e-04 -3.68e-04 3.12e-04 2.16e-04 -4.49e-04 5.89e-04 7.74e-05 2.18e-03 2.21e-05 4.40e-05 -8.13e-05 1.13e-04 -1.99e-04 9.28e-05 -2.84e-04 2.21e-05 7.61e-04
    1.18e-03 -3.61e-04 1.18e-04 -2.19e-04 -4.00e-05 2.83e-04 2.33e-05 -1.18e-04 -3.61e-04 2.06e-03 -2.36e-04 9.26e-05 3.95e-04 1.26e-03 7.85e-05 -2.70e-04 1.18e-04 -2.36e-04 5.50e-04 1.70e-04 -2.03e-04 -1.24e-04 -1.96e-04 -1.88e-05 -2.19e-04 9.26e-05 1.70e-04 6.23e-04 -1.91e-05 1.21e-05 -9.56e-05 1.35e-05 -4.00e-05 3.95e-04 -2.03e-04 -1.91e-05 5.38e-04 2.56e-04 2.87e-04 -1.19e-05 2.83e-04 1.26e-03 -1.24e-04 1.21e-05 2.56e-04 2.05e-03 4.36e-05 -3.22e-04 2.33e-05 7.85e-05 -1.96e-04 -9.56e-05 2.87e-04 4.36e-05 6.41e-04 -1.43e-05 -1.18e-04 -2.70e-04 -1.88e-05 1.35e-05 -1.19e-05 -3.22e-04 -1.43e-05 3.39e-04
    1.38e-03 -6.56e-04 1.74e-04 -1.27e-04 -1.59e-04 3.17e-04 -1.29e-04 -8.00e-05 -6.56e-04 2.12e-03 -2.54e-04 1.20e-04 4.53e-04 1.11e-03 1.15e-04 -1.96e-04 1.74e-04 -2.54e-04 6.67e-04 1.58e-04 -2.19e-04 -3.06e-04 -2.05e-04 -4.89e-05 -1.27e-04 1.20e-04 1.58e-04 6.63e-04 -1.49e-04 -2.71e-04 -9.90e-05 -7.46e-05 -1.59e-04 4.53e-04 -2.19e-04 -1.49e-04 5.93e-04 4.15e-04 2.88e-04 -3.49e-05 3.17e-04 1.11e-03 -3.06e-04 -2.71e-04 4.15e-04 2.17e-03 -3.67e-06 -3.46e-04 -1.29e-04 1.15e-04 -2.05e-04 -9.90e-05 2.88e-04 -3.67e-06 6.28e-04 2.85e-05 -8.00e-05 -1.96e-04 -4.89e-05 -7.46e-05 -3.49e-05 -3.46e-04 2.85e-05 4.39e-04
    1.79e-03 -9.78e-05 6.26e-05 5.65e-04 8.93e-04 2.14e-04 3.99e-05 9.24e-05 -9.78e-05 3.33e-04 -1.74e-04 1.85e-05 -4.05e-05 7.44e-07 6.91e-05 4.27e-05 6.26e-05 -1.74e-04 3.91e-04 -4.11e-05 -8.07e-05 3.95e-05 8.16e-05 9.61e-05 5.65e-04 1.85e-05 -4.11e-05 8.76e-04 5.61e-04 2.37e-04 -7.91e-05 -1.15e-04 8.93e-04 -4.05e-05 -8.07e-05 5.61e-04 1.50e-03 3.00e-04 4.96e-05 -1.71e-04 2.14e-04 7.44e-07 3.95e-05 2.37e-04 3.00e-04 6.59e-04 1.78e-04 1.29e-04 3.99e-05 6.91e-05 8.16e-05 -7.91e-05 4.96e-05 1.78e-04 3.85e-04 2.68e-04 9.24e-05 4.27e-05 9.61e-05 -1.15e-04 -1.71e-04 1.29e-04 2.68e-04 8.68e-04
    2.07e-03 -3.68e-04 1.36e-04 -2.75e-04 1.91e-04 1.21e-03 -2.76e-04 -6.61e-05 -3.68e-04 1.15e-03 2.28e-05 4.63e-05 5.24e-06 2.50e-04 1.63e-04 -1.88e-04 1.36e-04 2.28e-05 6.79e-04 1.07e-04 -5.10e-05 -5.00e-05 -9.40e-05 -5.73e-05 -2.75e-04 4.63e-05 1.07e-04 6.97e-04 -1.32e-04 -3.19e-04 1.83e-06 -2.57e-05 1.91e-04 5.24e-06 -5.10e-05 -1.32e-04 5.86e-04 2.84e-04 1.63e-04 -1.14e-05 1.21e-03 2.50e-04 -5.00e-05 -3.19e-04 2.84e-04 2.06e-03 -8.66e-05 -2.46e-04 -2.76e-04 1.63e-04 -9.40e-05 1.83e-06 1.63e-04 -8.66e-05 5.40e-04 -1.18e-04 -6.61e-05 -1.88e-04 -5.73e-05 -2.57e-05 -1.14e-05 -2.46e-04 -1.18e-04 4.91e-04
    6.72e-04 -1.70e-04 -3.93e-05 -2.97e-04 2.83e-05 2.12e-04 2.74e-04 -2.02e-04 -1.70e-04 1.02e-03 1.56e-04 -8.40e-05 2.18e-04 -3.11e-04 -4.14e-04 -7.99e-05 -3.93e-05 1.56e-04 9.73e-04 -1.93e-04 1.64e-04 1.42e-04 -5.32e-04 -1.81e-04 -2.97e-04 -8.40e-05 -1.93e-04 6.53e-04 6.36e-05 -1.26e-04 7.19e-05 2.27e-04 2.83e-05 2.18e-04 1.64e-04 6.36e-05 6.36e-04 -7.34e-05 -7.18e-05 -1.44e-04 2.12e-04 -3.11e-04 1.42e-04 -1.26e-04 -7.34e-05 9.08e-04 -1.02e-04 -7.84e-05 2.74e-04 -4.14e-04 -5.32e-04 7.19e-05 -7.18e-05 -1.02e-04 2.88e-03 -3.41e-04 -2.02e-04 -7.99e-05 -1.81e-04 2.27e-04 -1.44e-04 -7.84e-05 -3.41e-04 7.58e-04
    8.27e-04 -4.44e-04 -2.42e-04 -1.53e-04 1.39e-04 2.05e-04 -2.64e-04 8.16e-05 -4.44e-04 2.29e-03 1.04e-04 3.15e-05 4.56e-04 1.19e-04 -1.53e-04 8.29e-05 -2.42e-04 1.04e-04 5.87e-04 1.59e-05 -8.38e-05 -1.43e-04 2.46e-04 -1.14e-04 -1.53e-04 3.15e-05 1.59e-05 5.79e-04 -1.37e-04 -1.27e-04 1.84e-04 -4.54e-07 1.39e-04 4.56e-04 -8.38e-05 -1.37e-04 5.58e-04 2.23e-05 -1.54e-04 1.01e-04 2.05e-04 1.19e-04 -1.43e-04 -1.27e-04 2.23e-05 1.81e-03 -1.38e-03 2.02e-04 -2.64e-04 -1.53e-04 2.46e-04 1.84e-04 -1.54e-04 -1.38e-03 2.30e-03 -2.93e-04 8.16e-05 8.29e-05 -1.14e-04 -4.54e-07 1.01e-04 2.02e-04 -2.93e-04 4.84e-04
    4.81e-04 -1.43e-04 1.14e-04 -5.66e-05 -1.18e-04 -6.90e-05 8.85e-05 -2.12e-04 -1.43e-04 6.61e-04 -7.37e-06 -1.99e-05 -1.22e-05 -6.45e-05 -1.36e-04 1.16e-04 1.14e-04 -7.37e-06 2.63e-03 -4.02e-04 4.05e-04 -1.26e-04 1.12e-03 6.34e-04 -5.66e-05 -1.99e-05 -4.02e-04 6.09e-04 -4.22e-05 1.08e-04 -3.88e-04 -3.78e-05 -1.18e-04 -1.22e-05 4.05e-04 -4.22e-05 4.10e-04 7.96e-05 2.32e-04 1.24e-04 -6.90e-05 -6.45e-05 -1.26e-04 1.08e-04 7.96e-05 4.78e-04 -2.68e-04 -3.05e-05 8.85e-05 -1.36e-04 1.12e-03 -3.88e-04 2.32e-04 -2.68e-04 1.94e-03 -8.54e-05 -2.12e-04 1.16e-04 6.34e-04 -3.78e-05 1.24e-04 -3.05e-05 -8.54e-05 1.01e-03
    2.31e-03 -1.92e-06 2.14e-04 -2.39e-04 1.98e-04 1.28e-03 5.87e-05 -6.14e-05 -1.92e-06 5.46e-04 6.70e-05 -2.39e-05 -2.17e-04 -5.15e-05 -1.38e-06 -1.13e-05 2.14e-04 6.70e-05 5.26e-04 6.71e-05 -1.61e-04 9.07e-05 8.55e-05 -9.30e-05 -2.39e-04 -2.39e-05 6.71e-05 4.75e-04 -5.87e-05 -9.62e-05 3.34e-05 -7.49e-07 1.98e-04 -2.17e-04 -1.61e-04 -5.87e-05 1.49e-03 5.79e-04 1.64e-05 5.34e-05 1.28e-03 -5.15e-05 9.07e-05 -9.62e-05 5.79e-04 1.76e-03 3.44e-05 3.44e-05 5.87e-05 -1.38e-06 8.55e-05 3.34e-05 1.64e-05 3.44e-05 3.46e-04 -1.12e-04 -6.14e-05 -1.13e-05 -9.30e-05 -7.49e-07 5.34e-05 3.44e-05 -1.12e-04 3.77e-04
    4.90e-04 -8.02e-05 1.47e-04 -1.22e-04 1.11e-05 -2.23e-05 3.20e-04 -1.22e-04 -8.02e-05 7.06e-04 -3.60e-05 -2.23e-04 -1.21e-04 -2.34e-04 3.23e-04 -2.44e-04 1.47e-04 -3.60e-05 2.26e-03 -3.56e-04 3.84e-04 -1.01e-04 -2.02e-04 1.47e-04 -1.22e-04 -2.23e-04 -3.56e-04 6.67e-04 1.25e-04 1.16e-04 -1.27e-04 2.53e-04 1.11e-05 -1.21e-04 3.84e-04 1.25e-04 6.01e-04 6.86e-05 -4.90e-06 3.90e-05 -2.23e-05 -2.34e-04 -1.01e-04 1.16e-04 6.86e-05 5.64e-04 -3.50e-04 2.54e-04 3.20e-04 3.23e-04 -2.02e-04 -1.27e-04 -4.90e-06 -3.50e-04 2.50e-03 -1.05e-03 -1.22e-04 -2.44e-04 1.47e-04 2.53e-04 3.90e-05 2.54e-04 -1.05e-03 1.13e-03
    5.40e-04 -4.67e-04 2.37e-04 -9.20e-05 -1.93e-04 -1.60e-04 -2.79e-04 -1.31e-04 -4.67e-04 1.65e-03 -6.05e-04 2.22e-04 1.94e-04 1.65e-04 5.01e-04 1.75e-04 2.37e-04 -6.05e-04 1.78e-03 -4.37e-04 1.57e-04 -1.05e-04 9.71e-04 4.18e-05 -9.20e-05 2.22e-04 -4.37e-04 4.90e-04 -1.19e-04 1.55e-04 -5.42e-04 -1.56e-04 -1.93e-04 1.94e-04 1.57e-04 -1.19e-04 3.15e-04 9.00e-05 6.34e-04 1.05e-04 -1.60e-04 1.65e-04 -1.05e-04 1.55e-04 9.00e-05 4.57e-04 4.52e-06 -1.97e-04 -2.79e-04 5.01e-04 9.71e-04 -5.42e-04 6.34e-04 4.52e-06 2.46e-03 1.42e-04 -1.31e-04 1.75e-04 4.18e-05 -1.56e-04 1.05e-04 -1.97e-04 1.42e-04 5.84e-04
    3.87e-04 -8.75e-05 -2.61e-05 -1.76e-04 -1.49e-04 -2.62e-04 5.96e-05 -2.24e-04 -8.75e-05 2.49e-04 -1.11e-05 5.97e-05 -3.15e-04 1.19e-04 -5.69e-05 -1.92e-04 -2.61e-05 -1.11e-05 2.42e-04 1.51e-04 5.20e-05 -1.24e-04 -1.90e-04 6.51e-05 -1.76e-04 5.97e-05 1.51e-04 5.27e-04 5.49e-04 1.73e-05 -1.36e-04 7.95e-05 -1.49e-04 -3.15e-04 5.20e-05 5.49e-04 2.66e-03 1.41e-04 -6.09e-05 4.84e-04 -2.62e-04 1.19e-04 -1.24e-04 1.73e-05 1.41e-04 5.88e-04 6.66e-05 1.47e-04 5.96e-05 -5.69e-05 -1.90e-04 -1.36e-04 -6.09e-05 6.66e-05 3.58e-04 5.02e-05 -2.24e-04 -1.92e-04 6.51e-05 7.95e-05 4.84e-04 1.47e-04 5.02e-05 1.66e-03
    6.49e-04 -2.51e-04 6.11e-04 4.25e-05 -2.36e-04 -2.04e-04 -2.69e-04 1.26e-04 -2.51e-04 1.97e-03 -9.04e-04 5.25e-04 -8.85e-05 3.14e-04 -9.35e-05 -1.70e-04 6.11e-04 -9.04e-04 2.10e-03 -2.72e-04 -1.85e-04 -1.62e-04 -1.21e-04 1.94e-04 4.25e-05 5.25e-04 -2.72e-04 8.00e-04 -2.04e-04 6.47e-05 -2.57e-04 -2.87e-04 -2.36e-04 -8.85e-05 -1.85e-04 -2.04e-04 3.02e-04 1.28e-04 3.12e-04 6.36e-05 -2.04e-04 3.14e-04 -1.62e-04 6.47e-05 1.28e-04 6.60e-04 2.28e-04 -2.05e-04 -2.69e-04 -9.35e-05 -1.21e-04 -2.57e-04 3.12e-04 2.28e-04 1.22e-03 1.00e-04 1.26e-04 -1.70e-04 1.94e-04 -2.87e-04 6.36e-05 -2.05e-04 1.00e-04 8.46e-04
    4.10e-04 1.28e-05 2.21e-05 2.60e-05 -8.38e-05 -5.97e-05 -4.90e-05 -1.05e-04 1.28e-05 5.55e-04 -1.95e-04 -1.67e-05 -4.99e-05 3.80e-06 -1.62e-05 -1.27e-04 2.21e-05 -1.95e-04 2.62e-03 -2.59e-04 3.48e-04 -2.59e-04 1.99e-04 1.31e-03 2.60e-05 -1.67e-05 -2.59e-04 6.06e-04 2.17e-05 3.44e-05 -1.78e-05 -1.07e-04 -8.38e-05 -4.99e-05 3.48e-04 2.17e-05 4.25e-04 6.24e-05 7.69e-05 2.24e-04 -5.97e-05 3.80e-06 -2.59e-04 3.44e-05 6.24e-05 4.47e-04 -2.94e-04 -1.25e-04 -4.90e-05 -1.62e-05 1.99e-04 -1.78e-05 7.69e-05 -2.94e-04 1.67e-03 -2.80e-04 -1.05e-04 -1.27e-04 1.31e-03 -1.07e-04 2.24e-04 -1.25e-04 -2.80e-04 1.55e-03
    4.52e-04 -1.05e-04 8.77e-05 8.01e-06 -1.28e-04 -6.44e-05 7.76e-05 -1.96e-04 -1.05e-04 5.71e-04 -6.33e-05 -4.44e-05 -4.27e-05 -8.64e-05 -9.22e-05 -4.77e-06 8.77e-05 -6.33e-05 2.72e-03 -3.73e-04 4.21e-04 -1.81e-04 8.43e-04 8.89e-04 8.01e-06 -4.44e-05 -3.73e-04 5.37e-04 -1.27e-05 8.46e-05 -3.41e-04 -4.54e-05 -1.28e-04 -4.27e-05 4.21e-04 -1.27e-05 3.66e-04 7.18e-05 1.78e-04 1.36e-04 -6.44e-05 -8.64e-05 -1.81e-04 8.46e-05 7.18e-05 4.55e-04 -3.06e-04 -2.13e-05 7.76e-05 -9.22e-05 8.43e-04 -3.41e-04 1.78e-04 -3.06e-04 1.99e-03 -2.30e-04 -1.96e-04 -4.77e-06 8.89e-04 -4.54e-05 1.36e-04 -2.13e-05 -2.30e-04 1.20e-03
    5.30e-04 -2.06e-06 7.67e-05 3.92e-05 -1.57e-04 -1.04e-04 -1.53e-04 2.57e-05 -2.06e-06 2.70e-03 -2.79e-05 6.55e-04 -2.69e-04 2.93e-04 1.51e-04 -1.48e-04 7.67e-05 -2.79e-05 1.30e-03 -1.14e-04 2.05e-05 5.21e-05 5.14e-05 1.49e-04 3.92e-05 6.55e-04 -1.14e-04 8.37e-04 -1.61e-04 1.70e-05 -4.62e-05 -1.30e-04 -1.57e-04 -2.69e-04 2.05e-05 -1.61e-04 3.17e-04 7.00e-05 1.56e-04 8.45e-05 -1.04e-04 2.93e-04 5.21e-05 1.70e-05 7.00e-05 6.88e-04 2.63e-04 -1.08e-04 -1.53e-04 1.51e-04 5.14e-05 -4.62e-05 1.56e-04 2.63e-04 7.98e-04 5.15e-05 2.57e-05 -1.48e-04 1.49e-04 -1.30e-04 8.45e-05 -1.08e-04 5.15e-05 6.85e-04
    4.10e-04 -5.15e-04 2.81e-04 -1.01e-04 -9.29e-05 -1.00e-04 -3.80e-04 1.62e-04 -5.15e-04 1.86e-03 -7.29e-04 2.69e-04 2.04e-04 -2.62e-06 9.48e-04 6.50e-06 2.81e-04 -7.29e-04 1.58e-03 -4.63e-04 9.92e-05 5.03e-05 3.43e-04 1.71e-04 -1.01e-04 2.69e-04 -4.63e-04 6.70e-04 -2.54e-04 -1.05e-04 -4.45e-04 -2.27e-04 -9.29e-05 2.04e-04 9.92e-05 -2.54e-04 2.80e-04 1.24e-04 6.29e-04 7.85e-05 -1.00e-04 -2.62e-06 5.03e-05 -1.05e-04 1.24e-04 5.47e-04 1.07e-04 -1.79e-04 -3.80e-04 9.48e-04 3.43e-04 -4.45e-04 6.29e-04 1.07e-04 2.37e-03 3.81e-05 1.62e-04 6.50e-06 1.71e-04 -2.27e-04 7.85e-05 -1.79e-04 3.81e-05 5.34e-04
    5.67e-04 -5.70e-04 2.82e-04 -1.11e-04 -1.38e-04 -1.01e-04 -3.52e-04 1.05e-04 -5.70e-04 1.91e-03 -1.00e-03 2.30e-04 1.85e-04 1.14e-04 5.28e-04 -7.51e-06 2.82e-04 -1.00e-03 1.84e-03 -5.69e-04 7.52e-05 -1.35e-04 4.25e-04 1.15e-04 -1.11e-04 2.30e-04 -5.69e-04 7.37e-04 -2.05e-04 -2.63e-05 -4.22e-04 -1.81e-04 -1.38e-04 1.85e-04 7.52e-05 -2.05e-04 2.84e-04 1.18e-04 6.21e-04 8.96e-05 -1.01e-04 1.14e-04 -1.35e-04 -2.63e-05 1.18e-04 6.59e-04 5.74e-05 -2.70e-04 -3.52e-04 5.28e-04 4.25e-04 -4.22e-04 6.21e-04 5.74e-05 2.34e-03 3.88e-05 1.05e-04 -7.51e-06 1.15e-04 -1.81e-04 8.96e-05 -2.70e-04 3.88e-05 6.48e-04
    6.69e-04 -3.25e-04 -2.38e-04 -1.87e-04 7.33e-05 4.81e-04 -1.40e-04 1.00e-05 -3.25e-04 1.66e-03 1.16e-04 4.39e-05 4.19e-04 -3.84e-04 -4.91e-04 2.21e-05 -2.38e-04 1.16e-04 6.28e-04 -4.01e-05 8.95e-06 -1.30e-04 5.30e-05 -8.67e-05 -1.87e-04 4.39e-05 -4.01e-05 5.10e-04 4.41e-05 -2.43e-04 1.42e-04 1.36e-04 7.33e-05 4.19e-04 8.95e-06 4.41e-05 4.91e-04 -1.84e-04 -1.31e-04 -3.59e-05 4.81e-04 -3.84e-04 -1.30e-04 -2.43e-04 -1.84e-04 1.76e-03 -1.07e-03 6.42e-05 -1.40e-04 -4.91e-04 5.30e-05 1.42e-04 -1.31e-04 -1.07e-03 2.50e-03 -2.38e-04 1.00e-05 2.21e-05 -8.67e-05 1.36e-04 -3.59e-05 6.42e-05 -2.38e-04 4.80e-04
    5.11e-04 -3.88e-04 -6.17e-05 -6.02e-05 4.21e-05 7.78e-05 -3.13e-04 7.01e-05 -3.88e-04 2.74e-03 2.81e-05 1.44e-05 4.51e-04 7.84e-04 7.01e-04 9.62e-05 -6.17e-05 2.81e-05 4.13e-04 -3.49e-05 -5.29e-05 2.92e-05 3.84e-05 -5.44e-05 -6.02e-05 1.44e-05 -3.49e-05 4.74e-04 -1.88e-04 -1.26e-04 7.24e-05 -1.07e-04 4.21e-05 4.51e-04 -5.29e-05 -1.88e-04 3.10e-04 1.00e-04 1.30e-04 1.34e-04 7.78e-05 7.84e-04 2.92e-05 -1.26e-04 1.00e-04 1.63e-03 -8.77e-04 5.19e-05 -3.13e-04 7.01e-04 3.84e-05 7.24e-05 1.30e-04 -8.77e-04 2.01e-03 -1.36e-04 7.01e-05 9.62e-05 -5.44e-05 -1.07e-04 1.34e-04 5.19e-05 -1.36e-04 3.79e-04
    5.90e-04 2.24e-04 -1.74e-05 -8.21e-05 3.06e-04 2.25e-04 -6.42e-05 -2.82e-04 2.24e-04 6.94e-04 -2.01e-05 1.27e-05 2.64e-04 -9.29e-05 -2.34e-04 -3.36e-04 -1.74e-05 -2.01e-05 1.86e-04 5.41e-05 -2.90e-05 -7.46e-05 9.67e-05 -1.19e-05 -8.21e-05 1.27e-05 5.41e-05 1.57e-04 -3.61e-05 -1.50e-04 2.21e-05 1.53e-05 3.06e-04 2.64e-04 -2.90e-05 -3.61e-05 6.46e-04 -9.76e-05 -1.20e-04 -1.51e-04 2.25e-04 -9.29e-05 -7.46e-05 -1.50e-04 -9.76e-05 2.62e-03 1.87e-04 -1.55e-04 -6.42e-05 -2.34e-04 9.67e-05 2.21e-05 -1.20e-04 1.87e-04 7.77e-04 1.93e-04 -2.82e-04 -3.36e-04 -1.19e-05 1.53e-05 -1.51e-04 -1.55e-04 1.93e-04 4.43e-04
    2.36e-03 2.63e-05 2.52e-04 -2.02e-04 -7.52e-06 1.18e-03 5.22e-05 1.40e-05 2.63e-05 4.86e-04 5.71e-05 -4.59e-05 -1.61e-04 -6.19e-05 -7.50e-06 3.20e-05 2.52e-04 5.71e-05 5.34e-04 3.53e-05 -2.37e-04 2.11e-05 1.21e-04 -9.91e-05 -2.02e-04 -4.59e-05 3.53e-05 4.34e-04 -3.86e-05 -7.12e-05 2.62e-05 -2.79e-05 -7.52e-06 -1.61e-04 -2.37e-04 -3.86e-05 1.86e-03 7.24e-04 3.80e-05 3.40e-05 1.18e-03 -6.19e-05 2.11e-05 -7.12e-05 7.24e-04 1.71e-03 5.92e-05 6.00e-05 5.22e-05 -7.50e-06 1.21e-04 2.62e-05 3.80e-05 5.92e-05 2.31e-04 -1.08e-04 1.40e-05 3.20e-05 -9.91e-05 -2.79e-05 3.40e-05 6.00e-05 -1.08e-04 3.70e-04
    5.12e-04 -8.08e-05 -4.29e-05 5.04e-05 -1.50e-04 -1.86e-04 -8.56e-05 1.42e-04 -8.08e-05 2.40e-03 -9.24e-05 7.83e-04 -3.06e-04 1.60e-04 1.90e-04 -2.81e-04 -4.29e-05 -9.24e-05 1.30e-03 -2.07e-04 6.39e-05 1.11e-04 -3.95e-05 2.72e-04 5.04e-05 7.83e-04 -2.07e-04 8.82e-04 -2.68e-04 5.11e-06 -3.68e-05 -3.48e-04 -1.50e-04 -3.06e-04 6.39e-05 -2.68e-04 2.98e-04 1.32e-04 1.60e-04 1.17e-04 -1.86e-04 1.60e-04 1.11e-04 5.11e-06 1.32e-04 8.42e-04 4.03e-04 -1.16e-04 -8.56e-05 1.90e-04 -3.95e-05 -3.68e-05 1.60e-04 4.03e-04 1.08e-03 8.59e-05 1.42e-04 -2.81e-04 2.72e-04 -3.48e-04 1.17e-04 -1.16e-04 8.59e-05 8.24e-04
    7.03e-04 -4.33e-04 -2.51e-04 -2.22e-04 2.38e-05 3.64e-04 -3.10e-04 7.33e-05 -4.33e-04 1.79e-03 2.10e-04 4.06e-05 4.14e-04 -2.00e-04 -2.95e-04 1.77e-04 -2.51e-04 2.10e-04 6.27e-04 -4.99e-05 2.17e-05 -1.70e-05 7.24e-05 -1.10e-04 -2.22e-04 4.06e-05 -4.99e-05 5.51e-04 -7.94e-05 -2.17e-04 2.17e-04 7.73e-05 2.38e-05 4.14e-04 2.17e-05 -7.94e-05 3.75e-04 -1.36e-04 -1.49e-04 1.57e-05 3.64e-04 -2.00e-04 -1.70e-05 -2.17e-04 -1.36e-04 1.56e-03 -1.16e-03 8.35e-05 -3.10e-04 -2.95e-04 7.24e-05 2.17e-04 -1.49e-04 -1.16e-03 2.08e-03 -3.06e-04 7.33e-05 1.77e-04 -1.10e-04 7.73e-05 1.57e-05 8.35e-05 -3.06e-04 4.22e-04
    6.05e-04 -3.50e-04 -7.43e-05 -9.45e-05 7.60e-05 -1.80e-05 -1.53e-04 -4.29e-05 -3.50e-04 2.58e-03 -1.25e-04 -5.95e-06 4.43e-04 9.12e-04 4.00e-04 -2.63e-05 -7.43e-05 -1.25e-04 2.48e-04 1.65e-05 -1.33e-04 -1.06e-04 1.19e-04 -3.50e-05 -9.45e-05 -5.95e-06 1.65e-05 3.44e-04 -1.43e-04 5.62e-06 7.30e-05 -2.14e-05 7.60e-05 4.43e-04 -1.33e-04 -1.43e-04 3.83e-04 9.84e-05 4.92e-05 7.24e-05 -1.80e-05 9.12e-04 -1.06e-04 5.62e-06 9.84e-05 1.79e-03 -1.00e-03 6.89e-05 -1.53e-04 4.00e-04 1.19e-04 7.30e-05 4.92e-05 -1.00e-03 1.78e-03 -1.22e-04 -4.29e-05 -2.63e-05 -3.50e-05 -2.14e-05 7.24e-05 6.89e-05 -1.22e-04 2.85e-04
    5.88e-04 2.27e-05 8.28e-05 -2.95e-04 -3.07e-05 -4.06e-05 6.28e-05 -1.51e-04 2.27e-05 6.96e-04 2.66e-05 2.44e-04 2.49e-04 -6.85e-05 -5.50e-06 1.34e-04 8.28e-05 2.66e-05 1.01e-03 -3.19e-04 4.37e-05 -2.23e-04 2.14e-05 1.54e-04 -2.95e-04 2.44e-04 -3.19e-04 2.11e-03 3.57e-04 1.27e-05 -3.19e-04 1.14e-03 -3.07e-05 2.49e-04 4.37e-05 3.57e-04 8.69e-04 -9.00e-06 4.27e-05 2.84e-04 -4.06e-05 -6.85e-05 -2.23e-04 1.27e-05 -9.00e-06 5.13e-04 -8.98e-05 -1.50e-04 6.28e-05 -5.50e-06 2.14e-05 -3.19e-04 4.27e-05 -8.98e-05 4.85e-04 -2.50e-04 -1.51e-04 1.34e-04 1.54e-04 1.14e-03 2.84e-04 -1.50e-04 -2.50e-04 1.85e-03
    1.95e-03 -2.52e-04 2.58e-04 -9.53e-04 -1.67e-04 1.52e-04 8.55e-05 -3.21e-05 -2.52e-04 7.08e-04 -8.56e-05 2.25e-04 1.45e-04 -2.32e-04 3.22e-05 8.90e-05 2.58e-04 -8.56e-05 7.69e-04 -2.56e-04 1.38e-04 4.14e-05 5.31e-05 -2.18e-04 -9.53e-04 2.25e-04 -2.56e-04 1.87e-03 1.98e-04 -8.76e-05 -1.47e-05 3.23e-05 -1.67e-04 1.45e-04 1.38e-04 1.98e-04 1.01e-03 -2.15e-04 1.45e-04 -6.54e-05 1.52e-04 -2.32e-04 4.14e-05 -8.76e-05 -2.15e-04 8.25e-04 -8.10e-05 -1.41e-04 8.55e-05 3.22e-05 5.31e-05 -1.47e-05 1.45e-04 -8.10e-05 3.88e-04 1.87e-05 -3.21e-05 8.90e-05 -2.18e-04 3.23e-05 -6.54e-05 -1.41e-04 1.87e-05 8.19e-04
    7.68e-04 -1.59e-04 -2.37e-04 1.47e-04 -2.47e-04 -1.56e-04 -4.52e-05 -9.69e-05 -1.59e-04 2.35e-03 2.41e-04 7.33e-04 -3.08e-04 -4.48e-05 3.38e-04 -2.25e-04 -2.37e-04 2.41e-04 7.44e-04 -1.07e-04 7.06e-05 2.03e-04 4.49e-05 9.16e-05 1.47e-04 7.33e-04 -1.07e-04 1.01e-03 -2.95e-04 -1.34e-04 7.01e-05 -3.61e-04 -2.47e-04 -3.08e-04 7.06e-05 -2.95e-04 4.50e-04 1.42e-04 7.14e-05 2.36e-04 -1.56e-04 -4.48e-05 2.03e-04 -1.34e-04 1.42e-04 9.10e-04 3.25e-04 3.34e-05 -4.52e-05 3.38e-04 4.49e-05 7.01e-05 7.14e-05 3.25e-04 9.19e-04 7.31e-05 -9.69e-05 -2.25e-04 9.16e-05 -3.61e-04 2.36e-04 3.34e-05 7.31e-05 7.83e-04
    6.07e-04 -1.58e-04 5.64e-04 3.49e-05 -2.11e-04 -2.27e-04 -1.85e-04 2.14e-04 -1.58e-04 2.18e-03 -6.53e-04 6.04e-04 -1.44e-04 3.33e-04 -1.28e-05 -2.18e-04 5.64e-04 -6.53e-04 2.15e-03 -2.38e-04 -1.65e-04 -1.29e-04 -9.19e-05 2.73e-04 3.49e-05 6.04e-04 -2.38e-04 7.99e-04 -2.11e-04 6.64e-05 -1.86e-04 -3.01e-04 -2.11e-04 -1.44e-04 -1.65e-04 -2.11e-04 2.75e-04 1.18e-04 2.63e-04 4.09e-05 -2.27e-04 3.33e-04 -1.29e-04 6.64e-05 1.18e-04 6.98e-04 3.20e-04 -2.17e-04 -1.85e-04 -1.28e-05 -9.19e-05 -1.86e-04 2.63e-04 3.20e-04 1.11e-03 1.12e-04 2.14e-04 -2.18e-04 2.73e-04 -3.01e-04 4.09e-05 -2.17e-04 1.12e-04 8.34e-04
    4.97e-04 -5.25e-05 1.65e-04 2.22e-05 -1.32e-04 -1.86e-04 -9.85e-05 2.90e-04 -5.25e-05 2.39e-03 -2.42e-04 6.91e-04 -2.30e-04 3.18e-04 7.07e-05 -2.61e-04 1.65e-04 -2.42e-04 1.66e-03 -2.28e-04 2.20e-06 -3.65e-05 -1.33e-05 2.84e-04 2.22e-05 6.91e-04 -2.28e-04 8.43e-04 -2.02e-04 1.03e-04 -7.20e-05 -3.16e-04 -1.32e-04 -2.30e-04 2.20e-06 -2.02e-04 2.62e-04 1.13e-04 2.00e-04 4.05e-05 -1.86e-04 3.18e-04 -3.65e-05 1.03e-04 1.13e-04 8.13e-04 3.98e-04 -2.45e-04 -9.85e-05 7.07e-05 -1.33e-05 -7.20e-05 2.00e-04 3.98e-04 1.08e-03 1.24e-04 2.90e-04 -2.61e-04 2.84e-04 -3.16e-04 4.05e-05 -2.45e-04 1.24e-04 8.32e-04
    3.43e-04 -3.22e-05 7.96e-05 -2.83e-05 -7.50e-05 -6.43e-05 1.15e-04 -8.95e-05 -3.22e-05 5.01e-04 -1.38e-04 -6.04e-05 -7.19e-05 -6.61e-05 4.80e-05 -9.62e-05 7.96e-05 -1.38e-04 2.75e-03 -3.76e-04 3.97e-04 -1.91e-04 6.08e-04 9.81e-04 -2.83e-05 -6.04e-05 -3.76e-04 5.19e-04 3.85e-05 1.05e-04 -2.51e-04 -1.12e-05 -7.50e-05 -7.19e-05 3.97e-04 3.85e-05 3.45e-04 6.23e-05 1.10e-04 9.62e-05 -6.43e-05 -6.61e-05 -1.91e-04 1.05e-04 6.23e-05 4.31e-04 -3.27e-04 -3.42e-05 1.15e-04 4.80e-05 6.08e-04 -2.51e-04 1.10e-04 -3.27e-04 2.01e-03 -4.96e-04 -8.95e-05 -9.62e-05 9.81e-04 -1.12e-05 9.62e-05 -3.42e-05 -4.96e-04 1.29e-03
    1.95e-03 -4.96e-04 1.52e-04 -2.48e-04 7.70e-05 9.92e-04 -2.85e-04 -8.82e-05 -4.96e-04 1.38e-03 -4.61e-05 9.48e-05 1.29e-04 4.36e-04 2.03e-04 -2.26e-04 1.52e-04 -4.61e-05 6.44e-04 1.16e-04 -1.04e-04 -8.40e-05 -1.57e-04 -4.37e-05 -2.48e-04 9.48e-05 1.16e-04 6.98e-04 -1.39e-04 -2.96e-04 -1.45e-06 -3.02e-05 7.70e-05 1.29e-04 -1.04e-04 -1.39e-04 5.58e-04 2.72e-04 2.22e-04 -2.55e-05 9.92e-04 4.36e-04 -8.40e-05 -2.96e-04 2.72e-04 2.05e-03 -7.05e-05 -2.77e-04 -2.85e-04 2.03e-04 -1.57e-04 -1.45e-06 2.22e-04 -7.05e-05 5.79e-04 -8.68e-05 -8.82e-05 -2.26e-04 -4.37e-05 -3.02e-05 -2.55e-05 -2.77e-04 -8.68e-05 4.61e-04
    7.02e-04 -2.61e-04 6.98e-04 3.57e-05 -2.88e-04 -2.50e-04 -2.01e-04 1.17e-04 -2.61e-04 1.90e-03 -5.56e-04 3.83e-04 -2.06e-05 4.16e-04 -7.46e-05 -9.15e-05 6.98e-04 -5.56e-04 2.29e-03 -9.40e-05 -2.79e-04 -2.02e-04 -5.47e-05 1.64e-04 3.57e-05 3.83e-04 -9.40e-05 7.18e-04 -1.15e-04 1.74e-04 -9.58e-05 -3.00e-04 -2.88e-04 -2.06e-05 -2.79e-04 -1.15e-04 2.88e-04 1.26e-04 2.43e-04 3.66e-05 -2.50e-04 4.16e-04 -2.02e-04 1.74e-04 1.26e-04 7.09e-04 2.34e-04 -2.24e-04 -2.01e-04 -7.46e-05 -5.47e-05 -9.58e-05 2.43e-04 2.34e-04 1.08e-03 1.85e-04 1.17e-04 -9.15e-05 1.64e-04 -3.00e-04 3.66e-05 -2.24e-04 1.85e-04 9.35e-04
    1.01e-03 -1.18e-04 5.78e-05 2.22e-04 -3.42e-04 -1.27e-04 -1.06e-04 -1.18e-04 -1.18e-04 2.46e-03 3.12e-04 6.24e-04 -3.13e-04 7.93e-06 4.18e-04 -6.85e-05 5.78e-05 3.12e-04 4.66e-04 8.55e-05 -5.43e-05 4.06e-05 1.34e-04 -1.23e-05 2.22e-04 6.24e-04 8.55e-05 7.93e-04 -2.50e-04 -1.34e-04 9.91e-05 -1.66e-04 -3.42e-04 -3.13e-04 -5.43e-05 -2.50e-04 4.09e-04 1.12e-04 2.99e-05 9.34e-05 -1.27e-04 7.93e-06 4.06e-05 -1.34e-04 1.12e-04 7.82e-04 1.78e-04 4.02e-05 -1.06e-04 4.18e-04 1.34e-04 9.91e-05 2.99e-05 1.78e-04 4.45e-04 -8.47e-06 -1.18e-04 -6.85e-05 -1.23e-05 -1.66e-04 9.34e-05 4.02e-05 -8.47e-06 4.59e-04
    6.23e-04 -1.41e-04 1.82e-04 -2.68e-04 3.99e-05 6.25e-05 4.04e-04 -2.70e-04 -1.41e-04 7.47e-04 6.35e-05 -1.79e-04 1.90e-05 -2.95e-04 -5.27e-05 -9.72e-05 1.82e-04 6.35e-05 1.88e-03 -4.07e-04 3.97e-04 9.99e-05 -4.37e-04 -1.56e-04 -2.68e-04 -1.79e-04 -4.07e-04 6.01e-04 3.58e-05 9.75e-05 -1.14e-04 2.88e-04 3.99e-05 1.90e-05 3.97e-04 3.58e-05 5.05e-04 3.75e-05 -1.05e-05 -1.38e-04 6.25e-05 -2.95e-04 9.99e-05 9.75e-05 3.75e-05 6.02e-04 -1.68e-04 1.16e-04 4.04e-04 -5.27e-05 -4.37e-04 -1.14e-04 -1.05e-05 -1.68e-04 2.49e-03 -7.33e-04 -2.70e-04 -9.72e-05 -1.56e-04 2.88e-04 -1.38e-04 1.16e-04 -7.33e-04 9.00e-04
    1.93e-03 3.44e-05 1.61e-04 -2.82e-04 1.80e-04 1.25e-03 -3.16e-05 5.38e-06 3.44e-05 4.85e-04 6.57e-05 -1.23e-04 -2.16e-04 -1.79e-05 -2.29e-05 5.69e-06 1.61e-04 6.57e-05 6.39e-04 4.02e-05 -1.17e-04 1.26e-05 9.18e-05 -1.70e-04 -2.82e-04 -1.23e-04 4.02e-05 4.02e-04 1.34e-04 -1.44e-04 -3.00e-05 1.85e-06 1.80e-04 -2.16e-04 -1.17e-04 1.34e-04 1.04e-03 2.40e-04 -9.17e-06 9.82e-05 1.25e-03 -1.79e-05 1.26e-05 -1.44e-04 2.40e-04 1.63e-03 -4.16e-05 -5.49e-05 -3.16e-05 -2.29e-05 9.18e-05 -3.00e-05 -9.17e-06 -4.16e-05 4.26e-04 -1.61e-04 5.38e-06 5.69e-06 -1.70e-04 1.85e-06 9.82e-05 -5.49e-05 -1.61e-04 3.83e-04
    4.79e-04 -5.36e-04 2.40e-04 -1.20e-04 -1.05e-04 -1.24e-04 -3.24e-04 1.71e-04 -5.36e-04 1.89e-03 -7.67e-04 1.66e-04 2.00e-04 3.91e-05 6.50e-04 3.36e-05 2.40e-04 -7.67e-04 1.77e-03 -5.15e-04 8.05e-05 -1.28e-04 5.64e-04 5.15e-05 -1.20e-04 1.66e-04 -5.15e-04 7.55e-04 -2.40e-04 -7.45e-05 -4.19e-04 -2.24e-04 -1.05e-04 2.00e-04 8.05e-05 -2.40e-04 2.90e-04 8.94e-05 6.44e-04 6.44e-05 -1.24e-04 3.91e-05 -1.28e-04 -7.45e-05 8.94e-05 6.05e-04 9.05e-05 -2.75e-04 -3.24e-04 6.50e-04 5.64e-04 -4.19e-04 6.44e-04 9.05e-05 2.49e-03 1.06e-04 1.71e-04 3.36e-05 5.15e-05 -2.24e-04 6.44e-05 -2.75e-04 1.06e-04 6.22e-04
    1.02e-03 8.36e-05 -3.62e-04 3.48e-04 -1.22e-04 1.40e-04 1.39e-05 3.92e-04 8.36e-05 8.63e-04 8.18e-05 5.57e-04 -2.13e-04 -2.37e-04 1.14e-04 8.30e-05 -3.62e-04 8.18e-05 9.91e-04 -1.68e-04 2.97e-05 -2.74e-04 1.36e-04 -2.86e-04 3.48e-04 5.57e-04 -1.68e-04 2.64e-03 -3.25e-04 -2.16e-04 -8.59e-05 -1.37e-04 -1.22e-04 -2.13e-04 2.97e-05 -3.25e-04 6.70e-04 8.93e-05 -6.22e-06 -9.80e-05 1.40e-04 -2.37e-04 -2.74e-04 -2.16e-04 8.93e-05 9.05e-04 -7.93e-05 8.72e-05 1.39e-05 1.14e-04 1.36e-04 -8.59e-05 -6.22e-06 -7.93e-05 5.66e-04 1.11e-04 3.92e-04 8.30e-05 -2.86e-04 -1.37e-04 -9.80e-05 8.72e-05 1.11e-04 1.07e-03
    6.73e-04 -8.77e-05 1.19e-04 -2.77e-04 9.74e-05 1.49e-04 3.84e-04 -2.43e-04 -8.77e-05 8.39e-04 6.54e-05 -1.14e-04 6.23e-05 -2.94e-04 -1.44e-04 -6.70e-05 1.19e-04 6.54e-05 1.27e-03 -3.18e-04 2.62e-04 1.96e-04 -5.35e-04 -2.19e-04 -2.77e-04 -1.14e-04 -3.18e-04 6.81e-04 8.40e-05 -4.59e-05 1.98e-05 3.31e-04 9.74e-05 6.23e-05 2.62e-04 8.40e-05 6.23e-04 3.49e-05 2.88e-05 -1.73e-04 1.49e-04 -2.94e-04 1.96e-04 -4.59e-05 3.49e-05 6.40e-04 -5.29e-05 -2.76e-05 3.84e-04 -1.44e-04 -5.35e-04 1.98e-05 2.88e-05 -5.29e-05 2.67e-03 -4.63e-04 -2.43e-04 -6.70e-05 -2.19e-04 3.31e-04 -1.73e-04 -2.76e-05 -4.63e-04 8.27e-04
    5.20e-04 2.14e-05 -2.33e-05 -1.74e-04 -3.97e-06 -3.91e-05 1.08e-04 3.87e-05 2.14e-05 5.07e-04 6.32e-05 2.25e-04 2.47e-05 -1.03e-04 5.16e-05 1.87e-04 -2.33e-05 6.32e-05 5.30e-04 -3.57e-05 -4.51e-05 -6.89e-05 -8.02e-05 -1.53e-05 -1.74e-04 2.25e-04 -3.57e-05 1.97e-03 2.16e-04 2.93e-05 -1.70e-04 1.07e-03 -3.97e-06 2.47e-05 -4.51e-05 2.16e-04 1.50e-03 9.57e-05 7.93e-05 7.21e-04 -3.91e-05 -1.03e-04 -6.89e-05 2.93e-05 9.57e-05 3.83e-04 -1.03e-04 5.51e-05 1.08e-04 5.16e-05 -8.02e-05 -1.70e-04 7.93e-05 -1.03e-04 4.73e-04 7.28e-06 3.87e-05 1.87e-04 -1.53e-05 1.07e-03 7.21e-04 5.51e-05 7.28e-06 1.76e-03
    7.89e-04 1.08e-05 -8.02e-05 -4.56e-05 3.05e-04 3.94e-04 -2.99e-05 -5.44e-05 1.08e-05 1.27e-03 -3.28e-06 2.09e-05 4.33e-04 -4.08e-04 -3.76e-04 -2.75e-04 -8.02e-05 -3.28e-06 4.52e-04 4.35e-05 -1.22e-04 -4.05e-04 2.70e-04 -4.05e-05 -4.56e-05 2.09e-05 4.35e-05 3.63e-04 7.19e-05 -1.62e-04 3.67e-05 7.99e-05 3.05e-04 4.33e-04 -1.22e-04 7.19e-05 7.55e-04 -3.02e-06 -2.39e-04 -1.15e-04 3.94e-04 -4.08e-04 -4.05e-04 -1.62e-04 -3.02e-06 2.20e-03 -6.28e-04 8.39e-05 -2.99e-05 -3.76e-04 2.70e-04 3.67e-05 -2.39e-04 -6.28e-04 1.86e-03 2.10e-05 -5.44e-05 -2.75e-04 -4.05e-05 7.99e-05 -1.15e-04 8.39e-05 2.10e-05 4.22e-04
    1.87e-03 -2.38e-04 3.11e-04 -1.11e-03 -3.67e-04 1.38e-04 4.27e-05 7.57e-05 -2.38e-04 5.19e-04 -1.19e-04 1.76e-04 8.33e-05 -2.35e-04 5.67e-05 1.34e-04 3.11e-04 -1.19e-04 6.11e-04 -2.50e-04 7.38e-05 -3.04e-07 1.03e-04 -2.60e-04 -1.11e-03 1.76e-04 -2.50e-04 1.49e-03 1.01e-04 5.00e-05 -9.26e-05 1.58e-05 -3.67e-04 8.33e-05 7.38e-05 1.01e-04 1.21e-03 -2.45e-04 1.49e-04 -7.59e-05 1.38e-04 -2.35e-04 -3.04e-07 5.00e-05 -2.45e-04 7.39e-04 -1.16e-04 -1.36e-04 4.27e-05 5.67e-05 1.03e-04 -9.26e-05 1.49e-04 -1.16e-04 2.06e-04 -2.50e-05 7.57e-05 1.34e-04 -2.60e-04 1.58e-05 -7.59e-05 -1.36e-04 -2.50e-05 6.85e-04
    9.87e-04 -3.03e-04 -1.10e-04 2.39e-04 -3.63e-04 -2.07e-04 -7.80e-05 -1.84e-04 -3.03e-04 2.23e-03 3.12e-04 6.75e-04 -2.76e-04 -1.77e-04 3.69e-04 -1.15e-04 -1.10e-04 3.12e-04 6.09e-04 2.81e-05 -2.66e-05 1.06e-04 6.37e-05 6.59e-06 2.39e-04 6.75e-04 2.81e-05 1.00e-03 -3.65e-04 -2.80e-04 6.49e-05 -3.41e-04 -3.63e-04 -2.76e-04 -2.66e-05 -3.65e-04 5.44e-04 2.01e-04 6.39e-05 2.31e-04 -2.07e-04 -1.77e-04 1.06e-04 -2.80e-04 2.01e-04 9.96e-04 2.38e-04 1.63e-04 -7.80e-05 3.69e-04 6.37e-05 6.49e-05 6.39e-05 2.38e-04 8.59e-04 6.22e-05 -1.84e-04 -1.15e-04 6.59e-06 -3.41e-04 2.31e-04 1.63e-04 6.22e-05 6.91e-04
    4.39e-04 1.04e-04 1.12e-05 -2.17e-05 -4.88e-05 -5.44e-06 -1.74e-04 -1.08e-04 1.04e-04 5.96e-04 -2.51e-04 1.43e-04 -5.52e-05 2.96e-05 -6.05e-05 -1.44e-04 1.12e-05 -2.51e-04 2.06e-03 -2.79e-04 3.13e-04 -2.09e-04 -3.08e-04 1.17e-03 -2.17e-05 1.43e-04 -2.79e-04 8.12e-04 -2.93e-06 -8.10e-05 1.75e-04 -7.90e-05 -4.88e-05 -5.52e-05 3.13e-04 -2.93e-06 4.73e-04 4.43e-05 2.14e-06 2.60e-04 -5.44e-06 2.96e-05 -2.09e-04 -8.10e-05 4.43e-05 4.18e-04 -1.07e-04 -2.52e-04 -1.74e-04 -6.05e-05 -3.08e-04 1.75e-04 2.14e-06 -1.07e-04 8.91e-04 -1.08e-04 -1.08e-04 -1.44e-04 1.17e-03 -7.90e-05 2.60e-04 -2.52e-04 -1.08e-04 1.60e-03
    1.84e-03 -6.79e-04 2.29e-04 -2.39e-04 8.45e-05 7.85e-04 -3.64e-04 2.79e-05 -6.79e-04 1.29e-03 2.09e-05 1.17e-04 3.68e-05 1.40e-04 2.29e-04 -2.13e-04 2.29e-04 2.09e-05 6.12e-04 1.59e-04 -8.60e-05 -3.66e-04 -9.62e-05 -6.75e-05 -2.39e-04 1.17e-04 1.59e-04 5.79e-04 -7.13e-05 -4.89e-04 -7.14e-06 -7.06e-06 8.45e-05 3.68e-05 -8.60e-05 -7.13e-05 5.90e-04 4.13e-04 8.72e-05 -7.36e-06 7.85e-04 1.40e-04 -3.66e-04 -4.89e-04 4.13e-04 1.83e-03 -9.85e-05 -1.81e-04 -3.64e-04 2.29e-04 -9.62e-05 -7.14e-06 8.72e-05 -9.85e-05 4.68e-04 -6.20e-05 2.79e-05 -2.13e-04 -6.75e-05 -7.06e-06 -7.36e-06 -1.81e-04 -6.20e-05 3.35e-04
    1.24e-03 -2.62e-05 -3.49e-04 2.45e-04 -6.57e-05 2.44e-04 9.46e-06 3.58e-04 -2.62e-05 8.03e-04 -3.51e-05 4.88e-04 -1.49e-04 -2.85e-04 1.77e-04 1.69e-04 -3.49e-04 -3.51e-05 7.73e-04 -6.69e-05 4.00e-05 -7.52e-05 7.47e-05 -2.53e-04 2.45e-04 4.88e-04 -6.69e-05 2.57e-03 -1.94e-04 -2.58e-04 -5.99e-05 -7.42e-05 -6.57e-05 -1.49e-04 4.00e-05 -1.94e-04 6.11e-04 2.92e-05 -8.95e-06 -1.37e-04 2.44e-04 -2.85e-04 -7.52e-05 -2.58e-04 2.92e-05 8.80e-04 -8.13e-05 -9.38e-05 9.46e-06 1.77e-04 7.47e-05 -5.99e-05 -8.95e-06 -8.13e-05 3.66e-04 1.55e-04 3.58e-04 1.69e-04 -2.53e-04 -7.42e-05 -1.37e-04 -9.38e-05 1.55e-04 1.04e-03
    2.11e-04 3.14e-05 2.98e-06 -1.12e-04 -7.90e-05 -3.93e-05 2.25e-05 -8.90e-05 3.14e-05 2.35e-04 2.32e-06 6.41e-05 -1.46e-05 -5.12e-06 -5.10e-05 -6.57e-05 2.98e-06 2.32e-06 5.43e-04 1.65e-04 1.81e-04 -2.69e-04 -3.34e-04 -9.20e-05 -1.12e-04 6.41e-05 1.65e-04 7.29e-04 4.21e-04 -2.52e-04 -1.25e-04 1.01e-04 -7.90e-05 -1.46e-05 1.81e-04 4.21e-04 1.08e-03 -1.78e-04 -1.80e-04 -3.74e-04 -3.93e-05 -5.12e-06 -2.69e-04 -2.52e-04 -1.78e-04 5.08e-04 1.53e-04 -1.63e-04 2.25e-05 -5.10e-05 -3.34e-04 -1.25e-04 -1.80e-04 1.53e-04 6.33e-04 1.33e-04 -8.90e-05 -6.57e-05 -9.20e-05 1.01e-04 -3.74e-04 -1.63e-04 1.33e-04 2.03e-03
    7.82e-04 -2.45e-04 7.40e-04 5.75e-05 -3.32e-04 -2.80e-04 -1.16e-04 8.61e-05 -2.45e-04 1.46e-03 -2.59e-04 4.05e-05 5.23e-05 3.56e-04 -6.04e-05 1.60e-05 7.40e-04 -2.59e-04 2.40e-03 2.94e-05 -3.52e-04 -2.45e-04 1.72e-06 7.26e-05 5.75e-05 4.05e-05 2.94e-05 6.29e-04 -6.82e-05 1.43e-04 1.14e-05 -2.69e-04 -3.32e-04 5.23e-05 -3.52e-04 -6.82e-05 3.37e-04 1.72e-04 1.58e-04 8.68e-06 -2.80e-04 3.56e-04 -2.45e-04 1.43e-04 1.72e-04 7.63e-04 1.75e-04 -1.57e-04 -1.16e-04 -6.04e-05 1.72e-06 1.14e-05 1.58e-04 1.75e-04 9.24e-04 2.02e-04 8.61e-05 1.60e-05 7.26e-05 -2.69e-04 8.68e-06 -1.57e-04 2.02e-04 1.01e-03
    6.47e-04 -2.24e-04 -3.26e-04 7.06e-05 -1.63e-04 -1.27e-04 -6.93e-05 1.95e-05 -2.24e-04 2.18e-03 1.07e-04 5.08e-04 -1.85e-04 1.04e-04 3.06e-04 -3.74e-04 -3.26e-04 1.07e-04 9.03e-04 -1.65e-04 1.28e-04 1.27e-04 4.20e-05 1.18e-04 7.06e-05 5.08e-04 -1.65e-04 7.64e-04 -1.85e-04 1.10e-04 1.42e-04 -3.58e-04 -1.63e-04 -1.85e-04 1.28e-04 -1.85e-04 2.81e-04 7.59e-05 9.43e-05 1.78e-04 -1.27e-04 1.04e-04 1.27e-04 1.10e-04 7.59e-05 8.58e-04 4.07e-04 -2.06e-04 -6.93e-05 3.06e-04 4.20e-05 1.42e-04 9.43e-05 4.07e-04 8.52e-04 3.21e-05 1.95e-05 -3.74e-04 1.18e-04 -3.58e-04 1.78e-04 -2.06e-04 3.21e-05 8.41e-04
    1.23e-03 -8.19e-04 1.80e-04 -3.84e-05 -2.90e-04 1.22e-05 -1.40e-04 -6.56e-05 -8.19e-04 1.99e-03 -2.14e-04 1.50e-04 3.97e-04 6.90e-04 -1.18e-05 -1.45e-05 1.80e-04 -2.14e-04 6.51e-04 1.62e-04 -2.19e-04 -4.23e-04 -1.89e-04 -6.57e-05 -3.84e-05 1.50e-04 1.62e-04 8.10e-04 -2.22e-04 -4.46e-04 -1.83e-04 -1.19e-04 -2.90e-04 3.97e-04 -2.19e-04 -2.22e-04 5.76e-04 4.92e-04 2.27e-04 4.85e-06 1.22e-05 6.90e-04 -4.23e-04 -4.46e-04 4.92e-04 1.82e-03 -6.31e-05 -1.53e-04 -1.40e-04 -1.18e-05 -1.89e-04 -1.83e-04 2.27e-04 -6.31e-05 7.80e-04 6.52e-05 -6.56e-05 -1.45e-05 -6.57e-05 -1.19e-04 4.85e-06 -1.53e-04 6.52e-05 4.71e-04
    7.20e-04 7.07e-05 -1.23e-04 3.29e-04 -2.35e-04 1.06e-04 1.33e-05 2.48e-04 7.07e-05 7.80e-04 1.32e-04 6.45e-04 -2.26e-04 -1.82e-04 1.74e-04 -1.17e-04 -1.23e-04 1.32e-04 1.24e-03 -3.54e-04 1.51e-05 -3.37e-04 1.59e-04 -3.81e-04 3.29e-04 6.45e-04 -3.54e-04 2.53e-03 -3.49e-04 -1.12e-04 -6.29e-05 -2.34e-04 -2.35e-04 -2.26e-04 1.51e-05 -3.49e-04 5.23e-04 5.69e-05 -3.14e-05 -2.93e-05 1.06e-04 -1.82e-04 -3.37e-04 -1.12e-04 5.69e-05 6.15e-04 -8.60e-05 2.45e-04 1.33e-05 1.74e-04 1.59e-04 -6.29e-05 -3.14e-05 -8.60e-05 5.39e-04 5.89e-05 2.48e-04 -1.17e-04 -3.81e-04 -2.34e-04 -2.93e-05 2.45e-04 5.89e-05 1.07e-03
    1.85e-03 -6.07e-04 1.48e-04 -2.37e-04 -1.22e-05 8.28e-04 -3.12e-04 -9.69e-05 -6.07e-04 1.67e-03 -1.18e-04 1.17e-04 2.38e-04 5.96e-04 2.17e-04 -2.19e-04 1.48e-04 -1.18e-04 7.09e-04 1.10e-04 -1.19e-04 -1.65e-04 -1.84e-04 -2.19e-05 -2.37e-04 1.17e-04 1.10e-04 6.85e-04 -1.30e-04 -3.10e-04 -1.81e-05 -5.14e-05 -1.22e-05 2.38e-04 -1.19e-04 -1.30e-04 5.21e-04 3.32e-04 2.65e-04 -3.66e-05 8.28e-04 5.96e-04 -1.65e-04 -3.10e-04 3.32e-04 2.09e-03 -6.22e-05 -3.17e-04 -3.12e-04 2.17e-04 -1.84e-04 -1.81e-05 2.65e-04 -6.22e-05 5.83e-04 -5.07e-05 -9.69e-05 -2.19e-04 -2.19e-05 -5.14e-05 -3.66e-05 -3.17e-04 -5.07e-05 4.67e-04
    1.06e-03 -3.45e-04 -2.12e-04 2.19e-04 -3.62e-04 -1.30e-04 -1.12e-04 -2.00e-04 -3.45e-04 2.14e-03 2.72e-04 4.45e-04 -1.84e-04 -5.20e-05 4.26e-04 -3.32e-04 -2.12e-04 2.72e-04 6.29e-04 -3.74e-05 1.79e-05 2.05e-04 8.27e-05 -5.57e-05 2.19e-04 4.45e-04 -3.74e-05 8.56e-04 -2.41e-04 -1.56e-07 1.82e-04 -3.62e-04 -3.62e-04 -1.84e-04 1.79e-05 -2.41e-04 4.56e-04 8.94e-05 5.23e-05 2.82e-04 -1.30e-04 -5.20e-05 2.05e-04 -1.56e-07 8.94e-05 9.40e-04 2.82e-04 2.23e-05 -1.12e-04 4.26e-04 8.27e-05 1.82e-04 5.23e-05 2.82e-04 7.78e-04 -1.25e-05 -2.00e-04 -3.32e-04 -5.57e-05 -3.62e-04 2.82e-04 2.23e-05 -1.25e-05 7.67e-04
    4.37e-04 -1.99e-06 -7.88e-05 -1.29e-04 1.03e-04 -2.09e-05 3.90e-05 -3.31e-06 -1.99e-06 5.31e-04 1.74e-04 1.17e-04 1.62e-04 -1.74e-04 9.21e-05 1.10e-04 -7.88e-05 1.74e-04 6.12e-04 -8.85e-05 2.37e-05 -1.22e-04 -5.18e-06 3.83e-05 -1.29e-04 1.17e-04 -8.85e-05 1.73e-03 3.42e-04 5.85e-05 -2.48e-04 1.09e-03 1.03e-04 1.62e-04 2.37e-05 3.42e-04 1.06e-03 3.91e-05 6.13e-05 2.69e-04 -2.09e-05 -1.74e-04 -1.22e-04 5.85e-05 3.91e-05 4.21e-04 -1.53e-04 -1.22e-05 3.90e-05 9.21e-05 -5.18e-06 -2.48e-04 6.13e-05 -1.53e-04 4.30e-04 -1.36e-04 -3.31e-06 1.10e-04 3.83e-05 1.09e-03 2.69e-04 -1.22e-05 -1.36e-04 1.60e-03
    5.53e-04 -1.54e-04 3.54e-04 -2.23e-05 -1.74e-04 -2.20e-04 -8.75e-05 3.03e-04 -1.54e-04 2.15e-03 -2.63e-04 4.52e-04 -1.13e-04 3.96e-04 1.20e-04 -2.56e-04 3.54e-04 -2.63e-04 2.01e-03 -1.49e-04 -1.29e-04 -1.61e-04 -1.83e-05 2.97e-04 -2.23e-05 4.52e-04 -1.49e-04 6.10e-04 -8.69e-05 2.37e-04 4.04e-05 -2.64e-04 -1.74e-04 -1.13e-04 -1.29e-04 -8.69e-05 2.09e-04 1.07e-04 1.39e-04 -8.03e-06 -2.20e-04 3.96e-04 -1.61e-04 2.37e-04 1.07e-04 7.85e-04 3.53e-04 -3.11e-04 -8.75e-05 1.20e-04 -1.83e-05 4.04e-05 1.39e-04 3.53e-04 8.74e-04 1.14e-04 3.03e-04 -2.56e-04 2.97e-04 -2.64e-04 -8.03e-06 -3.11e-04 1.14e-04 8.88e-04
    2.32e-04 -4.13e-05 5.35e-05 -5.42e-05 -4.39e-05 -3.37e-05 7.87e-05 -1.35e-05 -4.13e-05 4.58e-04 -1.59e-04 -6.04e-05 -1.17e-04 -3.85e-05 1.87e-04 -1.46e-04 5.35e-05 -1.59e-04 2.49e-03 -3.31e-04 3.29e-04 -2.43e-04 1.92e-04 9.86e-04 -5.42e-05 -6.04e-05 -3.31e-04 5.08e-04 3.77e-05 2.38e-05 -5.46e-05 -1.41e-05 -4.39e-05 -1.17e-04 3.29e-04 3.77e-05 2.59e-04 7.03e-05 5.56e-05 8.25e-05 -3.37e-05 -3.85e-05 -2.43e-04 2.38e-05 7.03e-05 3.31e-04 -2.97e-04 2.16e-05 7.87e-05 1.87e-04 1.92e-04 -5.46e-05 5.56e-05 -2.97e-04 1.84e-03 -6.73e-04 -1.35e-05 -1.46e-04 9.86e-04 -1.41e-05 8.25e-05 2.16e-05 -6.73e-04 1.23e-03
    1.83e-03 -1.16e-03 3.56e-04 -1.95e-05 -3.75e-04 -8.06e-05 -3.62e-04 6.77e-06 -1.16e-03 1.92e-03 -6.44e-05 2.76e-04 1.70e-04 -1.06e-04 2.56e-04 -4.59e-05 3.56e-04 -6.44e-05 7.51e-04 1.97e-04 -2.48e-04 -5.34e-04 -7.87e-05 -1.08e-04 -1.95e-05 2.76e-04 1.97e-04 8.57e-04 -2.38e-04 -5.64e-04 -2.27e-05 -1.04e-04 -3.75e-04 1.70e-04 -2.48e-04 -2.38e-04 6.82e-04 4.63e-04 1.33e-04 2.09e-05 -8.06e-05 -1.06e-04 -5.34e-04 -5.64e-04 4.63e-04 1.49e-03 -8.58e-05 2.11e-05 -3.62e-04 2.56e-04 -7.87e-05 -2.27e-05 1.33e-04 -8.58e-05 6.94e-04 3.53e-05 6.77e-06 -4.59e-05 -1.08e-04 -1.04e-04 2.09e-05 2.11e-05 3.53e-05 5.31e-04
    4.19e-04 -5.41e-05 7.06e-04 1.25e-04 -1.94e-04 -9.35e-05 -2.44e-05 -2.78e-05 -5.41e-05 1.17e-03 -6.95e-05 7.91e-05 -5.27e-05 2.61e-04 -1.98e-05 4.73e-05 7.06e-04 -6.95e-05 2.61e-03 -6.99e-06 -3.52e-04 -1.41e-04 1.07e-04 7.70e-05 1.25e-04 7.91e-05 -6.99e-06 5.72e-04 -1.16e-04 -4.39e-05 5.49e-06 -1.96e-04 -1.94e-04 -5.27e-05 -3.52e-04 -1.16e-04 1.51e-04 3.11e-05 4.65e-05 1.10e-05 -9.35e-05 2.61e-04 -1.41e-04 -4.39e-05 3.11e-05 2.43e-04 -5.00e-06 -1.95e-05 -2.44e-05 -1.98e-05 1.07e-04 5.49e-06 4.65e-05 -5.00e-06 3.87e-04 -3.75e-05 -2.78e-05 4.73e-05 7.70e-05 -1.96e-04 1.10e-05 -1.95e-05 -3.75e-05 6.31e-04
    4.92e-04 7.78e-06 1.20e-04 -3.14e-04 -3.62e-05 -2.86e-05 5.16e-05 -2.74e-04 7.78e-06 7.03e-04 -1.24e-05 3.08e-04 1.62e-04 -4.66e-05 -6.54e-05 -3.56e-06 1.20e-04 -1.24e-05 1.23e-03 -5.19e-04 1.56e-04 -2.68e-04 2.29e-05 2.41e-04 -3.14e-04 3.08e-04 -5.19e-04 2.04e-03 2.70e-04 -6.73e-06 -3.23e-04 9.26e-04 -3.62e-05 1.62e-04 1.56e-04 2.70e-04 6.13e-04 -6.16e-05 2.55e-05 2.57e-04 -2.86e-05 -4.66e-05 -2.68e-04 -6.73e-06 -6.16e-05 4.76e-04 -6.85e-05 -1.70e-04 5.16e-05 -6.54e-05 2.29e-05 -3.23e-04 2.55e-05 -6.85e-05 5.15e-04 -2.95e-04 -2.74e-04 -3.56e-06 2.41e-04 9.26e-04 2.57e-04 -1.70e-04 -2.95e-04 1.80e-03
    5.59e-04 -1.31e-04 1.09e-04 -1.08e-04 -5.26e-05 -6.67e-05 3.35e-04 -1.72e-04 -1.31e-04 6.91e-04 -5.36e-05 -2.20e-04 -1.35e-04 -2.32e-04 1.83e-04 -1.81e-04 1.09e-04 -5.36e-05 2.45e-03 -3.80e-04 4.18e-04 -3.24e-05 -9.64e-05 3.87e-04 -1.08e-04 -2.20e-04 -3.80e-04 6.95e-04 9.72e-05 1.68e-04 -1.48e-04 1.73e-04 -5.26e-05 -1.35e-04 4.18e-04 9.72e-05 6.14e-04 1.15e-04 -2.95e-06 6.47e-05 -6.67e-05 -2.32e-04 -3.24e-05 1.68e-04 1.15e-04 6.08e-04 -3.70e-04 2.11e-04 3.35e-04 1.83e-04 -9.64e-05 -1.48e-04 -2.95e-06 -3.70e-04 2.14e-03 -8.91e-04 -1.72e-04 -1.81e-04 3.87e-04 1.73e-04 6.47e-05 2.11e-04 -8.91e-04 1.07e-03
    2.39e-04 3.94e-05 1.40e-04 -6.81e-05 9.13e-05 -3.01e-05 1.46e-04 -6.71e-05 3.94e-05 4.69e-04 3.21e-05 -1.33e-04 -1.02e-04 -1.83e-04 4.24e-04 -3.81e-04 1.40e-04 3.21e-05 1.39e-03 -1.12e-04 4.23e-04 -3.36e-04 -3.95e-04 -2.97e-04 -6.81e-05 -1.33e-04 -1.12e-04 5.41e-04 2.39e-04 -1.36e-05 6.84e-05 3.38e-04 9.13e-05 -1.02e-04 4.23e-04 2.39e-04 6.16e-04 -1.20e-04 -1.28e-04 4.16e-05 -3.01e-05 -1.83e-04 -3.36e-04 -1.36e-05 -1.20e-04 4.48e-04 -1.23e-04 2.56e-04 1.46e-04 4.24e-04 -3.95e-04 6.84e-05 -1.28e-04 -1.23e-04 2.46e-03 -5.18e-04 -6.71e-05 -3.81e-04 -2.97e-04 3.38e-04 4.16e-05 2.56e-04 -5.18e-04 1.13e-03
    1.94e-03 -2.46e-04 2.24e-04 -6.58e-04 -8.99e-05 2.65e-04 9.67e-05 1.72e-05 -2.46e-04 6.41e-04 -6.38e-05 3.49e-04 7.33e-05 -2.56e-04 1.50e-04 6.91e-05 2.24e-04 -6.38e-05 5.57e-04 -9.08e-05 1.22e-04 1.11e-04 6.08e-05 -2.92e-04 -6.58e-04 3.49e-04 -9.08e-05 2.07e-03 9.75e-05 -1.55e-04 2.07e-05 1.02e-04 -8.99e-05 7.33e-05 1.22e-04 9.75e-05 6.93e-04 -1.65e-04 3.77e-05 -9.03e-05 2.65e-04 -2.56e-04 1.11e-04 -1.55e-04 -1.65e-04 6.83e-04 -2.52e-05 -1.75e-04 9.67e-05 1.50e-04 6.08e-05 2.07e-05 3.77e-05 -2.52e-05 1.74e-04 1.97e-05 1.72e-05 6.91e-05 -2.92e-04 1.02e-04 -9.03e-05 -1.75e-04 1.97e-05 7.23e-04
    8.73e-04 1.00e-04 -4.05e-04 3.85e-04 -2.15e-04 9.48e-05 -7.91e-05 3.32e-04 1.00e-04 6.47e-04 2.15e-05 1.16e-04 -2.03e-04 -2.15e-04 1.34e-04 2.34e-04 -4.05e-04 2.15e-05 9.95e-04 -3.13e-04 4.78e-05 -2.30e-04 1.51e-04 -2.31e-04 3.85e-04 1.16e-04 -3.13e-04 2.29e-03 -3.41e-04 -2.98e-04 -2.56e-04 -1.16e-04 -2.15e-04 -2.03e-04 4.78e-05 -3.41e-04 4.46e-04 1.27e-04 1.32e-05 -1.14e-04 9.48e-05 -2.15e-04 -2.30e-04 -2.98e-04 1.27e-04 8.43e-04 -1.06e-04 -1.21e-04 -7.91e-05 1.34e-04 1.51e-04 -2.56e-04 1.32e-05 -1.06e-04 3.31e-04 1.21e-04 3.32e-04 2.34e-04 -2.31e-04 -1.16e-04 -1.14e-04 -1.21e-04 1.21e-04 1.02e-03
    6.36e-04 -1.92e-04 1.67e-04 -2.40e-04 -4.69e-05 -2.19e-05 4.08e-04 -2.55e-04 -1.92e-04 7.23e-04 6.45e-05 -1.90e-04 -5.14e-05 -2.72e-04 1.17e-07 3.00e-06 1.67e-04 6.45e-05 2.08e-03 -4.38e-04 4.02e-04 1.37e-04 -1.38e-04 6.32e-05 -2.40e-04 -1.90e-04 -4.38e-04 5.91e-04 -1.00e-05 1.68e-04 -2.67e-04 2.11e-04 -4.69e-05 -5.14e-05 4.02e-04 -1.00e-05 4.24e-04 9.17e-05 -4.37e-06 -6.03e-05 -2.19e-05 -2.72e-04 1.37e-04 1.68e-04 9.17e-05 5.72e-04 -2.48e-04 1.53e-04 4.08e-04 1.17e-07 -1.38e-04 -2.67e-04 -4.37e-06 -2.48e-04 1.86e-03 -8.50e-04 -2.55e-04 3.00e-06 6.32e-05 2.11e-04 -6.03e-05 1.53e-04 -8.50e-04 8.70e-04
    7.97e-04 -1.99e-04 -4.77e-05 -6.88e-06 1.68e-04 -1.48e-04 2.95e-05 1.24e-04 -1.99e-04 3.32e-04 -4.97e-05 6.82e-05 -5.23e-05 1.20e-04 8.03e-05 -2.33e-04 -4.77e-05 -4.97e-05 2.36e-04 3.16e-05 9.64e-05 -6.15e-05 1.51e-05 1.06e-04 -6.88e-06 6.82e-05 3.16e-05 9.92e-04 8.41e-04 2.22e-04 -1.70e-04 -1.17e-04 1.68e-04 -5.23e-05 9.64e-05 8.41e-04 2.09e-03 3.19e-04 4.63e-05 3.23e-04 -1.48e-04 1.20e-04 -6.15e-05 2.22e-04 3.19e-04 5.24e-04 1.61e-05 -7.58e-06 2.95e-05 8.03e-05 1.51e-05 -1.70e-04 4.63e-05 1.61e-05 2.71e-04 1.50e-04 1.24e-04 -2.33e-04 1.06e-04 -1.17e-04 3.23e-04 -7.58e-06 1.50e-04 1.07e-03
    5.96e-04 -4.87e-04 -2.16e-05 -1.71e-04 -2.34e-05 2.05e-05 -3.92e-04 5.99e-05 -4.87e-04 2.24e-03 2.38e-04 3.84e-05 3.73e-04 1.59e-04 8.59e-04 1.91e-04 -2.16e-05 2.38e-04 7.02e-04 -1.69e-04 6.65e-05 2.78e-04 -2.45e-05 -2.13e-05 -1.71e-04 3.84e-05 -1.69e-04 6.48e-04 -1.93e-04 -1.94e-04 3.87e-05 -1.06e-04 -2.34e-05 3.73e-04 6.65e-05 -1.93e-04 3.02e-04 2.57e-05 1.95e-04 1.06e-04 2.05e-05 1.59e-04 2.78e-04 -1.94e-04 2.57e-05 8.40e-04 -2.17e-04 -1.41e-05 -3.92e-04 8.59e-04 -2.45e-05 3.87e-05 1.95e-04 -2.17e-04 1.49e-03 -1.48e-04 5.99e-05 1.91e-04 -2.13e-05 -1.06e-04 1.06e-04 -1.41e-05 -1.48e-04 4.52e-04
    9.93e-04 1.32e-04 -1.93e-07 -1.41e-04 5.53e-04 -2.49e-04 -1.94e-04 -3.25e-04 1.32e-04 5.40e-04 -2.70e-05 -5.56e-05 -6.21e-06 1.52e-04 -2.62e-04 -2.62e-04 -1.93e-07 -2.70e-05 2.52e-04 1.03e-04 -1.35e-04 -6.65e-05 -2.16e-05 -3.56e-05 -1.41e-04 -5.56e-05 1.03e-04 2.09e-04 -1.61e-04 -1.15e-04 7.98e-05 7.93e-05 5.53e-04 -6.21e-06 -1.35e-04 -1.61e-04 1.46e-03 -3.09e-04 -8.68e-05 -2.44e-04 -2.49e-04 1.52e-04 -6.65e-05 -1.15e-04 -3.09e-04 1.70e-03 -1.66e-05 -7.37e-05 -1.94e-04 -2.62e-04 -2.16e-05 7.98e-05 -8.68e-05 -1.66e-05 5.73e-04 1.70e-04 -3.25e-04 -2.62e-04 -3.56e-05 7.93e-05 -2.44e-04 -7.37e-05 1.70e-04 6.26e-04
    1.23e-03 -4.47e-04 8.55e-05 2.60e-04 -4.60e-04 -2.67e-04 -1.18e-04 -1.83e-04 -4.47e-04 2.15e-03 2.61e-04 6.29e-04 -2.15e-04 -2.49e-04 3.61e-04 -1.50e-05 8.55e-05 2.61e-04 5.52e-04 1.36e-04 -1.44e-04 -1.62e-05 6.77e-05 -2.86e-05 2.60e-04 6.29e-04 1.36e-04 9.82e-04 -3.90e-04 -4.04e-04 2.81e-05 -2.99e-04 -4.60e-04 -2.15e-04 -1.44e-04 -3.90e-04 5.45e-04 3.00e-04 6.25e-05 1.90e-04 -2.67e-04 -2.49e-04 -1.62e-05 -4.04e-04 3.00e-04 1.03e-03 1.57e-04 2.08e-04 -1.18e-04 3.61e-04 6.77e-05 2.81e-05 6.25e-05 1.57e-04 8.01e-04 4.51e-05 -1.83e-04 -1.50e-05 -2.86e-05 -2.99e-04 1.90e-04 2.08e-04 4.51e-05 5.90e-04
    8.07e-04 -2.73e-05 -1.05e-04 -6.38e-06 -5.16e-05 -9.16e-05 2.16e-04 1.38e-04 -2.73e-05 5.12e-04 -6.99e-05 2.50e-04 -2.63e-04 -2.32e-05 6.86e-05 -2.18e-04 -1.05e-04 -6.99e-05 3.65e-04 3.39e-05 1.39e-04 6.34e-05 -8.22e-05 1.10e-04 -6.38e-06 2.50e-04 3.39e-05 1.45e-03 4.68e-05 3.19e-04 -2.54e-04 -5.33e-05 -5.16e-05 -2.63e-04 1.39e-04 4.68e-05 1.64e-03 1.06e-04 2.77e-05 9.44e-04 -9.16e-05 -2.32e-05 6.34e-05 3.19e-04 1.06e-04 4.37e-04 -1.06e-04 -2.16e-05 2.16e-04 6.86e-05 -8.22e-05 -2.54e-04 2.77e-05 -1.06e-04 4.04e-04 1.29e-04 1.38e-04 -2.18e-04 1.10e-04 -5.33e-05 9.44e-04 -2.16e-05 1.29e-04 1.46e-03
    4.80e-04 -3.78e-04 -1.03e-04 -9.30e-05 1.03e-05 -4.61e-05 -3.24e-04 1.67e-04 -3.78e-04 1.91e-03 2.00e-05 5.66e-05 3.30e-04 3.49e-04 9.06e-04 5.70e-05 -1.03e-04 2.00e-05 4.34e-04 -1.57e-04 2.28e-05 2.51e-04 -1.18e-04 2.81e-05 -9.30e-05 5.66e-05 -1.57e-04 6.86e-04 -2.64e-04 -2.47e-04 -4.99e-05 -2.22e-04 1.03e-05 3.30e-04 2.28e-05 -2.64e-04 2.97e-04 1.33e-04 2.50e-04 1.64e-04 -4.61e-05 3.49e-04 2.51e-04 -2.47e-04 1.33e-04 8.19e-04 -1.59e-04 6.15e-05 -3.24e-04 9.06e-04 -1.18e-04 -4.99e-05 2.50e-04 -1.59e-04 1.68e-03 -1.69e-04 1.67e-04 5.70e-05 2.81e-05 -2.22e-04 1.64e-04 6.15e-05 -1.69e-04 4.52e-04
    1.58e-03 -7.93e-04 2.07e-04 -1.95e-04 -1.13e-04 4.69e-04 -2.46e-04 -1.03e-04 -7.93e-04 1.84e-03 -2.03e-04 1.39e-04 3.90e-04 7.65e-04 1.68e-04 -1.91e-04 2.07e-04 -2.03e-04 5.86e-04 2.02e-04 -2.31e-04 -2.55e-04 -2.43e-04 -4.13e-05 -1.95e-04 1.39e-04 2.02e-04 6.23e-04 -1.65e-04 -3.43e-04 -1.23e-05 -7.44e-05 -1.13e-04 3.90e-04 -2.31e-04 -1.65e-04 5.68e-04 3.53e-04 3.09e-04 -3.93e-05 4.69e-04 7.65e-04 -2.55e-04 -3.43e-04 3.53e-04 2.04e-03 -5.35e-05 -3.07e-04 -2.46e-04 1.68e-04 -2.43e-04 -1.23e-05 3.09e-04 -5.35e-05 5.81e-04 3.34e-05 -1.03e-04 -1.91e-04 -4.13e-05 -7.44e-05 -3.93e-05 -3.07e-04 3.34e-05 4.20e-04
    5.49e-04 -2.83e-05 -3.75e-05 -1.28e-04 -5.09e-05 -3.57e-05 1.13e-04 4.33e-05 -2.83e-05 4.60e-04 3.83e-05 2.18e-04 -9.88e-05 -7.09e-05 5.72e-05 9.56e-05 -3.75e-05 3.83e-05 4.72e-04 -1.74e-05 -2.95e-06 -2.63e-05 -1.58e-04 -3.88e-05 -1.28e-04 2.18e-04 -1.74e-05 1.93e-03 1.58e-04 1.14e-04 -1.40e-04 7.90e-04 -5.09e-05 -9.88e-05 -2.95e-06 1.58e-04 1.72e-03 1.03e-04 7.46e-05 9.17e-04 -3.57e-05 -7.09e-05 -2.63e-05 1.14e-04 1.03e-04 3.62e-04 -1.06e-04 1.54e-04 1.13e-04 5.72e-05 -1.58e-04 -1.40e-04 7.46e-05 -1.06e-04 4.83e-04 5.99e-05 4.33e-05 9.56e-05 -3.88e-05 7.90e-04 9.17e-04 1.54e-04 5.99e-05 1.65e-03
    1.99e-03 6.98e-05 3.70e-04 -2.29e-04 -9.60e-04 3.83e-04 -4.23e-05 2.33e-04 6.98e-05 4.34e-04 -7.28e-05 -6.39e-05 4.68e-05 3.18e-05 1.32e-05 1.23e-04 3.70e-04 -7.28e-05 6.29e-04 -1.14e-04 -2.70e-04 -1.41e-04 1.41e-04 -1.17e-04 -2.29e-04 -6.39e-05 -1.14e-04 6.55e-04 2.53e-04 1.74e-04 -5.48e-06 -7.14e-05 -9.60e-04 4.68e-05 -2.70e-04 2.53e-04 1.71e-03 1.72e-04 6.55e-05 -1.77e-05 3.83e-04 3.18e-05 -1.41e-04 1.74e-04 1.72e-04 9.90e-04 2.34e-05 -3.93e-05 -4.23e-05 1.32e-05 1.41e-04 -5.48e-06 6.55e-05 2.34e-05 2.54e-04 -9.07e-05 2.33e-04 1.23e-04 -1.17e-04 -7.14e-05 -1.77e-05 -3.93e-05 -9.07e-05 5.01e-04
    4.98e-04 -2.64e-04 -3.83e-05 4.43e-06 -3.80e-05 -2.12e-04 1.16e-04 1.43e-05 -2.64e-04 1.72e-03 -1.78e-04 -1.55e-05 3.67e-04 1.01e-03 -2.59e-05 -4.65e-05 -3.83e-05 -1.78e-04 4.64e-04 8.26e-05 -1.44e-04 -8.67e-05 -1.90e-04 6.61e-06 4.43e-06 -1.55e-05 8.26e-05 7.64e-04 -2.93e-04 -2.92e-04 -2.43e-04 -1.83e-04 -3.80e-05 3.67e-04 -1.44e-04 -2.93e-04 4.51e-04 4.01e-04 1.77e-04 1.36e-04 -2.12e-04 1.01e-03 -8.67e-05 -2.92e-04 4.01e-04 1.53e-03 -5.26e-05 1.59e-05 1.16e-04 -2.59e-05 -1.90e-04 -2.43e-04 1.77e-04 -5.26e-05 1.25e-03 -9.27e-05 1.43e-05 -4.65e-05 6.61e-06 -1.83e-04 1.36e-04 1.59e-05 -9.27e-05 4.97e-04
    6.73e-04 1.53e-04 -2.75e-04 3.90e-04 -2.64e-04 6.05e-07 -1.42e-05 3.19e-04 1.53e-04 6.56e-04 6.92e-05 2.56e-04 -2.48e-04 -2.47e-04 1.29e-04 9.95e-05 -2.75e-04 6.92e-05 1.20e-03 -4.37e-04 5.69e-05 -3.38e-04 2.05e-04 -3.00e-04 3.90e-04 2.56e-04 -4.37e-04 2.33e-03 -3.56e-04 -2.50e-04 -2.43e-04 -2.37e-04 -2.64e-04 -2.48e-04 5.69e-05 -3.56e-04 5.03e-04 1.42e-04 -7.01e-07 -8.28e-05 6.05e-07 -2.47e-04 -3.38e-04 -2.50e-04 1.42e-04 7.83e-04 -1.41e-04 1.36e-04 -1.42e-05 1.29e-04 2.05e-04 -2.43e-04 -7.01e-07 -1.41e-04 3.97e-04 1.35e-04 3.19e-04 9.95e-05 -3.00e-04 -2.37e-04 -8.28e-05 1.36e-04 1.35e-04 1.04e-03
    1.44e-03 -9.11e-05 -2.15e-04 5.58e-05 -1.48e-05 2.91e-04 1.21e-05 2.47e-04 -9.11e-05 7.70e-04 -9.90e-05 4.50e-04 -1.00e-04 -2.83e-04 1.57e-04 2.00e-04 -2.15e-04 -9.90e-05 6.18e-04 7.42e-06 4.76e-05 8.90e-05 4.50e-05 -2.62e-04 5.58e-05 4.50e-04 7.42e-06 2.37e-03 -7.89e-05 -2.83e-04 -3.41e-05 -1.52e-05 -1.48e-05 -1.00e-04 4.76e-05 -7.89e-05 5.80e-04 -3.47e-05 -2.62e-05 -1.58e-04 2.91e-04 -2.83e-04 8.90e-05 -2.83e-04 -3.47e-05 8.91e-04 -5.20e-05 -1.63e-04 1.21e-05 1.57e-04 4.50e-05 -3.41e-05 -2.62e-05 -5.20e-05 3.01e-04 1.31e-04 2.47e-04 2.00e-04 -2.62e-04 -1.52e-05 -1.58e-04 -1.63e-04 1.31e-04 9.86e-04
    6.40e-04 7.23e-05 6.55e-04 1.62e-04 -2.68e-04 -1.41e-04 2.00e-05 -1.13e-04 7.23e-05 5.03e-04 1.88e-04 -1.34e-04 -5.86e-05 2.92e-05 -8.18e-06 1.78e-04 6.55e-04 1.88e-04 2.31e-03 1.03e-05 -3.14e-04 -8.00e-05 1.98e-04 -1.04e-04 1.62e-04 -1.34e-04 1.03e-05 7.74e-04 -9.93e-05 -1.45e-04 2.33e-05 -2.48e-04 -2.68e-04 -5.86e-05 -3.14e-04 -9.93e-05 3.42e-04 8.03e-05 1.16e-05 3.84e-05 -1.41e-04 2.92e-05 -8.00e-05 -1.45e-04 8.03e-05 4.21e-04 -3.11e-05 -1.34e-06 2.00e-05 -8.18e-06 1.98e-04 2.33e-05 1.16e-05 -3.11e-05 3.88e-04 6.16e-05 -1.13e-04 1.78e-04 -1.04e-04 -2.48e-04 3.84e-05 -1.34e-06 6.16e-05 8.60e-04
    1.76e-03 -8.99e-04 2.42e-04 -2.27e-04 -6.80e-05 5.96e-04 -4.19e-04 -1.40e-05 -8.99e-04 1.56e-03 -7.19e-05 1.51e-04 1.98e-04 3.39e-04 2.34e-04 -2.37e-04 2.42e-04 -7.19e-05 6.10e-04 1.95e-04 -1.67e-04 -4.33e-04 -1.70e-04 -8.85e-05 -2.27e-04 1.51e-04 1.95e-04 5.53e-04 -1.21e-04 -5.00e-04 2.44e-05 -4.33e-05 -6.80e-05 1.98e-04 -1.67e-04 -1.21e-04 5.22e-04 4.31e-04 1.78e-04 -5.68e-05 5.96e-04 3.39e-04 -4.33e-04 -5.00e-04 4.31e-04 1.89e-03 -1.04e-04 -2.25e-04 -4.19e-04 2.34e-04 -1.70e-04 2.44e-05 1.78e-04 -1.04e-04 5.09e-04 -5.71e-05 -1.40e-05 -2.37e-04 -8.85e-05 -4.33e-05 -5.68e-05 -2.25e-04 -5.71e-05 3.95e-04
    1.33e-03 -4.06e-04 -2.45e-05 3.16e-04 -4.64e-04 -1.52e-04 -1.44e-04 -2.28e-04 -4.06e-04 2.06e-03 3.14e-04 4.35e-04 -1.60e-04 -1.40e-04 4.52e-04 -2.34e-04 -2.45e-05 3.14e-04 5.53e-04 7.24e-05 -1.03e-04 1.16e-04 9.74e-05 -1.09e-04 3.16e-04 4.35e-04 7.24e-05 8.27e-04 -3.10e-04 -1.31e-04 1.87e-04 -3.24e-04 -4.64e-04 -1.60e-04 -1.03e-04 -3.10e-04 5.12e-04 1.29e-04 4.62e-05 2.65e-04 -1.52e-04 -1.40e-04 1.16e-04 -1.31e-04 1.29e-04 9.59e-04 1.98e-04 1.82e-04 -1.44e-04 4.52e-04 9.74e-05 1.87e-04 4.62e-05 1.98e-04 7.21e-04 -1.39e-05 -2.28e-04 -2.34e-04 -1.09e-04 -3.24e-04 2.65e-04 1.82e-04 -1.39e-05 6.63e-04
    6.05e-04 -3.60e-04 -5.60e-05 -4.94e-06 -6.33e-05 -2.39e-04 2.49e-04 -3.51e-05 -3.60e-04 1.77e-03 -1.46e-04 -7.75e-06 3.93e-04 9.99e-04 -2.05e-04 -1.12e-07 -5.60e-05 -1.46e-04 4.18e-04 1.55e-04 -2.00e-04 -1.67e-04 -1.13e-04 -6.49e-05 -4.94e-06 -7.75e-06 1.55e-04 7.80e-04 -2.79e-04 -2.30e-04 -1.89e-04 -1.81e-04 -6.33e-05 3.93e-04 -2.00e-04 -2.79e-04 5.39e-04 3.52e-04 1.79e-04 1.51e-04 -2.39e-04 9.99e-04 -1.67e-04 -2.30e-04 3.52e-04 1.49e-03 -2.99e-05 1.86e-05 2.49e-04 -2.05e-04 -1.13e-04 -1.89e-04 1.79e-04 -2.99e-05 9.00e-04 4.26e-05 -3.51e-05 -1.12e-07 -6.49e-05 -1.81e-04 1.51e-04 1.86e-05 4.26e-05 4.88e-04
    8.19e-04 -3.92e-04 -2.53e-04 -1.17e-04 2.08e-04 -1.30e-04 2.28e-05 9.03e-05 -3.92e-04 2.22e-03 1.10e-05 -1.81e-05 4.38e-04 4.64e-04 -1.85e-04 -3.80e-05 -2.53e-04 1.10e-05 4.83e-04 6.93e-05 -1.14e-04 -7.36e-05 1.83e-04 -1.34e-04 -1.17e-04 -1.81e-05 6.93e-05 6.40e-04 -1.67e-04 -2.79e-05 2.27e-04 -6.29e-05 2.08e-04 4.38e-04 -1.14e-04 -1.67e-04 6.24e-04 1.26e-04 -1.29e-04 1.05e-04 -1.30e-04 4.64e-04 -7.36e-05 -2.79e-05 1.26e-04 1.44e-03 -8.46e-04 2.56e-04 2.28e-05 -1.85e-04 1.83e-04 2.27e-04 -1.29e-04 -8.46e-04 1.79e-03 -2.48e-04 9.03e-05 -3.80e-05 -1.34e-04 -6.29e-05 1.05e-04 2.56e-04 -2.48e-04 5.26e-04
    1.74e-03 -1.20e-03 2.95e-04 -1.34e-04 -2.46e-04 1.55e-04 -3.78e-04 -9.17e-06 -1.20e-03 1.83e-03 -1.62e-04 2.33e-04 2.68e-04 1.57e-04 2.00e-04 -1.05e-04 2.95e-04 -1.62e-04 6.54e-04 1.87e-04 -2.21e-04 -5.29e-04 -1.45e-04 -8.85e-05 -1.34e-04 2.33e-04 1.87e-04 7.40e-04 -2.20e-04 -6.02e-04 -5.58e-05 -9.14e-05 -2.46e-04 2.68e-04 -2.21e-04 -2.20e-04 6.22e-04 5.18e-04 1.92e-04 6.56e-06 1.55e-04 1.57e-04 -5.29e-04 -6.02e-04 5.18e-04 1.66e-03 -1.26e-04 -7.76e-05 -3.78e-04 2.00e-04 -1.45e-04 -5.58e-05 1.92e-04 -1.26e-04 6.83e-04 3.71e-05 -9.17e-06 -1.05e-04 -8.85e-05 -9.14e-05 6.56e-06 -7.76e-05 3.71e-05 4.53e-04
    5.53e-04 -1.92e-04 4.14e-04 -3.52e-05 -1.84e-04 -2.44e-04 -5.23e-05 3.22e-04 -1.92e-04 1.87e-03 -2.52e-04 2.79e-04 -3.48e-05 4.04e-04 3.22e-05 -1.65e-04 4.14e-04 -2.52e-04 2.05e-03 -1.30e-04 -1.66e-04 -2.53e-04 4.17e-05 2.48e-04 -3.52e-05 2.79e-04 -1.30e-04 5.29e-04 -5.37e-05 2.54e-04 4.21e-05 -2.28e-04 -1.84e-04 -3.48e-05 -1.66e-04 -5.37e-05 1.81e-04 1.27e-04 1.32e-04 -4.95e-05 -2.44e-04 4.04e-04 -2.53e-04 2.54e-04 1.27e-04 7.71e-04 3.11e-04 -3.23e-04 -5.23e-05 3.22e-05 4.17e-05 4.21e-05 1.32e-04 3.11e-04 8.06e-04 1.49e-04 3.22e-04 -1.65e-04 2.48e-04 -2.28e-04 -4.95e-05 -3.23e-04 1.49e-04 8.81e-04
    5.75e-04 -2.49e-04 -1.78e-04 -2.24e-04 3.35e-05 4.47e-04 4.32e-05 -5.95e-05 -2.49e-04 1.19e-03 7.12e-05 6.05e-05 3.49e-04 -3.33e-04 -5.09e-04 4.48e-05 -1.78e-04 7.12e-05 6.04e-04 -6.49e-05 4.29e-05 -3.90e-05 -7.15e-05 -7.59e-05 -2.24e-04 6.05e-05 -6.49e-05 5.00e-04 9.16e-05 -2.50e-04 1.19e-04 1.77e-04 3.35e-05 3.49e-04 4.29e-05 9.16e-05 4.72e-04 -1.70e-04 -9.14e-05 -6.24e-05 4.47e-04 -3.33e-04 -3.90e-05 -2.50e-04 -1.70e-04 1.50e-03 -5.79e-04 -4.71e-05 4.32e-05 -5.09e-04 -7.15e-05 1.19e-04 -9.14e-05 -5.79e-04 2.54e-03 -1.93e-04 -5.95e-05 4.48e-05 -7.59e-05 1.77e-04 -6.24e-05 -4.71e-05 -1.93e-04 4.69e-04
    4.26e-04 -2.63e-04 -1.37e-04 2.69e-05 2.41e-06 3.36e-05 -3.21e-04 8.54e-05 -2.63e-04 2.17e-03 -4.01e-05 -1.94e-05 3.02e-04 5.08e-04 2.44e-04 9.49e-05 -1.37e-04 -4.01e-05 2.21e-04 -1.07e-05 -2.70e-05 -6.86e-05 1.75e-04 -1.05e-04 2.69e-05 -1.94e-05 -1.07e-05 2.82e-04 -1.81e-05 -7.29e-05 5.02e-05 -4.11e-05 2.41e-06 3.02e-04 -2.70e-05 -1.81e-05 2.96e-04 9.66e-05 1.11e-04 -7.55e-05 3.36e-05 5.08e-04 -6.86e-05 -7.29e-05 9.66e-05 1.60e-03 -1.20e-03 4.13e-05 -3.21e-04 2.44e-04 1.75e-04 5.02e-05 1.11e-04 -1.20e-03 1.81e-03 -9.28e-05 8.54e-05 9.49e-05 -1.05e-04 -4.11e-05 -7.55e-05 4.13e-05 -9.28e-05 3.47e-04
    8.71e-04 -8.77e-05 6.15e-04 1.12e-04 -3.68e-04 -2.95e-04 -2.14e-05 4.28e-05 -8.77e-05 8.95e-04 1.56e-05 -2.90e-04 5.32e-05 1.58e-04 -7.92e-05 1.89e-04 6.15e-04 1.56e-05 2.17e-03 -3.00e-05 -2.65e-04 -2.52e-04 1.60e-04 -1.49e-04 1.12e-04 -2.90e-04 -3.00e-05 7.81e-04 -7.81e-05 -9.82e-05 3.04e-05 -2.32e-04 -3.68e-04 5.32e-05 -2.65e-04 -7.81e-05 4.67e-04 2.40e-04 7.82e-05 -3.23e-06 -2.95e-04 1.58e-04 -2.52e-04 -9.82e-05 2.40e-04 8.51e-04 8.34e-05 -2.29e-05 -2.14e-05 -7.92e-05 1.60e-04 3.04e-05 7.82e-05 8.34e-05 6.56e-04 1.71e-04 4.28e-05 1.89e-04 -1.49e-04 -2.32e-04 -3.23e-06 -2.29e-05 1.71e-04 1.09e-03
    1.96e-03 -1.22e-04 4.07e-04 -5.13e-04 -2.08e-04 1.78e-04 1.51e-04 -1.22e-04 -1.22e-04 5.69e-04 -3.95e-05 1.32e-04 7.55e-05 -2.54e-04 1.33e-04 1.66e-04 4.07e-04 -3.95e-05 5.52e-04 -1.66e-04 4.44e-05 7.68e-08 1.11e-04 -2.98e-04 -5.13e-04 1.32e-04 -1.66e-04 1.78e-03 1.52e-04 -1.30e-05 -5.20e-05 1.74e-04 -2.08e-04 7.55e-05 4.44e-05 1.52e-04 6.40e-04 -2.06e-04 2.42e-05 -1.01e-05 1.78e-04 -2.54e-04 7.68e-08 -1.30e-05 -2.06e-04 7.64e-04 -3.79e-05 -1.68e-04 1.51e-04 1.33e-04 1.11e-04 -5.20e-05 2.42e-05 -3.79e-05 1.59e-04 1.53e-05 -1.22e-04 1.66e-04 -2.98e-04 1.74e-04 -1.01e-05 -1.68e-04 1.53e-05 6.02e-04
    1.05e-03 -1.59e-04 -3.37e-04 7.19e-05 4.81e-05 2.86e-04 7.45e-06 3.95e-04 -1.59e-04 6.49e-04 7.74e-05 6.17e-04 -8.20e-05 -3.28e-04 1.84e-04 -1.64e-04 -3.37e-04 7.74e-05 6.61e-04 -7.23e-06 8.81e-05 -7.37e-05 5.13e-05 -3.53e-04 7.19e-05 6.17e-04 -7.23e-06 1.93e-03 -2.44e-05 -7.44e-05 8.90e-05 -2.48e-04 4.81e-05 -8.20e-05 8.81e-05 -2.44e-05 6.18e-04 1.04e-05 -6.93e-05 -1.83e-04 2.86e-04 -3.28e-04 -7.37e-05 -7.44e-05 1.04e-05 6.45e-04 -1.33e-04 1.51e-04 7.45e-06 1.84e-04 5.13e-05 8.90e-05 -6.93e-05 -1.33e-04 2.11e-04 4.09e-05 3.95e-04 -1.64e-04 -3.53e-04 -2.48e-04 -1.83e-04 1.51e-04 4.09e-05 8.36e-04
    6.97e-04 -3.18e-04 -5.88e-05 -1.17e-04 8.35e-05 -2.97e-04 2.70e-04 -4.31e-05 -3.18e-04 2.11e-03 -2.02e-04 3.84e-05 4.10e-04 1.10e-03 -9.25e-05 -1.15e-04 -5.88e-05 -2.02e-04 3.99e-04 1.39e-04 -1.93e-04 -7.04e-05 -1.98e-05 -8.63e-05 -1.17e-04 3.84e-05 1.39e-04 5.35e-04 -1.51e-04 1.01e-04 9.36e-06 -6.00e-05 8.35e-05 4.10e-04 -1.93e-04 -1.51e-04 5.89e-04 1.59e-04 1.76e-04 1.03e-04 -2.97e-04 1.10e-03 -7.04e-05 1.01e-04 1.59e-04 1.30e-03 -9.94e-05 6.04e-05 2.70e-04 -9.25e-05 -1.98e-05 9.36e-06 1.76e-04 -9.94e-05 8.94e-04 -2.99e-05 -4.31e-05 -1.15e-04 -8.63e-05 -6.00e-05 1.03e-04 6.04e-05 -2.99e-05 3.61e-04
    4.25e-04 7.89e-05 6.54e-05 -9.23e-05 -3.63e-05 7.00e-06 -1.22e-04 -1.58e-04 7.89e-05 6.17e-04 -2.81e-04 2.54e-04 -9.64e-05 1.59e-05 -1.25e-04 -1.36e-04 6.54e-05 -2.81e-04 1.90e-03 -4.13e-04 3.32e-04 -2.02e-04 -2.93e-04 9.89e-04 -9.23e-05 2.54e-04 -4.13e-04 1.05e-03 2.15e-06 -1.26e-04 1.27e-04 5.70e-05 -3.63e-05 -9.64e-05 3.32e-04 2.15e-06 4.69e-04 7.80e-06 1.97e-05 2.51e-04 7.00e-06 1.59e-05 -2.02e-04 -1.26e-04 7.80e-06 3.61e-04 -2.97e-07 -2.78e-04 -1.22e-04 -1.25e-04 -2.93e-04 1.27e-04 1.97e-05 -2.97e-07 5.96e-04 -1.22e-04 -1.58e-04 -1.36e-04 9.89e-04 5.70e-05 2.51e-04 -2.78e-04 -1.22e-04 1.58e-03
    2.01e-03 -1.38e-04 2.41e-04 -3.49e-04 -1.66e-04 3.17e-04 1.31e-04 -2.66e-05 -1.38e-04 6.45e-04 -8.80e-05 3.19e-04 1.12e-05 -2.36e-04 1.45e-04 1.25e-04 2.41e-04 -8.80e-05 5.46e-04 -2.78e-05 3.21e-05 1.17e-04 6.42e-05 -2.71e-04 -3.49e-04 3.19e-04 -2.78e-05 2.23e-03 1.17e-06 -1.57e-04 1.30e-05 1.84e-04 -1.66e-04 1.12e-05 3.21e-05 1.17e-06 4.96e-04 -1.53e-04 -2.08e-05 -6.65e-05 3.17e-04 -2.36e-04 1.17e-04 -1.57e-04 -1.53e-04 7.20e-04 7.56e-06 -2.18e-04 1.31e-04 1.45e-04 6.42e-05 1.30e-05 -2.08e-05 7.56e-06 1.73e-04 3.35e-05 -2.66e-05 1.25e-04 -2.71e-04 1.84e-04 -6.65e-05 -2.18e-04 3.35e-05 7.34e-04
    3.45e-04 -2.19e-04 6.81e-07 -6.28e-06 -4.24e-06 -5.39e-05 -1.79e-04 1.68e-04 -2.19e-04 1.73e-03 -1.13e-04 4.97e-05 2.90e-04 5.37e-04 7.12e-04 1.73e-06 6.81e-07 -1.13e-04 4.92e-04 -1.56e-04 8.00e-06 2.54e-04 -2.04e-04 6.15e-05 -6.28e-06 4.97e-05 -1.56e-04 6.50e-04 -2.85e-04 -2.90e-04 -1.94e-04 -2.34e-04 -4.24e-06 2.90e-04 8.00e-06 -2.85e-04 3.01e-04 2.47e-04 3.16e-04 1.26e-04 -5.39e-05 5.37e-04 2.54e-04 -2.90e-04 2.47e-04 9.72e-04 -3.56e-05 8.09e-05 -1.79e-04 7.12e-04 -2.04e-04 -1.94e-04 3.16e-04 -3.56e-05 1.66e-03 -1.40e-04 1.68e-04 1.73e-06 6.15e-05 -2.34e-04 1.26e-04 8.09e-05 -1.40e-04 4.11e-04
    1.71e-03 -9.13e-04 1.66e-04 -1.88e-04 -1.55e-04 4.23e-04 -3.93e-04 -5.03e-05 -9.13e-04 1.74e-03 -2.00e-04 8.35e-05 3.45e-04 4.42e-04 1.79e-04 -1.78e-04 1.66e-04 -2.00e-04 6.99e-04 9.68e-05 -1.96e-04 -3.93e-04 -1.83e-04 -4.51e-05 -1.88e-04 8.35e-05 9.68e-05 7.05e-04 -1.45e-04 -4.68e-04 -3.41e-05 -4.69e-05 -1.55e-04 3.45e-04 -1.96e-04 -1.45e-04 5.10e-04 4.29e-04 2.47e-04 -2.48e-05 4.23e-04 4.42e-04 -3.93e-04 -4.68e-04 4.29e-04 1.84e-03 -8.75e-05 -2.74e-04 -3.93e-04 1.79e-04 -1.83e-04 -3.41e-05 2.47e-04 -8.75e-05 5.39e-04 -1.20e-05 -5.03e-05 -1.78e-04 -4.51e-05 -4.69e-05 -2.48e-05 -2.74e-04 -1.20e-05 4.31e-04
    1.11e-03 5.33e-05 -4.35e-04 3.27e-04 -1.79e-04 1.72e-04 -9.51e-05 2.40e-04 5.33e-05 6.91e-04 -8.98e-05 -1.73e-05 -1.60e-04 -2.42e-04 1.14e-04 3.38e-04 -4.35e-04 -8.98e-05 8.23e-04 -2.41e-04 6.37e-05 -2.68e-05 1.12e-04 -2.18e-04 3.27e-04 -1.73e-05 -2.41e-04 2.17e-03 -2.61e-04 -3.00e-04 -2.95e-04 -8.63e-05 -1.79e-04 -1.60e-04 6.37e-05 -2.61e-04 3.86e-04 6.06e-05 2.40e-05 -1.41e-04 1.72e-04 -2.42e-04 -2.68e-05 -3.00e-04 6.06e-05 9.72e-04 -4.09e-05 -3.32e-04 -9.51e-05 1.14e-04 1.12e-04 -2.95e-04 2.40e-05 -4.09e-05 2.73e-04 9.16e-05 2.40e-04 3.38e-04 -2.18e-04 -8.63e-05 -1.41e-04 -3.32e-04 9.16e-05 1.07e-03
    1.77e-03 -9.99e-04 3.85e-04 1.49e-04 -4.86e-04 -2.58e-04 -2.93e-04 -2.57e-05 -9.99e-04 1.96e-03 8.12e-05 3.86e-04 2.21e-05 -3.18e-04 3.24e-04 -5.67e-06 3.85e-04 8.12e-05 6.71e-04 1.65e-04 -2.65e-04 -4.32e-04 -2.59e-06 -1.01e-04 1.49e-04 3.86e-04 1.65e-04 9.07e-04 -3.01e-04 -4.86e-04 -4.90e-06 -1.39e-04 -4.86e-04 2.21e-05 -2.65e-04 -3.01e-04 6.53e-04 3.67e-04 8.04e-05 8.02e-05 -2.58e-04 -3.18e-04 -4.32e-04 -4.86e-04 3.67e-04 1.18e-03 -2.72e-05 1.92e-04 -2.93e-04 3.24e-04 -2.59e-06 -4.90e-06 8.04e-05 -2.72e-05 6.90e-04 3.55e-05 -2.57e-05 -5.67e-06 -1.01e-04 -1.39e-04 8.02e-05 1.92e-04 3.55e-05 4.94e-04
    7.21e-04 -2.44e-04 2.86e-05 -3.74e-04 -3.46e-05 1.63e-04 2.24e-04 -2.38e-04 -2.44e-04 9.22e-04 2.14e-04 -8.61e-05 1.83e-04 -2.99e-04 -4.26e-04 1.15e-04 2.86e-05 2.14e-04 1.21e-03 -3.11e-04 2.68e-04 2.76e-04 -5.06e-04 -2.04e-04 -3.74e-04 -8.61e-05 -3.11e-04 6.74e-04 -5.59e-05 -1.12e-04 4.73e-05 2.35e-04 -3.46e-05 1.83e-04 2.68e-04 -5.59e-05 4.20e-04 -1.90e-05 -3.26e-05 -1.43e-04 1.63e-04 -2.99e-04 2.76e-04 -1.12e-04 -1.90e-05 7.42e-04 -1.01e-04 -1.06e-04 2.24e-04 -4.26e-04 -5.06e-04 4.73e-05 -3.26e-05 -1.01e-04 1.82e-03 -3.80e-04 -2.38e-04 1.15e-04 -2.04e-04 2.35e-04 -1.43e-04 -1.06e-04 -3.80e-04 6.37e-04
    9.72e-04 -3.29e-04 -3.59e-04 1.61e-04 -2.58e-04 -6.96e-05 -1.29e-04 -7.58e-05 -3.29e-04 1.99e-03 1.86e-04 2.48e-04 -9.69e-05 2.34e-05 3.60e-04 -4.07e-04 -3.59e-04 1.86e-04 7.92e-04 -1.32e-04 9.73e-05 1.75e-04 8.00e-05 -7.36e-06 1.61e-04 2.48e-04 -1.32e-04 8.03e-04 -1.67e-04 1.44e-04 1.82e-04 -3.43e-04 -2.58e-04 -9.69e-05 9.73e-05 -1.67e-04 3.52e-04 2.26e-05 9.82e-05 2.18e-04 -6.96e-05 2.34e-05 1.75e-04 1.44e-04 2.26e-05 9.62e-04 3.22e-04 -2.40e-04 -1.29e-04 3.60e-04 8.00e-05 1.82e-04 9.82e-05 3.22e-04 7.85e-04 -5.54e-05 -7.58e-05 -4.07e-04 -7.36e-06 -3.43e-04 2.18e-04 -2.40e-04 -5.54e-05 8.85e-04
    2.35e-04 5.16e-05 -2.89e-05 -1.16e-04 -1.27e-04 -5.66e-05 9.48e-05 -1.20e-04 5.16e-05 2.50e-04 -9.62e-06 5.19e-05 -5.55e-05 -3.88e-06 -3.20e-05 -6.41e-05 -2.89e-05 -9.62e-06 5.07e-04 1.21e-04 7.99e-05 -2.32e-04 -3.47e-04 -4.38e-05 -1.16e-04 5.19e-05 1.21e-04 8.18e-04 4.46e-04 -2.09e-04 -1.41e-04 -1.88e-05 -1.27e-04 -5.55e-05 7.99e-05 4.46e-04 1.43e-03 -7.76e-05 -1.54e-04 -3.70e-04 -5.66e-05 -3.88e-06 -2.32e-04 -2.09e-04 -7.76e-05 5.14e-04 9.88e-05 -1.05e-04 9.48e-05 -3.20e-05 -3.47e-04 -1.41e-04 -1.54e-04 9.88e-05 5.82e-04 9.97e-05 -1.20e-04 -6.41e-05 -4.38e-05 -1.88e-05 -3.70e-04 -1.05e-04 9.97e-05 1.80e-03
    5.41e-04 -2.19e-04 -2.54e-04 -3.13e-06 -1.31e-04 -1.65e-04 -5.98e-05 1.74e-04 -2.19e-04 1.99e-03 5.02e-06 3.00e-04 -8.82e-05 2.06e-04 2.52e-04 -3.69e-04 -2.54e-04 5.02e-06 1.22e-03 -1.98e-04 1.04e-04 -4.06e-05 4.72e-05 1.67e-04 -3.13e-06 3.00e-04 -1.98e-04 6.36e-04 -8.84e-05 2.20e-04 1.66e-04 -2.92e-04 -1.31e-04 -8.82e-05 1.04e-04 -8.84e-05 2.03e-04 6.55e-05 1.05e-04 7.62e-05 -1.65e-04 2.06e-04 -4.06e-05 2.20e-04 6.55e-05 8.56e-04 3.95e-04 -4.02e-04 -5.98e-05 2.52e-04 4.72e-05 1.66e-04 1.05e-04 3.95e-04 7.67e-04 3.37e-05 1.74e-04 -3.69e-04 1.67e-04 -2.92e-04 7.62e-05 -4.02e-04 3.37e-05 8.98e-04
    4.57e-04 -6.51e-05 1.48e-04 2.55e-05 -1.35e-04 -1.59e-04 -1.93e-04 2.20e-04 -6.51e-05 1.57e-03 -4.00e-04 6.60e-04 -2.68e-04 -3.59e-05 -3.71e-05 -2.53e-04 1.48e-04 -4.00e-04 1.28e-03 -3.64e-04 9.30e-05 1.31e-04 -1.96e-04 2.56e-04 2.55e-05 6.60e-04 -3.64e-04 1.06e-03 -2.51e-04 -1.73e-04 -2.56e-04 -3.14e-04 -1.35e-04 -2.68e-04 9.30e-05 -2.51e-04 4.72e-04 1.40e-04 3.49e-04 9.35e-05 -1.59e-04 -3.59e-05 1.31e-04 -1.73e-04 1.40e-04 7.46e-04 3.45e-04 -9.38e-05 -1.93e-04 -3.71e-05 -1.96e-04 -2.56e-04 3.49e-04 3.45e-04 1.32e-03 9.23e-05 2.20e-04 -2.53e-04 2.56e-04 -3.14e-04 9.35e-05 -9.38e-05 9.23e-05 7.64e-04
    1.60e-03 -1.49e-04 -2.26e-05 -1.46e-04 -2.19e-05 3.13e-04 4.97e-05 1.19e-04 -1.49e-04 6.88e-04 -1.57e-04 3.98e-04 -4.58e-05 -3.00e-04 1.22e-04 2.03e-04 -2.26e-05 -1.57e-04 4.88e-04 3.81e-06 7.56e-05 1.72e-04 3.90e-05 -2.99e-04 -1.46e-04 3.98e-04 3.81e-06 2.09e-03 -1.95e-05 -2.87e-04 -4.28e-05 -1.18e-05 -2.19e-05 -4.58e-05 7.56e-05 -1.95e-05 5.47e-04 -1.62e-04 -2.05e-05 -1.87e-04 3.13e-04 -3.00e-04 1.72e-04 -2.87e-04 -1.62e-04 8.64e-04 -2.66e-05 -2.37e-04 4.97e-05 1.22e-04 3.90e-05 -4.28e-05 -2.05e-05 -2.66e-05 2.16e-04 7.50e-05 1.19e-04 2.03e-04 -2.99e-04 -1.18e-05 -1.87e-04 -2.37e-04 7.50e-05 9.19e-04
    1.77e-03 -4.90e-05 1.12e-04 -1.26e-04 2.26e-04 5.05e-04 2.72e-05 -8.85e-06 -4.90e-05 4.65e-04 -5.37e-05 -4.29e-05 -2.23e-04 -1.15e-04 -1.53e-04 -2.71e-05 1.12e-04 -5.37e-05 4.15e-04 1.27e-04 -2.56e-04 6.30e-05 5.93e-05 -8.85e-05 -1.26e-04 -4.29e-05 1.27e-04 2.62e-04 -1.63e-04 -1.70e-05 8.47e-05 -2.56e-05 2.26e-04 -2.23e-04 -2.56e-04 -1.63e-04 1.87e-03 8.20e-04 7.11e-05 -3.87e-06 5.05e-04 -1.15e-04 6.30e-05 -1.70e-05 8.20e-04 1.16e-03 -2.40e-05 1.59e-04 2.72e-05 -1.53e-04 5.93e-05 8.47e-05 7.11e-05 -2.40e-05 4.03e-04 -5.65e-05 -8.85e-06 -2.71e-05 -8.85e-05 -2.56e-05 -3.87e-06 1.59e-04 -5.65e-05 3.64e-04
    6.62e-04 1.61e-04 6.51e-05 3.21e-04 -3.00e-04 -2.32e-05 6.47e-05 1.02e-04 1.61e-04 7.16e-04 1.02e-04 4.69e-04 -2.72e-04 -1.28e-04 1.40e-04 -5.04e-05 6.51e-05 1.02e-04 1.49e-03 -5.06e-04 9.33e-06 -3.27e-04 2.46e-04 -3.73e-04 3.21e-04 4.69e-04 -5.06e-04 2.37e-03 -3.96e-04 -1.60e-04 -1.24e-04 -2.63e-04 -3.00e-04 -2.72e-04 9.33e-06 -3.96e-04 5.41e-04 8.93e-05 2.15e-05 2.62e-05 -2.32e-05 -1.28e-04 -3.27e-04 -1.60e-04 8.93e-05 6.32e-04 -8.18e-05 2.75e-04 6.47e-05 1.40e-04 2.46e-04 -1.24e-04 2.15e-05 -8.18e-05 5.65e-04 6.69e-05 1.02e-04 -5.04e-05 -3.73e-04 -2.63e-04 2.62e-05 2.75e-04 6.69e-05 1.08e-03
    8.29e-04 2.63e-05 5.49e-04 2.00e-04 -3.85e-04 -2.60e-04 3.76e-05 -2.21e-05 2.63e-05 6.88e-04 1.54e-04 -2.69e-04 -5.45e-06 4.46e-05 -8.45e-06 2.12e-04 5.49e-04 1.54e-04 2.09e-03 -1.35e-04 -1.91e-04 -2.22e-04 2.48e-04 -2.38e-04 2.00e-04 -2.69e-04 -1.35e-04 1.10e-03 -1.53e-04 -2.14e-04 9.54e-06 -2.26e-04 -3.85e-04 -5.45e-06 -1.91e-04 -1.53e-04 5.14e-04 2.28e-04 4.44e-05 2.39e-05 -2.60e-04 4.46e-05 -2.22e-04 -2.14e-04 2.28e-04 7.72e-04 9.73e-06 9.36e-05 3.76e-05 -8.45e-06 2.48e-04 9.54e-06 4.44e-05 9.73e-06 5.75e-04 1.39e-04 -2.21e-05 2.12e-04 -2.38e-04 -2.26e-04 2.39e-05 9.36e-05 1.39e-04 1.08e-03
    6.95e-04 -1.81e-04 3.91e-04 1.87e-05 -2.82e-04 -3.17e-04 -3.19e-05 2.18e-04 -1.81e-04 1.06e-03 -4.09e-05 -3.09e-04 9.30e-05 2.46e-04 3.76e-06 2.64e-05 3.91e-04 -4.09e-05 1.97e-03 -1.28e-04 -1.78e-04 -3.44e-04 1.02e-04 -1.55e-05 1.87e-05 -3.09e-04 -1.28e-04 5.97e-04 -4.34e-05 7.15e-05 4.61e-05 -1.48e-04 -2.82e-04 9.30e-05 -1.78e-04 -4.34e-05 3.07e-04 2.32e-04 1.10e-04 -7.09e-05 -3.17e-04 2.46e-04 -3.44e-04 7.15e-05 2.32e-04 7.93e-04 2.09e-04 -2.54e-04 -3.19e-05 3.76e-06 1.02e-04 4.61e-05 1.10e-04 2.09e-04 6.36e-04 1.60e-04 2.18e-04 2.64e-05 -1.55e-05 -1.48e-04 -7.09e-05 -2.54e-04 1.60e-04 9.70e-04
    2.13e-03 -6.17e-05 4.18e-04 -2.31e-04 -2.84e-04 2.58e-04 1.67e-04 -1.32e-04 -6.17e-05 5.75e-04 -3.67e-05 7.21e-05 6.85e-06 -2.40e-04 1.29e-04 1.56e-04 4.18e-04 -3.67e-05 4.93e-04 -8.70e-05 -4.79e-05 7.40e-05 1.01e-04 -2.95e-04 -2.31e-04 7.21e-05 -8.70e-05 1.81e-03 7.01e-05 -1.39e-05 -5.70e-05 2.04e-04 -2.84e-04 6.85e-06 -4.79e-05 7.01e-05 5.00e-04 -1.68e-04 -3.95e-05 -2.04e-05 2.58e-04 -2.40e-04 7.40e-05 -1.39e-05 -1.68e-04 7.21e-04 6.24e-06 -1.94e-04 1.67e-04 1.29e-04 1.01e-04 -5.70e-05 -3.95e-05 6.24e-06 1.64e-04 5.26e-06 -1.32e-04 1.56e-04 -2.95e-04 2.04e-04 -2.04e-05 -1.94e-04 5.26e-06 6.21e-04
    1.56e-03 -6.70e-04 3.07e-04 2.98e-04 -5.12e-04 -2.47e-04 -1.96e-04 -1.53e-04 -6.70e-04 2.07e-03 2.36e-04 5.64e-04 -1.04e-04 -2.89e-04 4.49e-04 -6.05e-05 3.07e-04 2.36e-04 5.15e-04 1.61e-04 -2.54e-04 -1.53e-04 4.80e-05 -6.35e-05 2.98e-04 5.64e-04 1.61e-04 8.23e-04 -3.89e-04 -4.03e-04 1.10e-04 -1.79e-04 -5.12e-04 -1.04e-04 -2.54e-04 -3.89e-04 5.08e-04 2.76e-04 1.88e-05 1.73e-04 -2.47e-04 -2.89e-04 -1.53e-04 -4.03e-04 2.76e-04 1.00e-03 9.22e-05 2.34e-04 -1.96e-04 4.49e-04 4.80e-05 1.10e-04 1.88e-05 9.22e-05 5.51e-04 6.61e-05 -1.53e-04 -6.05e-05 -6.35e-05 -1.79e-04 1.73e-04 2.34e-04 6.61e-05 4.40e-04
    4.67e-04 2.47e-05 1.72e-04 -3.24e-04 -1.24e-06 -4.53e-05 4.42e-05 -2.72e-04 2.47e-05 6.66e-04 -1.28e-04 3.71e-04 5.78e-05 5.98e-06 -1.26e-04 -3.35e-05 1.72e-04 -1.28e-04 1.47e-03 -6.28e-04 2.03e-04 -2.63e-04 -2.64e-05 4.26e-04 -3.24e-04 3.71e-04 -6.28e-04 1.93e-03 1.55e-04 -4.92e-05 -2.63e-04 7.31e-04 -1.24e-06 5.78e-05 2.03e-04 1.55e-04 5.09e-04 -5.64e-05 2.98e-05 2.43e-04 -4.53e-05 5.98e-06 -2.63e-04 -4.92e-05 -5.64e-05 4.31e-04 -2.14e-05 -2.25e-04 4.42e-05 -1.26e-04 -2.64e-05 -2.63e-04 2.98e-05 -2.14e-05 5.24e-04 -2.86e-04 -2.72e-04 -3.35e-05 4.26e-04 7.31e-04 2.43e-04 -2.25e-04 -2.86e-04 1.76e-03
    1.05e-03 -6.63e-05 6.82e-05 5.23e-04 3.65e-04 6.13e-05 9.60e-06 1.46e-04 -6.63e-05 2.60e-04 -1.17e-04 1.10e-04 4.93e-05 7.60e-05 9.71e-05 -9.55e-06 6.82e-05 -1.17e-04 2.38e-04 -8.07e-05 -1.96e-05 -4.02e-06 4.56e-05 5.65e-05 5.23e-04 1.10e-04 -8.07e-05 1.18e-03 6.30e-04 2.70e-04 -1.19e-04 -2.53e-04 3.65e-04 4.93e-05 -1.96e-05 6.30e-04 1.14e-03 2.55e-04 1.00e-04 -2.89e-05 6.13e-05 7.60e-05 -4.02e-06 2.70e-04 2.55e-04 4.78e-04 1.36e-04 7.71e-05 9.60e-06 9.71e-05 4.56e-05 -1.19e-04 1.00e-04 1.36e-04 3.44e-04 2.83e-04 1.46e-04 -9.55e-06 5.65e-05 -2.53e-04 -2.89e-05 7.71e-05 2.83e-04 6.84e-04
    1.77e-03 -8.41e-04 4.38e-04 2.57e-04 -5.32e-04 -2.65e-04 -2.31e-04 -7.39e-05 -8.41e-04 1.86e-03 1.92e-04 4.00e-04 -4.71e-05 -3.40e-04 4.53e-04 -1.03e-04 4.38e-04 1.92e-04 5.52e-04 1.67e-04 -2.85e-04 -3.11e-04 5.36e-05 -9.36e-05 2.57e-04 4.00e-04 1.67e-04 7.32e-04 -2.99e-04 -3.26e-04 1.13e-04 -1.21e-04 -5.32e-04 -4.71e-05 -2.85e-04 -2.99e-04 5.14e-04 2.45e-04 -1.71e-05 1.31e-04 -2.65e-04 -3.40e-04 -3.11e-04 -3.26e-04 2.45e-04 9.01e-04 3.05e-05 2.80e-04 -2.31e-04 4.53e-04 5.36e-05 1.13e-04 -1.71e-05 3.05e-05 4.44e-04 2.47e-05 -7.39e-05 -1.03e-04 -9.36e-05 -1.21e-04 1.31e-04 2.80e-04 2.47e-05 3.57e-04
    1.29e-03 -9.34e-04 2.06e-04 -4.65e-05 -2.96e-04 7.47e-05 -2.22e-04 -6.28e-05 -9.34e-04 1.68e-03 -2.05e-04 1.89e-04 3.83e-04 4.45e-04 4.29e-05 -7.54e-05 2.06e-04 -2.05e-04 5.18e-04 2.09e-04 -2.66e-04 -4.13e-04 -2.24e-04 -7.93e-05 -4.65e-05 1.89e-04 2.09e-04 7.13e-04 -2.16e-04 -5.57e-04 -1.17e-04 -1.34e-04 -2.96e-04 3.83e-04 -2.66e-04 -2.16e-04 5.29e-04 4.48e-04 2.83e-04 3.43e-05 7.47e-05 4.45e-04 -4.13e-04 -5.57e-04 4.48e-04 1.56e-03 -4.25e-05 -8.32e-05 -2.22e-04 4.29e-05 -2.24e-04 -1.17e-04 2.83e-04 -4.25e-05 6.62e-04 1.08e-04 -6.28e-05 -7.54e-05 -7.93e-05 -1.34e-04 3.43e-05 -8.32e-05 1.08e-04 3.88e-04
    4.69e-04 7.35e-05 1.52e-04 -1.84e-04 -2.63e-05 -1.74e-06 -8.38e-05 -3.08e-04 7.35e-05 6.59e-04 -2.99e-04 3.82e-04 -5.57e-05 3.11e-05 -1.69e-04 -1.51e-04 1.52e-04 -2.99e-04 1.86e-03 -6.35e-04 3.50e-04 -2.06e-04 -1.99e-04 7.54e-04 -1.84e-04 3.82e-04 -6.35e-04 1.49e-03 -1.40e-05 -1.06e-04 -6.63e-05 2.50e-04 -2.63e-05 -5.57e-05 3.50e-04 -1.40e-05 4.72e-04 -3.03e-05 2.29e-05 2.68e-04 -1.74e-06 3.11e-05 -2.06e-04 -1.06e-04 -3.03e-05 3.82e-04 2.38e-05 -2.99e-04 -8.38e-05 -1.69e-04 -1.99e-04 -6.63e-05 2.29e-05 2.38e-05 5.78e-04 -1.90e-04 -3.08e-04 -1.51e-04 7.54e-04 2.50e-04 2.68e-04 -2.99e-04 -1.90e-04 1.70e-03
    4.46e-04 2.30e-05 1.22e-04 -2.50e-04 -4.16e-05 2.08e-05 6.93e-05 -3.47e-04 2.30e-05 6.07e-04 -8.36e-05 4.15e-04 9.91e-05 -2.00e-05 -3.65e-05 -2.59e-04 1.22e-04 -8.36e-05 1.35e-03 -7.54e-04 2.04e-04 -3.20e-04 5.02e-05 7.65e-05 -2.50e-04 4.15e-04 -7.54e-04 1.96e-03 2.23e-04 7.31e-05 -3.49e-04 5.27e-04 -4.16e-05 9.91e-05 2.04e-04 2.23e-04 5.46e-04 -1.30e-04 -1.70e-05 2.35e-04 2.08e-05 -2.00e-05 -3.20e-04 7.31e-05 -1.30e-04 3.88e-04 -2.48e-05 -9.67e-05 6.93e-05 -3.65e-05 5.02e-05 -3.49e-04 -1.70e-05 -2.48e-05 4.80e-04 -2.93e-04 -3.47e-04 -2.59e-04 7.65e-05 5.27e-04 2.35e-04 -9.67e-05 -2.93e-04 1.60e-03
    4.48e-04 -3.27e-05 -6.70e-05 -1.45e-04 1.44e-05 -4.04e-06 7.43e-06 -1.64e-04 -3.27e-05 6.44e-04 1.42e-04 1.09e-04 2.12e-04 -1.53e-04 8.55e-05 -1.56e-04 -6.70e-05 1.42e-04 8.75e-04 -2.97e-04 8.48e-05 -2.65e-04 2.72e-05 -8.43e-05 -1.45e-04 1.09e-04 -2.97e-04 1.57e-03 4.71e-04 1.22e-04 -3.12e-04 6.84e-04 1.44e-05 2.12e-04 8.48e-05 4.71e-04 8.33e-04 -1.10e-05 1.42e-05 1.77e-04 -4.04e-06 -1.53e-04 -2.65e-04 1.22e-04 -1.10e-05 5.19e-04 -1.16e-04 5.12e-05 7.43e-06 8.55e-05 2.72e-05 -3.12e-04 1.42e-05 -1.16e-04 4.99e-04 -2.47e-04 -1.64e-04 -1.56e-04 -8.43e-05 6.84e-04 1.77e-04 5.12e-05 -2.47e-04 1.51e-03
    7.79e-04 1.03e-04 -3.49e-04 3.85e-04 -2.72e-04 -4.51e-05 -7.17e-05 3.21e-04 1.03e-04 6.63e-04 5.58e-05 -1.18e-04 -1.65e-04 -1.99e-04 1.04e-04 2.42e-04 -3.49e-04 5.58e-05 1.29e-03 -4.50e-04 7.59e-05 -2.98e-04 2.27e-04 -2.38e-04 3.85e-04 -1.18e-04 -4.50e-04 2.05e-03 -3.26e-04 -2.82e-04 -3.09e-04 -1.91e-04 -2.72e-04 -1.65e-04 7.59e-05 -3.26e-04 5.06e-04 1.83e-04 6.40e-05 -1.17e-04 -4.51e-05 -1.99e-04 -2.98e-04 -2.82e-04 1.83e-04 9.11e-04 -1.03e-04 -1.12e-04 -7.17e-05 1.04e-04 2.27e-04 -3.09e-04 6.40e-05 -1.03e-04 3.95e-04 1.17e-04 3.21e-04 2.42e-04 -2.38e-04 -1.91e-04 -1.17e-04 -1.12e-04 1.17e-04 1.08e-03
    2.24e-03 -1.46e-04 4.48e-04 5.03e-04 -5.56e-04 -6.89e-06 -4.11e-05 -1.34e-04 -1.46e-04 1.02e-03 1.73e-04 -3.47e-04 -7.22e-05 -1.58e-04 2.22e-04 -1.88e-04 4.48e-04 1.73e-04 8.10e-04 3.62e-05 -2.14e-04 4.71e-05 1.53e-04 -2.80e-04 5.03e-04 -3.47e-04 3.62e-05 1.02e-03 -8.08e-05 2.04e-04 -1.08e-04 -3.11e-06 -5.56e-04 -7.22e-05 -2.14e-04 -8.08e-05 6.32e-04 -1.11e-04 -3.81e-05 1.21e-04 -6.89e-06 -1.58e-04 4.71e-05 2.04e-04 -1.11e-04 1.02e-03 7.04e-05 -1.79e-05 -4.11e-05 2.22e-04 1.53e-04 -1.08e-04 -3.81e-05 7.04e-05 5.29e-04 -1.65e-04 -1.34e-04 -1.88e-04 -2.80e-04 -3.11e-06 1.21e-04 -1.79e-05 -1.65e-04 8.55e-04
    2.74e-04 -5.22e-05 9.49e-05 -1.56e-05 -1.00e-04 -9.89e-05 1.26e-04 -8.49e-05 -5.22e-05 3.50e-04 -1.14e-04 -1.21e-04 -1.15e-04 -4.99e-05 1.59e-04 -1.38e-04 9.49e-05 -1.14e-04 2.41e-03 -3.80e-04 3.30e-04 -1.41e-04 3.32e-04 6.82e-04 -1.56e-05 -1.21e-04 -3.80e-04 4.28e-04 3.99e-05 4.44e-05 -2.05e-04 1.95e-05 -1.00e-04 -1.15e-04 3.30e-04 3.99e-05 3.05e-04 8.72e-05 1.76e-04 7.41e-05 -9.89e-05 -4.99e-05 -1.41e-04 4.44e-05 8.72e-05 3.31e-04 -2.23e-04 8.88e-05 1.26e-04 1.59e-04 3.32e-04 -2.05e-04 1.76e-04 -2.23e-04 2.01e-03 -8.66e-04 -8.49e-05 -1.38e-04 6.82e-04 1.95e-05 7.41e-05 8.88e-05 -8.66e-04 1.05e-03
    6.31e-04 9.16e-05 4.95e-05 1.89e-04 -2.37e-04 9.20e-05 8.37e-05 -2.75e-05 9.16e-05 7.08e-04 6.60e-05 6.07e-04 -2.01e-04 -1.15e-04 1.40e-04 -3.22e-04 4.95e-05 6.60e-05 1.41e-03 -6.07e-04 8.17e-05 -2.94e-04 1.80e-04 -4.74e-04 1.89e-04 6.07e-04 -6.07e-04 2.21e-03 -2.21e-04 -1.85e-06 -1.15e-04 -2.99e-04 -2.37e-04 -2.01e-04 8.17e-05 -2.21e-04 5.12e-04 -3.35e-05 -6.40e-05 4.59e-05 9.20e-05 -1.15e-04 -2.94e-04 -1.85e-06 -3.35e-05 5.10e-04 -7.47e-05 3.00e-04 8.37e-05 1.40e-04 1.80e-04 -1.15e-04 -6.40e-05 -7.47e-05 5.55e-04 -5.51e-05 -2.75e-05 -3.22e-04 -4.74e-04 -2.99e-04 4.59e-05 3.00e-04 -5.51e-05 1.05e-03
    7.07e-04 8.74e-05 1.58e-04 -6.34e-05 -1.83e-04 2.79e-05 9.47e-05 -3.75e-04 8.74e-05 7.95e-04 -1.27e-04 5.08e-04 -1.10e-04 -7.75e-05 3.03e-05 -4.03e-04 1.58e-04 -1.27e-04 1.65e-03 -9.38e-04 2.18e-04 -1.96e-04 8.06e-05 -2.79e-04 -6.34e-05 5.08e-04 -9.38e-04 2.12e-03 -6.31e-05 1.93e-05 -2.25e-04 -9.71e-05 -1.83e-04 -1.10e-04 2.18e-04 -6.31e-05 6.18e-04 -8.30e-05 -3.70e-05 2.02e-04 2.79e-05 -7.75e-05 -1.96e-04 1.93e-05 -8.30e-05 5.32e-04 -1.26e-05 1.02e-04 9.47e-05 3.03e-05 8.06e-05 -2.25e-04 -3.70e-05 -1.26e-05 6.79e-04 -1.67e-04 -3.75e-04 -4.03e-04 -2.79e-04 -9.71e-05 2.02e-04 1.02e-04 -1.67e-04 1.39e-03
    1.45e-03 3.58e-06 -2.04e-04 2.67e-04 -1.36e-04 2.53e-04 -2.71e-06 1.25e-04 3.58e-06 6.56e-04 -1.51e-04 -9.87e-06 -1.04e-04 -2.29e-04 1.10e-04 3.31e-04 -2.04e-04 -1.51e-04 5.64e-04 -1.04e-04 7.59e-06 1.83e-04 8.77e-05 -2.20e-04 2.67e-04 -9.87e-06 -1.04e-04 2.03e-03 -1.60e-04 -2.70e-04 -2.42e-04 1.94e-05 -1.36e-04 -1.04e-04 7.59e-06 -1.60e-04 5.47e-04 -8.70e-05 4.79e-05 -1.26e-04 2.53e-04 -2.29e-04 1.83e-04 -2.70e-04 -8.70e-05 9.56e-04 1.37e-05 -3.63e-04 -2.71e-06 1.10e-04 8.77e-05 -2.42e-04 4.79e-05 1.37e-05 2.99e-04 5.54e-05 1.25e-04 3.31e-04 -2.20e-04 1.94e-05 -1.26e-04 -3.63e-04 5.54e-05 1.01e-03
    1.45e-03 -2.73e-04 8.85e-05 -1.70e-04 1.90e-04 6.52e-04 -3.33e-04 1.66e-04 -2.73e-04 7.94e-04 1.35e-04 -1.04e-04 -1.86e-04 -3.20e-05 1.28e-04 -1.23e-04 8.85e-05 1.35e-04 6.97e-04 6.39e-05 -3.04e-05 -3.78e-04 3.69e-05 -1.89e-04 -1.70e-04 -1.04e-04 6.39e-05 5.11e-04 5.03e-05 -3.74e-04 -1.84e-05 -8.38e-06 1.90e-04 -1.86e-04 -3.04e-05 5.03e-05 7.46e-04 2.73e-04 -1.15e-05 5.20e-05 6.52e-04 -3.20e-05 -3.78e-04 -3.74e-04 2.73e-04 1.55e-03 -1.03e-04 -1.01e-05 -3.33e-04 1.28e-04 3.69e-05 -1.84e-05 -1.15e-05 -1.03e-04 4.42e-04 -1.82e-04 1.66e-04 -1.23e-04 -1.89e-04 -8.38e-06 5.20e-05 -1.01e-05 -1.82e-04 4.77e-04
    8.33e-04 -1.61e-05 -3.11e-05 1.51e-04 7.85e-04 -1.80e-05 -1.28e-05 -9.52e-05 -1.61e-05 4.08e-04 -2.54e-05 -1.52e-05 -2.26e-05 -1.60e-05 -8.12e-05 -1.11e-04 -3.11e-05 -2.54e-05 1.03e-04 -2.53e-05 -1.94e-04 -6.42e-05 2.80e-05 1.09e-04 1.51e-04 -1.52e-05 -2.53e-05 3.26e-04 4.02e-04 1.35e-04 -4.60e-05 -1.32e-04 7.85e-04 -2.26e-05 -1.94e-04 4.02e-04 1.56e-03 1.92e-04 1.46e-04 -2.90e-04 -1.80e-05 -1.60e-05 -6.42e-05 1.35e-04 1.92e-04 5.08e-04 1.34e-05 1.47e-04 -1.28e-05 -8.12e-05 2.80e-05 -4.60e-05 1.46e-04 1.34e-05 3.40e-04 2.33e-04 -9.52e-05 -1.11e-04 1.09e-04 -1.32e-04 -2.90e-04 1.47e-04 2.33e-04 6.06e-04
    2.09e-03 -3.68e-04 5.14e-04 4.77e-04 -5.97e-04 -1.68e-04 -1.27e-04 -3.86e-05 -3.68e-04 1.33e-03 2.61e-04 -5.85e-05 -1.25e-04 -2.50e-04 3.22e-04 -2.45e-04 5.14e-04 2.61e-04 7.29e-04 1.01e-04 -2.82e-04 -2.21e-04 1.75e-04 -2.29e-04 4.77e-04 -5.85e-05 1.01e-04 7.74e-04 -1.72e-04 5.14e-05 5.25e-05 -3.51e-05 -5.97e-04 -1.25e-04 -2.82e-04 -1.72e-04 6.67e-04 4.44e-05 -5.17e-05 1.16e-04 -1.68e-04 -2.50e-04 -2.21e-04 5.14e-05 4.44e-05 1.04e-03 5.41e-05 1.84e-04 -1.27e-04 3.22e-04 1.75e-04 5.25e-05 -5.17e-05 5.41e-05 4.91e-04 -1.98e-04 -3.86e-05 -2.45e-04 -2.29e-04 -3.51e-05 1.16e-04 1.84e-04 -1.98e-04 5.82e-04
    5.89e-04 1.73e-04 1.78e-06 3.84e-04 -2.88e-04 -1.37e-04 6.14e-05 2.33e-04 1.73e-04 5.26e-04 1.09e-04 4.96e-05 -1.97e-04 -1.09e-04 9.23e-05 1.33e-04 1.78e-06 1.09e-04 1.55e-03 -4.82e-04 2.27e-05 -3.66e-04 2.94e-04 -2.58e-04 3.84e-04 4.96e-05 -4.82e-04 1.99e-03 -4.37e-04 -2.80e-04 -2.75e-04 -2.05e-04 -2.88e-04 -1.97e-04 2.27e-05 -4.37e-04 4.36e-04 1.95e-04 -7.27e-06 -6.26e-05 -1.37e-04 -1.09e-04 -3.66e-04 -2.80e-04 1.95e-04 6.57e-04 -9.31e-05 1.76e-04 6.14e-05 9.23e-05 2.94e-04 -2.75e-04 -7.27e-06 -9.31e-05 3.48e-04 1.00e-04 2.33e-04 1.33e-04 -2.58e-04 -2.05e-04 -6.26e-05 1.76e-04 1.00e-04 1.00e-03
    6.46e-04 -2.75e-04 -1.65e-04 -3.19e-04 -3.31e-05 2.60e-04 6.74e-05 -7.21e-05 -2.75e-04 1.03e-03 1.69e-04 -1.79e-05 2.74e-04 -2.80e-04 -5.35e-04 8.39e-05 -1.65e-04 1.69e-04 7.27e-04 -1.29e-04 1.37e-04 1.80e-04 -3.48e-04 -1.27e-04 -3.19e-04 -1.79e-05 -1.29e-04 6.12e-04 -1.44e-05 -2.22e-04 1.33e-04 1.60e-04 -3.31e-05 2.74e-04 1.37e-04 -1.44e-05 4.17e-04 -1.20e-04 -2.14e-04 -9.58e-05 2.60e-04 -2.80e-04 1.80e-04 -2.22e-04 -1.20e-04 9.32e-04 -2.08e-04 -1.52e-04 6.74e-05 -5.35e-04 -3.48e-04 1.33e-04 -2.14e-04 -2.08e-04 1.97e-03 -2.15e-04 -7.21e-05 8.39e-05 -1.27e-04 1.60e-04 -9.58e-05 -1.52e-04 -2.15e-04 4.97e-04
    5.75e-04 1.63e-04 2.89e-05 -1.15e-04 2.78e-04 4.68e-04 7.91e-05 -2.12e-04 1.63e-04 5.11e-04 -3.00e-05 2.71e-05 1.96e-04 1.90e-04 -1.05e-04 -2.06e-04 2.89e-05 -3.00e-05 7.53e-05 2.14e-05 -3.88e-05 -3.25e-05 2.69e-05 2.56e-05 -1.15e-04 2.71e-05 2.14e-05 7.98e-05 -7.78e-05 -2.01e-06 -5.77e-06 3.30e-05 2.78e-04 1.96e-04 -3.88e-05 -7.78e-05 4.07e-04 5.89e-05 5.70e-05 -1.14e-04 4.68e-04 1.90e-04 -3.25e-05 -2.01e-06 5.89e-05 2.18e-03 3.06e-05 -3.03e-04 7.91e-05 -1.05e-04 2.69e-05 -5.77e-06 5.70e-05 3.06e-05 2.37e-04 3.89e-05 -2.12e-04 -2.06e-04 2.56e-05 3.30e-05 -1.14e-04 -3.03e-04 3.89e-05 1.60e-04
    6.66e-04 -4.65e-04 -2.74e-04 -1.96e-04 9.87e-05 -4.71e-06 -2.97e-04 1.70e-04 -4.65e-04 1.77e-03 1.94e-04 6.93e-06 3.43e-04 2.38e-04 6.64e-05 1.40e-04 -2.74e-04 1.94e-04 4.52e-04 2.07e-07 -3.08e-05 1.06e-04 4.54e-05 -1.49e-04 -1.96e-04 6.93e-06 2.07e-07 6.54e-04 -2.06e-04 -2.05e-04 2.66e-04 -1.27e-04 9.87e-05 3.43e-04 -3.08e-05 -2.06e-04 3.56e-04 6.48e-05 -9.54e-05 1.45e-04 -4.71e-06 2.38e-04 1.06e-04 -2.05e-04 6.48e-05 1.17e-03 -7.60e-04 1.37e-04 -2.97e-04 6.64e-05 4.54e-05 2.66e-04 -9.54e-05 -7.60e-04 1.62e-03 -2.82e-04 1.70e-04 1.40e-04 -1.49e-04 -1.27e-04 1.45e-04 1.37e-04 -2.82e-04 4.82e-04
    5.68e-04 -4.43e-04 5.65e-05 -2.53e-04 -7.01e-05 1.14e-05 -1.04e-04 -9.61e-05 -4.43e-04 1.59e-03 3.96e-04 6.57e-06 2.64e-04 -7.63e-05 5.14e-04 4.01e-04 5.65e-05 3.96e-04 1.76e-03 -4.53e-04 2.47e-04 2.33e-04 6.90e-04 8.93e-05 -2.53e-04 6.57e-06 -4.53e-04 5.74e-04 -1.31e-04 -3.12e-05 -1.83e-04 1.64e-05 -7.01e-05 2.64e-04 2.47e-04 -1.31e-04 2.96e-04 4.77e-06 2.48e-04 2.20e-05 1.14e-05 -7.63e-05 2.33e-04 -3.12e-05 4.77e-06 5.75e-04 -1.31e-06 -1.63e-04 -1.04e-04 5.14e-04 6.90e-04 -1.83e-04 2.48e-04 -1.31e-06 1.16e-03 -3.69e-05 -9.61e-05 4.01e-04 8.93e-05 1.64e-05 2.20e-05 -1.63e-04 -3.69e-05 5.26e-04
    1.67e-03 -4.86e-04 2.29e-04 4.00e-04 -5.62e-04 -1.63e-04 -1.87e-04 -1.64e-04 -4.86e-04 2.00e-03 2.97e-04 3.89e-04 -1.01e-04 -1.49e-04 4.77e-04 -2.27e-04 2.29e-04 2.97e-04 5.28e-04 1.47e-04 -2.63e-04 -4.27e-06 7.85e-05 -1.25e-04 4.00e-04 3.89e-04 1.47e-04 7.09e-04 -3.22e-04 -1.52e-04 1.74e-04 -2.38e-04 -5.62e-04 -1.01e-04 -2.63e-04 -3.22e-04 4.91e-04 1.37e-04 1.08e-05 2.37e-04 -1.63e-04 -1.49e-04 -4.27e-06 -1.52e-04 1.37e-04 9.72e-04 1.42e-04 2.62e-04 -1.87e-04 4.77e-04 7.85e-05 1.74e-04 1.08e-05 1.42e-04 6.20e-04 -2.53e-05 -1.64e-04 -2.27e-04 -1.25e-04 -2.38e-04 2.37e-04 2.62e-04 -2.53e-05 5.43e-04
    7.26e-04 -7.45e-05 2.45e-04 1.33e-04 -3.27e-04 -2.94e-04 -7.07e-06 2.02e-04 -7.45e-05 7.80e-04 7.13e-05 -4.39e-04 6.35e-05 8.12e-05 -9.28e-07 1.89e-04 2.45e-04 7.13e-05 1.86e-03 -2.47e-04 -8.13e-05 -3.82e-04 1.96e-04 -1.44e-04 1.33e-04 -4.39e-04 -2.47e-04 1.01e-03 -1.41e-04 -1.37e-04 -4.84e-05 -1.44e-04 -3.27e-04 6.35e-05 -8.13e-05 -1.41e-04 4.24e-04 2.62e-04 1.15e-04 -7.23e-05 -2.94e-04 8.12e-05 -3.82e-04 -1.37e-04 2.62e-04 8.81e-04 8.10e-05 -1.42e-04 -7.07e-06 -9.28e-07 1.96e-04 -4.84e-05 1.15e-04 8.10e-05 5.23e-04 1.33e-04 2.02e-04 1.89e-04 -1.44e-04 -1.44e-04 -7.23e-05 -1.42e-04 1.33e-04 1.04e-03
    1.75e-03 -3.00e-05 1.62e-04 -2.29e-06 -1.40e-04 2.71e-04 9.98e-05 -2.68e-05 -3.00e-05 6.21e-04 -1.78e-04 9.70e-05 -4.59e-05 -2.64e-04 9.61e-05 3.23e-04 1.62e-04 -1.78e-04 4.96e-04 -7.85e-05 2.76e-06 2.27e-04 9.33e-05 -3.21e-04 -2.29e-06 9.70e-05 -7.85e-05 1.90e-03 1.21e-06 -1.85e-04 -1.48e-04 1.02e-04 -1.40e-04 -4.59e-05 2.76e-06 1.21e-06 4.84e-04 -2.14e-04 1.11e-05 -9.87e-05 2.71e-04 -2.64e-04 2.27e-04 -1.85e-04 -2.14e-04 9.13e-04 4.64e-05 -2.88e-04 9.98e-05 9.61e-05 9.33e-05 -1.48e-04 1.11e-05 4.64e-05 2.42e-04 3.32e-05 -2.68e-05 3.23e-04 -3.21e-04 1.02e-04 -9.87e-05 -2.88e-04 3.32e-05 9.05e-04
    7.92e-04 1.23e-04 5.30e-04 2.66e-04 -3.57e-04 -2.05e-04 7.51e-05 -1.61e-04 1.23e-04 5.98e-04 1.32e-04 -5.25e-05 -9.71e-05 2.56e-05 4.13e-05 1.24e-04 5.30e-04 1.32e-04 2.01e-03 -2.59e-04 -1.55e-04 -1.36e-04 2.57e-04 -3.17e-04 2.66e-04 -5.25e-05 -2.59e-04 1.39e-03 -2.36e-04 -2.36e-04 -1.22e-05 -2.86e-04 -3.57e-04 -9.71e-05 -1.55e-04 -2.36e-04 6.06e-04 1.87e-04 2.08e-05 9.99e-05 -2.05e-04 2.56e-05 -1.36e-04 -2.36e-04 1.87e-04 6.11e-04 -3.42e-05 2.16e-04 7.51e-05 4.13e-05 2.57e-04 -1.22e-05 2.08e-05 -3.42e-05 5.99e-04 1.15e-04 -1.61e-04 1.24e-04 -3.17e-04 -2.86e-04 9.99e-05 2.16e-04 1.15e-04 1.11e-03
    2.16e-03 -3.43e-06 4.80e-04 4.07e-04 -4.42e-04 5.97e-05 3.86e-05 -1.24e-04 -3.43e-06 6.41e-04 5.70e-05 -4.36e-04 -8.69e-06 -1.75e-04 1.71e-04 3.28e-05 4.80e-04 5.70e-05 7.03e-04 -7.16e-05 -1.52e-04 2.67e-05 2.27e-04 -3.14e-04 4.07e-04 -4.36e-04 -7.16e-05 1.20e-03 2.97e-05 1.79e-04 -1.72e-04 1.15e-04 -4.42e-04 -8.69e-06 -1.52e-04 2.97e-05 5.73e-04 -1.93e-04 -5.38e-05 2.75e-05 5.97e-05 -1.75e-04 2.67e-05 1.79e-04 -1.93e-04 9.61e-04 8.57e-05 -1.23e-04 3.86e-05 1.71e-04 2.27e-04 -1.72e-04 -5.38e-05 8.57e-05 2.85e-04 -1.96e-04 -1.24e-04 3.28e-05 -3.14e-04 1.15e-04 2.75e-05 -1.23e-04 -1.96e-04 7.60e-04
    5.44e-04 4.86e-05 1.22e-04 -2.04e-04 -1.15e-04 1.85e-05 3.09e-05 -4.19e-04 4.86e-05 7.12e-04 -2.10e-04 4.47e-04 -2.15e-05 -5.01e-05 -6.23e-05 -3.59e-04 1.22e-04 -2.10e-04 1.63e-03 -9.67e-04 3.19e-04 -2.25e-04 1.76e-06 1.25e-05 -2.04e-04 4.47e-04 -9.67e-04 1.99e-03 7.31e-05 6.26e-06 -2.92e-04 1.82e-04 -1.15e-04 -2.15e-05 3.19e-04 7.31e-05 5.79e-04 -1.15e-04 -7.43e-06 2.47e-04 1.85e-05 -5.01e-05 -2.25e-04 6.26e-06 -1.15e-04 4.51e-04 -2.78e-05 -2.48e-05 3.09e-05 -6.23e-05 1.76e-06 -2.92e-04 -7.43e-06 -2.78e-05 6.47e-04 -2.43e-04 -4.19e-04 -3.59e-04 1.25e-05 1.82e-04 2.47e-04 -2.48e-05 -2.43e-04 1.55e-03
    6.61e-04 -2.18e-04 1.86e-04 -2.83e-05 -2.25e-04 -2.97e-04 -5.74e-05 2.88e-04 -2.18e-04 1.34e-03 -4.84e-05 -1.84e-04 7.58e-05 2.63e-04 6.83e-05 -1.18e-04 1.86e-04 -4.84e-05 1.81e-03 -1.57e-04 -8.69e-05 -3.17e-04 6.24e-05 5.12e-05 -2.83e-05 -1.84e-04 -1.57e-04 5.54e-04 -1.98e-05 1.82e-04 9.70e-05 -1.65e-04 -2.25e-04 7.58e-05 -8.69e-05 -1.98e-05 2.64e-04 1.97e-04 1.42e-04 -6.63e-05 -2.97e-04 2.63e-04 -3.17e-04 1.82e-04 1.97e-04 8.54e-04 2.96e-04 -3.91e-04 -5.74e-05 6.83e-05 6.24e-05 9.70e-05 1.42e-04 2.96e-04 7.08e-04 1.51e-04 2.88e-04 -1.18e-04 5.12e-05 -1.65e-04 -6.63e-05 -3.91e-04 1.51e-04 1.02e-03
    1.07e-03 2.67e-05 -4.64e-04 4.04e-04 -2.40e-04 2.85e-05 -1.17e-04 2.53e-04 2.67e-05 7.25e-04 -1.16e-06 -3.14e-04 -1.22e-04 -1.74e-04 8.73e-05 3.36e-04 -4.64e-04 -1.16e-06 1.08e-03 -3.84e-04 8.33e-05 -1.20e-04 1.82e-04 -1.98e-04 4.04e-04 -3.14e-04 -3.84e-04 1.89e-03 -3.00e-04 -2.14e-04 -3.18e-04 -1.29e-04 -2.40e-04 -1.22e-04 8.33e-05 -3.00e-04 5.93e-04 3.74e-05 1.18e-04 -1.54e-04 2.85e-05 -1.74e-04 -1.20e-04 -2.14e-04 3.74e-05 1.05e-03 -2.07e-05 -4.02e-04 -1.17e-04 8.73e-05 1.82e-04 -3.18e-04 1.18e-04 -2.07e-05 4.15e-04 5.64e-05 2.53e-04 3.36e-04 -1.98e-04 -1.29e-04 -1.54e-04 -4.02e-04 5.64e-05 1.07e-03
    3.54e-04 -2.12e-04 -1.08e-04 7.14e-06 -9.37e-06 -1.30e-04 -1.85e-05 1.40e-04 -2.12e-04 1.57e-03 -4.44e-05 1.40e-05 2.78e-04 6.80e-04 2.99e-04 1.53e-05 -1.08e-04 -4.44e-05 2.85e-04 -8.50e-05 -2.32e-05 1.08e-04 -1.98e-04 2.96e-06 7.14e-06 1.40e-05 -8.50e-05 7.34e-04 -2.97e-04 -2.82e-04 -1.89e-04 -2.79e-04 -9.37e-06 2.78e-04 -2.32e-05 -2.97e-04 3.36e-04 2.89e-04 1.99e-04 1.61e-04 -1.30e-04 6.80e-04 1.08e-04 -2.82e-04 2.89e-04 1.04e-03 -5.34e-05 8.88e-05 -1.85e-05 2.99e-04 -1.98e-04 -1.89e-04 1.99e-04 -5.34e-05 1.36e-03 -1.45e-04 1.40e-04 1.53e-05 2.96e-06 -2.79e-04 1.61e-04 8.88e-05 -1.45e-04 5.00e-04
    8.71e-04 1.18e-04 4.12e-05 -8.65e-05 4.53e-04 -2.68e-04 -2.25e-04 -2.92e-04 1.18e-04 2.55e-04 -2.44e-05 -5.97e-05 -9.74e-05 1.40e-04 -1.44e-04 -1.79e-04 4.12e-05 -2.44e-05 1.45e-04 6.42e-05 -2.07e-04 -7.13e-05 -1.01e-05 3.78e-05 -8.65e-05 -5.97e-05 6.42e-05 1.31e-04 -1.66e-04 -1.28e-04 9.29e-05 7.67e-05 4.53e-04 -9.74e-05 -2.07e-04 -1.66e-04 1.76e-03 -2.11e-04 -1.31e-04 -1.99e-04 -2.68e-04 1.40e-04 -7.13e-05 -1.28e-04 -2.11e-04 1.34e-03 -2.39e-05 6.20e-05 -2.25e-04 -1.44e-04 -1.01e-05 9.29e-05 -1.31e-04 -2.39e-05 3.61e-04 1.35e-04 -2.92e-04 -1.79e-04 3.78e-05 7.67e-05 -1.99e-04 6.20e-05 1.35e-04 5.56e-04
    7.23e-04 -1.59e-04 -1.15e-04 5.36e-07 9.70e-05 -6.84e-05 8.33e-05 1.64e-04 -1.59e-04 3.88e-04 -1.92e-05 1.42e-04 -3.17e-05 5.26e-05 9.31e-05 -2.91e-04 -1.15e-04 -1.92e-05 2.77e-04 3.73e-06 1.10e-04 -7.69e-05 3.92e-05 1.00e-04 5.36e-07 1.42e-04 3.73e-06 1.08e-03 6.45e-04 2.66e-04 -2.02e-04 -2.05e-04 9.70e-05 -3.17e-05 1.10e-04 6.45e-04 1.65e-03 2.13e-04 4.27e-05 4.28e-04 -6.84e-05 5.26e-05 -7.69e-05 2.66e-04 2.13e-04 4.82e-04 -7.12e-05 -4.16e-05 8.33e-05 9.31e-05 3.92e-05 -2.02e-04 4.27e-05 -7.12e-05 3.11e-04 1.57e-04 1.64e-04 -2.91e-04 1.00e-04 -2.05e-04 4.28e-04 -4.16e-05 1.57e-04 1.20e-03
    7.15e-04 5.10e-05 2.42e-04 2.39e-04 -3.59e-04 -3.01e-04 3.64e-05 9.98e-05 5.10e-05 5.77e-04 1.34e-04 -2.77e-04 -4.67e-05 -2.11e-05 1.96e-05 2.07e-04 2.42e-04 1.34e-04 1.81e-03 -3.53e-04 -3.67e-05 -3.27e-04 2.54e-04 -2.63e-04 2.39e-04 -2.77e-04 -3.53e-04 1.40e-03 -2.47e-04 -2.53e-04 -1.06e-04 -2.01e-04 -3.59e-04 -4.67e-05 -3.67e-05 -2.47e-04 5.32e-04 2.73e-04 8.23e-05 -5.34e-06 -3.01e-04 -2.11e-05 -3.27e-04 -2.53e-04 2.73e-04 7.60e-04 1.83e-07 1.02e-04 3.64e-05 1.96e-05 2.54e-04 -1.06e-04 8.23e-05 1.83e-07 5.05e-04 1.59e-04 9.98e-05 2.07e-04 -2.63e-04 -2.01e-04 -5.34e-06 1.02e-04 1.59e-04 1.07e-03
    1.86e-03 -3.93e-04 2.21e-04 5.01e-04 -5.49e-04 -8.32e-05 -2.03e-04 -1.85e-04 -3.93e-04 1.29e-03 2.10e-04 -2.06e-04 -5.11e-05 -1.37e-04 3.18e-04 -3.13e-04 2.21e-04 2.10e-04 6.12e-04 4.23e-05 -2.07e-04 5.18e-05 1.67e-04 -2.51e-04 5.01e-04 -2.06e-04 4.23e-05 7.76e-04 -1.89e-04 1.20e-04 1.03e-05 -1.33e-04 -5.49e-04 -5.11e-05 -2.07e-04 -1.89e-04 5.22e-04 -7.28e-05 2.94e-05 1.61e-04 -8.32e-05 -1.37e-04 5.18e-05 1.20e-04 -7.28e-05 1.01e-03 1.43e-04 2.50e-05 -2.03e-04 3.18e-04 1.67e-04 1.03e-05 2.94e-05 1.43e-04 5.21e-04 -2.58e-04 -1.85e-04 -3.13e-04 -2.51e-04 -1.33e-04 1.61e-04 2.50e-05 -2.58e-04 7.05e-04
    6.69e-04 7.32e-05 1.95e-04 7.02e-05 -2.76e-04 2.75e-05 1.11e-04 -3.15e-04 7.32e-05 7.10e-04 -5.92e-05 5.30e-04 -1.78e-04 -6.99e-05 7.60e-05 -3.82e-04 1.95e-04 -5.92e-05 1.65e-03 -8.31e-04 1.39e-04 -2.13e-04 1.83e-04 -4.53e-04 7.02e-05 5.30e-04 -8.31e-04 2.12e-03 -1.95e-04 -1.42e-05 -1.89e-04 -2.78e-04 -2.76e-04 -1.78e-04 1.39e-04 -1.95e-04 5.43e-04 -4.51e-05 -6.68e-05 1.45e-04 2.75e-05 -6.99e-05 -2.13e-04 -1.42e-05 -4.51e-05 4.86e-04 -1.96e-05 2.94e-04 1.11e-04 7.60e-05 1.83e-04 -1.89e-04 -6.68e-05 -1.96e-05 5.79e-04 -9.01e-05 -3.15e-04 -3.82e-04 -4.53e-04 -2.78e-04 1.45e-04 2.94e-04 -9.01e-05 1.16e-03
    7.75e-04 -4.27e-04 -5.86e-06 -2.72e-05 -9.09e-05 -1.41e-04 1.80e-04 -5.63e-05 -4.27e-04 1.55e-03 -1.20e-04 -8.81e-06 3.76e-04 7.98e-04 -1.83e-04 7.03e-06 -5.86e-06 -1.20e-04 4.35e-04 1.62e-04 -2.09e-04 -1.40e-04 -1.43e-04 -1.00e-04 -2.72e-05 -8.81e-06 1.62e-04 7.14e-04 -2.49e-04 -2.36e-04 -1.32e-04 -1.97e-04 -9.09e-05 3.76e-04 -2.09e-04 -2.49e-04 5.71e-04 2.64e-04 2.51e-04 1.38e-04 -1.41e-04 7.98e-04 -1.40e-04 -2.36e-04 2.64e-04 1.38e-03 -3.08e-05 3.27e-05 1.80e-04 -1.83e-04 -1.43e-04 -1.32e-04 2.51e-04 -3.08e-05 6.75e-04 9.75e-05 -5.63e-05 7.03e-06 -1.00e-04 -1.97e-04 1.38e-04 3.27e-05 9.75e-05 4.22e-04
    3.62e-04 6.10e-05 1.35e-04 -3.26e-04 -1.69e-05 -1.16e-05 7.25e-05 -3.51e-04 6.10e-05 5.40e-04 -2.53e-04 4.59e-04 1.55e-05 3.36e-05 -1.45e-04 -1.90e-04 1.35e-04 -2.53e-04 1.51e-03 -8.22e-04 3.23e-04 -2.83e-04 -7.59e-05 3.18e-04 -3.26e-04 4.59e-04 -8.22e-04 1.77e-03 1.15e-04 -7.74e-05 -2.56e-04 3.98e-04 -1.69e-05 1.55e-05 3.23e-04 1.15e-04 4.70e-04 -1.22e-04 4.21e-05 2.06e-04 -1.16e-05 3.36e-05 -2.83e-04 -7.74e-05 -1.22e-04 2.84e-04 5.49e-05 -1.88e-04 7.25e-05 -1.45e-04 -7.59e-05 -2.56e-04 4.21e-05 5.49e-05 4.28e-04 -2.37e-04 -3.51e-04 -1.90e-04 3.18e-04 3.98e-04 2.06e-04 -1.88e-04 -2.37e-04 1.57e-03
    2.01e-03 -5.41e-04 4.28e-04 4.53e-04 -6.32e-04 -1.90e-04 -2.05e-04 -7.49e-05 -5.41e-04 1.61e-03 2.59e-04 1.64e-04 -1.01e-04 -2.37e-04 4.09e-04 -2.69e-04 4.28e-04 2.59e-04 6.71e-04 1.64e-04 -3.13e-04 -1.86e-04 1.26e-04 -1.79e-04 4.53e-04 1.64e-04 1.64e-04 6.55e-04 -2.67e-04 -8.42e-05 1.49e-04 -1.38e-04 -6.32e-04 -1.01e-04 -3.13e-04 -2.67e-04 5.79e-04 1.08e-04 -9.34e-06 1.35e-04 -1.90e-04 -2.37e-04 -1.86e-04 -8.42e-05 1.08e-04 1.03e-03 8.11e-05 2.35e-04 -2.05e-04 4.09e-04 1.26e-04 1.49e-04 -9.34e-06 8.11e-05 5.33e-04 -1.42e-04 -7.49e-05 -2.69e-04 -1.79e-04 -1.38e-04 1.35e-04 2.35e-04 -1.42e-04 5.50e-04
    1.35e-03 -4.57e-04 -1.83e-04 3.55e-04 -4.47e-04 -7.28e-05 -1.90e-04 -2.13e-04 -4.57e-04 1.76e-03 2.54e-04 9.91e-05 -6.97e-05 -7.82e-05 3.99e-04 -3.82e-04 -1.83e-04 2.54e-04 6.05e-04 -2.74e-05 -3.92e-05 1.85e-04 1.17e-04 -1.35e-04 3.55e-04 9.91e-05 -2.74e-05 7.06e-04 -2.35e-04 9.13e-05 1.49e-04 -3.21e-04 -4.47e-04 -6.97e-05 -3.92e-05 -2.35e-04 4.28e-04 6.29e-07 9.01e-05 2.48e-04 -7.28e-05 -7.82e-05 1.85e-04 9.13e-05 6.29e-07 9.40e-04 2.41e-04 -5.37e-05 -1.90e-04 3.99e-04 1.17e-04 1.49e-04 9.01e-05 2.41e-04 6.87e-04 -1.35e-04 -2.13e-04 -3.82e-04 -1.35e-04 -3.21e-04 2.48e-04 -5.37e-05 -1.35e-04 7.46e-04
    5.64e-04 -2.45e-04 -1.08e-04 -1.18e-05 -1.47e-04 -1.99e-04 -3.39e-05 2.82e-04 -2.45e-04 1.68e-03 -2.74e-05 6.72e-05 1.13e-05 2.45e-04 1.70e-04 -2.82e-04 -1.08e-04 -2.74e-05 1.47e-03 -2.06e-04 2.66e-05 -1.95e-04 6.21e-05 1.47e-04 -1.18e-05 6.72e-05 -2.06e-04 5.44e-04 -3.55e-05 2.73e-04 1.41e-04 -2.15e-04 -1.47e-04 1.13e-05 2.66e-05 -3.55e-05 1.88e-04 7.34e-05 1.21e-04 -1.37e-05 -1.99e-04 2.45e-04 -1.95e-04 2.73e-04 7.34e-05 8.46e-04 3.46e-04 -4.85e-04 -3.39e-05 1.70e-04 6.21e-05 1.41e-04 1.21e-04 3.46e-04 7.19e-04 4.57e-05 2.82e-04 -2.82e-04 1.47e-04 -2.15e-04 -1.37e-05 -4.85e-04 4.57e-05 8.99e-04
    2.07e-03 3.10e-05 4.20e-04 2.25e-04 -3.16e-04 2.58e-04 1.81e-04 -9.84e-05 3.10e-05 3.67e-04 -1.86e-05 -2.30e-04 5.10e-05 -9.11e-05 1.07e-04 3.39e-05 4.20e-04 -1.86e-05 4.84e-04 -7.85e-05 -5.52e-05 9.91e-05 1.53e-04 -1.70e-04 2.25e-04 -2.30e-04 -7.85e-05 1.27e-03 -2.37e-05 1.23e-04 -1.41e-04 1.64e-04 -3.16e-04 5.10e-05 -5.52e-05 -2.37e-05 4.10e-04 -1.93e-04 2.84e-05 4.98e-05 2.58e-04 -9.11e-05 9.91e-05 1.23e-04 -1.93e-04 6.94e-04 8.02e-05 -2.23e-04 1.81e-04 1.07e-04 1.53e-04 -1.41e-04 2.84e-05 8.02e-05 1.84e-04 -1.65e-05 -9.84e-05 3.39e-05 -1.70e-04 1.64e-04 4.98e-05 -2.23e-04 -1.65e-05 5.61e-04
    7.50e-04 -4.85e-04 -2.18e-04 -2.63e-04 -2.95e-06 1.69e-04 -3.58e-04 9.84e-05 -4.85e-04 1.67e-03 2.54e-04 4.35e-05 3.33e-04 -5.74e-05 -4.82e-05 2.04e-04 -2.18e-04 2.54e-04 7.21e-04 -9.16e-05 3.71e-05 1.87e-04 -1.07e-05 -1.15e-04 -2.63e-04 4.35e-05 -9.16e-05 6.49e-04 -1.33e-04 -2.08e-04 2.96e-04 1.66e-05 -2.95e-06 3.33e-04 3.71e-05 -1.33e-04 3.39e-04 -6.05e-05 -9.81e-05 3.40e-05 1.69e-04 -5.74e-05 1.87e-04 -2.08e-04 -6.05e-05 1.12e-03 -6.88e-04 1.13e-06 -3.58e-04 -4.82e-05 -1.07e-05 2.96e-04 -9.81e-05 -6.88e-04 1.44e-03 -2.82e-04 9.84e-05 2.04e-04 -1.15e-04 1.66e-05 3.40e-05 1.13e-06 -2.82e-04 4.66e-04
    7.84e-04 -3.78e-04 5.78e-04 7.35e-05 -3.66e-04 -3.07e-04 -2.64e-04 -1.21e-04 -3.78e-04 1.06e-03 -4.92e-04 9.72e-05 1.01e-04 3.15e-04 -2.18e-04 6.43e-05 5.78e-04 -4.92e-04 1.36e-03 8.09e-07 -2.47e-04 -1.36e-04 -1.18e-05 -1.82e-04 7.35e-05 9.72e-05 8.09e-07 5.33e-04 -8.44e-05 2.00e-04 -1.95e-04 -2.88e-04 -3.66e-04 1.01e-04 -2.47e-04 -8.44e-05 3.59e-04 1.98e-04 3.16e-04 9.57e-05 -3.07e-04 3.15e-04 -1.36e-04 2.00e-04 1.98e-04 6.72e-04 1.24e-04 -1.26e-04 -2.64e-04 -2.18e-04 -1.18e-05 -1.95e-04 3.16e-04 1.24e-04 1.26e-03 2.70e-04 -1.21e-04 6.43e-05 -1.82e-04 -2.88e-04 9.57e-05 -1.26e-04 2.70e-04 8.91e-04
    7.54e-04 -7.17e-05 -3.46e-04 1.89e-04 -4.30e-05 1.94e-04 1.13e-05 2.81e-04 -7.17e-05 6.79e-04 1.48e-04 5.32e-04 -8.62e-05 -2.74e-04 2.05e-04 -3.31e-04 -3.46e-04 1.48e-04 8.08e-04 -2.07e-04 9.58e-05 -2.01e-04 1.16e-04 -3.65e-04 1.89e-04 5.32e-04 -2.07e-04 1.78e-03 1.37e-05 5.48e-05 -1.76e-05 -3.40e-04 -4.30e-05 -8.62e-05 9.58e-05 1.37e-05 6.36e-04 1.45e-05 -9.51e-05 -9.92e-05 1.94e-04 -2.74e-04 -2.01e-04 5.48e-05 1.45e-05 5.97e-04 -1.53e-04 2.61e-04 1.13e-05 2.05e-04 1.16e-04 -1.76e-05 -9.51e-05 -1.53e-04 2.69e-04 -2.07e-05 2.81e-04 -3.31e-04 -3.65e-04 -3.40e-04 -9.92e-05 2.61e-04 -2.07e-05 9.11e-04
    6.81e-04 -2.42e-04 1.36e-04 -3.67e-04 -5.06e-05 5.84e-05 3.32e-04 -2.93e-04 -2.42e-04 7.87e-04 1.79e-04 -9.60e-05 4.35e-05 -2.74e-04 -2.60e-04 1.44e-04 1.36e-04 1.79e-04 1.60e-03 -4.53e-04 3.38e-04 2.78e-04 -3.44e-04 -1.68e-04 -3.67e-04 -9.60e-05 -4.53e-04 6.57e-04 -3.74e-05 3.14e-05 -1.07e-04 2.20e-04 -5.06e-05 4.35e-05 3.38e-04 -3.74e-05 3.88e-04 5.14e-05 3.25e-05 -1.43e-04 5.84e-05 -2.74e-04 2.78e-04 3.14e-05 5.14e-05 6.14e-04 -4.84e-05 -3.44e-05 3.32e-04 -2.60e-04 -3.44e-04 -1.07e-04 3.25e-05 -4.84e-05 1.47e-03 -4.92e-04 -2.93e-04 1.44e-04 -1.68e-04 2.20e-04 -1.43e-04 -3.44e-05 -4.92e-04 6.60e-04
    6.78e-04 9.79e-05 -3.69e-05 1.83e-04 -2.91e-04 -5.26e-05 7.67e-05 -6.56e-05 9.79e-05 7.63e-04 1.06e-05 3.65e-04 -2.32e-04 -1.81e-04 1.14e-04 -2.56e-04 -3.69e-05 1.06e-05 1.49e-03 -7.32e-04 1.28e-04 -3.13e-04 1.89e-04 -4.15e-04 1.83e-04 3.65e-04 -7.32e-04 2.12e-03 -1.77e-04 -5.06e-05 -2.59e-04 -3.49e-04 -2.91e-04 -2.32e-04 1.28e-04 -1.77e-04 6.35e-04 7.04e-05 -2.66e-05 8.53e-05 -5.26e-05 -1.81e-04 -3.13e-04 -5.06e-05 7.04e-05 7.00e-04 -9.40e-05 3.74e-04 7.67e-05 1.14e-04 1.89e-04 -2.59e-04 -2.66e-05 -9.40e-05 6.11e-04 1.38e-05 -6.56e-05 -2.56e-04 -4.15e-04 -3.49e-04 8.53e-05 3.74e-04 1.38e-05 1.14e-03
    1.98e-03 1.79e-05 4.33e-04 2.12e-04 -3.53e-04 7.42e-05 7.14e-05 -1.57e-04 1.79e-05 5.98e-04 -6.62e-05 -3.27e-04 3.20e-05 -1.67e-04 1.13e-04 1.58e-04 4.33e-04 -6.62e-05 6.69e-04 -1.11e-04 -9.21e-05 7.50e-05 1.90e-04 -3.19e-04 2.12e-04 -3.27e-04 -1.11e-04 1.43e-03 7.27e-05 1.02e-04 -2.10e-04 1.44e-04 -3.53e-04 3.20e-05 -9.21e-05 7.27e-05 5.64e-04 -2.33e-04 -2.48e-05 -3.23e-05 7.42e-05 -1.67e-04 7.50e-05 1.02e-04 -2.33e-04 9.61e-04 6.85e-05 -2.24e-04 7.14e-05 1.13e-04 1.90e-04 -2.10e-04 -2.48e-05 6.85e-05 2.33e-04 -1.32e-04 -1.57e-04 1.58e-04 -3.19e-04 1.44e-04 -3.23e-05 -2.24e-04 -1.32e-04 8.25e-04
    5.61e-04 1.92e-04 3.07e-04 2.02e-04 -2.81e-04 -5.18e-05 1.51e-04 -2.29e-04 1.92e-04 4.86e-04 8.54e-06 5.17e-04 -2.55e-04 -6.40e-05 6.16e-05 -1.73e-04 3.07e-04 8.54e-06 1.56e-03 -6.80e-04 3.02e-05 -1.63e-04 2.69e-04 -4.42e-04 2.02e-04 5.17e-04 -6.80e-04 2.00e-03 -3.69e-04 -1.24e-04 -1.27e-04 -3.16e-04 -2.81e-04 -2.55e-04 3.02e-05 -3.69e-04 3.87e-04 3.95e-05 -7.47e-05 1.11e-04 -5.18e-05 -6.40e-05 -1.63e-04 -1.24e-04 3.95e-05 3.63e-04 5.78e-07 3.21e-04 1.51e-04 6.16e-05 2.69e-04 -1.27e-04 -7.47e-05 5.78e-07 3.72e-04 2.92e-06 -2.29e-04 -1.73e-04 -4.42e-04 -3.16e-04 1.11e-04 3.21e-04 2.92e-06 9.83e-04
    3.28e-04 1.45e-04 5.34e-04 2.04e-04 -1.72e-04 -4.70e-05 1.01e-04 3.41e-05 1.45e-04 2.24e-04 2.29e-04 1.78e-04 -9.97e-05 7.88e-06 5.77e-05 1.14e-04 5.34e-04 2.29e-04 1.75e-03 -8.49e-05 -2.91e-04 -1.51e-04 2.58e-04 -1.29e-04 2.04e-04 1.78e-04 -8.49e-05 1.11e-03 -2.24e-04 -1.58e-04 3.97e-05 -1.67e-04 -1.72e-04 -9.97e-05 -2.91e-04 -2.24e-04 1.24e-04 4.57e-05 -5.73e-05 -1.40e-05 -4.70e-05 7.88e-06 -1.51e-04 -1.58e-04 4.57e-05 9.81e-05 -2.64e-05 1.21e-04 1.01e-04 5.77e-05 2.58e-04 3.97e-05 -5.73e-05 -2.64e-05 7.09e-05 -1.49e-05 3.41e-05 1.14e-04 -1.29e-04 -1.67e-04 -1.40e-05 1.21e-04 -1.49e-05 5.00e-04
    1.73e-03 -8.58e-04 2.24e-04 1.50e-04 -3.47e-04 -2.12e-04 -3.38e-04 3.00e-05 -8.58e-04 1.53e-03 9.07e-05 8.65e-05 1.17e-04 -2.47e-04 2.51e-04 -1.50e-04 2.24e-04 9.07e-05 7.64e-04 1.70e-04 -2.97e-04 -4.24e-04 -3.42e-05 -2.00e-04 1.50e-04 8.65e-05 1.70e-04 8.42e-04 -2.68e-04 -4.50e-04 3.54e-05 -2.20e-04 -3.47e-04 1.17e-04 -2.97e-04 -2.68e-04 8.21e-04 2.98e-04 1.16e-04 1.33e-04 -2.12e-04 -2.47e-04 -4.24e-04 -4.50e-04 2.98e-04 1.32e-03 -7.04e-05 3.31e-04 -3.38e-04 2.51e-04 -3.42e-05 3.54e-05 1.16e-04 -7.04e-05 7.97e-04 -4.67e-05 3.00e-05 -1.50e-04 -2.00e-04 -2.20e-04 1.33e-04 3.31e-04 -4.67e-05 6.25e-04
    1.62e-03 -7.34e-04 1.93e-04 -7.34e-05 -7.94e-05 2.94e-04 -4.46e-04 1.71e-04 -7.34e-04 1.21e-03 6.51e-05 -1.33e-05 2.32e-05 1.61e-05 1.98e-04 -1.75e-04 1.93e-04 6.51e-05 6.95e-04 1.78e-04 -1.90e-04 -5.26e-04 -7.73e-05 -1.87e-04 -7.34e-05 -1.33e-05 1.78e-04 6.24e-04 -1.17e-04 -4.59e-04 1.42e-06 -8.07e-05 -7.94e-05 2.32e-05 -1.90e-04 -1.17e-04 6.39e-04 3.75e-04 1.21e-04 2.40e-05 2.94e-04 1.61e-05 -5.26e-04 -4.59e-04 3.75e-04 1.55e-03 -1.67e-04 1.07e-04 -4.46e-04 1.98e-04 -7.73e-05 1.42e-06 1.21e-04 -1.67e-04 6.14e-04 -9.98e-05 1.71e-04 -1.75e-04 -1.87e-04 -8.07e-05 2.40e-05 1.07e-04 -9.98e-05 4.49e-04
    6.90e-04 -3.48e-04 7.81e-05 -7.34e-05 -2.33e-04 -1.94e-04 -1.40e-04 -2.71e-04 -3.48e-04 8.79e-04 -9.56e-05 -1.95e-05 7.99e-05 8.55e-05 -1.19e-04 2.59e-04 7.81e-05 -9.56e-05 1.52e-03 -2.95e-04 2.19e-04 -3.15e-05 7.71e-04 4.80e-05 -7.34e-05 -1.95e-05 -2.95e-04 6.10e-04 -6.50e-05 1.63e-04 -4.02e-04 -1.39e-04 -2.33e-04 7.99e-05 2.19e-04 -6.50e-05 4.57e-04 1.72e-04 3.80e-04 1.33e-04 -1.94e-04 8.55e-05 -3.15e-05 1.63e-04 1.72e-04 6.06e-04 -1.11e-04 -7.19e-05 -1.40e-04 -1.19e-04 7.71e-04 -4.02e-04 3.80e-04 -1.11e-04 1.58e-03 1.70e-04 -2.71e-04 2.59e-04 4.80e-05 -1.39e-04 1.33e-04 -7.19e-05 1.70e-04 7.27e-04
    1.74e-03 -8.53e-04 3.41e-04 3.09e-04 -5.08e-04 -3.00e-04 -3.00e-04 -1.86e-05 -8.53e-04 1.57e-03 1.98e-04 1.85e-04 2.75e-05 -3.60e-04 4.02e-04 -2.35e-04 3.41e-04 1.98e-04 7.25e-04 1.09e-04 -2.82e-04 -3.49e-04 4.51e-05 -1.82e-04 3.09e-04 1.85e-04 1.09e-04 8.04e-04 -2.91e-04 -3.34e-04 7.76e-05 -2.20e-04 -5.08e-04 2.75e-05 -2.82e-04 -2.91e-04 6.44e-04 2.50e-04 2.96e-05 1.87e-04 -3.00e-04 -3.60e-04 -3.49e-04 -3.34e-04 2.50e-04 1.10e-03 2.33e-05 3.35e-04 -3.00e-04 4.02e-04 4.51e-05 7.76e-05 2.96e-05 2.33e-05 6.79e-04 -4.68e-05 -1.86e-05 -2.35e-04 -1.82e-04 -2.20e-04 1.87e-04 3.35e-04 -4.68e-05 5.79e-04
    4.34e-04 -2.14e-05 -6.74e-05 -1.60e-04 -4.87e-05 2.26e-05 2.59e-05 -2.14e-04 -2.14e-05 6.27e-04 1.02e-04 1.93e-04 1.42e-04 -1.12e-04 4.71e-05 -3.29e-04 -6.74e-05 1.02e-04 1.01e-03 -5.31e-04 1.12e-04 -3.30e-04 4.82e-05 -1.33e-04 -1.60e-04 1.93e-04 -5.31e-04 1.63e-03 4.13e-04 1.67e-04 -3.46e-04 3.73e-04 -4.87e-05 1.42e-04 1.12e-04 4.13e-04 6.58e-04 -6.62e-05 6.71e-06 1.56e-04 2.26e-05 -1.12e-04 -3.30e-04 1.67e-04 -6.62e-05 5.03e-04 -1.16e-04 1.03e-04 2.59e-05 4.71e-05 4.82e-05 -3.46e-04 6.71e-06 -1.16e-04 5.46e-04 -2.53e-04 -2.14e-04 -3.29e-04 -1.33e-04 3.73e-04 1.56e-04 1.03e-04 -2.53e-04 1.37e-03
    1.03e-03 -6.83e-04 6.90e-05 2.08e-04 -3.91e-04 -2.51e-04 -1.40e-04 -1.76e-04 -6.83e-04 1.62e-03 9.77e-05 4.43e-04 3.71e-05 -2.51e-04 1.19e-04 2.64e-05 6.90e-05 9.77e-05 5.28e-04 1.37e-04 -2.12e-04 -1.77e-04 -4.96e-05 -6.45e-05 2.08e-04 4.43e-04 1.37e-04 9.10e-04 -4.31e-04 -6.09e-04 -6.46e-05 -2.59e-04 -3.91e-04 3.71e-05 -2.12e-04 -4.31e-04 5.43e-04 3.93e-04 1.37e-04 2.17e-04 -2.51e-04 -2.51e-04 -1.77e-04 -6.09e-04 3.93e-04 1.10e-03 6.83e-05 2.55e-04 -1.40e-04 1.19e-04 -4.96e-05 -6.46e-05 1.37e-04 6.83e-05 8.45e-04 1.15e-04 -1.76e-04 2.64e-05 -6.45e-05 -2.59e-04 2.17e-04 2.55e-04 1.15e-04 5.50e-04
    1.80e-03 6.61e-05 3.06e-04 -1.60e-04 -3.09e-04 7.46e-04 -4.09e-05 1.45e-04 6.61e-05 4.18e-04 -4.13e-05 -1.18e-04 -4.58e-07 2.31e-05 2.62e-06 9.34e-05 3.06e-04 -4.13e-05 6.68e-04 -5.11e-05 -3.24e-04 -8.13e-05 1.29e-04 -1.57e-04 -1.60e-04 -1.18e-04 -5.11e-05 4.81e-04 1.96e-04 3.04e-05 -3.18e-06 -5.51e-05 -3.09e-04 -4.58e-07 -3.24e-04 1.96e-04 1.43e-03 3.67e-04 -3.01e-05 8.50e-05 7.46e-04 2.31e-05 -8.13e-05 3.04e-05 3.67e-04 1.18e-03 1.07e-06 2.44e-05 -4.09e-05 2.62e-06 1.29e-04 -3.18e-06 -3.01e-05 1.07e-06 2.22e-04 -1.36e-04 1.45e-04 9.34e-05 -1.57e-04 -5.51e-05 8.50e-05 2.44e-05 -1.36e-04 4.75e-04
    4.93e-04 1.20e-04 1.24e-04 -6.02e-05 -1.17e-04 -2.46e-05 -9.93e-05 -5.02e-04 1.20e-04 5.34e-04 -2.98e-04 3.26e-04 -1.39e-04 -3.63e-05 -1.21e-04 -3.24e-04 1.24e-04 -2.98e-04 1.65e-03 -7.17e-04 3.57e-04 -9.13e-05 -1.53e-04 3.86e-04 -6.02e-05 3.26e-04 -7.17e-04 1.35e-03 -7.80e-05 -6.59e-05 -8.09e-05 4.52e-07 -1.17e-04 -1.39e-04 3.57e-04 -7.80e-05 4.53e-04 -1.55e-05 1.86e-06 2.97e-04 -2.46e-05 -3.63e-05 -9.13e-05 -6.59e-05 -1.55e-05 3.57e-04 1.48e-05 -8.97e-05 -9.93e-05 -1.21e-04 -1.53e-04 -8.09e-05 1.86e-06 1.48e-05 6.11e-04 -5.39e-05 -5.02e-04 -3.24e-04 3.86e-04 4.52e-07 2.97e-04 -8.97e-05 -5.39e-05 1.35e-03
    4.36e-04 -9.22e-06 -1.53e-04 1.46e-04 -1.21e-04 -1.68e-04 2.21e-05 8.61e-05 -9.22e-06 1.23e-03 2.23e-05 5.88e-04 -2.58e-04 -2.01e-04 -7.72e-05 -9.91e-05 -1.53e-04 2.23e-05 6.78e-04 -2.12e-04 6.18e-05 2.10e-04 -1.36e-04 1.80e-04 1.46e-04 5.88e-04 -2.12e-04 9.39e-04 -4.03e-04 -3.34e-04 -1.99e-04 -3.80e-04 -1.21e-04 -2.58e-04 6.18e-05 -4.03e-04 3.52e-04 2.01e-04 1.69e-04 2.10e-04 -1.68e-04 -2.01e-04 2.10e-04 -3.34e-04 2.01e-04 7.62e-04 2.98e-04 7.31e-05 2.21e-05 -7.72e-05 -1.36e-04 -1.99e-04 1.69e-04 2.98e-04 1.09e-03 9.95e-05 8.61e-05 -9.91e-05 1.80e-04 -3.80e-04 2.10e-04 7.31e-05 9.95e-05 7.26e-04
    7.36e-04 4.79e-05 -1.45e-04 2.73e-04 -2.95e-04 -1.51e-04 -1.17e-05 6.71e-05 4.79e-05 6.98e-04 2.38e-05 7.58e-05 -2.05e-04 -2.07e-04 7.73e-05 -4.43e-05 -1.45e-04 2.38e-05 1.47e-03 -6.52e-04 1.24e-04 -2.52e-04 1.73e-04 -3.54e-04 2.73e-04 7.58e-05 -6.52e-04 1.93e-03 -1.61e-04 -1.31e-04 -2.75e-04 -3.30e-04 -2.95e-04 -2.05e-04 1.24e-04 -1.61e-04 6.36e-04 2.02e-04 4.78e-05 2.48e-05 -1.51e-04 -2.07e-04 -2.52e-04 -1.31e-04 2.02e-04 8.11e-04 -8.94e-05 2.67e-04 -1.17e-05 7.73e-05 1.73e-04 -2.75e-04 4.78e-05 -8.94e-05 5.86e-04 1.32e-04 6.71e-05 -4.43e-05 -3.54e-04 -3.30e-04 2.48e-05 2.67e-04 1.32e-04 1.09e-03
    1.37e-03 -1.97e-05 -3.42e-04 4.19e-04 -2.40e-04 1.31e-04 -1.23e-04 7.87e-05 -1.97e-05 6.95e-04 -7.64e-05 -3.42e-04 -7.98e-05 -1.74e-04 8.73e-05 3.34e-04 -3.42e-04 -7.64e-05 7.02e-04 -2.24e-04 2.04e-05 1.53e-04 1.48e-04 -1.99e-04 4.19e-04 -3.42e-04 -2.24e-04 1.74e-03 -2.33e-04 -1.86e-04 -3.15e-04 -6.21e-05 -2.40e-04 -7.98e-05 2.04e-05 -2.33e-04 5.24e-04 -9.77e-05 8.90e-05 -1.44e-04 1.31e-04 -1.74e-04 1.53e-04 -1.86e-04 -9.77e-05 9.97e-04 6.11e-05 -4.98e-04 -1.23e-04 8.73e-05 1.48e-04 -3.15e-04 8.90e-05 6.11e-05 3.13e-04 -4.53e-05 7.87e-05 3.34e-04 -1.99e-04 -6.21e-05 -1.44e-04 -4.98e-04 -4.53e-05 9.74e-04
    3.20e-04 8.45e-05 1.46e-04 -3.16e-04 3.22e-05 -8.34e-06 3.43e-05 -1.98e-04 8.45e-05 5.28e-04 -2.58e-04 4.24e-04 1.39e-06 3.39e-06 -1.63e-04 -6.42e-05 1.46e-04 -2.58e-04 1.51e-03 -6.51e-04 2.51e-04 -2.93e-04 -1.82e-04 4.97e-04 -3.16e-04 4.24e-04 -6.51e-04 1.51e-03 9.98e-05 -1.09e-04 -1.10e-04 3.95e-04 3.22e-05 1.39e-06 2.51e-04 9.98e-05 4.52e-04 -6.65e-05 7.42e-05 1.71e-04 -8.34e-06 3.39e-06 -2.93e-04 -1.09e-04 -6.65e-05 2.50e-04 1.84e-05 -2.37e-04 3.43e-05 -1.63e-04 -1.82e-04 -1.10e-04 7.42e-05 1.84e-05 4.64e-04 -1.51e-04 -1.98e-04 -6.42e-05 4.97e-04 3.95e-04 1.71e-04 -2.37e-04 -1.51e-04 1.51e-03
    1.72e-03 2.27e-05 3.78e-04 -9.53e-06 -2.21e-04 9.61e-05 1.03e-04 -1.39e-04 2.27e-05 4.77e-04 -1.20e-04 -9.15e-05 4.22e-05 -2.31e-04 5.84e-05 3.48e-04 3.78e-04 -1.20e-04 5.80e-04 -1.54e-04 -3.79e-05 1.24e-04 1.44e-04 -3.53e-04 -9.53e-06 -9.15e-05 -1.54e-04 1.58e-03 1.25e-04 -2.10e-05 -1.92e-04 1.09e-04 -2.21e-04 4.22e-05 -3.79e-05 1.25e-04 5.81e-04 -2.75e-04 -1.13e-05 -9.99e-05 9.61e-05 -2.31e-04 1.24e-04 -2.10e-05 -2.75e-04 9.53e-04 4.18e-05 -1.91e-04 1.03e-04 5.84e-05 1.44e-04 -1.92e-04 -1.13e-05 4.18e-05 2.01e-04 -5.34e-05 -1.39e-04 3.48e-04 -3.53e-04 1.09e-04 -9.99e-05 -1.91e-04 -5.34e-05 8.10e-04
    1.69e-03 5.37e-05 3.53e-04 -9.19e-05 -5.64e-04 3.94e-04 -5.68e-05 2.56e-04 5.37e-05 4.05e-04 -1.03e-04 -8.43e-05 4.63e-05 4.42e-05 2.60e-05 1.09e-04 3.53e-04 -1.03e-04 6.52e-04 -8.55e-05 -3.45e-04 -1.82e-04 1.18e-04 -1.39e-04 -9.19e-05 -8.43e-05 -8.55e-05 5.74e-04 2.25e-04 1.40e-04 -3.70e-05 -6.54e-05 -5.64e-04 4.63e-05 -3.45e-04 2.25e-04 1.45e-03 3.77e-04 1.38e-05 3.44e-05 3.94e-04 4.42e-05 -1.82e-04 1.40e-04 3.77e-04 9.78e-04 -6.65e-06 7.44e-06 -5.68e-05 2.60e-05 1.18e-04 -3.70e-05 1.38e-05 -6.65e-06 2.18e-04 -1.31e-04 2.56e-04 1.09e-04 -1.39e-04 -6.54e-05 3.44e-05 7.44e-06 -1.31e-04 5.16e-04
    8.08e-04 -2.97e-04 -4.30e-04 1.26e-04 -2.08e-04 -7.47e-05 -1.37e-04 7.35e-05 -2.97e-04 1.65e-03 1.01e-04 -2.68e-05 -8.80e-06 9.38e-05 2.68e-04 -3.71e-04 -4.30e-04 1.01e-04 9.53e-04 -2.00e-04 1.33e-04 2.89e-05 9.00e-05 2.89e-05 1.26e-04 -2.68e-05 -2.00e-04 6.59e-04 -9.27e-05 2.58e-04 1.49e-04 -2.51e-04 -2.08e-04 -8.80e-06 1.33e-04 -9.27e-05 2.47e-04 6.06e-07 1.31e-04 8.63e-05 -7.47e-05 9.38e-05 2.89e-05 2.58e-04 6.06e-07 9.01e-04 3.08e-04 -4.84e-04 -1.37e-04 2.68e-04 9.00e-05 1.49e-04 1.31e-04 3.08e-04 6.70e-04 -8.40e-05 7.35e-05 -3.71e-04 2.89e-05 -2.51e-04 8.63e-05 -4.84e-04 -8.40e-05 8.86e-04
    1.55e-03 -3.04e-05 3.36e-04 -3.90e-04 -7.65e-04 1.03e-04 9.40e-07 1.73e-04 -3.04e-05 4.56e-04 -6.14e-05 3.07e-06 3.51e-05 -1.20e-04 1.50e-05 1.18e-04 3.36e-04 -6.14e-05 6.51e-04 -1.97e-04 -8.71e-05 -1.13e-04 1.70e-04 -1.87e-04 -3.90e-04 3.07e-06 -1.97e-04 8.83e-04 2.10e-04 2.05e-04 -5.44e-05 -2.21e-05 -7.65e-04 3.51e-05 -8.71e-05 2.10e-04 1.25e-03 -1.42e-04 8.14e-05 -6.14e-06 1.03e-04 -1.20e-04 -1.13e-04 2.05e-04 -1.42e-04 7.26e-04 -6.20e-05 -4.01e-05 9.40e-07 1.50e-05 1.70e-04 -5.44e-05 8.14e-05 -6.20e-05 2.29e-04 -5.92e-05 1.73e-04 1.18e-04 -1.87e-04 -2.21e-05 -6.14e-06 -4.01e-05 -5.92e-05 5.46e-04
    7.00e-04 1.80e-04 3.29e-04 3.00e-04 -3.50e-04 -1.85e-04 1.28e-04 -4.05e-05 1.80e-04 5.42e-04 1.51e-04 1.13e-04 -1.84e-04 -2.93e-05 7.93e-05 9.70e-05 3.29e-04 1.51e-04 1.81e-03 -4.39e-04 -4.60e-05 -2.64e-04 3.43e-04 -3.55e-04 3.00e-04 1.13e-04 -4.39e-04 1.77e-03 -4.06e-04 -2.78e-04 -1.33e-04 -2.82e-04 -3.50e-04 -1.84e-04 -4.60e-05 -4.06e-04 4.53e-04 1.94e-04 -1.51e-05 4.18e-05 -1.85e-04 -2.93e-05 -2.64e-04 -2.78e-04 1.94e-04 5.83e-04 3.16e-06 2.35e-04 1.28e-04 7.93e-05 3.43e-04 -1.33e-04 -1.51e-05 3.16e-06 4.09e-04 7.35e-05 -4.05e-05 9.70e-05 -3.55e-04 -2.82e-04 4.18e-05 2.35e-04 7.35e-05 1.11e-03
    1.61e-03 -1.17e-05 -2.23e-05 2.95e-04 -1.78e-04 2.22e-04 -1.46e-05 8.46e-06 -1.17e-05 6.20e-04 -1.55e-04 -1.73e-04 -5.08e-05 -1.84e-04 1.12e-04 2.67e-04 -2.23e-05 -1.55e-04 5.36e-04 -8.63e-05 1.34e-05 2.46e-04 1.16e-04 -2.21e-04 2.95e-04 -1.73e-04 -8.63e-05 1.79e-03 -8.34e-05 -1.40e-04 -2.79e-04 8.30e-05 -1.78e-04 -5.08e-05 1.34e-05 -8.34e-05 4.55e-04 -1.15e-04 1.45e-05 -6.98e-05 2.22e-04 -1.84e-04 2.46e-04 -1.40e-04 -1.15e-04 9.61e-04 2.88e-05 -4.06e-04 -1.46e-05 1.12e-04 1.16e-04 -2.79e-04 1.45e-05 2.88e-05 2.50e-04 -4.35e-05 8.46e-06 2.67e-04 -2.21e-04 8.30e-05 -6.98e-05 -4.06e-04 -4.35e-05 9.68e-04
    1.61e-03 1.44e-04 2.71e-04 -1.86e-04 -6.55e-04 6.51e-04 -2.20e-05 1.12e-04 1.44e-04 1.99e-04 1.09e-04 -1.23e-04 -1.37e-05 -8.13e-06 3.16e-05 4.12e-05 2.71e-04 1.09e-04 4.45e-04 -1.16e-04 -1.89e-04 -1.79e-04 1.03e-04 -1.34e-04 -1.86e-04 -1.23e-04 -1.16e-04 3.53e-04 2.15e-04 3.54e-07 -2.55e-05 2.08e-05 -6.55e-04 -1.37e-05 -1.89e-04 2.15e-04 1.54e-03 2.42e-04 5.83e-05 1.21e-05 6.51e-04 -8.13e-06 -1.79e-04 3.54e-07 2.42e-04 1.04e-03 -1.38e-05 2.42e-05 -2.20e-05 3.16e-05 1.03e-04 -2.55e-05 5.83e-05 -1.38e-05 8.92e-05 -8.71e-05 1.12e-04 4.12e-05 -1.34e-04 2.08e-05 1.21e-05 2.42e-05 -8.71e-05 2.31e-04
    1.81e-03 -1.46e-04 1.05e-04 4.49e-04 -4.36e-04 2.78e-05 -1.39e-04 -1.73e-04 -1.46e-04 9.00e-04 4.40e-05 -4.97e-04 -5.29e-06 -1.18e-04 1.81e-04 -2.33e-05 1.05e-04 4.40e-05 6.57e-04 -9.45e-05 -1.30e-04 1.95e-04 1.72e-04 -2.77e-04 4.49e-04 -4.97e-04 -9.45e-05 1.19e-03 -1.10e-04 9.65e-05 -2.02e-04 -1.29e-05 -4.36e-04 -5.29e-06 -1.30e-04 -1.10e-04 4.59e-04 -1.62e-04 3.74e-05 4.46e-05 2.78e-05 -1.18e-04 1.95e-04 9.65e-05 -1.62e-04 1.03e-03 1.23e-04 -2.76e-04 -1.39e-04 1.81e-04 1.72e-04 -2.02e-04 3.74e-05 1.23e-04 4.11e-04 -2.31e-04 -1.73e-04 -2.33e-05 -2.77e-04 -1.29e-05 4.46e-05 -2.76e-04 -2.31e-04 9.17e-04
    6.11e-04 1.71e-04 2.87e-04 2.32e-04 -3.29e-04 -6.38e-05 1.39e-04 -1.77e-04 1.71e-04 4.91e-04 3.62e-05 3.34e-04 -2.60e-04 -3.75e-05 1.43e-04 -7.38e-05 2.87e-04 3.62e-05 1.63e-03 -5.34e-04 2.18e-05 -2.06e-04 2.60e-04 -4.08e-04 2.32e-04 3.34e-04 -5.34e-04 1.84e-03 -3.43e-04 -2.42e-04 -1.06e-04 -2.87e-04 -3.29e-04 -2.60e-04 2.18e-05 -3.43e-04 4.64e-04 6.29e-05 1.56e-05 1.23e-04 -6.38e-05 -3.75e-05 -2.06e-04 -2.42e-04 6.29e-05 4.05e-04 -6.90e-06 3.28e-04 1.39e-04 1.43e-04 2.60e-04 -1.06e-04 1.56e-05 -6.90e-06 4.77e-04 1.50e-05 -1.77e-04 -7.38e-05 -4.08e-04 -2.87e-04 1.23e-04 3.28e-04 1.50e-05 9.39e-04
    7.38e-04 2.79e-05 -1.53e-04 2.94e-04 -3.53e-04 -2.13e-04 -3.01e-05 1.77e-04 2.79e-05 5.48e-04 1.31e-04 -2.72e-04 -1.08e-04 -1.82e-04 2.93e-05 2.51e-04 -1.53e-04 1.31e-04 1.51e-03 -5.20e-04 9.46e-05 -2.75e-04 2.07e-04 -3.07e-04 2.94e-04 -2.72e-04 -5.20e-04 1.67e-03 -3.03e-04 -1.86e-04 -2.08e-04 -1.96e-04 -3.53e-04 -1.08e-04 9.46e-05 -3.03e-04 4.85e-04 2.98e-04 1.42e-04 -1.78e-05 -2.13e-04 -1.82e-04 -2.75e-04 -1.86e-04 2.98e-04 8.54e-04 -4.34e-05 8.11e-06 -3.01e-05 2.93e-05 2.07e-04 -2.08e-04 1.42e-04 -4.34e-05 5.03e-04 1.51e-04 1.77e-04 2.51e-04 -3.07e-04 -1.96e-04 -1.78e-05 8.11e-06 1.51e-04 1.04e-03
    6.97e-04 -1.62e-04 -2.19e-05 7.61e-05 -2.66e-04 -2.97e-04 -6.93e-05 2.84e-04 -1.62e-04 1.01e-03 4.98e-05 -4.54e-04 8.82e-05 1.28e-04 5.93e-05 1.18e-05 -2.19e-05 4.98e-05 1.67e-03 -2.54e-04 3.35e-06 -3.26e-04 1.29e-04 -6.96e-05 7.61e-05 -4.54e-04 -2.54e-04 8.66e-04 -8.84e-05 6.27e-05 -5.84e-06 -1.10e-04 -2.66e-04 8.82e-05 3.35e-06 -8.84e-05 3.51e-04 2.40e-04 1.76e-04 -7.88e-05 -2.97e-04 1.28e-04 -3.26e-04 6.27e-05 2.40e-04 9.37e-04 1.84e-04 -3.87e-04 -6.93e-05 5.93e-05 1.29e-04 -5.84e-06 1.76e-04 1.84e-04 6.04e-04 8.39e-05 2.84e-04 1.18e-05 -6.96e-05 -1.10e-04 -7.88e-05 -3.87e-04 8.39e-05 1.07e-03
    5.29e-04 -2.36e-04 -9.58e-05 -6.28e-06 -1.53e-04 -2.10e-04 -3.37e-05 3.15e-04 -2.36e-04 1.39e-03 -3.09e-07 -1.47e-04 8.35e-05 2.39e-04 1.27e-04 -1.92e-04 -9.58e-05 -3.09e-07 1.45e-03 -2.16e-04 2.90e-05 -2.13e-04 6.23e-05 1.18e-04 -6.28e-06 -1.47e-04 -2.16e-04 4.90e-04 -3.13e-05 2.47e-04 1.18e-04 -1.62e-04 -1.53e-04 8.35e-05 2.90e-05 -3.13e-05 1.53e-04 9.93e-05 1.19e-04 -4.95e-05 -2.10e-04 2.39e-04 -2.13e-04 2.47e-04 9.93e-05 7.99e-04 2.81e-04 -5.08e-04 -3.37e-05 1.27e-04 6.23e-05 1.18e-04 1.19e-04 2.81e-04 6.03e-04 3.00e-05 3.15e-04 -1.92e-04 1.18e-04 -1.62e-04 -4.95e-05 -5.08e-04 3.00e-05 8.69e-04
    3.35e-04 3.38e-05 1.77e-05 -4.63e-05 -3.46e-05 -1.07e-04 6.56e-05 1.55e-06 3.38e-05 6.36e-04 -5.17e-05 -1.42e-04 -1.66e-04 -1.41e-04 2.50e-04 -2.28e-04 1.77e-05 -5.17e-05 2.10e-03 -2.27e-04 4.35e-04 -1.78e-04 -4.05e-04 2.48e-04 -4.63e-05 -1.42e-04 -2.27e-04 6.99e-04 1.90e-04 1.05e-05 3.32e-04 6.28e-05 -3.46e-05 -1.66e-04 4.35e-04 1.90e-04 6.41e-04 1.65e-05 -8.90e-06 1.02e-04 -1.07e-04 -1.41e-04 -1.78e-04 1.05e-05 1.65e-05 5.10e-04 -2.23e-04 2.18e-04 6.56e-05 2.50e-04 -4.05e-04 3.32e-04 -8.90e-06 -2.23e-04 1.47e-03 -3.61e-04 1.55e-06 -2.28e-04 2.48e-04 6.28e-05 1.02e-04 2.18e-04 -3.61e-04 8.43e-04
    1.49e-03 -8.01e-04 1.53e-04 3.75e-05 -2.10e-04 -9.95e-06 -4.14e-04 1.76e-04 -8.01e-04 1.30e-03 7.95e-05 -4.08e-05 9.06e-05 -6.35e-05 1.71e-04 -1.06e-04 1.53e-04 7.95e-05 6.96e-04 1.22e-04 -2.05e-04 -5.38e-04 -6.87e-05 -2.15e-04 3.75e-05 -4.08e-05 1.22e-04 6.67e-04 -1.42e-04 -4.35e-04 -2.00e-05 -1.27e-04 -2.10e-04 9.06e-05 -2.05e-04 -1.42e-04 6.57e-04 3.35e-04 1.19e-04 6.75e-05 -9.95e-06 -6.35e-05 -5.38e-04 -4.35e-04 3.35e-04 1.36e-03 -1.66e-04 2.15e-04 -4.14e-04 1.71e-04 -6.87e-05 -2.00e-05 1.19e-04 -1.66e-04 6.86e-04 -8.72e-05 1.76e-04 -1.06e-04 -2.15e-04 -1.27e-04 6.75e-05 2.15e-04 -8.72e-05 5.15e-04
    6.38e-04 1.78e-05 -3.57e-04 2.75e-04 -7.75e-05 8.67e-05 4.70e-06 1.81e-04 1.78e-05 6.04e-04 4.56e-05 3.47e-04 -2.13e-04 -3.01e-04 1.57e-04 -1.99e-04 -3.57e-04 4.56e-05 8.28e-04 -3.37e-04 7.63e-05 -1.68e-04 9.91e-05 -3.30e-04 2.75e-04 3.47e-04 -3.37e-04 1.76e-03 2.01e-06 7.12e-06 -1.95e-04 -3.73e-04 -7.75e-05 -2.13e-04 7.63e-05 2.01e-06 5.43e-04 3.82e-05 -1.11e-04 -5.17e-05 8.67e-05 -3.01e-04 -1.68e-04 7.12e-06 3.82e-05 7.22e-04 -1.57e-04 2.45e-04 4.70e-06 1.57e-04 9.91e-05 -1.95e-04 -1.11e-04 -1.57e-04 2.79e-04 4.48e-05 1.81e-04 -1.99e-04 -3.30e-04 -3.73e-04 -5.17e-05 2.45e-04 4.48e-05 9.64e-04
    1.72e-03 -3.24e-04 8.58e-06 5.12e-04 -4.88e-04 -6.00e-05 -2.30e-04 -1.76e-04 -3.24e-04 1.15e-03 1.49e-04 -4.25e-04 1.10e-05 -9.53e-05 2.45e-04 -2.13e-04 8.58e-06 1.49e-04 6.43e-04 -8.11e-05 -1.32e-04 1.82e-04 1.63e-04 -2.52e-04 5.12e-04 -4.25e-04 -8.11e-05 9.05e-04 -1.59e-04 1.57e-04 -9.76e-05 -1.31e-04 -4.88e-04 1.10e-05 -1.32e-04 -1.59e-04 4.73e-04 -1.10e-04 7.78e-05 1.26e-04 -6.00e-05 -9.53e-05 1.82e-04 1.57e-04 -1.10e-04 1.02e-03 1.61e-04 -2.08e-04 -2.30e-04 2.45e-04 1.63e-04 -9.76e-05 7.78e-05 1.61e-04 5.08e-04 -2.78e-04 -1.76e-04 -2.13e-04 -2.52e-04 -1.31e-04 1.26e-04 -2.08e-04 -2.78e-04 8.42e-04
    7.13e-04 -4.25e-04 -1.68e-04 -3.09e-04 -6.08e-05 3.00e-04 -9.08e-05 -5.27e-05 -4.25e-04 1.41e-03 2.45e-04 6.38e-05 3.65e-04 -3.69e-04 -4.42e-04 1.82e-04 -1.68e-04 2.45e-04 7.89e-04 -1.64e-04 7.52e-05 2.11e-04 -1.86e-04 -9.07e-05 -3.09e-04 6.38e-05 -1.64e-04 5.96e-04 -5.82e-05 -2.69e-04 1.82e-04 1.34e-04 -6.08e-05 3.65e-04 7.52e-05 -5.82e-05 3.42e-04 -1.54e-04 -8.92e-05 -5.74e-05 3.00e-04 -3.69e-04 2.11e-04 -2.69e-04 -1.54e-04 1.05e-03 -4.44e-04 -9.85e-05 -9.08e-05 -4.42e-04 -1.86e-04 1.82e-04 -8.92e-05 -4.44e-04 1.64e-03 -2.64e-04 -5.27e-05 1.82e-04 -9.07e-05 1.34e-04 -5.74e-05 -9.85e-05 -2.64e-04 4.67e-04
    1.69e-03 -5.69e-04 3.54e-04 2.41e-04 -3.87e-04 -3.70e-04 -8.14e-05 8.05e-06 -5.69e-04 1.55e-03 1.73e-04 2.24e-04 -7.33e-05 -2.94e-04 2.89e-04 -2.25e-04 3.54e-04 1.73e-04 5.68e-04 3.08e-04 -2.98e-04 -2.79e-04 7.70e-05 -1.07e-04 2.41e-04 2.24e-04 3.08e-04 4.55e-04 -3.15e-04 -1.35e-04 1.61e-04 -8.43e-05 -3.87e-04 -7.33e-05 -2.98e-04 -3.15e-04 6.34e-04 1.26e-04 6.34e-05 1.17e-04 -3.70e-04 -2.94e-04 -2.79e-04 -1.35e-04 1.26e-04 7.69e-04 7.66e-06 2.25e-04 -8.14e-05 2.89e-04 7.70e-05 1.61e-04 6.34e-05 7.66e-06 5.98e-04 -1.15e-05 8.05e-06 -2.25e-04 -1.07e-04 -8.43e-05 1.17e-04 2.25e-04 -1.15e-05 2.25e-04
    1.58e-04 1.99e-05 -3.04e-05 -9.01e-05 -1.40e-04 -5.50e-05 8.87e-05 -2.03e-04 1.99e-05 1.10e-04 2.35e-06 4.45e-05 -1.18e-04 5.61e-05 -3.18e-05 -7.57e-05 -3.04e-05 2.35e-06 2.25e-04 1.77e-04 -7.95e-06 -1.66e-04 -2.16e-04 -2.42e-06 -9.01e-05 4.45e-05 1.77e-04 5.78e-04 3.25e-04 -1.21e-04 -2.61e-04 -6.61e-05 -1.40e-04 -1.18e-04 -7.95e-06 3.25e-04 1.58e-03 4.49e-05 -1.67e-04 -3.27e-04 -5.50e-05 5.61e-05 -1.66e-04 -1.21e-04 4.49e-05 3.54e-04 1.17e-04 -1.13e-05 8.87e-05 -3.18e-05 -2.16e-04 -2.61e-04 -1.67e-04 1.17e-04 3.39e-04 3.86e-05 -2.03e-04 -7.57e-05 -2.42e-06 -6.61e-05 -3.27e-04 -1.13e-05 3.86e-05 1.55e-03
    5.24e-04 -8.63e-05 -1.27e-04 -7.73e-05 1.62e-04 6.10e-06 5.75e-05 1.68e-04 -8.63e-05 4.87e-04 1.53e-04 1.12e-04 2.90e-05 -2.22e-04 1.56e-04 -3.84e-05 -1.27e-04 1.53e-04 4.71e-04 2.72e-05 1.18e-04 -4.41e-05 4.87e-05 1.71e-05 -7.73e-05 1.12e-04 2.72e-05 1.26e-03 1.46e-05 1.32e-04 -1.82e-04 7.11e-04 1.62e-04 2.90e-05 1.18e-04 1.46e-05 1.21e-03 4.62e-05 8.30e-05 1.88e-04 6.10e-06 -2.22e-04 -4.41e-05 1.32e-04 4.62e-05 3.87e-04 -1.65e-04 3.26e-05 5.75e-05 1.56e-04 4.87e-05 -1.82e-04 8.30e-05 -1.65e-04 2.16e-04 -7.73e-05 1.68e-04 -3.84e-05 1.71e-05 7.11e-04 1.88e-04 3.26e-05 -7.73e-05 1.25e-03
    5.78e-04 -3.29e-04 8.51e-05 -2.52e-04 -1.06e-04 -3.46e-05 3.30e-05 -1.82e-04 -3.29e-04 8.88e-04 3.16e-04 -1.02e-06 1.28e-04 -9.71e-05 -2.50e-06 3.66e-04 8.51e-05 3.16e-04 1.80e-03 -4.57e-04 2.62e-04 1.91e-04 6.44e-04 1.06e-04 -2.52e-04 -1.02e-06 -4.57e-04 5.91e-04 -9.13e-05 2.96e-05 -2.55e-04 2.38e-05 -1.06e-04 1.28e-04 2.62e-04 -9.13e-05 2.81e-04 7.58e-05 1.79e-04 -5.84e-06 -3.46e-05 -9.71e-05 1.91e-04 2.96e-05 7.58e-05 4.91e-04 -1.10e-05 -1.17e-04 3.30e-05 -2.50e-06 6.44e-04 -2.55e-04 1.79e-04 -1.10e-05 9.02e-04 -7.91e-05 -1.82e-04 3.66e-04 1.06e-04 2.38e-05 -5.84e-06 -1.17e-04 -7.91e-05 5.48e-04
    1.34e-03 -2.62e-04 -8.93e-05 -4.32e-04 1.72e-04 2.16e-04 -7.08e-06 2.70e-04 -2.62e-04 7.76e-04 -5.93e-05 4.00e-04 9.60e-06 -3.57e-04 2.18e-06 -1.46e-04 -8.93e-05 -5.93e-05 5.66e-04 1.92e-05 1.87e-04 6.11e-05 1.97e-05 -2.75e-04 -4.32e-04 4.00e-04 1.92e-05 1.26e-03 2.25e-05 -1.10e-04 7.88e-05 -1.65e-04 1.72e-04 9.60e-06 1.87e-04 2.25e-05 8.68e-04 -1.06e-04 1.44e-04 -2.42e-04 2.16e-04 -3.57e-04 6.11e-05 -1.10e-04 -1.06e-04 7.14e-04 -1.50e-04 -1.95e-06 -7.08e-06 2.18e-06 1.97e-05 7.88e-05 1.44e-04 -1.50e-04 2.90e-04 1.95e-05 2.70e-04 -1.46e-04 -2.75e-04 -1.65e-04 -2.42e-04 -1.95e-06 1.95e-05 7.59e-04
    1.68e-03 -5.25e-04 9.12e-05 4.67e-04 -5.48e-04 -1.36e-04 -2.49e-04 -2.15e-04 -5.25e-04 1.52e-03 2.49e-04 -3.43e-05 -2.07e-05 -1.37e-04 3.98e-04 -3.72e-04 9.12e-05 2.49e-04 5.78e-04 6.97e-05 -1.95e-04 5.96e-05 1.17e-04 -2.16e-04 4.67e-04 -3.43e-05 6.97e-05 6.81e-04 -2.36e-04 5.07e-05 1.19e-04 -2.40e-04 -5.48e-04 -2.07e-05 -1.95e-04 -2.36e-04 4.98e-04 -1.97e-05 5.12e-05 2.20e-04 -1.36e-04 -1.37e-04 5.96e-05 5.07e-05 -1.97e-05 9.97e-04 1.73e-04 7.93e-05 -2.49e-04 3.98e-04 1.17e-04 1.19e-04 5.12e-05 1.73e-04 5.89e-04 -1.94e-04 -2.15e-04 -3.72e-04 -2.16e-04 -2.40e-04 2.20e-04 7.93e-05 -1.94e-04 6.66e-04
    6.83e-04 -4.55e-04 -2.68e-05 8.39e-05 -2.03e-04 -1.53e-04 1.32e-04 -6.65e-05 -4.55e-04 1.26e-03 -1.03e-04 7.17e-05 3.10e-04 4.12e-04 -2.60e-04 1.07e-04 -2.68e-05 -1.03e-04 4.61e-04 9.70e-05 -1.61e-04 -1.78e-04 -2.08e-04 -6.19e-05 8.39e-05 7.17e-05 9.70e-05 8.29e-04 -3.19e-04 -4.88e-04 -2.57e-04 -2.29e-04 -2.03e-04 3.10e-04 -1.61e-04 -3.19e-04 4.66e-04 4.11e-04 1.97e-04 1.70e-04 -1.53e-04 4.12e-04 -1.78e-04 -4.88e-04 4.11e-04 1.26e-03 -3.17e-05 8.09e-05 1.32e-04 -2.60e-04 -2.08e-04 -2.57e-04 1.97e-04 -3.17e-05 9.29e-04 5.72e-05 -6.65e-05 1.07e-04 -6.19e-05 -2.29e-04 1.70e-04 8.09e-05 5.72e-05 5.46e-04
    4.88e-04 8.67e-05 -2.90e-05 2.62e-05 -1.60e-04 -5.16e-05 -1.94e-04 -3.60e-04 8.67e-05 4.87e-04 -2.35e-04 9.15e-05 -1.43e-04 -3.43e-05 -1.43e-04 -1.64e-04 -2.90e-05 -2.35e-04 1.64e-03 -3.43e-04 3.70e-04 2.00e-05 -2.03e-04 6.27e-04 2.62e-05 9.15e-05 -3.43e-04 9.14e-04 -4.53e-05 -9.88e-05 1.35e-04 -8.64e-05 -1.60e-04 -1.43e-04 3.70e-04 -4.53e-05 4.67e-04 7.85e-05 8.11e-05 2.91e-04 -5.16e-05 -3.43e-05 2.00e-05 -9.88e-05 7.85e-05 4.21e-04 -7.14e-05 -6.04e-05 -1.94e-04 -1.43e-04 -2.03e-04 1.35e-04 8.11e-05 -7.14e-05 8.40e-04 4.93e-05 -3.60e-04 -1.64e-04 6.27e-04 -8.64e-05 2.91e-04 -6.04e-05 4.93e-05 1.18e-03
    8.48e-04 -3.55e-05 -3.65e-04 3.28e-04 -3.11e-04 -1.30e-04 -1.15e-04 2.86e-04 -3.55e-05 7.45e-04 1.14e-04 -4.83e-04 -2.53e-05 -1.05e-04 7.81e-05 2.47e-04 -3.65e-04 1.14e-04 1.39e-03 -4.42e-04 8.90e-05 -2.50e-04 1.96e-04 -1.96e-04 3.28e-04 -4.83e-04 -4.42e-04 1.56e-03 -2.86e-04 -1.65e-04 -2.35e-04 -1.11e-04 -3.11e-04 -2.53e-05 8.90e-05 -2.86e-04 3.90e-04 1.97e-04 1.46e-04 -9.64e-05 -1.30e-04 -1.05e-04 -2.50e-04 -1.65e-04 1.97e-04 9.85e-04 -1.75e-05 -3.77e-04 -1.15e-04 7.81e-05 1.96e-04 -2.35e-04 1.46e-04 -1.75e-05 4.65e-04 4.95e-05 2.86e-04 2.47e-04 -1.96e-04 -1.11e-04 -9.64e-05 -3.77e-04 4.95e-05 1.03e-03
    6.59e-04 1.09e-04 6.58e-05 7.75e-05 -2.75e-04 -1.04e-04 4.74e-05 -3.93e-04 1.09e-04 6.27e-04 -1.13e-04 2.70e-04 -1.88e-04 -1.23e-04 -1.93e-05 -3.21e-04 6.58e-05 -1.13e-04 1.56e-03 -8.21e-04 2.38e-04 -1.52e-04 8.61e-05 -2.28e-04 7.75e-05 2.70e-04 -8.21e-04 1.83e-03 -5.43e-05 -2.63e-05 -2.40e-04 -2.10e-04 -2.75e-04 -1.88e-04 2.38e-04 -5.43e-05 6.49e-04 5.15e-05 -6.65e-06 2.44e-04 -1.04e-04 -1.23e-04 -1.52e-04 -2.63e-05 5.15e-05 5.44e-04 -4.04e-05 3.01e-04 4.74e-05 -1.93e-05 8.61e-05 -2.40e-04 -6.65e-06 -4.04e-05 6.49e-04 -2.71e-05 -3.93e-04 -3.21e-04 -2.28e-04 -2.10e-04 2.44e-04 3.01e-04 -2.71e-05 1.22e-03
    6.68e-04 1.32e-04 4.77e-04 2.39e-04 -3.41e-04 -1.02e-04 4.59e-05 -4.20e-04 1.32e-04 5.00e-04 3.39e-06 1.97e-04 -1.65e-04 2.10e-05 2.10e-05 -8.77e-05 4.77e-04 3.39e-06 1.66e-03 -4.90e-04 -6.98e-05 -1.72e-05 2.00e-04 -4.26e-04 2.39e-04 1.97e-04 -4.90e-04 1.43e-03 -3.06e-04 -1.99e-04 -7.04e-05 -3.52e-04 -3.41e-04 -1.65e-04 -6.98e-05 -3.06e-04 4.20e-04 7.76e-05 -3.51e-05 1.71e-04 -1.02e-04 2.10e-05 -1.72e-05 -1.99e-04 7.76e-05 3.83e-04 -2.39e-05 2.28e-04 4.59e-05 2.10e-05 2.00e-04 -7.04e-05 -3.51e-05 -2.39e-05 5.44e-04 7.75e-05 -4.20e-04 -8.77e-05 -4.26e-04 -3.52e-04 1.71e-04 2.28e-04 7.75e-05 1.07e-03
    1.18e-03 -6.56e-04 -1.01e-04 2.45e-04 -3.82e-04 -1.69e-04 -1.57e-04 -1.30e-04 -6.56e-04 1.57e-03 2.02e-04 2.91e-04 2.14e-05 -3.04e-04 1.76e-04 -1.21e-04 -1.01e-04 2.02e-04 5.31e-04 2.59e-05 -1.34e-04 -3.24e-06 -2.50e-05 -1.49e-04 2.45e-04 2.91e-04 2.59e-05 9.06e-04 -3.64e-04 -3.74e-04 4.53e-05 -3.68e-04 -3.82e-04 2.14e-05 -1.34e-04 -3.64e-04 5.11e-04 2.17e-04 1.26e-04 2.98e-04 -1.69e-04 -3.04e-04 -3.24e-06 -3.74e-04 2.17e-04 9.95e-04 1.33e-04 2.73e-04 -1.57e-04 1.76e-04 -2.50e-05 4.53e-05 1.26e-04 1.33e-04 8.98e-04 5.34e-05 -1.30e-04 -1.21e-04 -1.49e-04 -3.68e-04 2.98e-04 2.73e-04 5.34e-05 7.24e-04
    1.54e-03 -1.93e-04 3.23e-04 -6.18e-04 -4.14e-04 1.11e-04 3.89e-05 1.86e-04 -1.93e-04 5.15e-04 -1.35e-04 1.62e-04 5.01e-05 -1.98e-04 6.24e-05 8.20e-05 3.23e-04 -1.35e-04 6.28e-04 -2.32e-04 4.03e-05 -2.66e-05 1.30e-04 -2.30e-04 -6.18e-04 1.62e-04 -2.32e-04 1.14e-03 -6.59e-05 1.79e-04 -1.35e-04 -2.28e-05 -4.14e-04 5.01e-05 4.03e-05 -6.59e-05 1.19e-03 -2.55e-04 1.36e-04 -9.05e-05 1.11e-04 -1.98e-04 -2.66e-05 1.79e-04 -2.55e-04 6.69e-04 -1.41e-04 -9.70e-05 3.89e-05 6.24e-05 1.30e-04 -1.35e-04 1.36e-04 -1.41e-04 2.43e-04 -2.47e-05 1.86e-04 8.20e-05 -2.30e-04 -2.28e-05 -9.05e-05 -9.70e-05 -2.47e-05 6.83e-04
    1.21e-03 -7.82e-04 8.54e-05 6.66e-05 -3.37e-04 -5.85e-05 -1.59e-04 -5.37e-05 -7.82e-04 1.53e-03 -4.40e-05 1.93e-04 2.71e-04 1.57e-04 -1.32e-05 4.86e-05 8.54e-05 -4.40e-05 5.00e-04 2.01e-04 -2.64e-04 -3.15e-04 -1.62e-04 -1.02e-04 6.66e-05 1.93e-04 2.01e-04 8.43e-04 -3.13e-04 -6.47e-04 -1.23e-04 -2.38e-04 -3.37e-04 2.71e-04 -2.64e-04 -3.13e-04 5.84e-04 4.14e-04 2.24e-04 1.73e-04 -5.85e-05 1.57e-04 -3.15e-04 -6.47e-04 4.14e-04 1.33e-03 -8.36e-05 2.06e-04 -1.59e-04 -1.32e-05 -1.62e-04 -1.23e-04 2.24e-04 -8.36e-05 7.94e-04 1.11e-04 -5.37e-05 4.86e-05 -1.02e-04 -2.38e-04 1.73e-04 2.06e-04 1.11e-04 4.90e-04
    4.15e-04 -1.89e-05 -5.96e-05 -5.12e-05 -1.35e-04 -1.71e-04 8.96e-05 -1.03e-04 -1.89e-05 5.48e-04 -1.82e-04 -1.36e-04 -1.94e-04 -1.33e-04 5.25e-05 -7.90e-05 -5.96e-05 -1.82e-04 1.89e-03 -1.58e-04 3.87e-04 6.64e-05 -1.91e-04 6.28e-04 -5.12e-05 -1.36e-04 -1.58e-04 7.06e-04 1.22e-04 6.91e-05 2.07e-04 -5.59e-05 -1.35e-04 -1.94e-04 3.87e-04 1.22e-04 5.60e-04 1.86e-04 7.73e-05 1.89e-04 -1.71e-04 -1.33e-04 6.64e-05 6.91e-05 1.86e-04 5.40e-04 -2.50e-04 1.24e-04 8.96e-05 5.25e-05 -1.91e-04 2.07e-04 7.73e-05 -2.50e-04 1.18e-03 -3.15e-04 -1.03e-04 -7.90e-05 6.28e-04 -5.59e-05 1.89e-04 1.24e-04 -3.15e-04 8.11e-04
    1.41e-03 1.29e-04 1.98e-04 -2.31e-04 -1.04e-04 7.92e-04 -3.38e-05 6.78e-05 1.29e-04 3.14e-04 1.44e-04 -1.40e-04 -1.24e-04 -2.13e-05 4.93e-05 1.81e-05 1.98e-04 1.44e-04 5.46e-04 -1.02e-04 -1.10e-04 -2.05e-04 1.52e-04 -1.93e-04 -2.31e-04 -1.40e-04 -1.02e-04 4.39e-04 2.39e-04 -3.75e-05 -2.89e-05 1.61e-05 -1.04e-04 -1.24e-04 -1.10e-04 2.39e-04 9.60e-04 1.05e-04 -7.04e-05 1.46e-04 7.92e-04 -2.13e-05 -2.05e-04 -3.75e-05 1.05e-04 1.19e-03 -2.73e-05 2.70e-05 -3.38e-05 4.93e-05 1.52e-04 -2.89e-05 -7.04e-05 -2.73e-05 1.85e-04 -1.27e-04 6.78e-05 1.81e-05 -1.93e-04 1.61e-05 1.46e-04 2.70e-05 -1.27e-04 3.20e-04
    1.50e-03 -6.77e-04 2.02e-05 4.09e-04 -4.57e-04 -2.18e-04 -3.24e-04 -5.78e-05 -6.77e-04 1.41e-03 1.90e-04 -3.84e-05 3.50e-05 -2.74e-04 3.01e-04 -2.47e-04 2.02e-05 1.90e-04 6.46e-04 1.00e-05 -1.55e-04 -5.64e-05 6.35e-05 -2.10e-04 4.09e-04 -3.84e-05 1.00e-05 7.17e-04 -2.30e-04 -1.16e-04 1.24e-04 -3.15e-04 -4.57e-04 3.50e-05 -1.55e-04 -2.30e-04 5.24e-04 8.54e-05 1.03e-04 2.14e-04 -2.18e-04 -2.74e-04 -5.64e-05 -1.16e-04 8.54e-05 9.89e-04 1.17e-04 1.91e-04 -3.24e-04 3.01e-04 6.35e-05 1.24e-04 1.03e-04 1.17e-04 6.53e-04 -1.60e-04 -5.78e-05 -2.47e-04 -2.10e-04 -3.15e-04 2.14e-04 1.91e-04 -1.60e-04 7.37e-04
    7.01e-04 -3.89e-05 -6.13e-05 2.15e-04 -3.25e-04 -2.86e-04 -5.47e-05 2.33e-04 -3.89e-05 7.11e-04 1.21e-04 -4.48e-04 2.45e-05 -2.59e-05 4.16e-05 2.07e-04 -6.13e-05 1.21e-04 1.61e-03 -3.80e-04 4.59e-05 -3.38e-04 1.65e-04 -1.88e-04 2.15e-04 -4.48e-04 -3.80e-04 1.37e-03 -2.47e-04 -1.32e-04 -1.30e-04 -8.68e-05 -3.25e-04 2.45e-05 4.59e-05 -2.47e-04 3.99e-04 3.12e-04 1.64e-04 -5.46e-05 -2.86e-04 -2.59e-05 -3.38e-04 -1.32e-04 3.12e-04 8.77e-04 6.12e-05 -2.00e-04 -5.47e-05 4.16e-05 1.65e-04 -1.30e-04 1.64e-04 6.12e-05 5.75e-04 1.38e-04 2.33e-04 2.07e-04 -1.88e-04 -8.68e-05 -5.46e-05 -2.00e-04 1.38e-04 1.07e-03
    1.63e-03 -4.91e-04 2.06e-04 4.28e-04 -4.33e-04 -2.53e-04 -3.19e-04 1.39e-06 -4.91e-04 1.03e-03 1.94e-04 -3.08e-04 -2.43e-05 -2.69e-04 2.28e-04 -1.76e-04 2.06e-04 1.94e-04 7.23e-04 6.92e-06 -2.10e-04 -2.12e-04 1.53e-04 -2.83e-04 4.28e-04 -3.08e-04 6.92e-06 7.54e-04 -9.81e-05 -6.80e-05 -1.72e-05 -1.56e-04 -4.33e-04 -2.43e-05 -2.10e-04 -9.81e-05 6.41e-04 5.13e-05 1.96e-05 1.33e-04 -2.53e-04 -2.69e-04 -2.12e-04 -6.80e-05 5.13e-05 1.04e-03 2.56e-05 2.35e-04 -3.19e-04 2.28e-04 1.53e-04 -1.72e-05 1.96e-05 2.56e-05 5.27e-04 -2.56e-04 1.39e-06 -1.76e-04 -2.83e-04 -1.56e-04 1.33e-04 2.35e-04 -2.56e-04 7.26e-04
    6.28e-04 9.36e-05 -1.17e-05 2.57e-04 -3.42e-04 -2.63e-04 -1.64e-05 1.31e-04 9.36e-05 4.14e-04 1.37e-04 -2.89e-04 -8.42e-05 -7.86e-05 2.69e-05 2.25e-04 -1.17e-05 1.37e-04 1.51e-03 -4.28e-04 7.09e-05 -3.63e-04 2.46e-04 -2.51e-04 2.57e-04 -2.89e-04 -4.28e-04 1.47e-03 -3.01e-04 -2.61e-04 -1.85e-04 -1.12e-04 -3.42e-04 -8.42e-05 7.09e-05 -3.01e-04 3.86e-04 3.30e-04 8.32e-05 1.37e-05 -2.63e-04 -7.86e-05 -3.63e-04 -2.61e-04 3.30e-04 6.50e-04 1.80e-06 7.45e-05 -1.64e-05 2.69e-05 2.46e-04 -1.85e-04 8.32e-05 1.80e-06 4.09e-04 1.37e-04 1.31e-04 2.25e-04 -2.51e-04 -1.12e-04 1.37e-05 7.45e-05 1.37e-04 1.03e-03
    1.18e-03 -1.03e-04 -4.79e-04 4.19e-04 -3.01e-04 3.85e-05 -2.04e-04 1.15e-04 -1.03e-04 8.04e-04 3.06e-05 -5.55e-04 9.41e-06 -1.12e-04 1.00e-04 2.06e-04 -4.79e-04 3.06e-05 9.60e-04 -3.42e-04 6.60e-05 1.54e-06 1.68e-04 -1.59e-04 4.19e-04 -5.55e-04 -3.42e-04 1.42e-03 -2.14e-04 -6.25e-05 -2.90e-04 -9.99e-05 -3.01e-04 9.41e-06 6.60e-05 -2.14e-04 3.68e-04 8.43e-06 1.58e-04 -8.70e-05 3.85e-05 -1.12e-04 1.54e-06 -6.25e-05 8.43e-06 1.02e-03 6.32e-05 -5.92e-04 -2.04e-04 1.00e-04 1.68e-04 -2.90e-04 1.58e-04 6.32e-05 4.20e-04 -9.58e-05 1.15e-04 2.06e-04 -1.59e-04 -9.99e-05 -8.70e-05 -5.92e-04 -9.58e-05 9.79e-04
    1.56e-03 -4.55e-04 -9.75e-05 4.79e-04 -4.53e-04 -1.53e-04 -2.76e-04 -1.05e-04 -4.55e-04 1.33e-03 1.99e-04 -3.04e-04 3.11e-05 -1.38e-04 2.92e-04 -2.28e-04 -9.75e-05 1.99e-04 7.19e-04 -9.74e-05 -8.21e-05 9.48e-05 1.29e-04 -2.60e-04 4.79e-04 -3.04e-04 -9.74e-05 8.34e-04 -1.26e-04 1.25e-04 9.05e-06 -2.33e-04 -4.53e-04 3.11e-05 -8.21e-05 -1.26e-04 5.23e-04 -8.28e-05 9.66e-05 1.44e-04 -1.53e-04 -1.38e-04 9.48e-05 1.25e-04 -8.28e-05 1.05e-03 1.69e-04 -9.41e-05 -2.76e-04 2.92e-04 1.29e-04 9.05e-06 9.66e-05 1.69e-04 6.05e-04 -2.55e-04 -1.05e-04 -2.28e-04 -2.60e-04 -2.33e-04 1.44e-04 -9.41e-05 -2.55e-04 8.93e-04
    1.54e-03 -5.61e-04 4.44e-04 2.38e-04 -3.41e-04 -2.77e-04 -1.88e-04 1.47e-04 -5.61e-04 1.02e-03 2.12e-04 -6.89e-05 -7.58e-05 -3.06e-04 2.78e-04 -2.45e-04 4.44e-04 2.12e-04 6.08e-04 1.09e-04 -2.20e-04 -5.48e-04 1.26e-04 -1.83e-04 2.38e-04 -6.89e-05 1.09e-04 5.84e-04 -1.00e-04 -1.55e-04 3.83e-05 2.42e-05 -3.41e-04 -7.58e-05 -2.20e-04 -1.00e-04 6.26e-04 1.62e-04 -9.44e-05 6.83e-05 -2.77e-04 -3.06e-04 -5.48e-04 -1.55e-04 1.62e-04 9.70e-04 -6.94e-05 2.22e-04 -1.88e-04 2.78e-04 1.26e-04 3.83e-05 -9.44e-05 -6.94e-05 3.00e-04 -1.27e-04 1.47e-04 -2.45e-04 -1.83e-04 2.42e-05 6.83e-05 2.22e-04 -1.27e-04 3.51e-04
    7.18e-04 -2.23e-04 5.62e-04 1.68e-04 -3.88e-04 -2.72e-04 -1.84e-04 -3.06e-04 -2.23e-04 5.76e-04 -1.77e-04 -1.44e-04 8.42e-05 2.20e-04 -1.71e-04 1.49e-04 5.62e-04 -1.77e-04 1.25e-03 1.51e-04 -2.40e-04 -3.06e-05 -3.10e-05 -2.69e-04 1.68e-04 -1.44e-04 1.51e-04 4.18e-04 -8.78e-05 5.53e-06 -1.73e-05 -2.98e-04 -3.88e-04 8.42e-05 -2.40e-04 -8.78e-05 3.24e-04 2.14e-04 1.52e-04 1.59e-04 -2.72e-04 2.20e-04 -3.06e-05 5.53e-06 2.14e-04 4.83e-04 2.91e-05 3.49e-05 -1.84e-04 -1.71e-04 -3.10e-05 -1.73e-05 1.52e-04 2.91e-05 7.71e-04 2.90e-04 -3.06e-04 1.49e-04 -2.69e-04 -2.98e-04 1.59e-04 3.49e-05 2.90e-04 8.25e-04
    1.69e-03 -1.09e-04 3.96e-04 -1.66e-04 -6.16e-04 7.63e-05 -3.09e-05 3.13e-04 -1.09e-04 5.09e-04 -1.35e-04 9.80e-06 7.93e-05 -1.89e-05 3.34e-05 7.94e-05 3.96e-04 -1.35e-04 7.04e-04 -1.61e-04 -1.58e-04 -1.31e-04 1.58e-04 -1.34e-04 -1.66e-04 9.80e-06 -1.61e-04 8.83e-04 1.16e-04 2.66e-04 -7.62e-05 -5.25e-05 -6.16e-04 7.93e-05 -1.58e-04 1.16e-04 9.64e-04 7.41e-06 6.09e-05 -3.70e-05 7.63e-05 -1.89e-05 -1.31e-04 2.66e-04 7.41e-06 7.96e-04 -4.61e-05 -7.60e-05 -3.09e-05 3.34e-05 1.58e-04 -7.62e-05 6.09e-05 -4.61e-05 2.87e-04 -8.75e-05 3.13e-04 7.94e-05 -1.34e-04 -5.25e-05 -3.70e-05 -7.60e-05 -8.75e-05 6.31e-04
    1.23e-03 -4.08e-04 -3.77e-04 2.88e-04 -3.53e-04 -2.27e-05 -1.87e-04 -6.32e-05 -4.08e-04 1.54e-03 2.08e-04 -1.42e-04 1.69e-06 -1.09e-04 2.62e-04 -3.60e-04 -3.77e-04 2.08e-04 7.71e-04 -1.81e-04 6.28e-05 1.88e-04 8.77e-05 -1.34e-04 2.88e-04 -1.42e-04 -1.81e-04 7.76e-04 -1.70e-04 1.65e-04 7.90e-05 -2.98e-04 -3.53e-04 1.69e-06 6.28e-05 -1.70e-04 3.38e-04 -2.24e-05 1.34e-04 1.87e-04 -2.27e-05 -1.09e-04 1.88e-04 1.65e-04 -2.24e-05 9.70e-04 2.67e-04 -3.41e-04 -1.87e-04 2.62e-04 8.77e-05 7.90e-05 1.34e-04 2.67e-04 7.49e-04 -1.44e-04 -6.32e-05 -3.60e-04 -1.34e-04 -2.98e-04 1.87e-04 -3.41e-04 -1.44e-04 9.05e-04
    4.55e-04 -2.80e-04 5.18e-05 -2.69e-05 -5.69e-05 -7.84e-05 -2.76e-04 2.10e-04 -2.80e-04 1.13e-03 -2.39e-04 1.61e-04 1.11e-04 2.70e-05 4.63e-04 -6.13e-05 5.18e-05 -2.39e-04 7.42e-04 -2.81e-04 5.04e-05 2.84e-04 -2.24e-04 1.15e-04 -2.69e-05 1.61e-04 -2.81e-04 8.19e-04 -3.25e-04 -3.20e-04 -2.78e-04 -2.88e-04 -5.69e-05 1.11e-04 5.04e-05 -3.25e-04 3.36e-04 2.05e-04 3.62e-04 1.38e-04 -7.84e-05 2.70e-05 2.84e-04 -3.20e-04 2.05e-04 7.36e-04 4.02e-05 1.04e-05 -2.76e-04 4.63e-04 -2.24e-04 -2.78e-04 3.62e-04 4.02e-05 1.60e-03 -9.18e-05 2.10e-04 -6.13e-05 1.15e-04 -2.88e-04 1.38e-04 1.04e-05 -9.18e-05 5.71e-04
    6.96e-04 1.20e-04 1.79e-04 1.53e-04 -3.66e-04 -1.65e-04 3.95e-05 -3.73e-04 1.20e-04 5.09e-04 -8.73e-05 1.40e-04 -2.15e-04 -5.86e-05 2.74e-05 -1.13e-04 1.79e-04 -8.73e-05 1.56e-03 -6.10e-04 2.01e-04 -1.38e-04 1.43e-04 -3.77e-04 1.53e-04 1.40e-04 -6.10e-04 1.61e-03 -2.37e-04 -2.43e-04 -1.44e-04 -2.95e-04 -3.66e-04 -2.15e-04 2.01e-04 -2.37e-04 5.63e-04 1.49e-04 2.48e-05 2.13e-04 -1.65e-04 -5.86e-05 -1.38e-04 -2.43e-04 1.49e-04 5.29e-04 2.04e-05 3.11e-04 3.95e-05 2.74e-05 1.43e-04 -1.44e-04 2.48e-05 2.04e-05 5.74e-04 1.08e-04 -3.73e-04 -1.13e-04 -3.77e-04 -2.95e-04 2.13e-04 3.11e-04 1.08e-04 1.06e-03
    1.78e-03 -3.23e-05 1.21e-04 2.80e-04 -2.40e-04 2.03e-04 4.13e-05 -7.90e-05 -3.23e-05 5.27e-04 -1.04e-04 -3.06e-04 -1.84e-05 -1.71e-04 9.10e-05 2.30e-04 1.21e-04 -1.04e-04 4.91e-04 -1.11e-04 -5.94e-05 2.57e-04 1.39e-04 -2.49e-04 2.80e-04 -3.06e-04 -1.11e-04 1.50e-03 -4.16e-05 -1.53e-05 -2.42e-04 7.66e-05 -2.40e-04 -1.84e-05 -5.94e-05 -4.16e-05 4.28e-04 -1.79e-04 7.34e-05 -1.30e-05 2.03e-04 -1.71e-04 2.57e-04 -1.53e-05 -1.79e-04 8.92e-04 1.08e-04 -3.89e-04 4.13e-05 9.10e-05 1.39e-04 -2.42e-04 7.34e-05 1.08e-04 2.72e-04 -6.34e-05 -7.90e-05 2.30e-04 -2.49e-04 7.66e-05 -1.30e-05 -3.89e-04 -6.34e-05 8.97e-04
    1.07e-03 -7.84e-04 6.54e-05 1.39e-04 -3.67e-04 -1.86e-04 -1.15e-04 -1.38e-04 -7.84e-04 1.37e-03 -4.09e-06 2.67e-04 1.91e-04 -1.38e-04 8.82e-06 8.43e-05 6.54e-05 -4.09e-06 4.69e-04 9.28e-05 -1.96e-04 -1.91e-04 -1.23e-04 -5.32e-05 1.39e-04 2.67e-04 9.28e-05 8.38e-04 -3.90e-04 -6.33e-04 -1.36e-04 -2.71e-04 -3.67e-04 1.91e-04 -1.96e-04 -3.90e-04 5.08e-04 3.91e-04 1.93e-04 2.28e-04 -1.86e-04 -1.38e-04 -1.91e-04 -6.33e-04 3.91e-04 9.99e-04 -1.27e-05 2.59e-04 -1.15e-04 8.82e-06 -1.23e-04 -1.36e-04 1.93e-04 -1.27e-05 8.46e-04 1.32e-04 -1.38e-04 8.43e-05 -5.32e-05 -2.71e-04 2.28e-04 2.59e-04 1.32e-04 4.93e-04
    1.45e-03 -4.55e-05 4.94e-04 2.21e-04 -3.16e-04 -2.43e-04 2.89e-05 4.62e-05 -4.55e-05 5.99e-04 1.58e-04 -3.22e-04 -5.65e-05 -2.44e-04 1.68e-04 -3.16e-05 4.94e-04 1.58e-04 6.49e-04 -1.12e-04 -1.22e-04 -3.37e-04 2.05e-04 -2.89e-04 2.21e-04 -3.22e-04 -1.12e-04 9.50e-04 1.28e-04 1.54e-04 -1.14e-04 1.40e-04 -3.16e-04 -5.65e-05 -1.22e-04 1.28e-04 6.66e-04 -1.17e-04 -9.40e-05 7.00e-05 -2.43e-04 -2.44e-04 -3.37e-04 1.54e-04 -1.17e-04 8.64e-04 -2.19e-05 8.89e-05 2.89e-05 1.68e-04 2.05e-04 -1.14e-04 -9.40e-05 -2.19e-05 2.48e-04 -1.69e-04 4.62e-05 -3.16e-05 -2.89e-04 1.40e-04 7.00e-05 8.89e-05 -1.69e-04 5.22e-04
    8.05e-04 -2.83e-04 -1.41e-04 -1.74e-04 2.94e-04 -1.35e-04 2.56e-04 -2.83e-05 -2.83e-04 1.85e-03 -1.50e-05 3.98e-06 4.31e-04 1.49e-06 -3.94e-04 -1.98e-04 -1.41e-04 -1.50e-05 4.30e-04 1.30e-04 -1.51e-04 -1.81e-04 1.93e-04 -1.14e-04 -1.74e-04 3.98e-06 1.30e-04 4.06e-04 -8.70e-05 6.47e-05 1.15e-04 1.88e-05 2.94e-04 4.31e-04 -1.51e-04 -8.70e-05 8.18e-04 6.51e-06 -6.14e-05 -1.69e-06 -1.35e-04 1.49e-06 -1.81e-04 6.47e-05 6.51e-06 1.09e-03 -4.57e-04 1.89e-04 2.56e-04 -3.94e-04 1.93e-04 1.15e-04 -6.14e-05 -4.57e-04 1.28e-03 -9.16e-06 -2.83e-05 -1.98e-04 -1.14e-04 1.88e-05 -1.69e-06 1.89e-04 -9.16e-06 3.61e-04
    4.78e-04 8.29e-06 -4.28e-05 1.50e-04 -2.05e-04 -9.60e-05 -1.62e-04 -3.69e-04 8.29e-06 4.16e-04 -2.02e-04 -3.95e-05 -1.07e-04 -2.74e-06 -1.87e-04 -5.63e-05 -4.28e-05 -2.02e-04 1.48e-03 -7.30e-05 3.49e-04 -2.91e-06 8.20e-05 6.35e-04 1.50e-04 -3.95e-05 -7.30e-05 4.69e-04 -4.71e-05 -6.28e-06 -3.05e-05 -1.75e-04 -2.05e-04 -1.07e-04 3.49e-04 -4.71e-05 3.70e-04 1.22e-04 1.12e-04 3.04e-04 -9.60e-05 -2.74e-06 -2.91e-06 -6.28e-06 1.22e-04 4.03e-04 -2.56e-04 -4.12e-05 -1.62e-04 -1.87e-04 8.20e-05 -3.05e-05 1.12e-04 -2.56e-04 1.22e-03 1.43e-04 -3.69e-04 -5.63e-05 6.35e-04 -1.75e-04 3.04e-04 -4.12e-05 1.43e-04 1.04e-03
    4.84e-04 1.61e-05 -8.45e-05 -1.08e-04 -9.47e-05 1.09e-05 3.75e-05 -2.21e-04 1.61e-05 6.40e-04 4.59e-05 2.51e-04 2.24e-05 -1.25e-04 3.14e-05 -4.05e-04 -8.45e-05 4.59e-05 1.13e-03 -6.55e-04 1.43e-04 -3.14e-04 7.25e-05 -1.49e-04 -1.08e-04 2.51e-04 -6.55e-04 1.72e-03 3.10e-04 1.61e-04 -3.48e-04 1.21e-04 -9.47e-05 2.24e-05 1.43e-04 3.10e-04 6.02e-04 -6.88e-05 -3.32e-06 1.72e-04 1.09e-05 -1.25e-04 -3.14e-04 1.61e-04 -6.88e-05 5.03e-04 -1.02e-04 1.81e-04 3.75e-05 3.14e-05 7.25e-05 -3.48e-04 -3.32e-06 -1.02e-04 5.54e-04 -2.19e-04 -2.21e-04 -4.05e-04 -1.49e-04 1.21e-04 1.72e-04 1.81e-04 -2.19e-04 1.29e-03
    7.12e-04 -2.38e-04 1.20e-04 -3.32e-04 -9.12e-05 -5.51e-06 3.25e-04 -3.01e-04 -2.38e-04 8.13e-04 1.67e-04 -1.32e-04 1.79e-05 -2.59e-04 -2.35e-04 2.12e-04 1.20e-04 1.67e-04 1.77e-03 -4.39e-04 3.59e-04 2.66e-04 -4.29e-05 -3.19e-06 -3.32e-04 -1.32e-04 -4.39e-04 7.07e-04 -4.75e-05 8.80e-05 -1.80e-04 1.80e-04 -9.12e-05 1.79e-05 3.59e-04 -4.75e-05 4.40e-04 1.17e-04 2.18e-05 -8.08e-05 -5.51e-06 -2.59e-04 2.66e-04 8.80e-05 1.17e-04 6.04e-04 -7.85e-05 -2.45e-05 3.25e-04 -2.35e-04 -4.29e-05 -1.80e-04 2.18e-05 -7.85e-05 1.08e-03 -4.39e-04 -3.01e-04 2.12e-04 -3.19e-06 1.80e-04 -8.08e-05 -2.45e-05 -4.39e-04 6.50e-04
    2.79e-04 6.00e-05 2.70e-07 -1.77e-04 -1.39e-04 -7.21e-05 1.32e-04 -9.13e-05 6.00e-05 2.93e-04 -1.32e-05 1.15e-04 1.10e-05 7.29e-07 -7.29e-06 9.91e-06 2.70e-07 -1.32e-05 4.75e-04 -7.18e-05 -8.71e-05 -1.41e-04 -2.08e-04 -1.81e-04 -1.77e-04 1.15e-04 -7.18e-05 1.26e-03 5.34e-04 2.52e-05 -1.06e-04 2.26e-04 -1.39e-04 1.10e-05 -8.71e-05 5.34e-04 1.21e-03 1.15e-04 -1.92e-05 6.61e-04 -7.21e-05 7.29e-07 -1.41e-04 2.52e-05 1.15e-04 2.82e-04 -5.29e-05 2.55e-04 1.32e-04 -7.29e-06 -2.08e-04 -1.06e-04 -1.92e-05 -5.29e-05 5.10e-04 5.29e-05 -9.13e-05 9.91e-06 -1.81e-04 2.26e-04 6.61e-04 2.55e-04 5.29e-05 9.94e-04
    8.49e-04 -2.04e-05 2.53e-05 3.49e-04 1.30e-04 4.36e-05 9.47e-05 2.23e-04 -2.04e-05 3.33e-04 -1.56e-04 2.09e-04 5.95e-05 5.07e-05 9.96e-05 -7.35e-05 2.53e-05 -1.56e-04 2.49e-04 -1.12e-04 -4.53e-06 2.47e-06 7.74e-06 8.93e-05 3.49e-04 2.09e-04 -1.12e-04 1.36e-03 4.05e-04 3.19e-04 -1.76e-04 -3.54e-04 1.30e-04 5.95e-05 -4.53e-06 4.05e-04 8.49e-04 2.09e-04 6.13e-05 -2.68e-05 4.36e-05 5.07e-05 2.47e-06 3.19e-04 2.09e-04 3.99e-04 4.22e-05 -3.68e-05 9.47e-05 9.96e-05 7.74e-06 -1.76e-04 6.13e-05 4.22e-05 3.44e-04 2.79e-04 2.23e-04 -7.35e-05 8.93e-05 -3.54e-04 -2.68e-05 -3.68e-05 2.79e-04 8.23e-04
    7.16e-04 -4.69e-04 -1.58e-04 -3.11e-04 -5.38e-05 2.01e-04 -2.81e-04 -4.10e-06 -4.69e-04 1.45e-03 3.00e-04 6.11e-05 3.18e-04 -1.88e-04 -2.54e-04 2.71e-04 -1.58e-04 3.00e-04 8.13e-04 -1.92e-04 6.10e-05 2.66e-04 -6.35e-05 -9.08e-05 -3.11e-04 6.11e-05 -1.92e-04 6.28e-04 -1.15e-04 -2.50e-04 2.38e-04 6.93e-05 -5.38e-05 3.18e-04 6.10e-05 -1.15e-04 2.84e-04 -8.13e-05 -4.95e-05 -1.90e-05 2.01e-04 -1.88e-04 2.66e-04 -2.50e-04 -8.13e-05 9.67e-04 -5.01e-04 -8.49e-05 -2.81e-04 -2.54e-04 -6.35e-05 2.38e-04 -4.95e-05 -5.01e-04 1.18e-03 -2.19e-04 -4.10e-06 2.71e-04 -9.08e-05 6.93e-05 -1.90e-05 -8.49e-05 -2.19e-04 4.69e-04
    2.18e-04 7.09e-05 -1.42e-05 -1.47e-05 1.54e-04 3.41e-04 5.82e-05 -8.55e-05 7.09e-05 2.49e-04 -1.23e-05 3.18e-05 1.48e-04 -1.31e-04 -2.91e-05 -1.04e-04 -1.42e-05 -1.23e-05 8.43e-05 -7.81e-06 2.24e-06 -1.23e-04 9.31e-05 -2.33e-05 -1.47e-05 3.18e-05 -7.81e-06 7.84e-05 3.43e-05 -1.02e-04 1.55e-04 7.76e-05 1.54e-04 1.48e-04 2.24e-06 3.43e-05 2.53e-04 -4.81e-06 1.61e-05 -6.30e-05 3.41e-04 -1.31e-04 -1.23e-04 -1.02e-04 -4.81e-06 1.57e-03 -3.27e-04 -3.46e-05 5.82e-05 -2.91e-05 9.31e-05 1.55e-04 1.61e-05 -3.27e-04 1.98e-03 6.85e-05 -8.55e-05 -1.04e-04 -2.33e-05 7.76e-05 -6.30e-05 -3.46e-05 6.85e-05 2.98e-04
    1.47e-03 -2.31e-04 -2.80e-04 4.89e-04 -3.53e-04 -7.07e-06 -2.56e-04 -5.23e-05 -2.31e-04 1.06e-03 1.03e-04 -5.42e-04 4.52e-05 -6.73e-05 2.16e-04 -1.31e-05 -2.80e-04 1.03e-04 7.30e-04 -2.27e-04 -1.58e-05 1.84e-04 1.46e-04 -2.30e-04 4.89e-04 -5.42e-04 -2.27e-04 1.17e-03 -1.36e-04 7.88e-05 -1.96e-04 -1.18e-04 -3.53e-04 4.52e-05 -1.58e-05 -1.36e-04 3.97e-04 -1.02e-04 1.18e-04 2.73e-05 -7.07e-06 -6.73e-05 1.84e-04 7.88e-05 -1.02e-04 1.07e-03 1.37e-04 -4.14e-04 -2.56e-04 2.16e-04 1.46e-04 -1.96e-04 1.18e-04 1.37e-04 4.81e-04 -2.12e-04 -5.23e-05 -1.31e-05 -2.30e-04 -1.18e-04 2.73e-05 -4.14e-04 -2.12e-04 9.75e-04
    1.46e-03 -1.58e-04 -2.97e-04 4.84e-04 -3.37e-04 5.24e-05 -1.95e-04 -2.78e-05 -1.58e-04 8.24e-04 5.11e-06 -5.32e-04 5.48e-05 -9.48e-05 1.69e-04 1.54e-04 -2.97e-04 5.11e-06 8.05e-04 -2.94e-04 2.18e-05 1.85e-04 1.67e-04 -1.99e-04 4.84e-04 -5.32e-04 -2.94e-04 1.44e-03 -1.65e-04 1.93e-05 -2.95e-04 -5.39e-05 -3.37e-04 5.48e-05 2.18e-05 -1.65e-04 3.44e-04 -1.12e-04 1.20e-04 -6.67e-05 5.24e-05 -9.48e-05 1.85e-04 1.93e-05 -1.12e-04 9.80e-04 1.26e-04 -5.74e-04 -1.95e-04 1.69e-04 1.67e-04 -2.95e-04 1.20e-04 1.26e-04 3.64e-04 -1.17e-04 -2.78e-05 1.54e-04 -1.99e-04 -5.39e-05 -6.67e-05 -5.74e-04 -1.17e-04 9.34e-04
    1.21e-03 -4.63e-04 2.40e-05 -1.44e-04 1.04e-04 3.08e-04 -2.52e-04 5.97e-05 -4.63e-04 9.82e-04 -4.79e-05 -5.21e-05 6.99e-05 1.55e-04 1.64e-05 -1.41e-04 2.40e-05 -4.79e-05 4.02e-04 2.53e-04 -1.55e-04 -1.48e-04 -7.55e-05 -1.45e-04 -1.44e-04 -5.21e-05 2.53e-04 3.37e-04 -9.97e-05 -2.20e-04 4.59e-05 -9.34e-05 1.04e-04 6.99e-05 -1.55e-04 -9.97e-05 5.75e-04 9.57e-05 3.57e-04 -4.93e-05 3.08e-04 1.55e-04 -1.48e-04 -2.20e-04 9.57e-05 1.18e-03 -8.14e-05 6.37e-05 -2.52e-04 1.64e-05 -7.55e-05 4.59e-05 3.57e-04 -8.14e-05 5.85e-04 -3.73e-05 5.97e-05 -1.41e-04 -1.45e-04 -9.34e-05 -4.93e-05 6.37e-05 -3.73e-05 2.07e-04
    5.67e-04 8.14e-06 -1.92e-04 8.49e-05 -9.64e-05 7.77e-05 6.85e-05 2.30e-05 8.14e-06 6.46e-04 1.28e-04 4.44e-04 -4.72e-05 -2.01e-04 1.42e-04 -4.87e-04 -1.92e-04 1.28e-04 1.02e-03 -5.63e-04 1.17e-04 -2.95e-04 1.84e-04 -3.58e-04 8.49e-05 4.44e-04 -5.63e-04 1.78e-03 9.20e-05 1.82e-04 -2.12e-04 -2.44e-04 -9.64e-05 -4.72e-05 1.17e-04 9.20e-05 6.36e-04 -1.62e-05 -1.05e-04 2.59e-05 7.77e-05 -2.01e-04 -2.95e-04 1.82e-04 -1.62e-05 5.34e-04 -1.49e-04 2.96e-04 6.85e-05 1.42e-04 1.84e-04 -2.12e-04 -1.05e-04 -1.49e-04 3.63e-04 -1.24e-04 2.30e-05 -4.87e-04 -3.58e-04 -2.44e-04 2.59e-05 2.96e-04 -1.24e-04 1.01e-03
    4.93e-04 -3.02e-04 1.33e-04 -2.95e-06 -1.02e-04 -1.15e-04 -3.10e-04 2.46e-04 -3.02e-04 1.12e-03 -4.64e-04 3.39e-04 -1.32e-05 -9.04e-05 3.02e-04 -1.64e-04 1.33e-04 -4.64e-04 1.04e-03 -4.31e-04 7.01e-05 2.00e-04 -1.81e-04 2.05e-04 -2.95e-06 3.39e-04 -4.31e-04 8.89e-04 -3.62e-04 -2.34e-04 -4.15e-04 -3.33e-04 -1.02e-04 -1.32e-05 7.01e-05 -3.62e-04 3.40e-04 1.77e-04 4.46e-04 1.21e-04 -1.15e-04 -9.04e-05 2.00e-04 -2.34e-04 1.77e-04 6.97e-04 1.84e-04 -1.03e-04 -3.10e-04 3.02e-04 -1.81e-04 -4.15e-04 4.46e-04 1.84e-04 1.66e-03 -6.10e-06 2.46e-04 -1.64e-04 2.05e-04 -3.33e-04 1.21e-04 -1.03e-04 -6.10e-06 6.48e-04
    4.02e-04 -1.74e-04 1.67e-04 4.48e-05 -1.06e-04 -1.33e-04 -2.28e-04 2.31e-04 -1.74e-04 8.28e-04 -5.64e-04 4.52e-04 -1.14e-04 -1.37e-04 -2.22e-05 -1.96e-04 1.67e-04 -5.64e-04 9.93e-04 -4.07e-04 5.62e-05 1.55e-04 -2.83e-04 2.55e-04 4.48e-05 4.52e-04 -4.07e-04 8.18e-04 -3.46e-04 -2.74e-04 -3.53e-04 -3.49e-04 -1.06e-04 -1.14e-04 5.62e-05 -3.46e-04 2.61e-04 2.01e-04 3.28e-04 1.25e-04 -1.33e-04 -1.37e-04 1.55e-04 -2.74e-04 2.01e-04 6.49e-04 2.16e-04 -8.06e-05 -2.28e-04 -2.22e-05 -2.83e-04 -3.53e-04 3.28e-04 2.16e-04 1.27e-03 5.50e-06 2.31e-04 -1.96e-04 2.55e-04 -3.49e-04 1.25e-04 -8.06e-05 5.50e-06 6.09e-04
    1.19e-03 3.31e-05 3.99e-04 -2.78e-04 -2.11e-04 -1.29e-04 9.69e-05 -6.05e-05 3.31e-05 4.61e-04 -1.03e-05 -1.41e-04 8.01e-05 -2.43e-04 7.37e-05 2.19e-04 3.99e-04 -1.03e-05 6.09e-04 -2.79e-04 2.29e-05 -1.96e-04 1.80e-04 -3.46e-04 -2.78e-04 -1.41e-04 -2.79e-04 1.15e-03 3.25e-04 1.30e-04 -1.50e-04 1.45e-04 -2.11e-04 8.01e-05 2.29e-05 3.25e-04 7.80e-04 -2.64e-04 -3.18e-06 2.95e-06 -1.29e-04 -2.43e-04 -1.96e-04 1.30e-04 -2.64e-04 8.18e-04 -4.41e-05 -5.17e-05 9.69e-05 7.37e-05 1.80e-04 -1.50e-04 -3.18e-06 -4.41e-05 1.75e-04 -1.07e-04 -6.05e-05 2.19e-04 -3.46e-04 1.45e-04 2.95e-06 -5.17e-05 -1.07e-04 6.08e-04
    7.35e-04 -2.53e-04 -3.80e-04 1.18e-04 -2.09e-04 -1.51e-04 -1.08e-04 2.04e-04 -2.53e-04 1.32e-03 6.28e-05 -3.18e-04 8.25e-05 1.20e-04 1.82e-04 -2.06e-04 -3.80e-04 6.28e-05 1.22e-03 -2.59e-04 1.03e-04 -1.56e-04 1.15e-04 9.83e-06 1.18e-04 -3.18e-04 -2.59e-04 6.60e-04 -6.88e-05 2.32e-04 7.55e-05 -1.73e-04 -2.09e-04 8.25e-05 1.03e-04 -6.88e-05 2.34e-04 6.95e-05 1.58e-04 -1.30e-05 -1.51e-04 1.20e-04 -1.56e-04 2.32e-04 6.95e-05 9.09e-04 2.37e-04 -6.10e-04 -1.08e-04 1.82e-04 1.15e-04 7.55e-05 1.58e-04 2.37e-04 5.90e-04 -6.28e-05 2.04e-04 -2.06e-04 9.83e-06 -1.73e-04 -1.30e-05 -6.10e-04 -6.28e-05 9.19e-04
    5.96e-04 -2.17e-04 -2.99e-04 1.21e-04 -1.91e-04 -1.85e-04 -1.00e-04 3.07e-04 -2.17e-04 1.01e-03 6.48e-05 -4.84e-04 1.19e-04 1.24e-04 1.11e-04 -4.29e-05 -2.99e-04 6.48e-05 1.31e-03 -3.05e-04 9.14e-05 -2.28e-04 1.29e-04 5.36e-06 1.21e-04 -4.84e-04 -3.05e-04 7.01e-04 -8.43e-05 1.50e-04 -8.28e-06 -1.12e-04 -1.91e-04 1.19e-04 9.14e-05 -8.43e-05 2.21e-04 1.21e-04 1.76e-04 -1.07e-04 -1.85e-04 1.24e-04 -2.28e-04 1.50e-04 1.21e-04 8.12e-04 1.75e-04 -5.84e-04 -1.00e-04 1.11e-04 1.29e-04 -8.28e-06 1.76e-04 1.75e-04 4.64e-04 -2.61e-05 3.07e-04 -4.29e-05 5.36e-06 -1.12e-04 -1.07e-04 -5.84e-04 -2.61e-05 8.35e-04
    7.03e-04 -4.18e-04 2.18e-04 2.06e-05 -3.07e-04 -2.18e-04 -2.56e-04 -2.39e-04 -4.18e-04 9.60e-04 -4.38e-04 8.94e-05 9.81e-05 1.66e-04 -1.77e-04 1.67e-04 2.18e-04 -4.38e-04 1.08e-03 -1.83e-04 4.53e-05 -9.95e-05 4.53e-04 -1.54e-04 2.06e-05 8.94e-05 -1.83e-04 5.47e-04 -7.46e-05 2.06e-04 -3.87e-04 -2.17e-04 -3.07e-04 9.81e-05 4.53e-05 -7.46e-05 3.85e-04 1.80e-04 4.33e-04 1.61e-04 -2.18e-04 1.66e-04 -9.95e-05 2.06e-04 1.80e-04 5.90e-04 -2.84e-05 -9.79e-05 -2.56e-04 -1.77e-04 4.53e-04 -3.87e-04 4.33e-04 -2.84e-05 1.64e-03 2.11e-04 -2.39e-04 1.67e-04 -1.54e-04 -2.17e-04 1.61e-04 -9.79e-05 2.11e-04 7.33e-04
    2.66e-04 1.70e-05 -1.27e-05 -8.92e-05 -1.04e-04 -3.79e-05 -8.24e-05 -4.77e-05 1.70e-05 4.39e-04 -2.17e-04 9.25e-05 -1.41e-04 -3.54e-05 1.46e-06 -8.33e-05 -1.27e-05 -2.17e-04 1.77e-03 -2.22e-04 2.86e-04 -1.39e-04 -2.66e-04 7.74e-04 -8.92e-05 9.25e-05 -2.22e-04 6.96e-04 9.56e-05 -9.71e-05 2.66e-04 -9.49e-05 -1.04e-04 -1.41e-04 2.86e-04 9.56e-05 3.55e-04 2.67e-05 4.05e-05 1.90e-04 -3.79e-05 -3.54e-05 -1.39e-04 -9.71e-05 2.67e-05 2.89e-04 -1.32e-04 7.96e-05 -8.24e-05 1.46e-06 -2.66e-04 2.66e-04 4.05e-05 -1.32e-04 7.15e-04 -8.49e-05 -4.77e-05 -8.33e-05 7.74e-04 -9.49e-05 1.90e-04 7.96e-05 -8.49e-05 8.71e-04
    7.10e-04 -4.30e-04 3.52e-04 4.66e-05 -3.40e-04 -2.61e-04 -2.55e-04 -2.01e-04 -4.30e-04 9.51e-04 -5.58e-04 1.08e-04 1.10e-04 2.26e-04 -2.21e-04 8.07e-05 3.52e-04 -5.58e-04 8.98e-04 -1.13e-04 -6.96e-05 -1.27e-04 1.34e-04 -2.50e-04 4.66e-05 1.08e-04 -1.13e-04 5.13e-04 -7.29e-05 1.98e-04 -2.79e-04 -2.62e-04 -3.40e-04 1.10e-04 -6.96e-05 -7.29e-05 3.39e-04 2.03e-04 3.57e-04 1.36e-04 -2.61e-04 2.26e-04 -1.27e-04 1.98e-04 2.03e-04 6.02e-04 5.39e-05 -9.72e-05 -2.55e-04 -2.21e-04 1.34e-04 -2.79e-04 3.57e-04 5.39e-05 1.37e-03 2.21e-04 -2.01e-04 8.07e-05 -2.50e-04 -2.62e-04 1.36e-04 -9.72e-05 2.21e-04 7.35e-04
    9.09e-04 -4.84e-04 1.03e-05 -5.99e-05 -5.53e-05 -5.64e-05 7.80e-05 -8.67e-05 -4.84e-04 1.28e-03 -1.18e-04 -5.36e-05 3.48e-04 4.91e-04 -1.78e-04 3.43e-05 1.03e-05 -1.18e-04 3.48e-04 2.28e-04 -1.86e-04 -1.31e-04 -1.38e-04 -1.24e-04 -5.99e-05 -5.36e-05 2.28e-04 5.27e-04 -2.01e-04 -2.46e-04 -7.79e-05 -1.57e-04 -5.53e-05 3.48e-04 -1.86e-04 -2.01e-04 6.13e-04 1.80e-04 3.09e-04 1.09e-04 -5.64e-05 4.91e-04 -1.31e-04 -2.46e-04 1.80e-04 1.20e-03 -4.33e-05 8.86e-05 7.80e-05 -1.78e-04 -1.38e-04 -7.79e-05 3.09e-04 -4.33e-05 5.75e-04 9.72e-05 -8.67e-05 3.43e-05 -1.24e-04 -1.57e-04 1.09e-04 8.86e-05 9.72e-05 2.42e-04
    6.70e-04 -1.60e-04 2.12e-05 -1.85e-04 -1.67e-04 -1.79e-04 2.78e-04 -2.84e-04 -1.60e-04 6.12e-04 5.89e-05 -1.80e-04 -1.09e-04 -2.02e-04 -1.09e-04 1.41e-04 2.12e-05 5.89e-05 1.90e-03 -2.97e-04 3.95e-04 2.15e-04 2.26e-04 2.36e-04 -1.85e-04 -1.80e-04 -2.97e-04 6.50e-04 5.32e-05 1.87e-04 -2.19e-04 1.33e-05 -1.67e-04 -1.09e-04 3.95e-04 5.32e-05 5.24e-04 2.07e-04 6.48e-05 8.80e-05 -1.79e-04 -2.02e-04 2.15e-04 1.87e-04 2.07e-04 6.14e-04 -1.77e-04 8.67e-05 2.78e-04 -1.09e-04 2.26e-04 -2.19e-04 6.48e-05 -1.77e-04 1.21e-03 -3.38e-04 -2.84e-04 1.41e-04 2.36e-04 1.33e-05 8.80e-05 8.67e-05 -3.38e-04 7.32e-04
    5.76e-04 -1.79e-04 1.07e-05 5.15e-05 -2.44e-04 -1.91e-04 -6.12e-05 -3.31e-04 -1.79e-04 4.61e-04 -5.28e-05 -1.11e-04 -4.08e-05 2.37e-06 -1.86e-04 1.26e-04 1.07e-05 -5.28e-05 1.54e-03 -1.96e-04 2.96e-04 -1.33e-05 5.65e-04 3.12e-04 5.15e-05 -1.11e-04 -1.96e-04 4.90e-04 -1.01e-05 1.44e-04 -2.88e-04 -1.31e-04 -2.44e-04 -4.08e-05 2.96e-04 -1.01e-05 4.34e-04 1.87e-04 2.50e-04 2.01e-04 -1.91e-04 2.37e-06 -1.33e-05 1.44e-04 1.87e-04 5.07e-04 -2.33e-04 2.72e-05 -6.12e-05 -1.86e-04 5.65e-04 -2.88e-04 2.50e-04 -2.33e-04 1.39e-03 1.10e-04 -3.31e-04 1.26e-04 3.12e-04 -1.31e-04 2.01e-04 2.72e-05 1.10e-04 7.59e-04
    1.52e-03 -4.55e-05 1.86e-05 4.70e-04 -2.30e-04 1.13e-04 -1.08e-04 -7.39e-05 -4.55e-05 5.82e-04 -1.16e-04 -4.23e-04 -3.53e-08 -1.11e-04 1.20e-04 1.16e-04 1.86e-05 -1.16e-04 5.01e-04 -5.14e-05 -2.56e-05 1.97e-04 1.36e-04 -2.64e-04 4.70e-04 -4.23e-04 -5.14e-05 1.26e-03 -2.28e-05 5.09e-05 -2.77e-04 2.37e-05 -2.30e-04 -3.53e-08 -2.56e-05 -2.28e-05 4.19e-04 -9.38e-05 2.95e-05 2.35e-05 1.13e-04 -1.11e-04 1.97e-04 5.09e-05 -9.38e-05 9.76e-04 6.41e-05 -2.93e-04 -1.08e-04 1.20e-04 1.36e-04 -2.77e-04 2.95e-05 6.41e-05 2.36e-04 -1.87e-04 -7.39e-05 1.16e-04 -2.64e-04 2.37e-05 2.35e-05 -2.93e-04 -1.87e-04 8.74e-04
    1.09e-03 -2.61e-04 -5.01e-04 3.52e-04 -3.35e-04 -6.36e-05 -2.31e-04 5.58e-05 -2.61e-04 1.01e-03 1.33e-04 -5.59e-04 7.65e-05 -7.39e-05 1.64e-04 -3.69e-05 -5.01e-04 1.33e-04 8.81e-04 -3.20e-04 7.35e-05 8.79e-05 1.25e-04 -1.69e-04 3.52e-04 -5.59e-04 -3.20e-04 1.03e-03 -1.76e-04 1.06e-04 -1.23e-04 -1.43e-04 -3.35e-04 7.65e-05 7.35e-05 -1.76e-04 3.19e-04 1.59e-05 1.73e-04 2.73e-05 -6.36e-05 -7.39e-05 8.79e-05 1.06e-04 1.59e-05 9.65e-04 1.41e-04 -5.94e-04 -2.31e-04 1.64e-04 1.25e-04 -1.23e-04 1.73e-04 1.41e-04 5.13e-04 -1.62e-04 5.58e-05 -3.69e-05 -1.69e-04 -1.43e-04 2.73e-05 -5.94e-04 -1.62e-04 9.33e-04
    1.19e-03 -3.65e-04 -3.98e-04 3.98e-04 -3.70e-04 -3.75e-05 -2.36e-04 -7.98e-05 -3.65e-04 1.09e-03 1.91e-04 -4.47e-04 1.09e-04 -3.16e-05 2.60e-04 -2.28e-04 -3.98e-04 1.91e-04 6.92e-04 -2.29e-04 3.33e-05 1.78e-04 1.29e-04 -1.28e-04 3.98e-04 -4.47e-04 -2.29e-04 6.92e-04 -1.60e-04 2.17e-04 -1.13e-05 -1.83e-04 -3.70e-04 1.09e-04 3.33e-05 -1.60e-04 2.59e-04 -7.13e-05 1.35e-04 1.18e-04 -3.75e-05 -3.16e-05 1.78e-04 2.17e-04 -7.13e-05 8.25e-04 1.90e-04 -4.78e-04 -2.36e-04 2.60e-04 1.29e-04 -1.13e-05 1.35e-04 1.90e-04 4.47e-04 -2.02e-04 -7.98e-05 -2.28e-04 -1.28e-04 -1.83e-04 1.18e-04 -4.78e-04 -2.02e-04 7.63e-04
    1.33e-03 -1.13e-04 3.18e-04 -8.97e-06 -4.57e-04 -5.21e-05 -1.02e-04 3.32e-04 -1.13e-04 4.47e-04 -1.84e-04 -4.11e-05 7.15e-05 3.08e-06 -3.17e-06 1.26e-04 3.18e-04 -1.84e-04 6.25e-04 -6.45e-05 -2.49e-04 -1.70e-04 1.28e-04 -1.55e-04 -8.97e-06 -4.11e-05 -6.45e-05 6.48e-04 2.71e-04 1.81e-04 -4.15e-05 -8.50e-05 -4.57e-04 7.15e-05 -2.49e-04 2.71e-04 7.64e-04 2.53e-04 3.14e-05 -2.30e-05 -5.21e-05 3.08e-06 -1.70e-04 1.81e-04 2.53e-04 7.01e-04 -2.44e-05 7.12e-06 -1.02e-04 -3.17e-06 1.28e-04 -4.15e-05 3.14e-05 -2.44e-05 2.33e-04 -1.29e-04 3.32e-04 1.26e-04 -1.55e-04 -8.50e-05 -2.30e-05 7.12e-06 -1.29e-04 5.20e-04
    1.18e-03 -3.55e-04 5.54e-05 -5.42e-04 9.27e-05 1.73e-04 4.43e-05 3.06e-04 -3.55e-04 5.12e-04 -2.64e-05 2.58e-04 2.97e-05 -2.98e-04 1.29e-04 -1.13e-04 5.54e-05 -2.64e-05 4.73e-04 -5.19e-05 2.11e-04 8.54e-05 1.81e-05 -2.15e-04 -5.42e-04 2.58e-04 -5.19e-05 8.87e-04 -1.97e-04 7.33e-07 -4.20e-05 -1.11e-04 9.27e-05 2.97e-05 2.11e-04 -1.97e-04 9.88e-04 -1.22e-04 1.64e-04 -1.74e-04 1.73e-04 -2.98e-04 8.54e-05 7.33e-07 -1.22e-04 5.17e-04 -1.74e-04 -2.42e-05 4.43e-05 1.29e-04 1.81e-05 -4.20e-05 1.64e-04 -1.74e-04 1.61e-04 -5.13e-05 3.06e-04 -1.13e-04 -2.15e-04 -1.11e-04 -1.74e-04 -2.42e-05 -5.13e-05 6.54e-04
    1.51e-03 -1.49e-04 3.08e-04 -2.64e-04 -4.56e-04 1.69e-04 2.96e-05 3.00e-04 -1.49e-04 3.91e-04 -8.44e-05 1.81e-04 2.14e-05 -1.04e-04 1.15e-04 -1.07e-05 3.08e-04 -8.44e-05 4.86e-04 -1.99e-04 8.41e-05 2.17e-05 9.02e-05 -7.43e-05 -2.64e-04 1.81e-04 -1.99e-04 1.01e-03 -1.63e-04 3.08e-04 -1.53e-04 6.76e-06 -4.56e-04 2.14e-05 8.41e-05 -1.63e-04 9.55e-04 -2.35e-04 1.29e-04 -1.14e-04 1.69e-04 -1.04e-04 2.17e-05 3.08e-04 -2.35e-04 5.04e-04 -1.35e-04 -7.47e-05 2.96e-05 1.15e-04 9.02e-05 -1.53e-04 1.29e-04 -1.35e-04 2.27e-04 1.95e-05 3.00e-04 -1.07e-05 -7.43e-05 6.76e-06 -1.14e-04 -7.47e-05 1.95e-05 5.69e-04
    6.44e-04 1.11e-04 3.41e-05 4.26e-05 -2.29e-04 -1.14e-04 -6.87e-05 -4.73e-04 1.11e-04 4.91e-04 -2.17e-04 1.38e-04 -1.34e-04 -4.99e-05 -1.36e-04 -1.67e-04 3.41e-05 -2.17e-04 1.46e-03 -5.76e-04 3.32e-04 5.14e-05 -1.12e-04 -2.27e-06 4.26e-05 1.38e-04 -5.76e-04 1.39e-03 -3.16e-05 -1.04e-04 -3.65e-05 -1.24e-04 -2.29e-04 -1.34e-04 3.32e-04 -3.16e-05 6.23e-04 5.38e-05 8.17e-05 2.99e-04 -1.14e-04 -4.99e-05 5.14e-05 -1.04e-04 5.38e-05 4.46e-04 1.69e-06 1.11e-04 -6.87e-05 -1.36e-04 -1.12e-04 -3.65e-05 8.17e-05 1.69e-06 6.77e-04 9.11e-05 -4.73e-04 -1.67e-04 -2.27e-06 -1.24e-04 2.99e-04 1.11e-04 9.11e-05 1.13e-03
    4.36e-04 2.61e-05 -1.08e-04 1.57e-05 -1.44e-04 2.67e-05 3.32e-05 -1.29e-04 2.61e-05 6.44e-04 5.63e-05 2.84e-04 -5.60e-05 -1.65e-04 8.87e-05 -4.16e-04 -1.08e-04 5.63e-05 1.13e-03 -6.31e-04 1.33e-04 -2.85e-04 5.36e-05 -2.11e-04 1.57e-05 2.84e-04 -6.31e-04 1.71e-03 2.55e-04 1.85e-04 -3.04e-04 -8.05e-05 -1.44e-04 -5.60e-05 1.33e-04 2.55e-04 5.97e-04 -5.41e-05 -1.38e-05 1.32e-04 2.67e-05 -1.65e-04 -2.85e-04 1.85e-04 -5.41e-05 5.36e-04 -1.17e-04 2.56e-04 3.32e-05 8.87e-05 5.36e-05 -3.04e-04 -1.38e-05 -1.17e-04 5.88e-04 -1.45e-04 -1.29e-04 -4.16e-04 -2.11e-04 -8.05e-05 1.32e-04 2.56e-04 -1.45e-04 1.14e-03
    7.69e-04 -1.71e-04 -4.41e-04 2.40e-04 -2.77e-04 -1.50e-04 -1.48e-04 2.71e-04 -1.71e-04 8.58e-04 1.35e-04 -5.78e-04 7.39e-05 -2.56e-05 1.15e-04 6.53e-05 -4.41e-04 1.35e-04 1.21e-03 -3.92e-04 1.10e-04 -1.96e-04 1.59e-04 -1.38e-04 2.40e-04 -5.78e-04 -3.92e-04 1.09e-03 -1.81e-04 3.14e-05 -1.29e-04 -9.72e-05 -2.77e-04 7.39e-05 1.10e-04 -1.81e-04 2.85e-04 1.49e-04 1.72e-04 -9.40e-05 -1.50e-04 -2.56e-05 -1.96e-04 3.14e-05 1.49e-04 9.08e-04 9.76e-05 -5.84e-04 -1.48e-04 1.15e-04 1.59e-04 -1.29e-04 1.72e-04 9.76e-05 4.57e-04 -4.91e-05 2.71e-04 6.53e-05 -1.38e-04 -9.72e-05 -9.40e-05 -5.84e-04 -4.91e-05 9.35e-04
    5.36e-04 -7.37e-05 -2.05e-05 2.29e-04 1.63e-04 -8.78e-05 9.28e-05 2.23e-04 -7.37e-05 1.76e-04 -8.69e-05 1.39e-04 5.86e-05 1.29e-04 5.34e-05 -1.67e-04 -2.05e-05 -8.69e-05 1.07e-04 -5.10e-05 1.66e-05 -8.22e-06 -2.27e-05 1.45e-04 2.29e-04 1.39e-04 -5.10e-05 1.08e-03 6.55e-04 2.51e-04 -1.62e-04 -3.50e-04 1.63e-04 5.86e-05 1.66e-05 6.55e-04 1.02e-03 2.98e-04 4.42e-05 4.97e-05 -8.78e-05 1.29e-04 -8.22e-06 2.51e-04 2.98e-04 2.81e-04 6.50e-05 -4.64e-05 9.28e-05 5.34e-05 -2.27e-05 -1.62e-04 4.42e-05 6.50e-05 2.42e-04 2.60e-04 2.23e-04 -1.67e-04 1.45e-04 -3.50e-04 4.97e-05 -4.64e-05 2.60e-04 8.46e-04
    8.42e-04 -1.26e-04 -4.25e-05 -1.57e-04 9.62e-05 -2.07e-05 1.01e-04 3.98e-04 -1.26e-04 5.31e-04 -1.33e-06 2.34e-04 -1.62e-04 -1.39e-04 1.28e-04 -1.07e-04 -4.25e-05 -1.33e-06 4.51e-04 3.56e-05 1.64e-04 4.66e-05 1.29e-05 4.11e-05 -1.57e-04 2.34e-04 3.56e-05 1.32e-03 -4.74e-04 2.46e-04 -1.83e-04 2.21e-04 9.62e-05 -1.62e-04 1.64e-04 -4.74e-04 1.40e-03 -2.24e-06 1.22e-04 1.47e-04 -2.07e-05 -1.39e-04 4.66e-05 2.46e-04 -2.24e-06 4.49e-04 -1.88e-04 1.05e-05 1.01e-04 1.28e-04 1.29e-05 -1.83e-04 1.22e-04 -1.88e-04 2.50e-04 1.69e-05 3.98e-04 -1.07e-04 4.11e-05 2.21e-04 1.47e-04 1.05e-05 1.69e-05 1.00e-03
    8.07e-04 -2.88e-04 -5.11e-04 2.21e-04 -2.39e-04 -8.63e-05 -1.78e-04 1.01e-04 -2.88e-04 1.19e-03 1.62e-04 -3.96e-04 8.52e-05 2.17e-05 2.38e-04 -2.42e-04 -5.11e-04 1.62e-04 8.72e-04 -2.79e-04 1.21e-04 7.84e-06 1.25e-04 -7.94e-05 2.21e-04 -3.96e-04 -2.79e-04 6.86e-04 -1.15e-04 2.18e-04 3.94e-05 -1.71e-04 -2.39e-04 8.52e-05 1.21e-04 -1.15e-04 2.01e-04 1.97e-05 1.43e-04 3.68e-05 -8.63e-05 2.17e-05 7.84e-06 2.18e-04 1.97e-05 8.39e-04 2.17e-04 -6.02e-04 -1.78e-04 2.38e-04 1.25e-04 3.94e-05 1.43e-04 2.17e-04 4.90e-04 -1.44e-04 1.01e-04 -2.42e-04 -7.94e-05 -1.71e-04 3.68e-05 -6.02e-04 -1.44e-04 8.46e-04
    7.28e-04 -1.20e-04 -9.23e-05 -1.61e-04 3.40e-04 -1.37e-05 2.12e-04 -5.30e-05 -1.20e-04 1.23e-03 -3.44e-05 2.63e-05 3.28e-04 -3.05e-04 -4.52e-04 -2.66e-04 -9.23e-05 -3.44e-05 3.31e-04 8.45e-05 -1.16e-04 -2.82e-04 1.22e-04 -7.99e-05 -1.61e-04 2.63e-05 8.45e-05 2.40e-04 -3.76e-05 -4.94e-05 5.57e-05 5.14e-05 3.40e-04 3.28e-04 -1.16e-04 -3.76e-05 8.12e-04 -7.50e-05 -4.49e-05 -6.49e-05 -1.37e-05 -3.05e-04 -2.82e-04 -4.94e-05 -7.50e-05 9.61e-04 -3.49e-04 1.61e-04 2.12e-04 -4.52e-04 1.22e-04 5.57e-05 -4.49e-05 -3.49e-04 1.08e-03 6.48e-05 -5.30e-05 -2.66e-04 -7.99e-05 5.14e-05 -6.49e-05 1.61e-04 6.48e-05 2.97e-04
    5.64e-04 -1.78e-04 -3.60e-04 1.93e-04 -2.20e-04 -2.19e-04 -1.07e-04 2.93e-04 -1.78e-04 8.22e-04 1.27e-04 -5.25e-04 1.16e-04 4.02e-05 1.71e-04 -5.39e-07 -3.60e-04 1.27e-04 1.14e-03 -3.38e-04 8.03e-05 -2.48e-04 1.63e-04 -1.02e-04 1.93e-04 -5.25e-04 -3.38e-04 7.37e-04 -1.66e-04 6.05e-05 -6.57e-05 -1.54e-04 -2.20e-04 1.16e-04 8.03e-05 -1.66e-04 1.85e-04 1.56e-04 1.33e-04 -8.18e-05 -2.19e-04 4.02e-05 -2.48e-04 6.05e-05 1.56e-04 7.80e-04 1.26e-04 -5.42e-04 -1.07e-04 1.71e-04 1.63e-04 -6.57e-05 1.33e-04 1.26e-04 2.89e-04 -5.65e-05 2.93e-04 -5.39e-07 -1.02e-04 -1.54e-04 -8.18e-05 -5.42e-04 -5.65e-05 8.28e-04
    8.27e-04 -1.07e-04 -1.24e-05 -6.32e-05 -5.90e-05 -8.15e-05 1.47e-04 2.81e-04 -1.07e-04 4.83e-04 -1.49e-04 2.89e-04 -1.46e-04 -5.29e-05 7.34e-05 -2.63e-04 -1.24e-05 -1.49e-04 3.26e-04 -3.03e-05 1.57e-04 4.01e-05 -8.94e-06 8.05e-05 -6.32e-05 2.89e-04 -3.03e-05 1.06e-03 -6.91e-05 2.86e-04 -2.29e-04 -2.24e-04 -5.90e-05 -1.46e-04 1.57e-04 -6.91e-05 7.74e-04 9.02e-05 -5.07e-06 2.42e-04 -8.15e-05 -5.29e-05 4.01e-05 2.86e-04 9.02e-05 4.34e-04 -1.58e-04 -8.36e-05 1.47e-04 7.34e-05 -8.94e-06 -2.29e-04 -5.07e-06 -1.58e-04 3.01e-04 1.27e-04 2.81e-04 -2.63e-04 8.05e-05 -2.24e-04 2.42e-04 -8.36e-05 1.27e-04 9.34e-04
    1.53e-03 -3.48e-05 2.46e-04 2.63e-04 -2.77e-04 -3.94e-05 8.02e-05 -1.17e-04 -3.48e-05 1.32e-03 1.48e-04 2.11e-04 -1.57e-04 -2.81e-05 1.76e-04 -1.66e-04 2.46e-04 1.48e-04 2.89e-04 2.27e-04 -2.30e-04 4.16e-05 9.00e-05 -4.48e-05 2.63e-04 2.11e-04 2.27e-04 3.59e-04 -2.64e-04 8.74e-05 1.28e-04 -5.47e-05 -2.77e-04 -1.57e-04 -2.30e-04 -2.64e-04 5.21e-04 -7.96e-05 1.25e-04 2.36e-05 -3.94e-05 -2.81e-05 4.16e-05 8.74e-05 -7.96e-05 4.59e-04 9.46e-06 9.83e-05 8.02e-05 1.76e-04 9.00e-05 1.28e-04 1.25e-04 9.46e-06 4.15e-04 -3.79e-06 -1.17e-04 -1.66e-04 -4.48e-05 -5.47e-05 2.36e-05 9.83e-05 -3.79e-06 1.45e-04
    1.08e-03 -1.92e-04 4.19e-05 -2.72e-04 3.07e-05 5.34e-05 7.50e-05 4.11e-04 -1.92e-04 5.54e-04 -2.83e-05 2.34e-04 -1.02e-04 -1.89e-04 1.48e-04 -1.14e-04 4.19e-05 -2.83e-05 5.19e-04 -5.13e-05 1.71e-04 3.82e-05 4.55e-05 -4.74e-05 -2.72e-04 2.34e-04 -5.13e-05 1.17e-03 -5.09e-04 2.19e-04 -1.57e-04 -4.10e-05 3.07e-05 -1.02e-04 1.71e-04 -5.09e-04 1.22e-03 -8.94e-05 1.61e-04 -1.04e-04 5.34e-05 -1.89e-04 3.82e-05 2.19e-04 -8.94e-05 5.18e-04 -1.92e-04 -3.41e-05 7.50e-05 1.48e-04 4.55e-05 -1.57e-04 1.61e-04 -1.92e-04 2.59e-04 1.04e-05 4.11e-04 -1.14e-04 -4.74e-05 -4.10e-05 -1.04e-04 -3.41e-05 1.04e-05 7.51e-04
    8.47e-04 -4.75e-05 2.46e-05 1.86e-04 2.31e-05 2.22e-06 1.33e-04 2.97e-04 -4.75e-05 3.61e-04 -1.26e-04 2.82e-04 -2.43e-06 2.59e-05 8.92e-05 -1.42e-04 2.46e-05 -1.26e-04 2.57e-04 -1.06e-04 7.60e-05 -3.21e-06 1.96e-05 8.82e-05 1.86e-04 2.82e-04 -1.06e-04 1.43e-03 1.19e-04 3.27e-04 -2.15e-04 -3.68e-04 2.31e-05 -2.43e-06 7.60e-05 1.19e-04 6.66e-04 1.33e-04 5.54e-05 8.51e-05 2.22e-06 2.59e-05 -3.21e-06 3.27e-04 1.33e-04 3.79e-04 -4.20e-05 -7.41e-05 1.33e-04 8.92e-05 1.96e-05 -2.15e-04 5.54e-05 -4.20e-05 3.27e-04 2.36e-04 2.97e-04 -1.42e-04 8.82e-05 -3.68e-04 8.51e-05 -7.41e-05 2.36e-04 9.03e-04
    1.31e-03 -1.58e-04 3.39e-04 3.14e-04 -2.60e-04 -2.99e-04 -1.52e-04 1.46e-04 -1.58e-04 7.11e-04 1.77e-04 -3.38e-04 -1.12e-04 -2.42e-04 1.73e-04 -8.40e-05 3.39e-04 1.77e-04 6.82e-04 -3.05e-05 -1.62e-04 -3.94e-04 2.01e-04 -3.05e-04 3.14e-04 -3.38e-04 -3.05e-05 7.50e-04 3.66e-05 4.68e-06 -8.47e-05 4.09e-05 -2.60e-04 -1.12e-04 -1.62e-04 3.66e-05 6.31e-04 -8.06e-06 -8.37e-05 6.20e-05 -2.99e-04 -2.42e-04 -3.94e-04 4.68e-06 -8.06e-06 9.48e-04 -2.51e-05 2.16e-04 -1.52e-04 1.73e-04 2.01e-04 -8.47e-05 -8.37e-05 -2.51e-05 3.27e-04 -2.42e-04 1.46e-04 -8.40e-05 -3.05e-04 4.09e-05 6.20e-05 2.16e-04 -2.42e-04 6.19e-04
    6.66e-04 2.96e-05 1.78e-04 1.90e-04 -3.94e-04 -2.92e-04 -2.53e-05 -2.72e-04 2.96e-05 4.67e-04 -2.24e-05 -1.68e-04 -1.11e-04 -4.19e-05 -6.44e-05 2.67e-05 1.78e-04 -2.24e-05 1.43e-03 -4.45e-04 1.48e-04 -5.11e-05 9.05e-05 -3.33e-04 1.90e-04 -1.68e-04 -4.45e-04 1.22e-03 -2.11e-04 -1.61e-04 -7.97e-05 -2.27e-04 -3.94e-04 -1.11e-04 1.48e-04 -2.11e-04 5.27e-04 2.94e-04 9.51e-05 1.64e-04 -2.92e-04 -4.19e-05 -5.11e-05 -1.61e-04 2.94e-04 5.63e-04 5.71e-05 2.63e-04 -2.53e-05 -6.44e-05 9.05e-05 -7.97e-05 9.51e-05 5.71e-05 5.67e-04 1.91e-04 -2.72e-04 2.67e-05 -3.33e-04 -2.27e-04 1.64e-04 2.63e-04 1.91e-04 9.74e-04
    5.77e-04 -2.44e-04 -2.44e-04 1.56e-04 -1.56e-04 -1.20e-04 8.90e-05 1.63e-05 -2.44e-04 1.10e-03 1.58e-04 3.82e-04 -8.33e-05 -2.81e-04 -9.55e-05 -4.28e-05 -2.44e-04 1.58e-04 4.96e-04 -9.39e-05 3.98e-06 1.03e-04 -1.37e-04 1.69e-05 1.56e-04 3.82e-04 -9.39e-05 8.90e-04 -3.83e-04 -4.04e-04 -4.45e-05 -3.97e-04 -1.56e-04 -8.33e-05 3.98e-06 -3.83e-04 3.77e-04 2.03e-04 1.17e-04 2.78e-04 -1.20e-04 -2.81e-04 1.03e-04 -4.04e-04 2.03e-04 8.10e-04 1.94e-04 1.64e-04 8.90e-05 -9.55e-05 -1.37e-04 -4.45e-05 1.17e-04 1.94e-04 9.91e-04 1.28e-04 1.63e-05 -4.28e-05 1.69e-05 -3.97e-04 2.78e-04 1.64e-04 1.28e-04 7.32e-04
    5.91e-04 -2.51e-04 -5.16e-04 2.22e-04 -1.84e-04 -1.01e-04 -1.73e-04 2.05e-04 -2.51e-04 8.63e-04 1.59e-04 -5.20e-04 1.45e-04 3.98e-05 1.75e-04 -9.11e-05 -5.16e-04 1.59e-04 8.51e-04 -3.35e-04 1.43e-04 -6.40e-05 1.20e-04 -5.05e-05 2.22e-04 -5.20e-04 -3.35e-04 6.49e-04 -1.04e-04 1.74e-04 -6.71e-05 -9.43e-05 -1.84e-04 1.45e-04 1.43e-04 -1.04e-04 1.35e-04 3.83e-05 1.68e-04 -5.07e-05 -1.01e-04 3.98e-05 -6.40e-05 1.74e-04 3.83e-05 6.81e-04 1.67e-04 -6.22e-04 -1.73e-04 1.75e-04 1.20e-04 -6.71e-05 1.68e-04 1.67e-04 3.36e-04 -9.40e-05 2.05e-04 -9.11e-05 -5.05e-05 -9.43e-05 -5.07e-05 -6.22e-04 -9.40e-05 6.90e-04
    8.27e-04 -2.24e-04 2.56e-05 -1.71e-04 9.65e-05 -1.25e-04 1.30e-04 -8.39e-05 -2.24e-04 1.40e-03 -1.31e-04 3.65e-05 2.39e-04 5.87e-04 -9.28e-05 -1.52e-04 2.56e-05 -1.31e-04 3.09e-04 1.54e-04 -1.73e-04 3.91e-05 -4.30e-05 -6.16e-05 -1.71e-04 3.65e-05 1.54e-04 3.57e-04 -1.01e-04 1.65e-04 1.67e-05 -1.52e-05 9.65e-05 2.39e-04 -1.73e-04 -1.01e-04 5.84e-04 1.66e-05 2.32e-04 1.71e-05 -1.25e-04 5.87e-04 3.91e-05 1.65e-04 1.66e-05 9.11e-04 -3.56e-05 3.56e-05 1.30e-04 -9.28e-05 -4.30e-05 1.67e-05 2.32e-04 -3.56e-05 6.24e-04 -7.29e-06 -8.39e-05 -1.52e-04 -6.16e-05 -1.52e-05 1.71e-05 3.56e-05 -7.29e-06 2.56e-04 ];

locQmati = [0.00e+00 1.66e-04 -4.23e-05 5.21e-05 4.49e-05 -1.93e-04 2.36e-04 -3.59e-05 -1.66e-04 -7.11e-22 1.66e-05 -3.69e-05 -5.35e-05 -6.98e-05 -7.88e-04 6.04e-05 4.23e-05 -1.66e-05 -8.93e-22 8.56e-06 -1.65e-05 3.02e-05 3.00e-04 8.62e-06 -5.21e-05 3.69e-05 -8.56e-06 -1.95e-22 -1.74e-05 -9.91e-05 -1.05e-04 -7.11e-05 -4.49e-05 5.35e-05 1.65e-05 1.74e-05 -7.80e-22 -1.18e-04 -3.17e-05 3.99e-06 1.93e-04 6.98e-05 -3.02e-05 9.91e-05 1.18e-04 1.08e-21 1.14e-03 1.32e-04 -2.36e-04 7.88e-04 -3.00e-04 1.05e-04 3.17e-05 -1.14e-03 1.18e-22 -8.37e-04 3.59e-05 -6.04e-05 -8.62e-06 7.11e-05 -3.99e-06 -1.32e-04 8.37e-04 4.43e-22
    0.00e+00 -1.43e-04 2.19e-04 1.30e-04 -4.73e-04 -2.94e-04 1.03e-04 1.36e-04 1.43e-04 1.50e-21 1.49e-04 -5.32e-06 -2.38e-04 -2.81e-04 -4.74e-05 1.21e-05 -2.19e-04 -1.49e-04 -5.93e-22 1.44e-04 -2.67e-04 -2.04e-04 9.14e-05 -3.53e-04 -1.30e-04 5.32e-06 -1.44e-04 -1.70e-21 2.29e-04 2.52e-04 7.88e-05 5.54e-06 4.73e-04 2.38e-04 2.67e-04 -2.29e-04 -3.35e-23 5.09e-04 -9.56e-05 1.38e-03 2.94e-04 2.81e-04 2.04e-04 -2.52e-04 -5.09e-04 9.68e-22 -3.06e-04 6.13e-04 -1.03e-04 4.74e-05 -9.14e-05 -7.88e-05 9.56e-05 3.06e-04 -4.08e-22 1.29e-04 -1.36e-04 -1.21e-05 3.53e-04 -5.54e-06 -1.38e-03 -6.13e-04 -1.29e-04 -3.95e-22
    0.00e+00 -2.37e-04 2.78e-04 2.35e-04 -3.84e-04 -2.43e-04 2.57e-04 9.93e-05 2.37e-04 -4.35e-22 1.33e-04 -3.65e-05 -2.96e-04 -4.51e-04 -3.59e-05 2.52e-05 -2.78e-04 -1.33e-04 -1.30e-21 1.23e-04 -2.06e-04 -1.63e-04 1.44e-04 -3.08e-04 -2.35e-04 3.65e-05 -1.23e-04 9.04e-22 2.03e-04 2.62e-04 1.49e-04 2.35e-05 3.84e-04 2.96e-04 2.06e-04 -2.03e-04 -7.37e-22 5.77e-04 -2.31e-04 1.20e-03 2.43e-04 4.51e-04 1.63e-04 -2.62e-04 -5.77e-04 1.59e-21 -4.07e-04 4.43e-04 -2.57e-04 3.59e-05 -1.44e-04 -1.49e-04 2.31e-04 4.07e-04 5.13e-22 1.32e-04 -9.93e-05 -2.52e-05 3.08e-04 -2.35e-05 -1.20e-03 -4.43e-04 -1.32e-04 -7.52e-21
    0.00e+00 -2.29e-04 3.53e-04 2.49e-04 -4.33e-04 6.97e-05 3.33e-04 1.88e-04 2.29e-04 2.96e-21 1.39e-04 -2.51e-06 -4.02e-04 -3.28e-04 1.07e-05 -1.55e-05 -3.53e-04 -1.39e-04 1.38e-21 1.29e-04 -1.80e-04 -1.10e-04 1.48e-04 -2.02e-04 -2.49e-04 2.51e-06 -1.29e-04 2.31e-22 2.15e-04 2.12e-04 1.10e-04 -1.77e-06 4.33e-04 4.02e-04 1.80e-04 -2.15e-04 -9.22e-22 5.63e-04 -2.77e-04 9.33e-04 -6.97e-05 3.28e-04 1.10e-04 -2.12e-04 -5.63e-04 -1.19e-21 -3.07e-04 4.92e-04 -3.33e-04 -1.07e-05 -1.48e-04 -1.10e-04 2.77e-04 3.07e-04 -4.37e-22 8.64e-05 -1.88e-04 1.55e-05 2.02e-04 1.77e-06 -9.33e-04 -4.92e-04 -8.64e-05 -1.65e-22
    0.00e+00 8.65e-05 -8.76e-05 -1.15e-04 -4.60e-05 -5.89e-05 3.07e-05 -8.87e-06 -8.65e-05 -2.09e-21 -3.51e-04 8.68e-05 -1.02e-04 1.97e-04 -4.60e-04 -3.54e-05 8.76e-05 3.51e-04 -1.10e-21 -8.49e-05 -2.15e-04 -1.24e-05 -5.19e-04 1.07e-04 1.15e-04 -8.68e-05 8.49e-05 -4.63e-23 3.06e-05 -1.69e-06 2.20e-04 -5.86e-05 4.60e-05 1.02e-04 2.15e-04 -3.06e-05 -2.14e-21 -3.48e-05 1.06e-04 -2.37e-05 5.89e-05 -1.97e-04 1.24e-05 1.69e-06 3.48e-05 3.56e-22 9.15e-05 -1.19e-04 -3.07e-05 4.60e-04 5.19e-04 -2.20e-04 -1.06e-04 -9.15e-05 -9.93e-22 2.78e-04 8.87e-06 3.54e-05 -1.07e-04 5.86e-05 2.37e-05 1.19e-04 -2.78e-04 2.26e-21
    0.00e+00 4.06e-04 -2.13e-05 -2.63e-05 1.56e-05 -3.58e-05 4.19e-04 7.78e-05 -4.06e-04 -1.42e-21 -1.12e-04 3.09e-04 -2.18e-04 2.46e-04 -5.75e-04 -2.41e-04 2.13e-05 1.12e-04 2.64e-21 -6.19e-06 -2.98e-05 -2.57e-05 -6.35e-05 2.09e-04 2.63e-05 -3.09e-04 6.19e-06 2.25e-21 -5.95e-05 -1.50e-05 -1.53e-04 -7.03e-05 -1.56e-05 2.18e-04 2.98e-05 5.95e-05 1.35e-22 1.46e-06 9.97e-05 4.46e-06 3.58e-05 -2.46e-04 2.57e-05 1.50e-05 -1.46e-06 -4.29e-22 -3.28e-04 -1.95e-05 -4.19e-04 5.75e-04 6.35e-05 1.53e-04 -9.97e-05 3.28e-04 1.16e-22 -1.39e-04 -7.78e-05 2.41e-04 -2.09e-04 7.03e-05 -4.46e-06 1.95e-05 1.39e-04 -1.00e-21
    0.00e+00 5.13e-04 -1.87e-05 5.64e-05 7.46e-05 -5.18e-05 5.11e-04 6.38e-05 -5.13e-04 1.41e-21 5.52e-05 3.97e-04 -2.52e-04 4.25e-05 -3.05e-04 -2.72e-05 1.87e-05 -5.52e-05 -1.52e-21 -1.01e-05 -3.87e-05 5.86e-05 -2.21e-04 1.18e-04 -5.64e-05 -3.97e-04 1.01e-05 -2.85e-21 -7.06e-05 -1.05e-04 -1.49e-04 -1.11e-04 -7.46e-05 2.52e-04 3.87e-05 7.06e-05 -1.51e-21 -4.74e-05 1.44e-04 2.14e-05 5.18e-05 -4.25e-05 -5.86e-05 1.05e-04 4.74e-05 -2.66e-22 -3.12e-04 5.80e-05 -5.11e-04 3.05e-04 2.21e-04 1.49e-04 -1.44e-04 3.12e-04 1.55e-21 -1.45e-04 -6.38e-05 2.72e-05 -1.18e-04 1.11e-04 -2.14e-05 -5.80e-05 1.45e-04 -2.29e-21
    0.00e+00 -1.37e-04 2.94e-05 1.26e-04 8.06e-04 1.62e-04 1.89e-04 1.65e-04 1.37e-04 2.34e-21 -1.05e-04 -1.51e-04 -2.42e-04 -2.16e-04 -9.10e-07 -6.75e-05 -2.94e-05 1.05e-04 1.57e-21 2.95e-05 1.89e-04 -1.68e-05 -7.91e-06 -2.20e-04 -1.26e-04 1.51e-04 -2.95e-05 5.02e-22 2.22e-04 3.11e-04 1.54e-04 2.98e-05 -8.06e-04 2.42e-04 -1.89e-04 -2.22e-04 1.11e-23 9.80e-04 5.03e-04 8.90e-04 -1.62e-04 2.16e-04 1.68e-05 -3.11e-04 -9.80e-04 9.77e-22 -5.51e-05 7.72e-05 -1.89e-04 9.10e-07 7.91e-06 -1.54e-04 -5.03e-04 5.51e-05 -6.40e-22 2.29e-04 -1.65e-04 6.75e-05 2.20e-04 -2.98e-05 -8.90e-04 -7.72e-05 -2.29e-04 -6.01e-22
    0.00e+00 1.14e-04 7.14e-06 -6.18e-05 -7.59e-05 -1.17e-04 3.14e-04 -5.70e-05 -1.14e-04 5.05e-23 3.89e-04 6.15e-04 -1.87e-04 1.84e-04 -4.67e-04 -2.62e-05 -7.14e-06 -3.89e-04 -2.29e-21 -6.03e-05 -7.61e-05 -2.22e-04 -1.11e-04 1.25e-04 6.18e-05 -6.15e-04 6.03e-05 2.31e-21 -1.63e-04 -3.86e-04 -1.45e-04 5.25e-05 7.59e-05 1.87e-04 7.61e-05 1.63e-04 -4.34e-22 2.34e-05 1.33e-04 -2.81e-05 1.17e-04 -1.84e-04 2.22e-04 3.86e-04 -2.34e-05 6.11e-22 -3.23e-04 -2.40e-05 -3.14e-04 4.67e-04 1.11e-04 1.45e-04 -1.33e-04 3.23e-04 -6.14e-22 -1.47e-04 5.70e-05 2.62e-05 -1.25e-04 -5.25e-05 2.81e-05 2.40e-05 1.47e-04 1.54e-21
    0.00e+00 -3.06e-04 1.82e-04 2.86e-04 -2.46e-04 -1.83e-04 2.43e-04 1.09e-04 3.06e-04 -6.02e-22 1.12e-04 -3.26e-05 -2.94e-04 -4.57e-04 -7.48e-06 -2.89e-05 -1.82e-04 -1.12e-04 -3.20e-21 1.47e-04 -1.76e-04 -5.71e-05 1.36e-04 -2.30e-04 -2.86e-04 3.26e-05 -1.47e-04 7.97e-21 1.70e-04 2.45e-04 7.88e-05 3.10e-05 2.46e-04 2.94e-04 1.76e-04 -1.70e-04 -4.49e-22 5.61e-04 -2.61e-04 9.18e-04 1.83e-04 4.57e-04 5.71e-05 -2.45e-04 -5.61e-04 7.83e-21 -4.24e-04 3.34e-04 -2.43e-04 7.48e-06 -1.36e-04 -7.88e-05 2.61e-04 4.24e-04 1.77e-21 7.13e-05 -1.09e-04 2.89e-05 2.30e-04 -3.10e-05 -9.18e-04 -3.34e-04 -7.13e-05 -5.13e-21
    0.00e+00 -2.30e-04 5.47e-05 1.37e-04 8.10e-04 7.12e-05 2.43e-05 5.43e-04 2.30e-04 2.11e-22 -1.52e-04 -3.19e-04 -6.00e-04 -2.43e-04 6.62e-05 -2.30e-04 -5.47e-05 1.52e-04 -2.90e-21 9.72e-05 1.05e-05 4.86e-05 -1.29e-05 -2.64e-04 -1.37e-04 3.19e-04 -9.72e-05 1.06e-21 -3.04e-04 2.78e-04 6.93e-05 -4.14e-04 -8.10e-04 6.00e-04 -1.05e-05 3.04e-04 1.35e-21 1.01e-03 -6.49e-06 -8.69e-05 -7.12e-05 2.43e-04 -4.86e-05 -2.78e-04 -1.01e-03 2.37e-21 7.70e-05 -3.96e-04 -2.43e-05 -6.62e-05 1.29e-05 -6.93e-05 6.49e-06 -7.70e-05 -9.83e-22 2.76e-04 -5.43e-04 2.30e-04 2.64e-04 4.14e-04 8.69e-05 3.96e-04 -2.76e-04 3.26e-21
    0.00e+00 -1.70e-04 4.72e-04 2.02e-04 -4.25e-04 1.50e-04 3.14e-04 1.48e-04 1.70e-04 5.48e-21 1.53e-04 2.38e-05 -4.07e-04 -2.48e-04 1.25e-05 2.22e-05 -4.72e-04 -1.53e-04 -6.81e-22 1.50e-04 -1.63e-04 -1.60e-04 1.20e-04 -2.02e-04 -2.02e-04 -2.38e-05 -1.50e-04 -4.63e-22 2.70e-04 1.71e-04 1.45e-04 1.30e-04 4.25e-04 4.07e-04 1.63e-04 -2.70e-04 1.36e-21 4.05e-04 -2.80e-04 8.64e-04 -1.50e-04 2.48e-04 1.60e-04 -1.71e-04 -4.05e-04 -2.45e-21 -2.11e-04 4.75e-04 -3.14e-04 -1.25e-05 -1.20e-04 -1.45e-04 2.80e-04 2.11e-04 7.40e-22 6.02e-05 -1.48e-04 -2.22e-05 2.02e-04 -1.30e-04 -8.64e-04 -4.75e-04 -6.02e-05 -1.71e-21
    0.00e+00 -7.51e-06 1.05e-04 3.30e-05 -1.96e-04 6.96e-05 -8.85e-05 -4.81e-05 7.51e-06 2.13e-21 4.91e-04 5.07e-04 -2.56e-05 -6.38e-05 -2.17e-04 3.65e-06 -1.05e-04 -4.91e-04 -2.56e-21 -4.94e-05 -9.88e-05 -4.37e-04 -6.88e-05 1.91e-04 -3.30e-05 -5.07e-04 4.94e-05 3.53e-21 -1.27e-04 -4.33e-04 -2.76e-05 6.81e-05 1.96e-04 2.56e-05 9.88e-05 1.27e-04 7.55e-22 4.64e-06 8.25e-05 6.15e-06 -6.96e-05 6.38e-05 4.37e-04 4.33e-04 -4.64e-06 9.27e-22 -3.70e-04 1.88e-06 8.85e-05 2.17e-04 6.88e-05 2.76e-05 -8.25e-05 3.70e-04 1.57e-22 -1.62e-04 4.81e-05 -3.65e-06 -1.91e-04 -6.81e-05 -6.15e-06 -1.88e-06 1.62e-04 -2.25e-22
    0.00e+00 -2.49e-04 1.66e-05 -3.26e-05 6.47e-04 5.29e-05 8.96e-05 6.21e-04 2.49e-04 -1.02e-21 -8.02e-05 -3.59e-04 -4.61e-04 -1.53e-04 1.51e-04 -2.39e-04 -1.66e-05 8.02e-05 -3.03e-21 8.19e-05 -8.71e-05 5.21e-05 5.42e-06 -3.88e-04 3.26e-05 3.59e-04 -8.19e-05 3.06e-22 -3.56e-04 1.78e-04 7.91e-05 -4.79e-04 -6.47e-04 4.61e-04 8.71e-05 3.56e-04 -3.25e-21 7.99e-04 -1.13e-04 5.85e-05 -5.29e-05 1.53e-04 -5.21e-05 -1.78e-04 -7.99e-04 1.05e-22 7.48e-05 -2.53e-04 -8.96e-05 -1.51e-04 -5.42e-06 -7.91e-05 1.13e-04 -7.48e-05 8.61e-22 4.40e-04 -6.21e-04 2.39e-04 3.88e-04 4.79e-04 -5.85e-05 2.53e-04 -4.40e-04 -8.75e-22
    0.00e+00 -4.77e-05 1.34e-04 4.10e-05 -1.18e-05 -2.08e-04 1.09e-04 1.55e-04 4.77e-05 1.14e-21 8.75e-05 -1.39e-05 -2.32e-05 -2.70e-04 -5.80e-05 -2.35e-05 -1.34e-04 -8.75e-05 -1.25e-22 7.22e-05 -1.47e-04 -1.47e-04 1.75e-06 -2.86e-04 -4.10e-05 1.39e-05 -7.22e-05 -6.45e-22 1.01e-04 2.06e-04 1.12e-04 -3.27e-05 1.18e-05 2.32e-05 1.47e-04 -1.01e-04 5.94e-22 4.36e-04 1.87e-04 1.40e-03 2.08e-04 2.70e-04 1.47e-04 -2.06e-04 -4.36e-04 -2.03e-22 -2.40e-04 2.92e-04 -1.09e-04 5.80e-05 -1.75e-06 -1.12e-04 -1.87e-04 2.40e-04 3.16e-21 1.42e-04 -1.55e-04 2.35e-05 2.86e-04 3.27e-05 -1.40e-03 -2.92e-04 -1.42e-04 -9.78e-22
    0.00e+00 -8.68e-05 1.17e-04 1.32e-04 -4.88e-05 -3.46e-04 3.36e-05 1.13e-04 8.68e-05 1.40e-21 1.14e-04 4.53e-05 -7.81e-05 -3.24e-04 -8.14e-05 -2.58e-05 -1.17e-04 -1.14e-04 3.31e-22 7.35e-05 -1.97e-04 -1.63e-04 3.47e-05 -2.52e-04 -1.32e-04 -4.53e-05 -7.35e-05 -4.76e-22 1.06e-04 1.61e-04 8.30e-05 -3.72e-05 4.88e-05 7.81e-05 1.97e-04 -1.06e-04 3.11e-21 3.78e-04 -3.23e-05 1.29e-03 3.46e-04 3.24e-04 1.63e-04 -1.61e-04 -3.78e-04 3.67e-21 -3.07e-04 2.28e-04 -3.36e-05 8.14e-05 -3.47e-05 -8.30e-05 3.23e-05 3.07e-04 -5.79e-23 9.83e-05 -1.13e-04 2.58e-05 2.52e-04 3.72e-05 -1.29e-03 -2.28e-04 -9.83e-05 1.05e-22
    0.00e+00 5.24e-04 9.93e-05 8.84e-05 1.03e-04 -7.07e-04 4.13e-04 -6.09e-05 -5.24e-04 2.52e-22 1.97e-04 2.83e-04 -2.82e-04 -1.25e-04 -1.16e-04 2.77e-04 -9.93e-05 -1.97e-04 1.90e-22 -3.03e-06 -1.83e-04 1.40e-04 -8.72e-05 -3.01e-05 -8.84e-05 -2.83e-04 3.03e-06 8.60e-23 -1.38e-04 -5.94e-05 2.08e-04 -4.23e-05 -1.03e-04 2.82e-04 1.83e-04 1.38e-04 1.07e-21 -4.83e-04 2.38e-04 1.31e-04 7.07e-04 1.25e-04 -1.40e-04 5.94e-05 4.83e-04 2.04e-21 1.87e-04 3.26e-05 -4.13e-04 1.16e-04 8.72e-05 -2.08e-04 -2.38e-04 -1.87e-04 1.63e-21 -3.29e-04 6.09e-05 -2.77e-04 3.01e-05 4.23e-05 -1.31e-04 -3.26e-05 3.29e-04 5.68e-22
    0.00e+00 2.06e-04 2.97e-06 -1.55e-05 -1.78e-05 -2.40e-04 4.27e-04 -7.29e-05 -2.06e-04 8.27e-22 2.75e-04 6.22e-04 -2.33e-04 2.95e-04 -3.42e-04 1.29e-06 -2.97e-06 -2.75e-04 -1.01e-21 -8.27e-05 -4.70e-05 -8.90e-05 -1.31e-04 1.25e-04 1.55e-05 -6.22e-04 8.27e-05 1.48e-21 -1.12e-04 -3.02e-04 -1.01e-04 -2.56e-05 1.78e-05 2.33e-04 4.70e-05 1.12e-04 -1.45e-21 -8.34e-06 1.72e-04 -2.51e-05 2.40e-04 -2.95e-04 8.90e-05 3.02e-04 8.34e-06 -1.73e-21 -2.20e-04 -4.56e-06 -4.27e-04 3.42e-04 1.31e-04 1.01e-04 -1.72e-04 2.20e-04 -7.49e-22 -1.49e-04 7.29e-05 -1.29e-06 -1.25e-04 2.56e-05 2.51e-05 4.56e-06 1.49e-04 2.94e-22
    0.00e+00 -1.46e-04 2.48e-05 4.64e-05 6.62e-04 2.20e-05 8.81e-05 4.01e-05 1.46e-04 -1.09e-21 -1.15e-04 -1.23e-04 -1.90e-05 -1.66e-04 1.67e-05 -1.15e-04 -2.48e-05 1.15e-04 2.67e-23 3.53e-06 1.12e-04 -8.73e-05 1.34e-05 -2.49e-04 -4.64e-05 1.23e-04 -3.53e-06 4.48e-21 7.22e-05 1.22e-04 1.54e-04 -3.04e-05 -6.62e-04 1.90e-05 -1.12e-04 -7.22e-05 -1.45e-21 6.07e-04 3.25e-04 1.08e-03 -2.20e-05 1.66e-04 8.73e-05 -1.22e-04 -6.07e-04 -2.29e-21 -1.56e-04 -5.13e-06 -8.81e-05 -1.67e-05 -1.34e-05 -1.54e-04 -3.25e-04 1.56e-04 2.60e-21 2.09e-04 -4.01e-05 1.15e-04 2.49e-04 3.04e-05 -1.08e-03 5.13e-06 -2.09e-04 -5.56e-21
    0.00e+00 -6.41e-05 5.65e-05 1.30e-04 2.27e-04 -1.72e-04 2.75e-05 1.34e-04 6.41e-05 3.27e-23 8.54e-05 2.03e-05 7.42e-05 -2.65e-04 -1.46e-05 -6.21e-05 -5.65e-05 -8.54e-05 -5.37e-22 4.53e-05 -5.93e-05 -1.25e-04 -3.25e-05 -2.36e-04 -1.30e-04 -2.03e-05 -4.53e-05 9.18e-22 5.20e-05 8.97e-05 5.14e-05 -2.88e-06 -2.27e-04 -7.42e-05 5.93e-05 -5.20e-05 1.08e-21 4.35e-04 2.38e-04 1.36e-03 1.72e-04 2.65e-04 1.25e-04 -8.97e-05 -4.35e-04 -3.84e-22 -3.36e-04 6.05e-05 -2.75e-05 1.46e-05 3.25e-05 -5.14e-05 -2.38e-04 3.36e-04 6.27e-22 1.36e-04 -1.34e-04 6.21e-05 2.36e-04 2.88e-06 -1.36e-03 -6.05e-05 -1.36e-04 -1.02e-21
    0.00e+00 -5.97e-05 -6.40e-05 -7.95e-06 -1.50e-05 7.06e-06 4.75e-05 1.03e-04 5.97e-05 3.06e-22 -1.05e-04 -1.15e-04 -1.08e-04 -1.17e-05 -1.48e-05 -4.59e-06 6.40e-05 1.05e-04 4.36e-22 -3.00e-04 -8.90e-05 2.28e-04 -1.17e-04 -1.02e-03 7.95e-06 1.15e-04 3.00e-04 -2.17e-21 1.36e-04 -3.41e-05 -2.05e-04 -6.03e-04 1.50e-05 1.08e-04 8.90e-05 -1.36e-04 -4.11e-22 1.70e-04 -6.97e-05 -7.32e-04 -7.06e-06 1.17e-05 -2.28e-04 3.41e-05 -1.70e-04 -1.70e-22 3.19e-04 5.84e-04 -4.75e-05 1.48e-05 1.17e-04 2.05e-04 6.97e-05 -3.19e-04 -9.44e-22 5.31e-04 -1.03e-04 4.59e-06 1.02e-03 6.03e-04 7.32e-04 -5.84e-04 -5.31e-04 2.64e-22
    0.00e+00 6.61e-05 -1.41e-04 -1.45e-04 -6.34e-05 -5.33e-05 4.45e-05 -1.09e-04 -6.61e-05 2.82e-21 -3.15e-04 4.87e-05 -1.44e-04 1.60e-04 -4.28e-04 4.79e-05 1.41e-04 3.15e-04 1.51e-22 -8.20e-05 -1.85e-04 7.57e-05 -2.22e-04 5.18e-05 1.45e-04 -4.87e-05 8.20e-05 -5.50e-22 4.91e-05 -3.17e-05 1.90e-04 -8.32e-05 6.34e-05 1.44e-04 1.85e-04 -4.91e-05 -4.19e-21 -3.52e-05 1.47e-04 -7.59e-05 5.33e-05 -1.60e-04 -7.57e-05 3.17e-05 3.52e-05 -2.77e-21 1.73e-04 -1.56e-04 -4.45e-05 4.28e-04 2.22e-04 -1.90e-04 -1.47e-04 -1.73e-04 -4.69e-22 2.00e-04 1.09e-04 -4.79e-05 -5.18e-05 8.32e-05 7.59e-05 1.56e-04 -2.00e-04 2.12e-21
    0.00e+00 1.15e-04 -4.83e-05 -7.69e-05 -9.05e-06 -1.63e-04 8.26e-05 1.72e-04 -1.15e-04 4.11e-22 -1.92e-04 2.85e-04 -2.42e-04 4.54e-04 -3.56e-04 -3.29e-04 4.83e-05 1.92e-04 8.35e-22 1.31e-04 -2.90e-04 -2.14e-04 -5.53e-04 3.78e-04 7.69e-05 -2.85e-04 -1.31e-04 -1.31e-21 -2.96e-05 8.68e-05 -5.73e-05 -1.95e-04 9.05e-06 2.42e-04 2.90e-04 2.96e-05 1.27e-23 3.08e-05 2.35e-04 -4.06e-06 1.63e-04 -4.54e-04 2.14e-04 -8.68e-05 -3.08e-05 8.58e-23 -1.68e-04 -5.60e-05 -8.26e-05 3.56e-04 5.53e-04 5.73e-05 -2.35e-04 1.68e-04 -3.17e-22 8.63e-05 -1.72e-04 3.29e-04 -3.78e-04 1.95e-04 4.06e-06 5.60e-05 -8.63e-05 -2.38e-21
    0.00e+00 1.54e-06 1.34e-04 7.72e-05 -2.50e-04 1.19e-04 -2.19e-04 -5.06e-05 -1.54e-06 1.86e-21 4.55e-04 4.49e-04 1.86e-05 -2.46e-04 -2.38e-05 -3.33e-05 -1.34e-04 -4.55e-04 -3.27e-21 -3.39e-05 -9.58e-05 -4.72e-04 -5.62e-05 2.05e-04 -7.72e-05 -4.49e-04 3.39e-05 -4.62e-22 -1.24e-04 -5.21e-04 4.85e-06 1.54e-04 2.50e-04 -1.86e-05 9.58e-05 1.24e-04 1.34e-21 2.02e-05 6.42e-05 -1.73e-05 -1.19e-04 2.46e-04 4.72e-04 5.21e-04 -2.02e-05 1.13e-21 -3.54e-04 -2.29e-05 2.19e-04 2.38e-05 5.62e-05 -4.85e-06 -6.42e-05 3.54e-04 -6.79e-22 -1.38e-04 5.06e-05 3.33e-05 -2.05e-04 -1.54e-04 1.73e-05 2.29e-05 1.38e-04 3.71e-22
    0.00e+00 5.47e-05 1.84e-04 -1.10e-05 -9.05e-06 2.01e-04 1.17e-04 3.46e-04 -5.47e-05 1.74e-21 6.56e-05 1.06e-06 2.72e-05 -1.17e-04 1.14e-04 3.69e-06 -1.84e-04 -6.56e-05 -1.81e-21 2.88e-05 -6.77e-05 -1.03e-05 -8.17e-05 -2.01e-04 1.10e-05 -1.06e-06 -2.88e-05 5.35e-21 1.06e-04 1.21e-04 1.87e-04 1.98e-04 9.05e-06 -2.72e-05 6.77e-05 -1.06e-04 -4.41e-22 4.92e-04 4.27e-04 1.11e-03 -2.01e-04 1.17e-04 1.03e-05 -1.21e-04 -4.92e-04 1.24e-21 -8.09e-07 1.66e-04 -1.17e-04 -1.14e-04 8.17e-05 -1.87e-04 -4.27e-04 8.09e-07 1.04e-21 9.70e-05 -3.46e-04 -3.69e-06 2.01e-04 -1.98e-04 -1.11e-03 -1.66e-04 -9.70e-05 -8.72e-22
    0.00e+00 1.77e-04 -5.75e-05 -4.30e-05 8.19e-05 -1.88e-04 1.61e-04 2.46e-04 -1.77e-04 -2.45e-21 -1.99e-04 3.49e-04 -2.60e-04 4.21e-04 -4.22e-04 -3.67e-04 5.75e-05 1.99e-04 1.36e-21 1.45e-04 -2.21e-04 -2.81e-04 -4.42e-04 4.52e-04 4.30e-05 -3.49e-04 -1.45e-04 -1.48e-21 -2.23e-05 1.40e-04 -1.24e-04 -2.28e-04 -8.19e-05 2.60e-04 2.21e-04 2.23e-05 -1.65e-22 2.07e-05 2.05e-04 3.77e-05 1.88e-04 -4.21e-04 2.81e-04 -1.40e-04 -2.07e-05 -4.60e-22 -2.13e-04 1.75e-05 -1.61e-04 4.22e-04 4.42e-04 1.24e-04 -2.05e-04 2.13e-04 -6.45e-22 -4.56e-05 -2.46e-04 3.67e-04 -4.52e-04 2.28e-04 -3.77e-05 -1.75e-05 4.56e-05 1.96e-22
    0.00e+00 2.80e-04 -8.69e-05 1.41e-05 5.81e-05 -1.69e-04 2.80e-04 1.39e-04 -2.80e-04 6.60e-23 -1.46e-04 5.46e-04 -4.30e-04 1.31e-04 -3.50e-04 -2.18e-04 8.69e-05 1.46e-04 -3.14e-21 1.00e-04 -3.01e-05 4.26e-05 -1.83e-04 2.26e-04 -1.41e-05 -5.46e-04 -1.00e-04 -1.91e-21 -1.14e-04 -1.05e-04 -1.86e-04 -2.17e-04 -5.81e-05 4.30e-04 3.01e-05 1.14e-04 -2.23e-21 4.59e-05 2.11e-04 4.61e-05 1.69e-04 -1.31e-04 -4.26e-05 1.05e-04 -4.59e-05 -1.35e-21 -4.21e-04 1.03e-04 -2.80e-04 3.50e-04 1.83e-04 1.86e-04 -2.11e-04 4.21e-04 -1.80e-21 -2.11e-04 -1.39e-04 2.18e-04 -2.26e-04 2.17e-04 -4.61e-05 -1.03e-04 2.11e-04 2.52e-22
    0.00e+00 1.84e-04 -1.08e-04 -1.49e-05 7.08e-05 -2.10e-04 2.30e-04 2.03e-04 -1.84e-04 -1.66e-22 -1.95e-04 4.80e-04 -3.92e-04 3.07e-04 -4.57e-04 -3.35e-04 1.08e-04 1.95e-04 -2.30e-21 1.04e-04 -1.16e-04 -1.61e-04 -2.47e-04 3.13e-04 1.49e-05 -4.80e-04 -1.04e-04 -2.23e-21 -5.09e-05 2.82e-05 -2.14e-04 -2.60e-04 -7.08e-05 3.92e-04 1.16e-04 5.09e-05 -3.96e-22 4.48e-05 1.95e-04 3.07e-05 2.10e-04 -3.07e-04 1.61e-04 -2.82e-05 -4.48e-05 4.72e-21 -3.43e-04 8.08e-05 -2.30e-04 4.57e-04 2.47e-04 2.14e-04 -1.95e-04 3.43e-04 -1.37e-21 -1.63e-04 -2.03e-04 3.35e-04 -3.13e-04 2.60e-04 -3.07e-05 -8.08e-05 1.63e-04 -4.30e-22
    0.00e+00 2.57e-04 -1.30e-06 7.12e-06 5.94e-05 -2.64e-04 4.71e-04 -6.90e-05 -2.57e-04 -1.06e-21 1.64e-04 6.04e-04 -2.83e-04 1.94e-04 -2.96e-04 -2.46e-05 1.30e-06 -1.64e-04 4.28e-22 -3.89e-05 -5.40e-05 -1.33e-05 -1.85e-04 8.77e-05 -7.12e-06 -6.04e-04 3.89e-05 -4.28e-21 -1.36e-04 -2.91e-04 -1.30e-04 -9.27e-05 -5.94e-05 2.83e-04 5.40e-05 1.36e-04 -9.23e-23 -1.80e-05 2.02e-04 -8.87e-06 2.64e-04 -1.94e-04 1.33e-05 2.91e-04 1.80e-05 -7.11e-22 -2.28e-04 3.42e-06 -4.71e-04 2.96e-04 1.85e-04 1.30e-04 -2.02e-04 2.28e-04 -4.29e-21 -1.57e-04 6.90e-05 2.46e-05 -8.77e-05 9.27e-05 8.87e-06 -3.42e-06 1.57e-04 -2.55e-23
    0.00e+00 2.00e-04 -1.55e-04 -1.18e-04 -1.23e-04 -1.83e-04 3.82e-04 -2.21e-04 -2.00e-04 -1.42e-21 -2.03e-04 -1.89e-06 -1.59e-04 6.19e-05 -6.48e-04 2.59e-04 1.55e-04 2.03e-04 1.42e-21 -3.21e-04 -3.62e-04 1.85e-04 1.05e-04 -1.64e-04 1.18e-04 1.89e-06 3.21e-04 1.61e-22 1.41e-04 -7.62e-05 3.09e-04 -1.94e-04 1.23e-04 1.59e-04 3.62e-04 -1.41e-04 9.87e-22 -1.00e-04 1.14e-04 -2.43e-04 1.83e-04 -6.19e-05 -1.85e-04 7.62e-05 1.00e-04 1.02e-22 6.09e-04 -1.33e-04 -3.82e-04 6.48e-04 -1.05e-04 -3.09e-04 -1.14e-04 -6.09e-04 4.58e-22 -5.20e-04 2.21e-04 -2.59e-04 1.64e-04 1.94e-04 2.43e-04 1.33e-04 5.20e-04 -2.05e-21
    0.00e+00 -4.50e-06 -1.41e-04 -8.85e-05 -7.89e-05 -1.78e-04 -7.45e-05 1.42e-04 4.50e-06 -1.09e-22 -1.20e-04 2.38e-04 -9.49e-05 3.74e-04 -2.46e-04 -2.28e-04 1.41e-04 1.20e-04 2.03e-21 9.03e-05 -3.58e-04 -1.89e-04 -6.55e-04 2.15e-04 8.85e-05 -2.38e-04 -9.03e-05 -1.03e-21 4.40e-05 1.17e-04 8.44e-05 -2.05e-04 7.89e-05 9.49e-05 3.58e-04 -4.40e-05 -1.64e-21 1.83e-05 2.10e-04 6.85e-07 1.78e-04 -3.74e-04 1.89e-04 -1.17e-04 -1.83e-05 2.81e-23 6.58e-05 -1.54e-04 7.45e-05 2.46e-04 6.55e-04 -8.44e-05 -2.10e-04 -6.58e-05 9.91e-22 1.53e-04 -1.42e-04 2.28e-04 -2.15e-04 2.05e-04 -6.85e-07 1.54e-04 -1.53e-04 1.18e-21
    0.00e+00 4.75e-04 1.14e-04 6.29e-05 3.94e-05 -6.63e-04 4.84e-04 -6.33e-05 -4.75e-04 4.81e-22 1.90e-04 4.05e-04 -3.13e-04 2.14e-04 2.64e-05 2.49e-04 -1.14e-04 -1.90e-04 -3.47e-21 -2.23e-05 -1.81e-04 2.07e-04 -2.03e-04 -2.34e-05 -6.29e-05 -4.05e-04 2.23e-05 -2.02e-21 -1.93e-04 -9.19e-05 2.67e-04 -7.62e-05 -3.94e-05 3.13e-04 1.81e-04 1.93e-04 2.14e-21 -3.12e-04 3.20e-04 1.25e-04 6.63e-04 -2.14e-04 -2.07e-04 9.19e-05 3.12e-04 -3.30e-22 1.33e-04 1.54e-04 -4.84e-04 -2.64e-05 2.03e-04 -2.67e-04 -3.20e-04 -1.33e-04 2.47e-21 -3.02e-04 6.33e-05 -2.49e-04 2.34e-05 7.62e-05 -1.25e-04 -1.54e-04 3.02e-04 -1.27e-21
    0.00e+00 7.93e-06 -2.95e-04 -8.58e-05 -1.22e-04 -1.74e-04 -6.42e-05 -5.73e-05 -7.93e-06 -8.04e-22 -1.29e-04 2.01e-04 -9.17e-05 3.02e-04 -3.28e-04 4.66e-05 2.95e-04 1.29e-04 -1.44e-21 7.17e-05 -4.25e-04 -9.21e-05 -4.45e-04 -3.80e-05 8.58e-05 -2.01e-04 -7.17e-05 -1.47e-21 7.40e-05 4.91e-05 1.32e-04 -2.06e-04 1.22e-04 9.17e-05 4.25e-04 -7.40e-05 4.08e-22 -3.88e-05 1.94e-04 7.70e-06 1.74e-04 -3.02e-04 9.21e-05 -4.91e-05 3.88e-05 2.21e-21 1.67e-04 -1.58e-04 6.42e-05 3.28e-04 4.45e-04 -1.32e-04 -1.94e-04 -1.67e-04 -8.56e-22 2.02e-04 5.73e-05 -4.66e-05 3.80e-05 2.06e-04 -7.70e-06 1.58e-04 -2.02e-04 -8.45e-23
    0.00e+00 1.22e-05 -9.90e-05 -1.02e-04 -4.22e-05 -2.07e-04 -4.69e-05 2.37e-04 -1.22e-05 -8.32e-22 -1.57e-04 2.81e-04 -1.41e-04 4.54e-04 -2.64e-04 -3.82e-04 9.90e-05 1.57e-04 7.88e-22 9.68e-05 -3.39e-04 -2.81e-04 -6.67e-04 3.76e-04 1.02e-04 -2.81e-04 -9.68e-05 4.92e-22 3.04e-05 1.58e-04 -1.00e-05 -2.09e-04 4.22e-05 1.41e-04 3.39e-04 -3.04e-05 1.07e-21 1.51e-05 2.22e-04 -2.92e-05 2.07e-04 -4.54e-04 2.81e-04 -1.58e-04 -1.51e-05 -1.27e-21 4.27e-05 -1.35e-04 4.69e-05 2.64e-04 6.67e-04 1.00e-05 -2.22e-04 -4.27e-05 -9.25e-22 9.36e-05 -2.37e-04 3.82e-04 -3.76e-04 2.09e-04 2.92e-05 1.35e-04 -9.36e-05 8.35e-22
    0.00e+00 -1.96e-04 -9.34e-06 -1.62e-04 5.78e-04 6.69e-05 1.64e-04 5.93e-04 1.96e-04 -6.50e-21 -1.43e-04 -3.35e-04 -3.09e-04 -9.04e-05 2.13e-04 -1.24e-04 9.34e-06 1.43e-04 -2.57e-21 7.99e-05 -2.61e-04 -1.38e-05 2.36e-05 -3.21e-04 1.62e-04 3.35e-04 -7.99e-05 7.69e-22 -3.66e-04 1.81e-04 1.85e-04 -2.26e-04 -5.78e-04 3.09e-04 2.61e-04 3.66e-04 -1.17e-21 4.65e-04 -2.22e-04 1.53e-04 -6.69e-05 9.04e-05 1.38e-05 -1.81e-04 -4.65e-04 -2.18e-22 7.79e-05 -2.32e-04 -1.64e-04 -2.13e-04 -2.36e-05 -1.85e-04 2.22e-04 -7.79e-05 -8.99e-22 3.93e-04 -5.93e-04 1.24e-04 3.21e-04 2.26e-04 -1.53e-04 2.32e-04 -3.93e-04 -4.14e-22
    0.00e+00 2.09e-04 -1.04e-04 -1.12e-04 -3.44e-06 -7.41e-05 9.20e-05 2.70e-04 -2.09e-04 -1.07e-21 6.77e-05 4.56e-05 1.64e-04 5.50e-04 5.27e-04 -4.68e-04 1.04e-04 -6.77e-05 -6.30e-22 -2.62e-04 3.29e-05 -4.78e-04 -4.73e-05 6.79e-04 1.12e-04 -4.56e-05 2.62e-04 -3.96e-21 -4.43e-05 2.13e-04 1.58e-04 -5.53e-05 3.44e-06 -1.64e-04 -3.29e-05 4.43e-05 -1.99e-22 -5.92e-05 -3.13e-05 -1.51e-04 7.41e-05 -5.50e-04 4.78e-04 -2.13e-04 5.92e-05 1.92e-22 1.45e-04 1.53e-06 -9.20e-05 -5.27e-04 4.73e-05 -1.58e-04 3.13e-05 -1.45e-04 -3.44e-22 -1.63e-04 -2.70e-04 4.68e-04 -6.79e-04 5.53e-05 1.51e-04 -1.53e-06 1.63e-04 -1.89e-21
    0.00e+00 -1.53e-04 9.23e-05 1.20e-04 6.71e-04 1.15e-04 7.04e-05 1.34e-04 1.53e-04 4.99e-22 -3.12e-05 -2.33e-04 -2.83e-04 -2.37e-04 -7.50e-06 -1.74e-04 -9.23e-05 3.12e-05 6.84e-21 5.43e-05 2.67e-04 9.33e-07 -7.32e-06 -3.37e-04 -1.20e-04 2.33e-04 -5.43e-05 1.32e-21 1.01e-04 1.67e-04 9.96e-05 -1.80e-04 -6.71e-04 2.83e-04 -2.67e-04 -1.01e-04 -1.29e-21 1.07e-03 3.11e-04 3.96e-04 -1.15e-04 2.37e-04 -9.33e-07 -1.67e-04 -1.07e-03 9.43e-22 4.20e-06 3.58e-05 -7.04e-05 7.50e-06 7.32e-06 -9.96e-05 -3.11e-04 -4.20e-06 4.19e-21 3.33e-04 -1.34e-04 1.74e-04 3.37e-04 1.80e-04 -3.96e-04 -3.58e-05 -3.33e-04 -1.58e-21
    0.00e+00 9.81e-05 -1.20e-04 -6.11e-06 1.72e-05 -1.91e-04 2.12e-04 2.21e-04 -9.81e-05 1.96e-21 -1.50e-04 3.89e-04 -2.92e-04 3.62e-04 -4.34e-04 -4.33e-04 1.20e-04 1.50e-04 -1.32e-21 1.13e-04 -1.06e-04 -3.43e-04 -3.50e-04 4.16e-04 6.11e-06 -3.89e-04 -1.13e-04 1.11e-21 -6.71e-05 1.37e-04 -1.91e-04 -2.25e-04 -1.72e-05 2.92e-04 1.06e-04 6.71e-05 3.90e-22 2.04e-05 2.01e-04 1.35e-05 1.91e-04 -3.62e-04 3.43e-04 -1.37e-04 -2.04e-05 -1.85e-21 -2.70e-04 3.80e-05 -2.12e-04 4.34e-04 3.50e-04 1.91e-04 -2.01e-04 2.70e-04 -1.88e-21 -1.54e-04 -2.21e-04 4.33e-04 -4.16e-04 2.25e-04 -1.35e-05 -3.80e-05 1.54e-04 -5.37e-22
    0.00e+00 5.60e-05 -1.14e-04 -7.47e-05 1.77e-06 -2.08e-04 7.81e-05 2.80e-04 -5.60e-05 -4.97e-21 -1.67e-04 3.43e-04 -2.12e-04 5.03e-04 -3.56e-04 -4.91e-04 1.14e-04 1.67e-04 -8.10e-21 1.19e-04 -2.51e-04 -3.82e-04 -5.86e-04 4.50e-04 7.47e-05 -3.43e-04 -1.19e-04 -7.72e-22 2.66e-07 2.03e-04 -8.90e-05 -2.21e-04 -1.77e-06 2.12e-04 2.51e-04 -2.66e-07 2.95e-22 2.13e-05 2.35e-04 -2.64e-05 2.08e-04 -5.03e-04 3.82e-04 -2.03e-04 -2.13e-05 -3.77e-22 -9.82e-05 -5.61e-05 -7.81e-05 3.56e-04 5.86e-04 8.90e-05 -2.35e-04 9.82e-05 8.47e-22 -2.63e-05 -2.80e-04 4.91e-04 -4.50e-04 2.21e-04 2.64e-05 5.61e-05 2.63e-05 -9.76e-22
    0.00e+00 4.26e-04 -6.08e-05 -3.31e-05 1.63e-04 -1.59e-04 4.01e-04 6.44e-05 -4.26e-04 7.81e-22 8.91e-05 3.81e-04 -2.34e-04 8.17e-05 -6.06e-05 1.44e-05 6.08e-05 -8.91e-05 3.61e-21 4.90e-05 -7.98e-05 8.44e-05 -2.12e-04 6.17e-05 3.31e-05 -3.81e-04 -4.90e-05 1.30e-21 -7.73e-05 -1.29e-04 -1.08e-04 -7.08e-05 -1.63e-04 2.34e-04 7.98e-05 7.73e-05 1.50e-21 -7.22e-05 1.89e-04 2.35e-05 1.59e-04 -8.17e-05 -8.44e-05 1.29e-04 7.22e-05 2.69e-21 -2.68e-04 7.95e-05 -4.01e-04 6.06e-05 2.12e-04 1.08e-04 -1.89e-04 2.68e-04 2.31e-21 -1.23e-04 -6.44e-05 -1.44e-05 -6.17e-05 7.08e-05 -2.35e-05 -7.95e-05 1.23e-04 -2.43e-21
    0.00e+00 1.11e-04 -1.81e-04 -7.11e-05 -3.76e-05 -1.64e-04 6.68e-05 2.84e-04 -1.11e-04 2.96e-22 -1.56e-05 1.84e-04 3.01e-05 5.28e-04 2.07e-04 -5.67e-04 1.81e-04 1.56e-05 1.19e-21 -1.30e-04 -6.83e-05 -5.42e-04 -3.18e-04 5.59e-04 7.11e-05 -1.84e-04 1.30e-04 -4.16e-21 -3.64e-05 2.33e-04 7.57e-05 -1.83e-04 3.76e-05 -3.01e-05 6.83e-05 3.64e-05 -6.88e-22 -4.56e-05 7.96e-05 -1.26e-04 1.64e-04 -5.28e-04 5.42e-04 -2.33e-04 4.56e-05 5.30e-21 1.17e-04 -5.92e-06 -6.68e-05 -2.07e-04 3.18e-04 -7.57e-05 -7.96e-05 -1.17e-04 9.28e-22 -1.91e-04 -2.84e-04 5.67e-04 -5.59e-04 1.83e-04 1.26e-04 5.92e-06 1.91e-04 9.20e-22
    0.00e+00 1.26e-04 -4.48e-05 9.96e-05 1.21e-04 -3.24e-04 1.04e-04 6.68e-05 -1.26e-04 -1.25e-21 6.97e-05 3.23e-05 1.05e-05 -2.03e-04 -4.09e-04 6.36e-05 4.48e-05 -6.97e-05 3.58e-22 -5.12e-05 -6.63e-05 3.89e-05 1.52e-04 -6.13e-05 -9.96e-05 -3.23e-05 5.12e-05 -1.36e-21 -3.85e-05 -1.46e-04 3.15e-05 -6.81e-05 -1.21e-04 -1.05e-05 6.63e-05 3.85e-05 2.10e-21 -2.72e-04 4.44e-05 1.01e-04 3.24e-04 2.03e-04 -3.89e-05 1.46e-04 2.72e-04 2.31e-21 3.42e-04 9.59e-05 -1.04e-04 4.09e-04 -1.52e-04 -3.15e-05 -4.44e-05 -3.42e-04 1.14e-21 -5.40e-04 -6.68e-05 -6.36e-05 6.13e-05 6.81e-05 -1.01e-04 -9.59e-05 5.40e-04 -5.07e-21
    0.00e+00 4.17e-04 6.64e-05 1.24e-04 1.48e-04 -6.43e-04 2.59e-04 1.10e-05 -4.17e-04 -6.06e-22 1.66e-04 1.77e-04 -1.80e-04 -4.03e-04 -3.43e-04 2.36e-04 -6.64e-05 -1.66e-04 -7.18e-22 -6.43e-05 -1.94e-04 1.77e-04 -7.08e-05 -5.70e-05 -1.24e-04 -1.77e-04 6.43e-05 -3.39e-22 -5.62e-05 -9.87e-05 2.40e-04 -9.70e-05 -1.48e-04 1.80e-04 1.94e-04 5.62e-05 3.44e-22 -5.20e-04 7.48e-05 1.10e-04 6.43e-04 4.03e-04 -1.77e-04 9.87e-05 5.20e-04 -5.57e-22 1.65e-04 1.17e-05 -2.59e-04 3.43e-04 7.08e-05 -2.40e-04 -7.48e-05 -1.65e-04 3.41e-21 -4.45e-04 -1.10e-05 -2.36e-04 5.70e-05 9.70e-05 -1.10e-04 -1.17e-05 4.45e-04 -5.76e-22
    0.00e+00 2.12e-04 -1.18e-04 -1.09e-04 -1.21e-04 -2.48e-04 2.84e-04 -2.20e-04 -2.12e-04 2.89e-21 -1.38e-04 4.68e-05 -1.98e-04 2.10e-06 -6.30e-04 2.52e-04 1.18e-04 1.38e-04 3.52e-23 -2.47e-04 -3.05e-04 1.04e-04 -4.83e-05 -3.14e-04 1.09e-04 -4.68e-05 2.47e-04 6.47e-22 1.40e-04 1.24e-05 1.46e-04 -2.04e-04 1.21e-04 1.98e-04 3.05e-04 -1.40e-04 -2.77e-21 -5.77e-05 -4.73e-05 -2.59e-04 2.48e-04 -2.10e-06 -1.04e-04 -1.24e-05 5.77e-05 -9.45e-22 7.21e-04 -1.43e-04 -2.84e-04 6.30e-04 4.83e-05 -1.46e-04 4.73e-05 -7.21e-04 2.60e-21 -4.81e-04 2.20e-04 -2.52e-04 3.14e-04 2.04e-04 2.59e-04 1.43e-04 4.81e-04 -2.41e-21
    0.00e+00 -6.08e-05 1.11e-04 1.31e-04 1.47e-04 -8.33e-04 -2.28e-05 2.10e-04 6.08e-05 -1.94e-21 9.71e-05 7.48e-05 2.08e-04 -6.19e-04 -8.10e-05 3.51e-06 -1.11e-04 -9.71e-05 9.96e-22 3.16e-05 -9.15e-05 9.14e-06 2.64e-05 -5.28e-05 -1.31e-04 -7.48e-05 -3.16e-05 -5.15e-22 -6.94e-05 1.76e-04 9.71e-05 1.54e-05 -1.47e-04 -2.08e-04 9.15e-05 6.94e-05 -6.73e-22 -7.42e-04 8.39e-05 4.44e-04 8.33e-04 6.19e-04 -9.14e-06 -1.76e-04 7.42e-04 -1.49e-21 -5.17e-04 -6.57e-04 2.28e-05 8.10e-05 -2.64e-05 -9.71e-05 -8.39e-05 5.17e-04 3.58e-22 -1.11e-04 -2.10e-04 -3.51e-06 5.28e-05 -1.54e-05 -4.44e-04 6.57e-04 1.11e-04 -3.24e-22
    0.00e+00 -2.09e-04 -3.18e-06 -1.36e-04 5.17e-04 6.06e-05 2.09e-04 5.47e-04 2.09e-04 -1.59e-21 -1.40e-04 -3.07e-04 -3.51e-04 -1.41e-04 1.34e-04 -1.57e-04 3.18e-06 1.40e-04 6.68e-21 4.55e-05 -1.34e-04 6.08e-05 -1.82e-05 -4.47e-04 1.36e-04 3.07e-04 -4.55e-05 5.63e-23 -3.09e-04 1.51e-04 7.72e-05 -5.04e-04 -5.17e-04 3.51e-04 1.34e-04 3.09e-04 -2.79e-22 5.15e-04 -1.91e-04 2.02e-04 -6.06e-05 1.41e-04 -6.08e-05 -1.51e-04 -5.15e-04 -1.09e-21 7.14e-05 -1.04e-04 -2.09e-04 -1.34e-04 1.82e-05 -7.72e-05 1.91e-04 -7.14e-05 -7.80e-21 5.35e-04 -5.47e-04 1.57e-04 4.47e-04 5.04e-04 -2.02e-04 1.04e-04 -5.35e-04 -3.31e-21
    0.00e+00 -1.71e-04 3.04e-05 1.09e-04 6.97e-04 7.78e-05 8.99e-05 3.69e-04 1.71e-04 -3.98e-22 -6.36e-05 -2.12e-04 -4.22e-04 -2.14e-04 2.54e-05 -2.22e-04 -3.04e-05 6.36e-05 1.44e-23 6.16e-05 1.42e-04 7.70e-05 -4.28e-05 -3.19e-04 -1.09e-04 2.12e-04 -6.16e-05 -1.09e-21 -8.47e-05 2.03e-04 2.70e-05 -4.00e-04 -6.97e-04 4.22e-04 -1.42e-04 8.47e-05 -1.44e-21 9.92e-04 1.60e-04 1.81e-04 -7.78e-05 2.14e-04 -7.70e-05 -2.03e-04 -9.92e-04 4.71e-22 4.38e-05 -1.08e-04 -8.99e-05 -2.54e-05 4.28e-05 -2.70e-05 -1.60e-04 -4.38e-05 -2.01e-21 3.99e-04 -3.69e-04 2.22e-04 3.19e-04 4.00e-04 -1.81e-04 1.08e-04 -3.99e-04 4.07e-21
    0.00e+00 2.93e-04 -4.37e-05 1.56e-05 8.14e-05 -1.95e-04 4.06e-04 -1.92e-05 -2.93e-04 1.63e-22 2.82e-05 6.17e-04 -3.57e-04 6.25e-05 -1.78e-04 -6.61e-05 4.37e-05 -2.82e-05 1.09e-21 -8.48e-06 -3.53e-05 1.21e-04 -2.18e-04 1.09e-04 -1.56e-05 -6.17e-04 8.48e-06 -6.27e-23 -1.38e-04 -2.65e-04 -1.21e-04 -1.53e-04 -8.14e-05 3.57e-04 3.53e-05 1.38e-04 3.17e-22 -2.26e-05 2.34e-04 2.33e-06 1.95e-04 -6.25e-05 -1.21e-04 2.65e-04 2.26e-05 -3.38e-21 -2.98e-04 5.29e-05 -4.06e-04 1.78e-04 2.18e-04 1.21e-04 -2.34e-04 2.98e-04 5.37e-22 -2.07e-04 1.92e-05 6.61e-05 -1.09e-04 1.53e-04 -2.33e-06 -5.29e-05 2.07e-04 -1.29e-21
    0.00e+00 2.86e-04 6.60e-04 3.30e-04 -2.90e-04 -6.10e-05 -1.39e-04 -5.32e-05 -2.86e-04 2.40e-21 7.47e-05 1.65e-04 1.29e-04 -1.56e-04 1.68e-04 -4.18e-05 -6.60e-04 -7.47e-05 8.54e-22 1.77e-04 -1.10e-04 -6.60e-04 5.58e-05 1.18e-04 -3.30e-04 -1.65e-04 -1.77e-04 -1.15e-21 -3.53e-05 -3.09e-04 2.45e-05 1.64e-04 2.90e-04 -1.29e-04 1.10e-04 3.53e-05 2.10e-21 5.88e-06 -1.56e-04 9.86e-05 6.10e-05 1.56e-04 6.60e-04 3.09e-04 -5.88e-06 -2.60e-23 -9.54e-05 -3.11e-05 1.39e-04 -1.68e-04 -5.58e-05 -2.45e-05 1.56e-04 9.54e-05 1.43e-22 6.35e-05 5.32e-05 4.18e-05 -1.18e-04 -1.64e-04 -9.86e-05 3.11e-05 -6.35e-05 -1.33e-22
    0.00e+00 2.58e-05 2.02e-04 1.68e-04 -3.08e-04 1.20e-04 -3.09e-04 -3.05e-05 -2.58e-05 2.03e-21 4.64e-04 3.82e-04 5.42e-05 -3.50e-04 1.09e-04 -2.09e-05 -2.02e-04 -4.64e-04 -4.81e-23 1.44e-06 -1.05e-04 -4.90e-04 -1.64e-05 1.80e-04 -1.68e-04 -3.82e-04 -1.44e-06 -1.29e-21 -1.36e-04 -5.58e-04 6.02e-05 1.89e-04 3.08e-04 -5.42e-05 1.05e-04 1.36e-04 1.75e-21 -1.79e-05 2.35e-05 -3.15e-06 -1.20e-04 3.50e-04 4.90e-04 5.58e-04 1.79e-05 1.23e-22 -3.73e-04 -6.19e-05 3.09e-04 -1.09e-04 1.64e-05 -6.02e-05 -2.35e-05 3.73e-04 -9.86e-22 -1.04e-04 3.05e-05 2.09e-05 -1.80e-04 -1.89e-04 3.15e-06 6.19e-05 1.04e-04 1.75e-21
    0.00e+00 4.45e-04 8.18e-05 1.49e-04 -1.52e-04 1.15e-04 7.93e-05 -6.36e-05 -4.45e-04 4.48e-22 3.32e-04 -2.88e-04 3.09e-04 7.31e-04 6.73e-04 1.84e-04 -8.18e-05 -3.32e-04 3.75e-23 -2.06e-04 1.56e-04 -9.25e-06 4.95e-05 2.20e-04 -1.49e-04 2.88e-04 2.06e-04 5.41e-22 -4.36e-05 3.12e-04 2.48e-04 1.21e-04 1.52e-04 -3.09e-04 -1.56e-04 4.36e-05 -1.45e-21 -1.53e-04 -1.16e-04 -7.03e-05 -1.15e-04 -7.31e-04 9.25e-06 -3.12e-04 1.53e-04 1.14e-22 -1.66e-04 1.27e-04 -7.93e-05 -6.73e-04 -4.95e-05 -2.48e-04 1.16e-04 1.66e-04 1.75e-21 5.43e-05 6.36e-05 -1.84e-04 -2.20e-04 -1.21e-04 7.03e-05 -1.27e-04 -5.43e-05 6.02e-22
    0.00e+00 -1.67e-04 8.14e-05 1.02e-04 6.28e-04 1.05e-04 6.35e-05 1.83e-04 1.67e-04 7.97e-22 -6.06e-05 -2.33e-04 -3.19e-04 -2.02e-04 9.26e-06 -2.01e-04 -8.14e-05 6.06e-05 -2.97e-23 3.24e-05 1.98e-04 -3.59e-06 -1.56e-05 -3.60e-04 -1.02e-04 2.33e-04 -3.24e-05 -3.72e-22 2.67e-05 1.62e-04 4.39e-05 -2.52e-04 -6.28e-04 3.19e-04 -1.98e-04 -2.67e-05 -4.97e-22 1.01e-03 2.13e-04 2.54e-04 -1.05e-04 2.02e-04 3.59e-06 -1.62e-04 -1.01e-03 -3.08e-22 3.91e-05 3.88e-05 -6.35e-05 -9.26e-06 1.56e-05 -4.39e-05 -2.13e-04 -3.91e-05 2.95e-21 3.71e-04 -1.83e-04 2.01e-04 3.60e-04 2.52e-04 -2.54e-04 -3.88e-05 -3.71e-04 -6.82e-21
    0.00e+00 9.41e-05 -1.98e-04 -9.84e-05 -3.56e-05 -1.67e-04 2.80e-05 2.15e-04 -9.41e-05 -2.61e-21 -5.72e-05 2.61e-04 -5.22e-05 5.49e-04 3.09e-05 -6.04e-04 1.98e-04 5.72e-05 -1.62e-21 -3.09e-05 -1.35e-04 -5.33e-04 -4.28e-04 5.21e-04 9.84e-05 -2.61e-04 3.09e-05 6.98e-21 -4.04e-05 1.77e-04 -8.24e-06 -1.73e-04 3.56e-05 5.22e-05 1.35e-04 4.04e-05 9.00e-23 -3.99e-05 1.53e-04 -1.09e-04 1.67e-04 -5.49e-04 5.33e-04 -1.77e-04 3.99e-05 4.91e-22 6.90e-05 -2.78e-05 -2.80e-05 -3.09e-05 4.28e-04 8.24e-06 -1.53e-04 -6.90e-05 1.02e-21 -1.51e-04 -2.15e-04 6.04e-04 -5.21e-04 1.73e-04 1.09e-04 2.78e-05 1.51e-04 2.87e-22
    0.00e+00 1.39e-04 2.24e-04 2.29e-04 -3.53e-04 2.81e-04 -3.51e-04 -5.97e-05 -1.39e-04 -1.00e-21 5.00e-04 3.38e-04 1.08e-04 -1.27e-04 1.15e-04 4.97e-05 -2.24e-04 -5.00e-04 -1.57e-21 5.88e-05 -1.65e-04 -5.24e-04 -7.78e-05 1.32e-04 -2.29e-04 -3.38e-04 -5.88e-05 1.26e-20 -1.02e-04 -5.42e-04 3.43e-05 2.13e-04 3.53e-04 -1.08e-04 1.65e-04 1.02e-04 -1.26e-21 2.06e-05 5.69e-05 7.12e-05 -2.81e-04 1.27e-04 5.24e-04 5.42e-04 -2.06e-05 -2.05e-21 -3.08e-04 1.91e-05 3.51e-04 -1.15e-04 7.78e-05 -3.43e-05 -5.69e-05 3.08e-04 2.92e-21 -8.52e-05 5.97e-05 -4.97e-05 -1.32e-04 -2.13e-04 -7.12e-05 -1.91e-05 8.52e-05 2.41e-22
    0.00e+00 8.77e-05 2.61e-04 2.54e-04 -3.69e-04 1.04e-04 -3.91e-04 -3.90e-05 -8.77e-05 -1.44e-21 4.51e-04 3.65e-04 1.40e-04 -3.84e-04 1.71e-04 4.11e-06 -2.61e-04 -4.51e-04 3.19e-21 4.18e-05 -1.43e-04 -5.56e-04 1.34e-05 1.95e-04 -2.54e-04 -3.65e-04 -4.18e-05 1.19e-22 -1.28e-04 -5.88e-04 8.99e-05 1.86e-04 3.69e-04 -1.40e-04 1.43e-04 1.28e-04 -1.67e-21 -3.07e-05 4.04e-06 -1.85e-05 -1.04e-04 3.84e-04 5.56e-04 5.88e-04 3.07e-05 -2.37e-21 -3.43e-04 -1.07e-04 3.91e-04 -1.71e-04 -1.34e-05 -8.99e-05 -4.04e-06 3.43e-04 5.90e-22 -5.90e-05 3.90e-05 -4.11e-06 -1.95e-04 -1.86e-04 1.85e-05 1.07e-04 5.90e-05 -1.02e-22
    0.00e+00 -2.83e-05 4.07e-04 -5.14e-05 5.08e-05 4.53e-04 5.36e-04 6.55e-04 2.83e-05 -6.19e-22 -3.62e-05 -1.91e-04 -9.11e-05 -1.36e-04 4.67e-05 2.03e-05 -4.07e-04 3.62e-05 2.17e-21 -6.53e-05 -1.38e-04 -4.16e-05 -1.33e-05 -1.46e-04 5.14e-05 1.91e-04 6.53e-05 2.34e-21 1.23e-04 2.74e-04 3.77e-04 3.58e-04 -5.08e-05 9.11e-05 1.38e-04 -1.23e-04 -1.61e-21 4.52e-04 4.12e-04 5.88e-04 -4.53e-04 1.36e-04 4.16e-05 -2.74e-04 -4.52e-04 9.92e-22 8.75e-05 1.85e-04 -5.36e-04 -4.67e-05 1.33e-05 -3.77e-04 -4.12e-04 -8.75e-05 -7.98e-22 8.82e-05 -6.55e-04 -2.03e-05 1.46e-04 -3.58e-04 -5.88e-04 -1.85e-04 -8.82e-05 -1.66e-22
    0.00e+00 3.27e-04 6.24e-04 3.76e-04 -4.57e-04 -7.86e-05 -2.51e-04 -1.23e-04 -3.27e-04 2.13e-22 1.56e-04 1.92e-04 2.49e-04 -1.64e-04 2.06e-04 -4.14e-06 -6.24e-04 -1.56e-04 -4.07e-21 1.60e-04 -1.87e-04 -7.04e-04 6.21e-05 1.39e-04 -3.76e-04 -1.92e-04 -1.60e-04 1.60e-21 -3.32e-05 -4.24e-04 4.29e-05 1.85e-04 4.57e-04 -2.49e-04 1.87e-04 3.32e-05 -2.03e-22 -1.33e-05 -1.52e-04 6.62e-05 7.86e-05 1.64e-04 7.04e-04 4.24e-04 1.33e-05 2.35e-21 -1.80e-04 -6.73e-05 2.51e-04 -2.06e-04 -6.21e-05 -4.29e-05 1.52e-04 1.80e-04 -3.64e-22 7.70e-05 1.23e-04 4.14e-06 -1.39e-04 -1.85e-04 -6.62e-05 6.73e-05 -7.70e-05 4.49e-21
    0.00e+00 2.32e-04 -1.89e-04 3.08e-05 -3.32e-05 -2.64e-04 1.14e-04 5.52e-05 -2.32e-04 8.94e-22 9.81e-05 1.62e-04 -1.16e-04 2.43e-05 -8.07e-04 3.26e-04 1.89e-04 -9.81e-05 -1.18e-22 -7.95e-05 -1.77e-04 2.83e-04 -3.62e-04 2.27e-05 -3.08e-05 -1.62e-04 7.95e-05 -1.30e-21 1.75e-05 6.09e-05 2.24e-04 -1.59e-04 3.32e-05 1.16e-04 1.77e-04 -1.75e-05 -6.83e-22 7.63e-05 -4.71e-04 6.64e-05 2.64e-04 -2.43e-05 -2.83e-04 -6.09e-05 -7.63e-05 -7.48e-22 5.14e-04 -8.81e-05 -1.14e-04 8.07e-04 3.62e-04 -2.24e-04 4.71e-04 -5.14e-04 -2.60e-21 -5.48e-04 -5.52e-05 -3.26e-04 -2.27e-05 1.59e-04 -6.64e-05 8.81e-05 5.48e-04 -1.59e-21
    0.00e+00 4.24e-04 4.76e-05 7.00e-05 8.40e-05 -5.85e-04 4.84e-04 -3.37e-05 -4.24e-04 1.22e-21 1.70e-04 4.98e-04 -3.59e-04 3.14e-04 6.63e-05 1.93e-04 -4.76e-05 -1.70e-04 -1.82e-22 -3.64e-05 -1.80e-04 2.32e-04 -2.87e-04 1.94e-05 -7.00e-05 -4.98e-04 3.64e-05 -2.68e-21 -1.63e-04 -1.60e-04 2.66e-04 -2.03e-04 -8.40e-05 3.59e-04 1.80e-04 1.63e-04 1.46e-22 -1.85e-04 3.10e-04 5.44e-05 5.85e-04 -3.14e-04 -2.32e-04 1.60e-04 1.85e-04 4.00e-22 1.46e-05 1.78e-04 -4.84e-04 -6.63e-05 2.87e-04 -2.66e-04 -3.10e-04 -1.46e-05 1.19e-21 -3.42e-04 3.37e-05 -1.93e-04 -1.94e-05 2.03e-04 -5.44e-05 -1.78e-04 3.42e-04 8.91e-22
    0.00e+00 9.17e-05 -3.43e-04 -1.46e-04 -1.11e-04 -1.12e-04 2.12e-05 -2.15e-04 -9.17e-05 -5.23e-22 -2.62e-04 1.13e-04 -1.01e-04 2.53e-04 -4.48e-04 8.15e-05 3.43e-04 2.62e-04 -2.56e-21 3.82e-05 -3.07e-04 4.72e-06 -2.15e-04 -9.06e-05 1.46e-04 -1.13e-04 -3.82e-05 1.37e-21 4.36e-05 -2.57e-05 1.78e-04 -1.39e-04 1.11e-04 1.01e-04 3.07e-04 -4.36e-05 9.68e-23 -7.33e-05 1.82e-04 -5.97e-05 1.12e-04 -2.53e-04 -4.72e-06 2.57e-05 7.33e-05 -4.24e-21 2.09e-04 -1.49e-04 -2.12e-05 4.48e-04 2.15e-04 -1.78e-04 -1.82e-04 -2.09e-04 -1.30e-21 1.87e-04 2.15e-04 -8.15e-05 9.06e-05 1.39e-04 5.97e-05 1.49e-04 -1.87e-04 1.70e-21
    0.00e+00 9.15e-05 6.10e-04 2.49e-04 -1.20e-04 2.44e-04 1.02e-05 1.43e-04 -9.15e-05 -8.28e-22 9.22e-05 4.68e-05 -1.68e-04 -1.90e-04 1.14e-04 -5.31e-05 -6.10e-04 -9.22e-05 -3.50e-21 1.49e-04 -8.44e-05 -3.90e-04 3.47e-05 -3.04e-05 -2.49e-04 -4.68e-05 -1.49e-04 -4.18e-21 7.53e-05 -1.30e-04 2.22e-05 1.15e-04 1.20e-04 1.68e-04 8.44e-05 -7.53e-05 -2.13e-22 2.19e-04 -1.82e-04 2.93e-04 -2.44e-04 1.90e-04 3.90e-04 1.30e-04 -2.19e-04 -1.41e-21 -8.40e-05 1.82e-04 -1.02e-05 -1.14e-04 -3.47e-05 -2.22e-05 1.82e-04 8.40e-05 1.21e-22 7.80e-05 -1.43e-04 5.31e-05 3.04e-05 -1.15e-04 -2.93e-04 -1.82e-04 -7.80e-05 -6.97e-22
    0.00e+00 1.91e-04 -2.96e-04 -1.35e-04 -1.53e-04 -2.02e-04 3.69e-04 -2.48e-04 -1.91e-04 1.06e-21 -2.17e-04 2.54e-05 -1.60e-04 1.53e-04 -6.74e-04 2.49e-04 2.96e-04 2.17e-04 8.83e-22 -1.93e-04 -3.68e-04 1.62e-04 1.87e-04 1.03e-05 1.35e-04 -2.54e-05 1.93e-04 -1.19e-21 9.90e-05 -9.80e-05 4.03e-04 -2.45e-04 1.53e-04 1.60e-04 3.68e-04 -9.90e-05 9.24e-22 -1.19e-04 2.45e-04 -1.22e-04 2.02e-04 -1.53e-04 -1.62e-04 9.80e-05 1.19e-04 4.03e-23 4.78e-04 -1.92e-04 -3.69e-04 6.74e-04 -1.87e-04 -4.03e-04 -2.45e-04 -4.78e-04 -9.57e-22 -3.97e-04 2.48e-04 -2.49e-04 -1.03e-05 2.45e-04 1.22e-04 1.92e-04 3.97e-04 -4.34e-22
    0.00e+00 5.67e-06 -1.24e-04 -1.09e-04 -6.02e-05 -1.93e-04 -4.03e-05 2.54e-04 -5.67e-06 1.53e-21 -1.41e-04 3.09e-04 -1.17e-04 4.98e-04 -1.95e-04 -4.86e-04 1.24e-04 1.41e-04 -2.14e-21 1.08e-04 -3.00e-04 -4.00e-04 -6.49e-04 3.91e-04 1.09e-04 -3.09e-04 -1.08e-04 1.49e-21 1.18e-05 1.71e-04 -1.12e-06 -2.35e-04 6.02e-05 1.17e-04 3.00e-04 -1.18e-05 -9.61e-22 1.02e-05 2.12e-04 -7.11e-05 1.93e-04 -4.98e-04 4.00e-04 -1.71e-04 -1.02e-05 -3.97e-22 7.80e-05 -1.04e-04 4.03e-05 1.95e-04 6.49e-04 1.12e-06 -2.12e-04 -7.80e-05 1.49e-21 6.01e-06 -2.54e-04 4.86e-04 -3.91e-04 2.35e-04 7.11e-05 1.04e-04 -6.01e-06 1.59e-21
    0.00e+00 -1.73e-04 -1.55e-05 7.85e-06 6.16e-04 2.28e-05 9.41e-05 3.70e-04 1.73e-04 3.59e-23 -7.64e-05 -2.24e-04 -3.92e-04 -1.78e-04 3.60e-05 -2.25e-04 1.55e-05 7.64e-05 2.05e-21 2.95e-05 5.88e-05 7.79e-05 -3.85e-05 -3.61e-04 -7.85e-06 2.24e-04 -2.95e-05 1.31e-22 -1.68e-04 1.76e-04 -4.66e-05 -4.57e-04 -6.16e-04 3.92e-04 -5.88e-05 1.68e-04 1.17e-21 8.62e-04 -2.65e-06 9.19e-05 -2.28e-05 1.78e-04 -7.79e-05 -1.76e-04 -8.62e-04 -1.19e-21 4.50e-05 -1.87e-05 -9.41e-05 -3.60e-05 3.85e-05 4.66e-05 2.65e-06 -4.50e-05 5.60e-22 4.49e-04 -3.70e-04 2.25e-04 3.61e-04 4.57e-04 -9.19e-05 1.87e-05 -4.49e-04 1.21e-21
    0.00e+00 2.66e-04 -1.01e-04 -9.58e-05 -4.51e-05 -2.72e-05 1.51e-04 2.66e-04 -2.66e-04 5.31e-22 6.63e-05 -4.78e-05 2.19e-04 6.18e-04 7.11e-04 -3.18e-04 1.01e-04 -6.63e-05 1.56e-21 -3.49e-04 3.86e-05 -3.78e-04 1.32e-04 7.88e-04 9.58e-05 4.78e-05 3.49e-04 1.80e-21 -1.75e-05 2.07e-04 2.55e-04 -7.15e-06 4.51e-05 -2.19e-04 -3.86e-05 1.75e-05 -1.06e-22 -8.99e-05 -6.64e-05 -1.86e-04 2.72e-05 -6.18e-04 3.78e-04 -2.07e-04 8.99e-05 8.11e-22 1.40e-04 -8.82e-06 -1.51e-04 -7.11e-04 -1.32e-04 -2.55e-04 6.64e-05 -1.40e-04 3.59e-22 -1.13e-04 -2.66e-04 3.18e-04 -7.88e-04 7.15e-06 1.86e-04 8.82e-06 1.13e-04 -3.54e-22
    0.00e+00 5.80e-05 -4.26e-04 -5.43e-05 -1.04e-04 1.61e-05 6.33e-05 -2.78e-04 -5.80e-05 -1.98e-22 -4.15e-04 -5.60e-05 -7.41e-05 2.06e-04 -3.53e-04 -1.84e-04 4.26e-04 4.15e-04 -5.82e-22 1.78e-04 -4.06e-05 1.49e-04 -3.71e-04 1.39e-05 5.43e-05 5.60e-05 -1.78e-04 2.47e-22 -6.10e-05 -1.21e-04 2.75e-04 -2.15e-04 1.04e-04 7.41e-05 4.06e-05 6.10e-05 2.22e-21 -7.66e-05 7.87e-05 -7.50e-05 -1.61e-05 -2.06e-04 -1.49e-04 1.21e-04 7.66e-05 -4.75e-22 1.97e-04 -8.58e-05 -6.33e-05 3.53e-04 3.71e-04 -2.75e-04 -7.87e-05 -1.97e-04 7.89e-22 2.45e-04 2.78e-04 1.84e-04 -1.39e-05 2.15e-04 7.50e-05 8.58e-05 -2.45e-04 1.64e-21
    0.00e+00 1.14e-04 -3.79e-04 -1.46e-04 -1.07e-04 -8.41e-05 6.55e-05 -2.82e-04 -1.14e-04 -4.76e-23 -3.45e-04 9.02e-05 -1.02e-04 2.40e-04 -5.12e-04 6.07e-05 3.79e-04 3.45e-04 -9.11e-22 6.65e-05 -2.63e-04 5.84e-05 -1.90e-04 -1.04e-04 1.46e-04 -9.02e-05 -6.65e-05 -9.62e-22 2.44e-05 -6.33e-05 2.39e-04 -1.92e-04 1.07e-04 1.02e-04 2.63e-04 -2.44e-05 6.83e-22 -9.55e-05 1.35e-04 -6.58e-05 8.41e-05 -2.40e-04 -5.84e-05 6.33e-05 9.55e-05 -6.43e-22 2.00e-04 -1.23e-04 -6.55e-05 5.12e-04 1.90e-04 -2.39e-04 -1.35e-04 -2.00e-04 3.84e-22 1.84e-04 2.82e-04 -6.07e-05 1.04e-04 1.92e-04 6.58e-05 1.23e-04 -1.84e-04 -9.43e-22
    0.00e+00 3.42e-04 -8.13e-05 5.28e-05 -8.71e-05 1.10e-04 1.12e-04 -4.43e-05 -3.42e-04 -1.14e-22 2.53e-04 -2.50e-04 2.69e-04 6.44e-04 8.23e-04 3.01e-05 8.13e-05 -2.53e-04 4.60e-22 -2.20e-04 6.51e-05 -1.57e-04 1.01e-04 3.34e-04 -5.28e-05 2.50e-04 2.20e-04 1.45e-21 -1.59e-05 2.28e-04 2.22e-04 6.91e-05 8.71e-05 -2.69e-04 -6.51e-05 1.59e-05 -1.61e-21 -1.30e-04 -9.38e-05 -7.86e-05 -1.10e-04 -6.44e-04 1.57e-04 -2.28e-04 1.30e-04 -1.57e-21 -1.35e-05 7.33e-05 -1.12e-04 -8.23e-04 -1.01e-04 -2.22e-04 9.38e-05 1.35e-05 2.18e-21 5.94e-05 4.43e-05 -3.01e-05 -3.34e-04 -6.91e-05 7.86e-05 -7.33e-05 -5.94e-05 -1.83e-21
    0.00e+00 5.72e-05 -1.51e-04 -9.17e-06 2.41e-05 -2.13e-04 1.82e-04 2.12e-04 -5.72e-05 2.04e-21 -1.27e-04 3.34e-04 -2.22e-04 4.52e-04 -3.33e-04 -5.17e-04 1.51e-04 1.27e-04 -1.11e-22 8.25e-05 -1.37e-04 -4.30e-04 -4.30e-04 4.94e-04 9.17e-06 -3.34e-04 -8.25e-05 -3.23e-21 -5.13e-05 1.82e-04 -1.46e-04 -2.05e-04 -2.41e-05 2.22e-04 1.37e-04 5.13e-05 -1.49e-22 -3.89e-06 1.84e-04 -3.98e-06 2.13e-04 -4.52e-04 4.30e-04 -1.82e-04 3.89e-06 4.39e-22 -1.75e-04 3.55e-05 -1.82e-04 3.33e-04 4.30e-04 1.46e-04 -1.84e-04 1.75e-04 4.01e-22 -1.53e-04 -2.12e-04 5.17e-04 -4.94e-04 2.05e-04 3.98e-06 -3.55e-05 1.53e-04 3.71e-21
    0.00e+00 8.27e-05 -1.55e-04 -5.38e-05 -1.70e-05 -1.53e-04 2.30e-05 1.11e-04 -8.27e-05 -2.31e-21 -1.10e-04 2.97e-04 -1.00e-04 5.44e-04 -9.94e-05 -6.18e-04 1.55e-04 1.10e-04 -5.77e-21 3.64e-05 -1.59e-04 -5.08e-04 -5.17e-04 5.56e-04 5.38e-05 -2.97e-04 -3.64e-05 -1.72e-22 -1.96e-05 1.56e-04 -6.40e-05 -1.51e-04 1.70e-05 1.00e-04 1.59e-04 1.96e-05 -8.36e-22 -3.22e-05 1.85e-04 -1.48e-05 1.53e-04 -5.44e-04 5.08e-04 -1.56e-04 3.22e-05 5.60e-21 -6.81e-06 -2.73e-06 -2.30e-05 9.94e-05 5.17e-04 6.40e-05 -1.85e-04 6.81e-06 1.66e-21 -9.28e-05 -1.11e-04 6.18e-04 -5.56e-04 1.51e-04 1.48e-05 2.73e-06 9.28e-05 -5.41e-22
    0.00e+00 4.49e-04 -7.62e-05 1.44e-04 1.63e-04 -4.82e-04 1.56e-04 9.03e-05 -4.49e-04 -3.03e-21 1.82e-04 2.54e-04 -3.03e-04 -4.25e-05 -2.84e-04 2.40e-04 7.62e-05 -1.82e-04 -1.41e-21 -9.21e-05 -1.43e-04 3.76e-04 -2.44e-04 -5.44e-06 -1.44e-04 -2.54e-04 9.21e-05 -1.39e-21 -1.18e-04 -1.76e-04 3.75e-04 -1.94e-04 -1.63e-04 3.03e-04 1.43e-04 1.18e-04 -5.50e-22 -2.17e-04 -2.08e-05 1.04e-04 4.82e-04 4.25e-05 -3.76e-04 1.76e-04 2.17e-04 -4.24e-21 7.40e-05 2.03e-04 -1.56e-04 2.84e-04 2.44e-04 -3.75e-04 2.08e-05 -7.40e-05 1.55e-21 -4.97e-04 -9.03e-05 -2.40e-04 5.44e-06 1.94e-04 -1.04e-04 -2.03e-04 4.97e-04 9.24e-23
    0.00e+00 3.10e-04 -1.45e-05 2.44e-05 6.94e-05 -2.46e-04 4.25e-04 -3.71e-05 -3.10e-04 5.69e-22 1.06e-04 5.05e-04 -2.78e-04 2.47e-04 -1.81e-05 1.83e-05 1.45e-05 -1.06e-04 -2.51e-21 -2.92e-06 -7.34e-05 3.46e-05 -2.52e-04 1.21e-04 -2.44e-05 -5.05e-04 2.92e-06 -7.66e-22 -1.08e-04 -2.16e-04 -2.49e-05 -1.36e-04 -6.94e-05 2.78e-04 7.34e-05 1.08e-04 -6.30e-22 -3.57e-05 2.65e-04 2.65e-05 2.46e-04 -2.47e-04 -3.46e-05 2.16e-04 3.57e-05 -1.21e-21 -2.19e-04 8.55e-05 -4.25e-04 1.81e-05 2.52e-04 2.49e-05 -2.65e-04 2.19e-04 3.28e-21 -1.55e-04 3.71e-05 -1.83e-05 -1.21e-04 1.36e-04 -2.65e-05 -8.55e-05 1.55e-04 3.28e-21
    0.00e+00 1.65e-04 1.10e-04 6.16e-05 7.43e-05 -7.19e-04 -6.04e-05 2.74e-05 -1.65e-04 -1.58e-21 5.19e-05 1.44e-04 2.69e-05 -6.96e-04 -8.43e-05 1.00e-04 -1.10e-04 -5.19e-05 4.62e-22 1.20e-05 -1.59e-04 4.49e-05 1.80e-05 -3.23e-05 -6.16e-05 -1.44e-04 -1.20e-05 6.41e-24 -3.89e-05 1.74e-04 6.98e-05 2.65e-05 -7.43e-05 -2.69e-05 1.59e-04 3.89e-05 -1.76e-21 -7.35e-04 1.35e-04 1.65e-04 7.19e-04 6.96e-04 -4.49e-05 -1.74e-04 7.35e-04 5.09e-21 -4.19e-04 -5.49e-04 6.04e-05 8.43e-05 -1.80e-05 -6.98e-05 -1.35e-04 4.19e-04 -3.53e-21 -2.35e-04 -2.74e-05 -1.00e-04 3.23e-05 -2.65e-05 -1.65e-04 5.49e-04 2.35e-04 5.28e-21
    0.00e+00 3.88e-05 5.95e-04 2.00e-04 -1.36e-04 1.73e-04 6.46e-05 1.20e-04 -3.88e-05 4.60e-22 1.21e-04 4.48e-05 -2.57e-04 -2.04e-04 4.32e-05 -3.55e-05 -5.95e-04 -1.21e-04 -2.46e-22 1.60e-04 -5.99e-05 -3.78e-04 4.57e-05 -8.25e-05 -2.00e-04 -4.48e-05 -1.60e-04 -1.67e-21 1.20e-04 -3.62e-05 9.52e-06 1.13e-04 1.36e-04 2.57e-04 5.99e-05 -1.20e-04 -1.34e-21 2.44e-04 -2.02e-04 3.61e-04 -1.73e-04 2.04e-04 3.78e-04 3.62e-05 -2.44e-04 -3.24e-21 -1.13e-04 2.44e-04 -6.46e-05 -4.32e-05 -4.57e-05 -9.52e-06 2.02e-04 1.13e-04 1.34e-21 5.23e-05 -1.20e-04 3.55e-05 8.25e-05 -1.13e-04 -3.61e-04 -2.44e-04 -5.23e-05 2.77e-21
    0.00e+00 2.68e-04 -1.92e-04 1.16e-04 -1.11e-04 4.17e-05 1.19e-04 -8.92e-05 -2.68e-04 7.94e-22 2.51e-04 -1.52e-04 2.67e-04 7.81e-04 9.63e-04 1.29e-04 1.92e-04 -2.51e-04 -8.71e-22 -2.52e-04 1.26e-04 -2.30e-04 1.30e-04 5.20e-04 -1.16e-04 1.52e-04 2.52e-04 -1.69e-21 -6.01e-06 3.93e-04 3.69e-04 3.32e-05 1.11e-04 -2.67e-04 -1.26e-04 6.01e-06 1.08e-22 -1.40e-04 -9.90e-05 -1.03e-05 -4.17e-05 -7.81e-04 2.30e-04 -3.93e-04 1.40e-04 4.39e-22 -1.02e-04 1.87e-04 -1.19e-04 -9.63e-04 -1.30e-04 -3.69e-04 9.90e-05 1.02e-04 3.25e-21 2.38e-05 8.92e-05 -1.29e-04 -5.20e-04 -3.32e-05 1.03e-05 -1.87e-04 -2.38e-05 -1.61e-21
    0.00e+00 4.00e-04 -8.71e-05 1.37e-04 1.74e-04 -4.93e-04 2.28e-04 9.36e-05 -4.00e-04 1.65e-21 1.82e-04 3.81e-04 -3.60e-04 1.77e-04 -3.88e-05 1.82e-04 8.71e-05 -1.82e-04 -4.69e-21 -7.10e-05 -1.65e-04 3.89e-04 -3.33e-04 6.34e-05 -1.37e-04 -3.81e-04 7.10e-05 -3.25e-22 -1.22e-04 -1.70e-04 4.43e-04 -2.62e-04 -1.74e-04 3.60e-04 1.65e-04 1.22e-04 5.64e-23 -1.36e-04 7.90e-05 7.21e-05 4.93e-04 -1.77e-04 -3.89e-04 1.70e-04 1.36e-04 -3.62e-22 -9.17e-05 2.33e-04 -2.28e-04 3.88e-05 3.33e-04 -4.43e-04 -7.90e-05 9.17e-05 8.35e-22 -4.54e-04 -9.36e-05 -1.82e-04 -6.34e-05 2.62e-04 -7.21e-05 -2.33e-04 4.54e-04 -1.66e-23
    0.00e+00 2.84e-04 8.49e-05 2.12e-05 -7.56e-05 -3.26e-04 4.44e-04 -1.02e-04 -2.84e-04 -9.15e-22 2.61e-04 4.55e-04 -1.99e-04 4.31e-04 -5.92e-05 1.04e-04 -8.49e-05 -2.61e-04 -1.07e-21 -1.60e-05 -1.04e-04 -8.49e-05 -1.88e-04 8.25e-05 -2.12e-05 -4.55e-04 1.60e-05 -3.18e-21 -1.08e-04 -1.31e-04 -7.68e-06 -6.37e-05 7.56e-05 1.99e-04 1.04e-04 1.08e-04 -5.52e-22 -6.69e-05 2.94e-04 1.69e-05 3.26e-04 -4.31e-04 8.49e-05 1.31e-04 6.69e-05 -5.12e-22 -7.48e-05 3.69e-05 -4.44e-04 5.92e-05 1.88e-04 7.68e-06 -2.94e-04 7.48e-05 -1.78e-21 -1.01e-04 1.02e-04 -1.04e-04 -8.25e-05 6.37e-05 -1.69e-05 -3.69e-05 1.01e-04 1.60e-21
    0.00e+00 -1.20e-04 -1.51e-04 -1.51e-04 -5.42e-05 1.80e-04 5.76e-05 -1.67e-04 1.20e-04 -1.19e-21 -9.87e-05 -5.17e-04 -5.31e-05 8.27e-05 1.58e-04 -5.37e-04 1.51e-04 9.87e-05 -3.19e-21 -2.57e-04 1.51e-04 1.29e-05 1.54e-04 -2.32e-04 1.51e-04 5.17e-04 2.57e-04 4.59e-22 2.34e-04 -1.24e-04 2.85e-07 -8.27e-05 5.42e-05 5.31e-05 -1.51e-04 -2.34e-04 5.95e-22 1.04e-04 5.37e-05 -5.56e-04 -1.80e-04 -8.27e-05 -1.29e-05 1.24e-04 -1.04e-04 -9.31e-23 -3.35e-05 1.07e-04 -5.76e-05 -1.58e-04 -1.54e-04 -2.85e-07 -5.37e-05 3.35e-05 6.96e-22 -3.79e-05 1.67e-04 5.37e-04 2.32e-04 8.27e-05 5.56e-04 -1.07e-04 3.79e-05 -7.80e-22
    0.00e+00 -3.49e-04 2.34e-04 6.32e-05 -2.52e-04 5.75e-04 -1.36e-04 -4.53e-04 3.49e-04 -1.79e-21 1.17e-05 -2.53e-04 -2.95e-06 -6.38e-05 1.70e-04 9.35e-05 -2.34e-04 -1.17e-05 7.63e-22 4.40e-04 7.55e-05 6.56e-05 4.57e-05 -7.77e-05 -6.32e-05 2.53e-04 -4.40e-04 -4.43e-22 -2.07e-04 -3.14e-04 5.57e-05 6.05e-04 2.52e-04 2.95e-06 -7.55e-05 2.07e-04 -1.74e-23 1.79e-05 2.17e-05 8.65e-05 -5.75e-04 6.38e-05 -6.56e-05 3.14e-04 -1.79e-05 -2.47e-21 -6.95e-05 -5.37e-05 1.36e-04 -1.70e-04 -4.57e-05 -5.57e-05 -2.17e-05 6.95e-05 -3.00e-21 -2.90e-05 4.53e-04 -9.35e-05 7.77e-05 -6.05e-04 -8.65e-05 5.37e-05 2.90e-05 -1.92e-22
    0.00e+00 1.96e-04 7.83e-05 1.78e-04 -1.94e-04 1.18e-04 1.47e-05 -1.80e-04 -1.96e-04 -3.57e-21 1.79e-04 -9.84e-05 2.50e-04 7.26e-04 8.83e-04 2.81e-04 -7.83e-05 -1.79e-04 1.95e-21 -2.20e-04 2.32e-04 -3.65e-05 1.11e-04 3.88e-04 -1.78e-04 9.84e-05 2.20e-04 -1.58e-21 -5.48e-05 4.35e-04 3.28e-04 9.09e-05 1.94e-04 -2.50e-04 -2.32e-04 5.48e-05 7.76e-22 -2.01e-04 -1.08e-04 -4.60e-05 -1.18e-04 -7.26e-04 3.65e-05 -4.35e-04 2.01e-04 -1.34e-21 -2.86e-04 2.82e-04 -1.47e-05 -8.83e-04 -1.11e-04 -3.28e-04 1.08e-04 2.86e-04 1.72e-21 1.20e-04 1.80e-04 -2.81e-04 -3.88e-04 -9.09e-05 4.60e-05 -2.82e-04 -1.20e-04 -3.56e-21
    0.00e+00 3.26e-04 -1.30e-04 -4.51e-05 -4.71e-05 1.58e-05 1.51e-04 2.16e-04 -3.26e-04 -1.31e-21 1.02e-04 -9.25e-05 2.49e-04 6.64e-04 8.17e-04 -2.08e-04 1.30e-04 -1.02e-04 -2.93e-21 -3.47e-04 5.16e-05 -3.56e-04 1.87e-04 8.19e-04 4.51e-05 9.25e-05 3.47e-04 -1.59e-22 -8.08e-06 2.36e-04 3.00e-04 4.17e-06 4.71e-05 -2.49e-04 -5.16e-05 8.08e-06 -4.23e-22 -9.59e-05 -1.07e-04 -1.39e-04 -1.58e-05 -6.64e-04 3.56e-04 -2.36e-04 9.59e-05 -2.69e-21 8.19e-05 2.53e-05 -1.51e-04 -8.17e-04 -1.87e-04 -3.00e-04 1.07e-04 -8.19e-05 3.51e-22 -7.03e-05 -2.16e-04 2.08e-04 -8.19e-04 -4.17e-06 1.39e-04 -2.53e-05 7.03e-05 -1.02e-21
    0.00e+00 3.19e-04 -2.41e-04 3.16e-05 -5.84e-05 4.26e-05 1.01e-04 2.41e-05 -3.19e-04 -5.14e-21 1.83e-04 -9.43e-05 2.61e-04 7.30e-04 9.44e-04 -4.00e-05 2.41e-04 -1.83e-04 -2.81e-21 -2.66e-04 5.93e-05 -3.08e-04 1.38e-04 6.34e-04 -3.16e-05 9.43e-05 2.66e-04 -1.41e-21 3.97e-06 3.22e-04 3.22e-04 -3.04e-06 5.84e-05 -2.61e-04 -5.93e-05 -3.97e-06 -4.00e-22 -1.03e-04 -8.26e-05 -3.42e-05 -4.26e-05 -7.30e-04 3.08e-04 -3.22e-04 1.03e-04 -6.37e-22 1.15e-05 1.00e-04 -1.01e-04 -9.44e-04 -1.38e-04 -3.22e-04 8.26e-05 -1.15e-05 7.38e-23 -1.18e-05 -2.41e-05 4.00e-05 -6.34e-04 3.04e-06 3.42e-05 -1.00e-04 1.18e-05 2.68e-22
    0.00e+00 8.52e-05 -3.27e-04 -1.39e-04 -1.44e-04 -6.34e-05 1.26e-04 -3.00e-04 -8.52e-05 2.04e-21 -4.01e-04 5.57e-05 -9.93e-05 2.27e-04 -5.35e-04 -3.32e-05 3.27e-04 4.01e-04 -2.76e-22 3.25e-05 -1.78e-04 1.25e-04 -1.27e-04 1.01e-06 1.39e-04 -5.57e-05 -3.25e-05 -6.99e-22 4.45e-05 -1.03e-04 2.86e-04 -2.08e-04 1.44e-04 9.93e-05 1.78e-04 -4.45e-05 2.81e-21 -8.25e-05 1.62e-04 -8.90e-05 6.34e-05 -2.27e-04 -1.25e-04 1.03e-04 8.25e-05 -1.13e-21 1.99e-04 -1.36e-04 -1.26e-04 5.35e-04 1.27e-04 -2.86e-04 -1.62e-04 -1.99e-04 -2.39e-21 1.22e-04 3.00e-04 3.32e-05 -1.01e-06 2.08e-04 8.90e-05 1.36e-04 -1.22e-04 -3.23e-22
    0.00e+00 2.75e-04 5.35e-04 3.89e-04 -4.75e-04 -4.09e-06 -3.15e-04 -1.26e-04 -2.75e-04 1.15e-21 2.44e-04 2.09e-04 2.82e-04 -1.65e-04 2.34e-04 3.11e-05 -5.35e-04 -2.44e-04 1.92e-21 1.44e-04 -1.85e-04 -6.64e-04 4.43e-05 1.55e-04 -3.89e-04 -2.09e-04 -1.44e-04 2.18e-21 -5.64e-05 -4.82e-04 6.73e-05 1.86e-04 4.75e-04 -2.82e-04 1.85e-04 5.64e-05 3.59e-21 -1.84e-05 -1.24e-04 6.71e-05 4.09e-06 1.65e-04 6.64e-04 4.82e-04 1.84e-05 4.49e-22 -2.44e-04 -7.34e-05 3.15e-04 -2.34e-04 -4.43e-05 -6.73e-05 1.24e-04 2.44e-04 2.93e-22 6.36e-05 1.26e-04 -3.11e-05 -1.55e-04 -1.86e-04 -6.71e-05 7.34e-05 -6.36e-05 -9.16e-23
    0.00e+00 2.39e-04 -6.33e-05 -1.49e-04 -5.68e-05 2.92e-05 1.42e-04 3.33e-04 -2.39e-04 4.04e-21 -1.19e-06 7.09e-05 2.11e-04 5.89e-04 6.89e-04 -2.76e-04 6.33e-05 1.19e-06 -5.66e-22 -4.14e-04 5.26e-05 -2.71e-04 2.75e-04 8.95e-04 1.49e-04 -7.09e-05 4.14e-04 -1.10e-21 -4.14e-05 1.20e-04 2.15e-04 3.03e-05 5.68e-05 -2.11e-04 -5.26e-05 4.14e-05 2.26e-22 -4.09e-05 -5.73e-05 -2.22e-04 -2.92e-05 -5.89e-04 2.71e-04 -1.20e-04 4.09e-05 -1.43e-21 1.75e-04 -6.98e-05 -1.42e-04 -6.89e-04 -2.75e-04 -2.15e-04 5.73e-05 -1.75e-04 -4.00e-22 -5.92e-05 -3.33e-04 2.76e-04 -8.95e-04 -3.03e-05 2.22e-04 6.98e-05 5.92e-05 6.40e-22
    0.00e+00 3.53e-04 1.79e-04 1.75e-04 -2.06e-04 2.27e-04 5.32e-05 2.10e-05 -3.53e-04 5.97e-22 2.05e-04 -1.41e-04 2.43e-04 6.28e-04 4.93e-04 2.26e-04 -1.79e-04 -2.05e-04 -1.36e-21 -1.15e-04 1.19e-04 1.21e-04 1.00e-04 1.91e-04 -1.75e-04 1.41e-04 1.15e-04 1.34e-21 -2.77e-05 3.13e-04 1.84e-04 1.63e-04 2.06e-04 -2.43e-04 -1.19e-04 2.77e-05 -1.23e-22 -2.23e-04 -1.23e-04 -1.12e-04 -2.27e-04 -6.28e-04 -1.21e-04 -3.13e-04 2.23e-04 -1.56e-21 -2.65e-04 1.58e-04 -5.32e-05 -4.93e-04 -1.00e-04 -1.84e-04 1.23e-04 2.65e-04 -2.06e-21 1.04e-04 -2.10e-05 -2.26e-04 -1.91e-04 -1.63e-04 1.12e-04 -1.58e-04 -1.04e-04 2.38e-21
    0.00e+00 2.44e-04 -2.25e-04 -1.39e-04 -1.97e-04 -3.24e-04 3.00e-04 -2.23e-04 -2.44e-04 3.25e-22 -8.29e-05 1.78e-04 -1.31e-04 1.08e-04 -7.16e-04 3.84e-04 2.25e-04 8.29e-05 -2.70e-21 -2.23e-04 -4.48e-04 1.18e-04 -4.44e-05 -4.18e-05 1.39e-04 -1.78e-04 2.23e-04 6.09e-22 1.92e-04 3.69e-05 1.94e-04 -1.84e-04 1.97e-04 1.31e-04 4.48e-04 -1.92e-04 -3.99e-23 -6.91e-05 -1.22e-04 -9.15e-05 3.24e-04 -1.08e-04 -1.18e-04 -3.69e-05 6.91e-05 -2.02e-21 6.45e-04 -2.84e-04 -3.00e-04 7.16e-04 4.44e-05 -1.94e-04 1.22e-04 -6.45e-04 -5.14e-21 -3.79e-04 2.23e-04 -3.84e-04 4.18e-05 1.84e-04 9.15e-05 2.84e-04 3.79e-04 8.09e-22
    0.00e+00 2.53e-04 6.89e-04 1.77e-04 -1.63e-04 -1.66e-04 -1.32e-05 -4.05e-05 -2.53e-04 -6.31e-22 1.51e-04 9.18e-05 -1.20e-04 -2.16e-04 1.01e-04 -1.17e-04 -6.89e-04 -1.51e-04 -3.97e-21 1.91e-04 -4.48e-05 -5.47e-04 7.90e-05 -9.94e-06 -1.77e-04 -9.18e-05 -1.91e-04 -1.66e-21 7.86e-05 -7.75e-05 -5.20e-05 1.48e-04 1.63e-04 1.20e-04 4.48e-05 -7.86e-05 -8.18e-22 -3.67e-05 -1.19e-04 2.20e-04 1.66e-04 2.16e-04 5.47e-04 7.75e-05 3.67e-05 -1.09e-21 1.81e-05 3.38e-06 1.32e-05 -1.01e-04 -7.90e-05 5.20e-05 1.19e-04 -1.81e-05 5.24e-21 2.98e-05 4.05e-05 1.17e-04 9.94e-06 -1.48e-04 -2.20e-04 -3.38e-06 -2.98e-05 1.09e-21
    0.00e+00 8.97e-06 -1.18e-04 -5.31e-05 1.66e-05 -1.19e-04 7.78e-05 1.33e-04 -8.97e-06 -2.79e-21 -1.38e-04 2.85e-04 -1.37e-04 4.96e-04 -2.10e-04 -5.52e-04 1.18e-04 1.38e-04 -5.39e-21 6.58e-05 -1.74e-04 -4.21e-04 -5.36e-04 5.14e-04 5.31e-05 -2.85e-04 -6.58e-05 7.88e-21 -2.17e-05 1.57e-04 -1.01e-04 -1.50e-04 -1.66e-05 1.37e-04 1.74e-04 2.17e-05 -5.54e-22 -2.97e-05 1.94e-04 1.48e-05 1.19e-04 -4.96e-04 4.21e-04 -1.57e-04 2.97e-05 5.29e-23 -1.01e-04 2.28e-06 -7.78e-05 2.10e-04 5.36e-04 1.01e-04 -1.94e-04 1.01e-04 -1.87e-21 -8.31e-05 -1.33e-04 5.52e-04 -5.14e-04 1.50e-04 -1.48e-05 -2.28e-06 8.31e-05 -5.89e-21
    0.00e+00 -2.46e-04 -5.39e-05 -4.51e-04 8.68e-05 3.22e-04 -3.90e-05 6.95e-05 2.46e-04 3.65e-21 -2.05e-04 -3.32e-04 -5.71e-05 1.06e-04 1.74e-04 3.78e-04 5.39e-05 2.05e-04 -5.50e-22 1.90e-04 -4.76e-05 -2.17e-04 1.03e-04 1.31e-04 4.51e-04 3.32e-04 -1.90e-04 -1.32e-21 -4.07e-04 1.66e-04 3.09e-04 7.65e-04 -8.68e-05 5.71e-05 4.76e-05 4.07e-04 1.49e-22 -7.89e-05 -9.35e-05 -2.86e-04 -3.22e-04 -1.06e-04 2.17e-04 -1.66e-04 7.89e-05 6.07e-22 -1.06e-04 -3.02e-04 3.90e-05 -1.74e-04 -1.03e-04 -3.09e-04 9.35e-05 1.06e-04 -1.44e-21 1.69e-04 -6.95e-05 -3.78e-04 -1.31e-04 -7.65e-04 2.86e-04 3.02e-04 -1.69e-04 8.12e-22
    0.00e+00 2.07e-04 -2.45e-04 -1.11e-04 -9.92e-05 -2.55e-04 1.69e-04 -5.37e-05 -2.07e-04 -2.57e-21 -4.13e-05 1.17e-04 -1.28e-04 1.12e-05 -7.20e-04 3.56e-04 2.45e-04 4.13e-05 6.61e-21 -1.23e-04 -2.50e-04 1.79e-04 -2.82e-04 2.42e-05 1.11e-04 -1.17e-04 1.23e-04 1.54e-21 1.15e-04 1.44e-04 1.59e-04 -1.47e-04 9.92e-05 1.28e-04 2.50e-04 -1.15e-04 -2.33e-22 1.04e-04 -3.18e-04 -6.73e-05 2.55e-04 -1.12e-05 -1.79e-04 -1.44e-04 -1.04e-04 -3.13e-21 6.46e-04 -2.84e-04 -1.69e-04 7.20e-04 2.82e-04 -1.59e-04 3.18e-04 -6.46e-04 -9.66e-22 -4.40e-04 5.37e-05 -3.56e-04 -2.42e-05 1.47e-04 6.73e-05 2.84e-04 4.40e-04 -9.29e-22
    0.00e+00 -1.34e-04 -5.40e-05 -9.63e-05 3.07e-04 1.37e-04 8.68e-05 6.53e-05 1.34e-04 -6.80e-22 -1.11e-04 -5.25e-04 -9.39e-06 -5.18e-05 1.41e-04 -2.62e-04 5.40e-05 1.11e-04 -4.12e-22 -2.65e-04 -1.46e-04 -4.47e-05 1.37e-04 -3.03e-04 9.63e-05 5.25e-04 2.65e-04 -3.32e-23 3.92e-04 5.94e-06 1.77e-04 3.92e-04 -3.07e-04 9.39e-06 1.46e-04 -3.92e-04 5.62e-22 1.73e-04 -1.56e-04 -3.08e-04 -1.37e-04 5.18e-05 4.47e-05 -5.94e-06 -1.73e-04 3.15e-23 -4.45e-05 -3.33e-06 -8.68e-05 -1.41e-04 -1.37e-04 -1.77e-04 1.56e-04 4.45e-05 1.44e-21 -1.01e-04 -6.53e-05 2.62e-04 3.03e-04 -3.92e-04 3.08e-04 3.33e-06 1.01e-04 9.47e-22
    0.00e+00 4.24e-04 1.03e-04 9.75e-05 2.75e-05 -6.69e-04 1.87e-04 -9.47e-05 -4.24e-04 -2.30e-21 1.30e-04 1.84e-04 -6.53e-05 -5.67e-04 -2.31e-04 2.14e-04 -1.03e-04 -1.30e-04 1.36e-21 9.35e-05 -2.21e-04 1.87e-05 4.46e-05 -3.59e-05 -9.75e-05 -1.84e-04 -9.35e-05 -2.99e-21 -9.87e-05 2.67e-06 6.51e-05 5.75e-05 -2.75e-05 6.53e-05 2.21e-04 9.87e-05 9.77e-22 -7.02e-04 1.54e-04 1.58e-04 6.69e-04 5.67e-04 -1.87e-05 -2.67e-06 7.02e-04 -3.23e-21 7.35e-05 -1.75e-04 -1.87e-04 2.31e-04 -4.46e-05 -6.51e-05 -1.54e-04 -7.35e-05 -1.87e-21 -3.66e-04 9.47e-05 -2.14e-04 3.59e-05 -5.75e-05 -1.58e-04 1.75e-04 3.66e-04 1.07e-21
    0.00e+00 -2.67e-04 3.68e-04 1.09e-04 3.44e-05 4.36e-04 -5.10e-05 -5.23e-04 2.67e-04 -1.98e-22 3.50e-05 -3.21e-04 1.15e-05 -2.05e-04 1.30e-04 8.89e-05 -3.68e-04 -3.50e-05 9.36e-22 2.72e-04 1.86e-04 3.78e-05 -1.28e-05 -2.38e-04 -1.09e-04 3.21e-04 -2.72e-04 8.39e-22 6.02e-07 -3.75e-04 1.14e-04 5.27e-04 -3.44e-05 -1.15e-05 -1.86e-04 -6.02e-07 -2.09e-22 -2.74e-05 -8.22e-05 9.97e-05 -4.36e-04 2.05e-04 -3.78e-05 3.75e-04 2.74e-05 -2.39e-22 -2.57e-05 -1.92e-05 5.10e-05 -1.30e-04 1.28e-05 -1.14e-04 8.22e-05 2.57e-05 2.26e-22 -4.02e-05 5.23e-04 -8.89e-05 2.38e-04 -5.27e-04 -9.97e-05 1.92e-05 4.02e-05 -4.59e-22
    0.00e+00 1.86e-04 2.25e-04 1.54e-04 -2.71e-04 1.40e-04 -1.18e-05 -1.07e-04 -1.86e-04 -1.86e-21 7.04e-05 -8.25e-05 2.48e-04 5.83e-04 7.79e-04 3.15e-04 -2.25e-04 -7.04e-05 -4.58e-22 -1.94e-04 2.15e-04 1.41e-04 1.34e-04 3.21e-04 -1.54e-04 8.25e-05 1.94e-04 -8.17e-22 -8.52e-05 4.26e-04 3.19e-04 1.56e-04 2.71e-04 -2.48e-04 -2.15e-04 8.52e-05 3.15e-21 -2.61e-04 -1.07e-04 -1.06e-04 -1.40e-04 -5.83e-04 -1.41e-04 -4.26e-04 2.61e-04 1.34e-21 -3.72e-04 2.30e-04 1.18e-05 -7.79e-04 -1.34e-04 -3.19e-04 1.07e-04 3.72e-04 -2.28e-21 1.37e-04 1.07e-04 -3.15e-04 -3.21e-04 -1.56e-04 1.06e-04 -2.30e-04 -1.37e-04 -4.47e-22
    0.00e+00 -3.73e-05 -3.84e-04 -4.34e-05 -9.03e-05 8.58e-05 9.71e-05 -2.70e-04 3.73e-05 8.39e-22 -3.93e-04 -1.92e-04 -6.92e-05 2.02e-04 -6.91e-05 -3.06e-04 3.84e-04 3.93e-04 1.59e-21 3.00e-04 9.21e-05 1.02e-04 -2.36e-04 -6.60e-05 4.34e-05 1.92e-04 -3.00e-04 5.30e-22 -1.30e-04 -5.22e-05 -5.69e-07 -1.76e-04 9.03e-05 6.92e-05 -9.21e-05 1.30e-04 1.57e-22 -1.77e-05 2.01e-05 -1.41e-04 -8.58e-05 -2.02e-04 -1.02e-04 5.22e-05 1.77e-05 1.09e-21 2.05e-04 -7.21e-05 -9.71e-05 6.91e-05 2.36e-04 5.69e-07 -2.01e-05 -2.05e-04 -8.32e-22 3.63e-04 2.70e-04 3.06e-04 6.60e-05 1.76e-04 1.41e-04 7.21e-05 -3.63e-04 2.63e-22
    0.00e+00 3.04e-04 4.82e-04 2.41e-04 -4.69e-04 -2.89e-04 -1.40e-04 -1.51e-04 -3.04e-04 6.43e-22 1.77e-04 1.72e-04 2.30e-04 -3.04e-04 2.10e-04 4.18e-05 -4.82e-04 -1.77e-04 -5.58e-22 1.69e-04 -2.32e-04 -6.52e-04 9.97e-05 1.35e-04 -2.41e-04 -1.72e-04 -1.69e-04 -1.26e-21 -4.09e-05 -2.50e-04 7.74e-05 2.28e-04 4.69e-04 -2.30e-04 2.32e-04 4.09e-05 2.30e-22 -1.34e-04 -1.75e-04 6.80e-05 2.89e-04 3.04e-04 6.52e-04 2.50e-04 1.34e-04 -2.22e-21 -1.74e-04 -2.48e-04 1.40e-04 -2.10e-04 -9.97e-05 -7.74e-05 1.75e-04 1.74e-04 1.03e-23 4.63e-05 1.51e-04 -4.18e-05 -1.35e-04 -2.28e-04 -6.80e-05 2.48e-04 -4.63e-05 4.55e-22
    0.00e+00 -2.53e-04 1.86e-04 -4.33e-04 1.60e-05 4.61e-04 1.18e-05 -1.27e-04 2.53e-04 9.49e-22 -2.30e-04 -3.56e-04 -4.38e-05 1.05e-04 2.00e-04 3.30e-04 -1.86e-04 2.30e-04 1.12e-21 2.13e-04 -1.56e-05 -2.97e-04 1.04e-04 6.51e-05 4.33e-04 3.56e-04 -2.13e-04 1.20e-21 -4.50e-04 1.15e-04 2.95e-04 8.19e-04 -1.60e-05 4.38e-05 1.56e-05 4.50e-04 -3.40e-21 -2.52e-05 -7.57e-05 -2.32e-04 -4.61e-04 -1.05e-04 2.97e-04 -1.15e-04 2.52e-05 -3.12e-21 -2.95e-05 -2.88e-04 -1.18e-05 -2.00e-04 -1.04e-04 -2.95e-04 7.57e-05 2.95e-05 3.59e-21 7.15e-05 1.27e-04 -3.30e-04 -6.51e-05 -8.19e-04 2.32e-04 2.88e-04 -7.15e-05 -1.61e-21
    0.00e+00 -7.38e-05 -1.33e-04 -6.76e-06 9.89e-05 -4.86e-06 1.49e-04 6.60e-05 7.38e-05 -6.78e-23 -1.20e-04 -1.66e-04 -1.22e-04 1.03e-05 -2.62e-05 -2.76e-07 1.33e-04 1.20e-04 6.31e-23 -3.44e-05 1.53e-04 8.17e-05 7.89e-06 -5.98e-04 6.76e-06 1.66e-04 3.44e-05 1.18e-21 1.72e-04 1.48e-04 -1.56e-04 -6.15e-04 -9.89e-05 1.22e-04 -1.53e-04 -1.72e-04 -3.12e-21 3.87e-04 4.80e-05 -6.12e-04 4.86e-06 -1.03e-05 -8.17e-05 -1.48e-04 -3.87e-04 2.42e-22 8.07e-05 5.76e-04 -1.49e-04 2.62e-05 -7.89e-06 1.56e-04 -4.80e-05 -8.07e-05 7.95e-23 5.83e-04 -6.60e-05 2.76e-07 5.98e-04 6.15e-04 6.12e-04 -5.76e-04 -5.83e-04 7.95e-23
    0.00e+00 2.74e-04 -7.37e-05 -2.06e-04 -6.08e-05 7.70e-05 1.20e-04 4.28e-04 -2.74e-04 7.09e-22 -1.04e-04 1.65e-04 1.95e-04 4.63e-04 4.41e-04 -2.36e-04 7.37e-05 1.04e-04 -1.66e-21 -3.77e-04 2.74e-05 -1.26e-04 3.54e-04 9.75e-04 2.06e-04 -1.65e-04 3.77e-04 8.37e-22 -1.09e-04 -5.23e-05 8.27e-05 5.94e-05 6.08e-05 -1.95e-04 -2.74e-05 1.09e-04 2.13e-21 -5.02e-06 -1.20e-05 -2.62e-04 -7.70e-05 -4.63e-04 1.26e-04 5.23e-05 5.02e-06 1.76e-21 1.38e-04 -1.83e-04 -1.20e-04 -4.41e-04 -3.54e-04 -8.27e-05 1.20e-05 -1.38e-04 6.07e-21 1.93e-06 -4.28e-04 2.36e-04 -9.75e-04 -5.94e-05 2.62e-04 1.83e-04 -1.93e-06 -1.19e-21
    0.00e+00 1.85e-04 -3.74e-05 1.21e-04 -1.59e-04 1.34e-04 -3.56e-05 -2.45e-04 -1.85e-04 1.44e-21 1.74e-04 1.58e-04 1.90e-04 7.16e-04 8.30e-04 1.29e-04 3.74e-05 -1.74e-04 -3.93e-22 -1.86e-04 1.71e-04 -2.13e-04 1.21e-04 4.83e-04 -1.21e-04 -1.58e-04 1.86e-04 -1.57e-22 -2.08e-05 3.47e-04 2.29e-04 4.21e-05 1.59e-04 -1.90e-04 -1.71e-04 2.08e-05 1.46e-21 -1.41e-04 -2.63e-05 3.57e-05 -1.34e-04 -7.16e-04 2.13e-04 -3.47e-04 1.41e-04 2.97e-21 -2.10e-04 3.38e-04 3.56e-05 -8.30e-04 -1.21e-04 -2.29e-04 2.63e-05 2.10e-04 2.05e-21 1.56e-04 2.45e-04 -1.29e-04 -4.83e-04 -4.21e-05 -3.57e-05 -3.38e-04 -1.56e-04 -2.16e-22
    0.00e+00 -2.30e-05 1.70e-04 1.84e-04 -3.49e-04 -6.47e-05 -2.79e-04 -5.11e-05 2.30e-05 1.15e-21 2.11e-04 3.38e-04 3.64e-05 -4.93e-04 2.24e-04 -3.50e-05 -1.70e-04 -2.11e-04 -4.72e-22 6.26e-06 -9.44e-05 -2.66e-04 8.54e-05 1.70e-04 -1.84e-04 -3.38e-04 -6.26e-06 -1.49e-21 -1.81e-04 -4.73e-04 1.58e-04 2.19e-04 3.49e-04 -3.64e-05 9.44e-05 1.81e-04 -7.17e-23 -1.67e-05 1.02e-05 -1.14e-04 6.47e-05 4.93e-04 2.66e-04 4.73e-04 1.67e-05 2.61e-22 -4.11e-04 -2.25e-04 2.79e-04 -2.24e-04 -8.54e-05 -1.58e-04 -1.02e-05 4.11e-04 -5.80e-22 -8.72e-05 5.11e-05 3.50e-05 -1.70e-04 -2.19e-04 1.14e-04 2.25e-04 8.72e-05 2.96e-21
    0.00e+00 -1.93e-04 -2.39e-04 -4.12e-04 7.96e-05 3.70e-04 -3.70e-05 2.89e-04 1.93e-04 -8.68e-22 -1.65e-04 -3.30e-04 -8.82e-05 2.34e-04 1.57e-04 3.49e-04 2.39e-04 1.65e-04 2.88e-22 1.76e-04 -7.54e-05 -6.46e-05 1.30e-04 1.78e-04 4.12e-04 3.30e-04 -1.76e-04 4.64e-22 -4.32e-04 1.62e-04 2.67e-04 6.11e-04 -7.96e-05 8.82e-05 7.54e-05 4.32e-04 -4.31e-22 -1.08e-04 -3.47e-05 -4.34e-04 -3.70e-04 -2.34e-04 6.46e-05 -1.62e-04 1.08e-04 -2.53e-22 -1.09e-04 -2.30e-04 3.70e-05 -1.57e-04 -1.30e-04 -2.67e-04 3.47e-05 1.09e-04 7.61e-22 2.40e-04 -2.89e-04 -3.49e-04 -1.78e-04 -6.11e-04 4.34e-04 2.30e-04 -2.40e-04 1.32e-21
    0.00e+00 2.49e-04 4.23e-04 3.93e-04 -4.80e-04 -7.80e-06 -3.59e-04 -1.48e-04 -2.49e-04 4.61e-21 2.48e-04 2.33e-04 2.73e-04 -2.11e-04 2.18e-04 2.53e-05 -4.23e-04 -2.48e-04 -7.36e-22 1.07e-04 -1.87e-04 -6.69e-04 2.70e-05 1.94e-04 -3.93e-04 -2.33e-04 -1.07e-04 3.93e-21 -5.78e-05 -5.28e-04 7.33e-05 2.12e-04 4.80e-04 -2.73e-04 1.87e-04 5.78e-05 -9.48e-23 -2.68e-05 -8.67e-05 2.15e-05 7.80e-06 2.11e-04 6.69e-04 5.28e-04 2.68e-05 -2.03e-21 -2.88e-04 -7.31e-05 3.59e-04 -2.18e-04 -2.70e-05 -7.33e-05 8.67e-05 2.88e-04 3.91e-21 7.19e-05 1.48e-04 -2.53e-05 -1.94e-04 -2.12e-04 -2.15e-05 7.31e-05 -7.19e-05 -1.39e-22
    0.00e+00 1.71e-04 2.52e-04 1.40e-04 -2.58e-04 2.42e-04 -5.40e-05 -2.04e-04 -1.71e-04 -1.31e-22 1.00e-04 1.85e-04 2.01e-04 6.22e-04 7.29e-04 2.17e-04 -2.52e-04 -1.00e-04 -1.73e-22 -1.59e-04 2.49e-04 3.44e-05 1.40e-04 3.54e-04 -1.40e-04 -1.85e-04 1.59e-04 -5.67e-22 -6.16e-05 4.23e-04 1.68e-04 1.09e-04 2.58e-04 -2.01e-04 -2.49e-04 6.16e-05 1.29e-21 -2.39e-04 -5.59e-05 -6.23e-05 -2.42e-04 -6.22e-04 -3.44e-05 -4.23e-04 2.39e-04 1.18e-21 -3.30e-04 4.14e-04 5.40e-05 -7.29e-04 -1.40e-04 -1.68e-04 5.59e-05 3.30e-04 -3.18e-21 2.11e-04 2.04e-04 -2.17e-04 -3.54e-04 -1.09e-04 6.23e-05 -4.14e-04 -2.11e-04 6.22e-23
    0.00e+00 -1.20e-04 -1.40e-04 1.06e-05 1.51e-04 1.26e-04 4.84e-05 5.78e-06 1.20e-04 6.00e-22 -1.32e-04 -5.47e-04 2.60e-06 2.89e-06 1.75e-04 -4.29e-04 1.40e-04 1.32e-04 7.00e-22 -3.99e-04 -2.98e-05 -1.45e-04 2.15e-04 -3.33e-04 -1.06e-05 5.47e-04 3.99e-04 -1.74e-22 3.67e-04 -1.01e-04 1.05e-04 -7.50e-06 -1.51e-04 -2.60e-06 2.98e-05 -3.67e-04 -1.02e-21 1.54e-04 -4.81e-05 -6.18e-04 -1.26e-04 -2.89e-06 1.45e-04 1.01e-04 -1.54e-04 -9.37e-22 -6.80e-05 1.38e-04 -4.84e-05 -1.75e-04 -2.15e-04 -1.05e-04 4.81e-05 6.80e-05 -2.32e-21 -1.50e-04 -5.78e-06 4.29e-04 3.33e-04 7.50e-06 6.18e-04 -1.38e-04 1.50e-04 7.14e-22
    0.00e+00 3.00e-04 -2.45e-04 -1.43e-05 -5.23e-05 6.51e-05 1.06e-04 1.04e-04 -3.00e-04 -8.90e-22 8.38e-05 1.04e-04 2.22e-04 6.86e-04 8.19e-04 -1.46e-04 2.45e-04 -8.38e-05 -3.96e-21 -2.82e-04 4.32e-05 -2.82e-04 2.56e-04 8.25e-04 1.43e-05 -1.04e-04 2.82e-04 -1.01e-21 -1.34e-05 1.62e-04 1.94e-04 3.02e-05 5.23e-05 -2.22e-04 -4.32e-05 1.34e-05 -3.58e-22 -3.63e-05 -2.33e-05 -6.28e-05 -6.51e-05 -6.86e-04 2.82e-04 -1.62e-04 3.63e-05 3.58e-23 3.36e-05 2.01e-05 -1.06e-04 -8.19e-04 -2.56e-04 -1.94e-04 2.33e-05 -3.36e-05 1.41e-21 2.01e-05 -1.04e-04 1.46e-04 -8.25e-04 -3.02e-05 6.28e-05 -2.01e-05 -2.01e-05 -1.76e-21
    0.00e+00 7.50e-05 -3.27e-04 -8.16e-05 -6.21e-05 -5.56e-05 1.28e-04 -2.10e-04 -7.50e-05 2.45e-22 -3.93e-04 -8.05e-05 -8.69e-05 1.91e-04 -4.62e-04 -7.83e-05 3.27e-04 3.93e-04 3.67e-21 -8.27e-06 -1.19e-04 2.04e-04 -1.37e-04 1.65e-04 8.16e-05 8.05e-05 8.27e-06 -8.84e-22 -6.59e-05 -1.40e-04 3.55e-04 -2.40e-04 6.21e-05 8.69e-05 1.19e-04 6.59e-05 -1.34e-22 -4.13e-05 1.92e-04 -7.03e-05 5.56e-05 -1.91e-04 -2.04e-04 1.40e-04 4.13e-05 -8.64e-22 1.54e-04 -1.24e-04 -1.28e-04 4.62e-04 1.37e-04 -3.55e-04 -1.92e-04 -1.54e-04 4.88e-23 -4.85e-06 2.10e-04 7.83e-05 -1.65e-04 2.40e-04 7.03e-05 1.24e-04 4.85e-06 -1.04e-22
    0.00e+00 1.40e-04 1.89e-04 1.82e-04 -4.57e-04 -7.18e-05 -1.48e-04 -5.73e-05 -1.40e-04 -1.29e-21 6.04e-05 9.91e-05 1.80e-04 -2.81e-04 3.72e-04 1.07e-04 -1.89e-04 -6.04e-05 -4.62e-22 8.95e-05 -1.76e-04 -1.19e-04 1.81e-04 1.94e-04 -1.82e-04 -9.91e-05 -8.95e-05 1.51e-21 -1.03e-04 -3.55e-05 2.00e-04 2.60e-04 4.57e-04 -1.80e-04 1.76e-04 1.03e-04 6.22e-22 -2.30e-04 -1.09e-04 -1.15e-04 7.18e-05 2.81e-04 1.19e-04 3.55e-05 2.30e-04 1.63e-22 -3.36e-04 -2.26e-04 1.48e-04 -3.72e-04 -1.81e-04 -2.00e-04 1.09e-04 3.36e-04 -3.05e-22 6.49e-05 5.73e-05 -1.07e-04 -1.94e-04 -2.60e-04 1.15e-04 2.26e-04 -6.49e-05 -1.15e-21
    0.00e+00 2.70e-04 7.29e-05 -2.12e-04 -3.22e-05 1.00e-04 -1.41e-05 2.95e-04 -2.70e-04 -2.27e-21 -1.91e-04 1.09e-05 1.49e-04 3.06e-04 3.54e-04 -4.70e-05 -7.29e-05 1.91e-04 1.55e-21 -3.85e-04 2.79e-05 -4.63e-05 1.60e-04 8.93e-04 2.12e-04 -1.09e-05 3.85e-04 8.99e-22 -8.44e-05 -2.99e-05 1.03e-04 1.07e-04 3.22e-05 -1.49e-04 -2.79e-05 8.44e-05 -3.32e-22 -5.68e-05 4.59e-06 -2.30e-04 -1.00e-04 -3.06e-04 4.63e-05 2.99e-05 5.68e-05 -2.30e-21 6.82e-05 -6.50e-05 1.41e-05 -3.54e-04 -1.60e-04 -1.03e-04 -4.59e-06 -6.82e-05 5.85e-22 6.90e-05 -2.95e-04 4.70e-05 -8.93e-04 -1.07e-04 2.30e-04 6.50e-05 -6.90e-05 2.38e-22
    0.00e+00 -1.02e-04 -1.92e-04 -2.40e-04 -7.79e-05 1.29e-04 7.43e-05 -1.73e-04 1.02e-04 -4.39e-22 -1.62e-04 -4.70e-04 -8.76e-05 1.17e-04 1.82e-04 -5.51e-04 1.92e-04 1.62e-04 4.71e-22 -1.54e-04 1.68e-04 5.66e-05 1.44e-04 -2.10e-04 2.40e-04 4.70e-04 1.54e-04 -3.25e-22 4.02e-05 -1.32e-04 -4.54e-05 -2.13e-04 7.79e-05 8.76e-05 -1.68e-04 -4.02e-05 -2.83e-22 1.19e-04 5.26e-05 -5.64e-04 -1.29e-04 -1.17e-04 -5.66e-05 1.32e-04 -1.19e-04 -9.96e-22 -1.35e-05 1.32e-04 -7.43e-05 -1.82e-04 -1.44e-04 4.54e-05 -5.26e-05 1.35e-05 -8.28e-23 9.68e-05 1.73e-04 5.51e-04 2.10e-04 2.13e-04 5.64e-04 -1.32e-04 -9.68e-05 6.34e-22
    0.00e+00 2.02e-04 -3.81e-04 -1.49e-04 -1.82e-04 -2.31e-04 3.01e-04 -2.53e-04 -2.02e-04 -2.00e-22 -2.11e-04 5.62e-05 -1.41e-04 2.29e-04 -6.52e-04 2.42e-04 3.81e-04 2.11e-04 5.99e-22 -8.38e-05 -4.05e-04 1.02e-04 1.93e-04 5.33e-05 1.49e-04 -5.62e-05 8.38e-05 -3.60e-22 8.06e-05 -1.06e-04 4.33e-04 -2.64e-04 1.82e-04 1.41e-04 4.05e-04 -8.06e-05 -1.39e-21 -1.11e-04 2.85e-04 -6.14e-05 2.31e-04 -2.29e-04 -1.02e-04 1.06e-04 1.11e-04 -8.03e-22 4.09e-04 -1.95e-04 -3.01e-04 6.52e-04 -1.93e-04 -4.33e-04 -2.85e-04 -4.09e-04 3.32e-22 -2.63e-04 2.53e-04 -2.42e-04 -5.33e-05 2.64e-04 6.14e-05 1.95e-04 2.63e-04 -3.17e-22
    0.00e+00 4.42e-05 -1.13e-04 -7.44e-05 -3.59e-05 -4.90e-05 2.08e-04 -1.40e-04 -4.42e-05 -1.90e-22 -1.72e-04 -1.83e-04 -1.85e-04 -1.84e-06 -4.48e-04 2.86e-05 1.13e-04 1.72e-04 -1.94e-21 -2.61e-04 -2.14e-04 1.92e-04 -5.50e-05 -3.35e-04 7.44e-05 1.83e-04 2.61e-04 -1.15e-21 5.51e-05 -9.31e-05 2.14e-04 -2.76e-04 3.59e-05 1.85e-04 2.14e-04 -5.51e-05 8.08e-22 -5.79e-05 1.14e-04 -3.73e-04 4.90e-05 1.84e-06 -1.92e-04 9.31e-05 5.79e-05 -1.48e-21 5.68e-04 7.68e-05 -2.08e-04 4.48e-04 5.50e-05 -2.14e-04 -1.14e-04 -5.68e-04 2.84e-22 -4.03e-04 1.40e-04 -2.86e-05 3.35e-04 2.76e-04 3.73e-04 -7.68e-05 4.03e-04 6.33e-22
    0.00e+00 -3.85e-04 2.29e-04 -5.19e-05 -2.90e-04 6.68e-04 -1.22e-04 -4.03e-04 3.85e-04 8.14e-22 -6.72e-05 -2.04e-04 -7.21e-05 3.93e-05 1.31e-04 1.51e-04 -2.29e-04 6.72e-05 1.63e-21 4.91e-04 -4.52e-05 1.69e-05 2.57e-06 -6.33e-05 5.19e-05 2.04e-04 -4.91e-04 4.19e-22 -2.54e-04 -2.16e-04 1.16e-04 7.20e-04 2.90e-04 7.21e-05 4.52e-05 2.54e-04 2.33e-21 4.20e-05 3.30e-05 -7.71e-06 -6.68e-04 -3.93e-05 -1.69e-05 2.16e-04 -4.20e-05 6.22e-22 -1.14e-04 -6.68e-05 1.22e-04 -1.31e-04 -2.57e-06 -1.16e-04 -3.30e-05 1.14e-04 -1.63e-21 -9.26e-07 4.03e-04 -1.51e-04 6.33e-05 -7.20e-04 7.71e-06 6.68e-05 9.26e-07 -8.49e-22
    0.00e+00 -7.14e-05 -6.17e-05 -3.88e-04 1.62e-05 5.22e-04 -2.58e-05 -3.51e-05 7.14e-05 6.21e-22 -1.76e-04 -5.35e-04 1.18e-04 2.39e-04 2.04e-04 2.38e-04 6.17e-05 1.76e-04 -1.42e-21 9.05e-05 -4.74e-06 -3.24e-04 1.07e-04 2.75e-04 3.88e-04 5.35e-04 -9.05e-05 -1.49e-21 -4.53e-04 2.05e-05 1.64e-04 7.02e-04 -1.62e-05 -1.18e-04 4.74e-06 4.53e-04 6.33e-22 -3.63e-05 -1.03e-04 -3.52e-04 -5.22e-04 -2.39e-04 3.24e-04 -2.05e-05 3.63e-05 2.20e-23 -6.99e-05 -4.87e-04 2.58e-05 -2.04e-04 -1.07e-04 -1.64e-04 1.03e-04 6.99e-05 2.85e-22 3.24e-05 3.51e-05 -2.38e-04 -2.75e-04 -7.02e-04 3.52e-04 4.87e-04 -3.24e-05 -2.14e-21
    0.00e+00 2.42e-04 -3.41e-04 -1.65e-04 -2.05e-04 -3.45e-04 3.29e-04 -2.40e-04 -2.42e-04 1.86e-22 -9.61e-05 2.13e-04 -1.24e-04 1.99e-04 -6.22e-04 3.98e-04 3.41e-04 9.61e-05 2.68e-21 -1.55e-04 -4.95e-04 6.25e-05 1.03e-04 1.34e-05 1.65e-04 -2.13e-04 1.55e-04 -2.45e-21 2.33e-04 4.86e-05 2.97e-04 -2.37e-04 2.05e-04 1.24e-04 4.95e-04 -2.33e-04 6.55e-22 -5.79e-05 1.01e-04 -5.20e-05 3.45e-04 -1.99e-04 -6.25e-05 -4.86e-05 5.79e-05 2.10e-21 5.21e-04 -2.92e-04 -3.29e-04 6.22e-04 -1.03e-04 -2.97e-04 -1.01e-04 -5.21e-04 3.38e-21 -2.06e-04 2.40e-04 -3.98e-04 -1.34e-05 2.37e-04 5.20e-05 2.92e-04 2.06e-04 8.19e-22
    0.00e+00 -1.41e-04 1.88e-04 1.42e-04 5.42e-04 2.22e-04 1.28e-04 4.24e-04 1.41e-04 -7.41e-22 -7.55e-05 -3.44e-04 -3.75e-04 -2.28e-04 1.17e-04 -6.11e-05 -1.88e-04 7.55e-05 1.77e-21 1.24e-04 1.02e-04 5.37e-05 -3.35e-05 -2.48e-04 -1.42e-04 3.44e-04 -1.24e-04 4.22e-21 3.00e-05 3.32e-04 3.12e-04 -2.81e-05 -5.42e-04 3.75e-04 -1.02e-04 -3.00e-05 4.92e-21 7.03e-04 2.16e-04 -2.11e-05 -2.22e-04 2.28e-04 -5.37e-05 -3.32e-04 -7.03e-04 -2.28e-21 7.10e-05 -2.19e-04 -1.28e-04 -1.17e-04 3.35e-05 -3.12e-04 -2.16e-04 -7.10e-05 2.37e-21 2.01e-04 -4.24e-04 6.11e-05 2.48e-04 2.81e-05 2.11e-05 2.19e-04 -2.01e-04 1.91e-21
    0.00e+00 2.53e-04 -5.21e-05 3.67e-05 6.58e-05 -2.31e-04 2.04e-04 1.26e-04 -2.53e-04 2.27e-21 -4.51e-05 5.49e-04 -4.16e-04 1.75e-04 -7.07e-05 -8.03e-05 5.21e-05 4.51e-05 1.06e-21 3.07e-05 -8.89e-05 1.59e-04 -2.58e-04 1.84e-04 -3.67e-05 -5.49e-04 -3.07e-05 -2.36e-22 -6.63e-05 -1.59e-04 5.15e-06 -2.50e-04 -6.58e-05 4.16e-04 8.89e-05 6.63e-05 3.64e-22 5.50e-05 1.72e-04 2.70e-05 2.31e-04 -1.75e-04 -1.59e-04 1.59e-04 -5.50e-05 1.12e-21 -4.15e-04 9.27e-05 -2.04e-04 7.07e-05 2.58e-04 -5.15e-06 -1.72e-04 4.15e-04 3.12e-21 -2.68e-04 -1.26e-04 8.03e-05 -1.84e-04 2.50e-04 -2.70e-05 -9.27e-05 2.68e-04 -2.04e-21
    0.00e+00 -1.27e-04 1.39e-04 1.53e-04 1.25e-04 -5.97e-04 -5.11e-05 2.70e-04 1.27e-04 -3.85e-21 1.41e-04 1.05e-04 2.03e-04 -4.45e-04 -1.39e-05 -8.40e-05 -1.39e-04 -1.41e-04 7.68e-23 6.64e-05 -9.40e-05 1.97e-05 4.19e-05 -4.86e-05 -1.53e-04 -1.05e-04 -6.64e-05 4.67e-22 -7.41e-05 1.56e-04 1.21e-04 3.52e-05 -1.25e-04 -2.03e-04 9.40e-05 7.41e-05 1.45e-21 -4.32e-04 1.49e-05 5.82e-04 5.97e-04 4.45e-04 -1.97e-05 -1.56e-04 4.32e-04 -3.64e-22 -4.90e-04 -5.19e-04 5.11e-05 1.39e-05 -4.19e-05 -1.21e-04 -1.49e-05 4.90e-04 -7.51e-22 -2.96e-05 -2.70e-04 8.40e-05 4.86e-05 -3.52e-05 -5.82e-04 5.19e-04 2.96e-05 -9.69e-22
    0.00e+00 1.88e-04 2.24e-04 1.53e-04 -2.74e-04 1.66e-04 -1.44e-05 -2.28e-05 -1.88e-04 -1.99e-21 -8.97e-06 -8.38e-05 2.18e-04 4.58e-04 6.79e-04 2.90e-04 -2.24e-04 8.97e-06 -1.36e-21 -1.25e-04 1.76e-04 2.71e-04 1.61e-04 2.60e-04 -1.53e-04 8.38e-05 1.25e-04 4.91e-21 -9.08e-05 3.36e-04 3.07e-04 2.27e-04 2.74e-04 -2.18e-04 -1.76e-04 9.08e-05 -9.94e-22 -2.78e-04 -1.45e-04 -1.87e-04 -1.66e-04 -4.58e-04 -2.71e-04 -3.36e-04 2.78e-04 1.67e-22 -3.99e-04 1.25e-04 1.44e-05 -6.79e-04 -1.61e-04 -3.07e-04 1.45e-04 3.99e-04 -1.63e-21 1.04e-04 2.28e-05 -2.90e-04 -2.60e-04 -2.27e-04 1.87e-04 -1.25e-04 -1.04e-04 1.28e-23
    0.00e+00 -2.11e-04 1.54e-04 -1.66e-04 6.06e-04 1.04e-04 1.55e-04 5.42e-04 2.11e-04 2.15e-21 -1.99e-04 -4.07e-04 -2.12e-04 -1.85e-04 2.41e-04 -1.06e-04 -1.54e-04 1.99e-04 7.19e-22 2.12e-05 -1.82e-04 -1.38e-05 -1.24e-05 -2.93e-04 1.66e-04 4.07e-04 -2.12e-05 -2.26e-21 -2.31e-04 7.80e-05 3.65e-04 -1.39e-04 -6.06e-04 2.12e-04 1.82e-04 2.31e-04 -5.76e-22 3.37e-04 -2.57e-04 -1.66e-04 -1.04e-04 1.85e-04 1.38e-05 -7.80e-05 -3.37e-04 7.19e-22 4.56e-05 -1.89e-04 -1.55e-04 -2.41e-04 1.24e-05 -3.65e-04 2.57e-04 -4.56e-05 -2.46e-22 2.57e-04 -5.42e-04 1.06e-04 2.93e-04 1.39e-04 1.66e-04 1.89e-04 -2.57e-04 3.72e-21
    0.00e+00 1.54e-04 -1.16e-04 2.86e-05 5.58e-05 -1.91e-04 3.55e-04 4.04e-05 -1.54e-04 1.15e-21 -7.03e-05 6.71e-04 -3.94e-04 -5.03e-05 -3.27e-04 -1.85e-04 1.16e-04 7.03e-05 6.83e-22 -4.06e-06 -4.12e-06 1.00e-04 -2.12e-04 1.17e-04 -2.86e-05 -6.71e-04 4.06e-06 -3.49e-21 -1.44e-04 -1.82e-04 -8.85e-05 -2.19e-04 -5.58e-05 3.94e-04 4.12e-06 1.44e-04 -4.00e-22 4.12e-05 1.61e-04 3.07e-05 1.91e-04 5.03e-05 -1.00e-04 1.82e-04 -4.12e-05 -1.02e-21 -4.06e-04 9.01e-05 -3.55e-04 3.27e-04 2.12e-04 8.85e-05 -1.61e-04 4.06e-04 8.19e-22 -2.65e-04 -4.04e-05 1.85e-04 -1.17e-04 2.19e-04 -3.07e-05 -9.01e-05 2.65e-04 3.88e-22
    0.00e+00 1.63e-04 3.38e-04 3.37e-04 -4.42e-04 2.94e-05 -4.26e-04 -7.28e-05 -1.63e-04 1.96e-21 3.71e-04 2.82e-04 2.01e-04 -3.64e-04 2.35e-04 3.16e-05 -3.38e-04 -3.71e-04 4.85e-22 7.76e-05 -1.74e-04 -5.57e-04 5.06e-05 1.70e-04 -3.37e-04 -2.82e-04 -7.76e-05 2.35e-21 -1.08e-04 -5.55e-04 1.08e-04 1.70e-04 4.42e-04 -2.01e-04 1.74e-04 1.08e-04 2.08e-22 -4.00e-05 -4.83e-05 1.36e-05 -2.94e-05 3.64e-04 5.57e-04 5.55e-04 4.00e-05 -5.82e-22 -3.16e-04 -1.76e-04 4.26e-04 -2.35e-04 -5.06e-05 -1.08e-04 4.83e-05 3.16e-04 -5.89e-22 1.98e-05 7.28e-05 -3.16e-05 -1.70e-04 -1.70e-04 -1.36e-05 1.76e-04 -1.98e-05 6.12e-22
    0.00e+00 -1.44e-04 -9.85e-06 -8.67e-05 4.27e-04 1.46e-04 1.28e-04 2.05e-04 1.44e-04 1.12e-21 -1.17e-04 -5.05e-04 -8.55e-05 -9.60e-05 1.47e-04 -1.69e-04 9.85e-06 1.17e-04 2.34e-21 -2.04e-04 -2.06e-04 -5.54e-05 1.15e-04 -3.29e-04 8.67e-05 5.05e-04 2.04e-04 -5.46e-22 2.16e-04 6.54e-05 2.22e-04 3.58e-04 -4.27e-04 8.55e-05 2.06e-04 -2.16e-04 3.51e-21 1.92e-04 -2.28e-04 -1.37e-04 -1.46e-04 9.60e-05 5.54e-05 -6.54e-05 -1.92e-04 -1.26e-21 2.05e-06 -4.78e-05 -1.28e-04 -1.47e-04 -1.15e-04 -2.22e-04 2.28e-04 -2.05e-06 -6.32e-21 -4.46e-06 -2.05e-04 1.69e-04 3.29e-04 -3.58e-04 1.37e-04 4.78e-05 4.46e-06 7.35e-22
    0.00e+00 -1.01e-04 5.02e-04 -8.57e-06 -1.42e-05 -1.73e-04 1.73e-04 -2.43e-04 1.01e-04 -1.92e-22 1.44e-04 -1.51e-04 -1.88e-04 -2.13e-04 6.13e-05 1.84e-05 -5.02e-04 -1.44e-04 -2.20e-21 1.25e-04 1.37e-04 -2.02e-04 7.29e-05 -2.17e-04 8.57e-06 1.51e-04 -1.25e-04 -1.69e-21 1.72e-04 1.22e-04 8.84e-05 2.15e-04 1.42e-05 1.88e-04 -1.37e-04 -1.72e-04 1.67e-21 -1.40e-04 -1.55e-04 3.87e-04 1.73e-04 2.13e-04 2.02e-04 -1.22e-04 1.40e-04 -5.97e-22 2.46e-05 1.33e-04 -1.73e-04 -6.13e-05 -7.29e-05 -8.84e-05 1.55e-04 -2.46e-05 7.90e-22 -1.89e-05 2.43e-04 -1.84e-05 2.17e-04 -2.15e-04 -3.87e-04 -1.33e-04 1.89e-05 1.06e-21
    0.00e+00 -7.67e-05 3.79e-05 -1.57e-06 -1.44e-04 -1.17e-04 1.94e-04 -1.06e-04 7.67e-05 1.92e-22 1.26e-04 6.78e-04 -2.44e-04 -3.45e-04 -2.67e-04 -1.14e-04 -3.79e-05 -1.26e-04 1.67e-21 -1.06e-04 2.61e-06 -3.80e-05 -5.12e-05 1.27e-04 1.57e-06 -6.78e-04 1.06e-04 1.56e-21 -2.06e-04 -4.72e-04 -7.14e-05 5.11e-06 1.44e-04 2.44e-04 -2.61e-06 2.06e-04 1.42e-21 5.13e-05 1.17e-04 -5.81e-05 1.17e-04 3.45e-04 3.80e-05 4.72e-04 -5.13e-05 5.03e-22 -3.36e-04 -1.21e-04 -1.94e-04 2.67e-04 5.12e-05 7.14e-05 -1.17e-04 3.36e-04 4.58e-21 -2.78e-04 1.06e-04 1.14e-04 -1.27e-04 -5.11e-06 5.81e-05 1.21e-04 2.78e-04 -1.47e-23
    0.00e+00 -2.96e-05 -2.14e-04 -3.45e-04 5.45e-05 4.28e-04 -5.94e-05 2.42e-04 2.96e-05 -6.16e-22 -1.67e-04 -5.15e-04 8.26e-05 2.06e-04 2.23e-04 3.36e-04 2.14e-04 1.67e-04 1.55e-21 8.33e-05 -3.18e-05 -2.40e-04 1.68e-04 2.41e-04 3.45e-04 5.15e-04 -8.33e-05 4.90e-22 -4.80e-04 1.00e-04 1.67e-04 5.40e-04 -5.45e-05 -8.26e-05 3.18e-05 4.80e-04 -4.38e-22 -9.19e-05 -9.52e-05 -4.56e-04 -4.28e-04 -2.06e-04 2.40e-04 -1.00e-04 9.19e-05 1.09e-21 -1.02e-04 -4.38e-04 5.94e-05 -2.23e-04 -1.68e-04 -1.67e-04 9.52e-05 1.02e-04 6.25e-22 1.60e-04 -2.42e-04 -3.36e-04 -2.41e-04 -5.40e-04 4.56e-04 4.38e-04 -1.60e-04 -2.38e-22
    0.00e+00 -3.01e-04 3.57e-04 -3.69e-04 -1.09e-04 5.39e-04 1.32e-07 -2.46e-04 3.01e-04 1.38e-21 -1.68e-04 -3.40e-04 -7.37e-05 3.89e-05 1.90e-04 2.96e-04 -3.57e-04 1.68e-04 1.13e-21 3.07e-04 -2.53e-05 -2.19e-04 8.22e-05 -2.75e-06 3.69e-04 3.40e-04 -3.07e-04 -2.83e-22 -4.45e-04 4.72e-05 2.55e-04 8.38e-04 1.09e-04 7.37e-05 2.53e-05 4.45e-04 2.91e-23 2.18e-05 -4.83e-05 -1.71e-04 -5.39e-04 -3.89e-05 2.19e-04 -4.72e-05 -2.18e-05 -1.63e-21 -5.15e-05 -1.95e-04 -1.32e-07 -1.90e-04 -8.22e-05 -2.55e-04 4.83e-05 5.15e-05 -7.09e-22 -1.11e-05 2.46e-04 -2.96e-04 2.75e-06 -8.38e-04 1.71e-04 1.95e-04 1.11e-05 1.10e-22
    0.00e+00 1.86e-04 4.99e-05 -2.28e-04 -7.92e-05 1.33e-04 6.00e-05 3.53e-04 -1.86e-04 1.63e-21 -2.08e-04 -1.08e-04 1.71e-04 1.34e-04 3.18e-05 8.36e-05 -4.99e-05 2.08e-04 -2.93e-23 -2.90e-04 -3.85e-05 1.38e-04 1.14e-04 8.14e-04 2.28e-04 1.08e-04 2.90e-04 -2.33e-21 -1.69e-04 2.68e-05 -3.45e-06 6.29e-05 7.92e-05 -1.71e-04 3.85e-05 1.69e-04 -5.06e-21 1.32e-06 -3.56e-05 -2.97e-04 -1.33e-04 -1.34e-04 -1.38e-04 -2.68e-05 -1.32e-06 -1.62e-21 7.05e-06 -8.15e-05 -6.00e-05 -3.18e-05 -1.14e-04 3.45e-06 3.56e-05 -7.05e-06 -6.99e-22 1.28e-04 -3.53e-04 -8.36e-05 -8.14e-04 -6.29e-05 2.97e-04 8.15e-05 -1.28e-04 6.25e-22
    0.00e+00 2.78e-04 3.75e-04 3.30e-04 -5.13e-04 -2.27e-04 -2.70e-04 -1.47e-04 -2.78e-04 -1.04e-21 2.27e-04 2.22e-04 2.45e-04 -3.57e-04 2.00e-04 6.56e-05 -3.75e-04 -2.27e-04 3.28e-22 1.65e-04 -2.49e-04 -6.12e-04 1.29e-04 1.74e-04 -3.30e-04 -2.22e-04 -1.65e-04 -1.65e-22 -8.15e-05 -3.56e-04 6.74e-05 2.25e-04 5.13e-04 -2.45e-04 2.49e-04 8.15e-05 -4.24e-22 -1.31e-04 -1.32e-04 2.16e-05 2.27e-04 3.57e-04 6.12e-04 3.56e-04 1.31e-04 -1.43e-21 -2.74e-04 -2.30e-04 2.70e-04 -2.00e-04 -1.29e-04 -6.74e-05 1.32e-04 2.74e-04 -2.67e-21 7.19e-05 1.47e-04 -6.56e-05 -1.74e-04 -2.25e-04 -2.16e-05 2.30e-04 -7.19e-05 3.14e-22
    0.00e+00 1.75e-04 3.34e-04 1.12e-04 -3.01e-04 3.11e-04 -7.42e-05 -1.33e-04 -1.75e-04 2.12e-21 8.55e-06 1.79e-04 2.11e-04 5.23e-04 6.33e-04 2.17e-04 -3.34e-04 -8.55e-06 3.57e-22 -1.00e-04 2.47e-04 2.06e-04 1.58e-04 2.91e-04 -1.12e-04 -1.79e-04 1.00e-04 1.87e-21 -6.87e-05 4.03e-04 1.63e-04 1.94e-04 3.01e-04 -2.11e-04 -2.47e-04 6.87e-05 -1.11e-21 -3.21e-04 -7.64e-05 -1.06e-04 -3.11e-04 -5.23e-04 -2.06e-04 -4.03e-04 3.21e-04 -4.18e-22 -3.75e-04 3.62e-04 7.42e-05 -6.33e-04 -1.58e-04 -1.63e-04 7.64e-05 3.75e-04 4.69e-21 2.10e-04 1.33e-04 -2.17e-04 -2.91e-04 -1.94e-04 1.06e-04 -3.62e-04 -2.10e-04 -1.19e-21
    0.00e+00 -6.20e-05 1.69e-04 2.42e-05 -2.61e-04 2.78e-05 -9.21e-05 -1.17e-04 6.20e-05 -4.99e-22 1.76e-04 6.63e-04 -1.31e-04 -2.77e-04 -2.22e-05 -7.59e-05 -1.69e-04 -1.76e-04 -3.53e-22 -6.85e-05 -2.45e-05 -1.21e-04 -6.64e-05 1.28e-04 -2.42e-05 -6.63e-04 6.85e-05 2.32e-21 -1.98e-04 -5.11e-04 2.01e-05 3.53e-05 2.61e-04 1.31e-04 2.45e-05 1.98e-04 -8.77e-22 6.15e-05 1.88e-04 -6.71e-05 -2.78e-05 2.77e-04 1.21e-04 5.11e-04 -6.15e-05 -2.64e-21 -3.44e-04 -1.06e-04 9.21e-05 2.22e-05 6.64e-05 -2.01e-05 -1.88e-04 3.44e-04 -8.50e-22 -2.38e-04 1.17e-04 7.59e-05 -1.28e-04 -3.53e-05 6.71e-05 1.06e-04 2.38e-04 -2.75e-22
    0.00e+00 3.06e-04 9.58e-05 6.87e-05 -6.50e-05 -4.54e-04 5.05e-04 -1.07e-04 -3.06e-04 -1.18e-21 1.81e-04 6.56e-04 -2.81e-04 4.00e-04 8.75e-06 1.75e-04 -9.58e-05 -1.81e-04 -2.89e-22 -3.80e-05 -1.43e-04 1.52e-04 -3.25e-04 5.86e-05 -6.87e-05 -6.56e-04 3.80e-05 -8.60e-22 -1.91e-04 -1.65e-04 1.07e-04 -1.69e-04 6.50e-05 2.81e-04 1.43e-04 1.91e-04 -6.12e-22 -1.09e-04 4.16e-04 6.08e-05 4.54e-04 -4.00e-04 -1.52e-04 1.65e-04 1.09e-04 4.63e-22 2.85e-05 2.06e-04 -5.05e-04 -8.75e-06 3.25e-04 -1.07e-04 -4.16e-04 -2.85e-05 -5.04e-22 -3.60e-04 1.07e-04 -1.75e-04 -5.86e-05 1.69e-04 -6.08e-05 -2.06e-04 3.60e-04 -1.24e-23
    0.00e+00 1.59e-04 3.21e-04 1.90e-04 -4.80e-04 -2.06e-04 -2.11e-04 -1.10e-04 -1.59e-04 -1.06e-21 1.43e-04 1.89e-04 1.70e-04 -4.49e-04 3.26e-04 5.69e-05 -3.21e-04 -1.43e-04 2.63e-22 8.01e-05 -1.70e-04 -3.60e-04 1.46e-04 1.78e-04 -1.90e-04 -1.89e-04 -8.01e-05 3.55e-22 -1.08e-04 -2.81e-04 1.71e-04 2.52e-04 4.80e-04 -1.70e-04 1.70e-04 1.08e-04 -7.96e-22 -1.49e-04 -9.53e-05 -6.57e-05 2.06e-04 4.49e-04 3.60e-04 2.81e-04 1.49e-04 -2.07e-21 -3.45e-04 -2.93e-04 2.11e-04 -3.26e-04 -1.46e-04 -1.71e-04 9.53e-05 3.45e-04 8.14e-22 4.43e-05 1.10e-04 -5.69e-05 -1.78e-04 -2.52e-04 6.57e-05 2.93e-04 -4.43e-05 -1.80e-21
    0.00e+00 2.57e-04 -2.74e-04 -7.28e-05 -2.93e-05 7.46e-05 7.72e-05 1.77e-04 -2.57e-04 2.04e-21 2.06e-05 1.83e-04 2.03e-04 6.23e-04 7.62e-04 -1.87e-04 2.74e-04 -2.06e-05 -4.30e-23 -2.68e-04 1.78e-05 -2.32e-04 3.03e-04 9.24e-04 7.28e-05 -1.83e-04 2.68e-04 -2.40e-21 -4.60e-05 3.91e-05 1.38e-04 4.89e-05 2.93e-05 -2.03e-04 -1.78e-05 4.60e-05 8.93e-23 -3.18e-06 6.24e-06 -9.43e-05 -7.46e-05 -6.23e-04 2.32e-04 -3.91e-05 3.18e-06 1.01e-21 7.08e-05 -9.25e-05 -7.72e-05 -7.62e-04 -3.03e-04 -1.38e-04 -6.24e-06 -7.08e-05 6.07e-22 1.93e-05 -1.77e-04 1.87e-04 -9.24e-04 -4.89e-05 9.43e-05 9.25e-05 -1.93e-05 -1.57e-21
    0.00e+00 3.27e-04 -1.42e-04 1.52e-04 9.42e-05 -3.76e-04 3.64e-05 1.36e-04 -3.27e-04 3.87e-21 1.47e-04 1.89e-04 -1.66e-04 -1.35e-04 -5.53e-04 2.40e-04 1.42e-04 -1.47e-04 -6.62e-22 -1.06e-04 -1.29e-04 3.42e-04 -2.52e-04 -3.38e-05 -1.52e-04 -1.89e-04 1.06e-04 -8.60e-22 -6.76e-05 -1.35e-04 3.07e-04 -1.79e-04 -9.42e-05 1.66e-04 1.29e-04 6.76e-05 7.28e-22 -1.76e-04 -2.26e-04 1.21e-04 3.76e-04 1.35e-04 -3.42e-04 1.35e-04 1.76e-04 -3.61e-21 2.11e-04 1.91e-04 -3.64e-05 5.53e-04 2.52e-04 -3.07e-04 2.26e-04 -2.11e-04 -2.21e-21 -5.84e-04 -1.36e-04 -2.40e-04 3.38e-05 1.79e-04 -1.21e-04 -1.91e-04 5.84e-04 4.75e-22
    0.00e+00 3.53e-04 2.85e-05 -4.95e-05 4.45e-06 -3.84e-04 3.67e-04 -1.63e-05 -3.53e-04 -9.83e-21 2.36e-04 3.05e-04 -2.38e-04 3.42e-04 5.21e-05 5.47e-05 -2.85e-05 -2.36e-04 -5.58e-22 3.35e-05 7.12e-06 4.25e-05 -7.77e-05 -2.95e-05 4.95e-05 -3.05e-04 -3.35e-05 -5.53e-22 -1.64e-04 -1.25e-04 -3.66e-05 5.12e-05 -4.45e-06 2.38e-04 -7.12e-06 1.64e-04 1.88e-22 -1.21e-04 2.45e-04 7.23e-05 3.84e-04 -3.42e-04 -4.25e-05 1.25e-04 1.21e-04 -1.89e-21 -6.72e-06 1.52e-04 -3.67e-04 -5.21e-05 7.77e-05 3.66e-05 -2.45e-04 6.72e-06 3.29e-21 1.33e-05 1.63e-05 -5.47e-05 2.95e-05 -5.12e-05 -7.23e-05 -1.52e-04 -1.33e-05 6.86e-21
    0.00e+00 2.88e-04 -1.69e-04 -2.44e-04 -6.55e-05 1.31e-04 9.90e-05 4.71e-04 -2.88e-04 -1.53e-21 -1.38e-04 2.00e-05 2.06e-04 3.03e-04 2.08e-04 -1.13e-04 1.69e-04 1.38e-04 -1.75e-21 -1.59e-04 -4.64e-05 8.04e-05 4.01e-04 9.19e-04 2.44e-04 -2.00e-05 1.59e-04 -1.25e-21 -2.19e-04 -1.80e-04 -6.09e-05 5.16e-05 6.55e-05 -2.06e-04 4.64e-05 2.19e-04 5.76e-22 -1.05e-05 9.60e-06 -3.15e-04 -1.31e-04 -3.03e-04 -8.04e-05 1.80e-04 1.05e-05 -4.67e-22 2.53e-05 -2.56e-04 -9.90e-05 -2.08e-04 -4.01e-04 6.09e-05 -9.60e-06 -2.53e-05 -2.33e-22 8.75e-05 -4.71e-04 1.13e-04 -9.19e-04 -5.16e-05 3.15e-04 2.56e-04 -8.75e-05 -4.76e-23
    0.00e+00 -3.90e-04 1.17e-04 1.52e-04 -3.98e-04 6.72e-04 -9.48e-05 -3.27e-04 3.90e-04 1.45e-21 8.09e-05 -2.74e-04 -7.53e-05 -5.06e-06 1.16e-04 4.56e-05 -1.17e-04 -8.09e-05 -2.02e-21 4.33e-04 -3.83e-05 1.28e-04 -1.04e-05 -4.06e-05 -1.52e-04 2.74e-04 -4.33e-04 2.68e-21 -1.52e-04 -2.28e-04 7.82e-05 6.66e-04 3.98e-04 7.53e-05 3.83e-05 1.52e-04 -1.91e-21 -7.53e-05 3.66e-05 1.07e-04 -6.72e-04 5.06e-06 -1.28e-04 2.28e-04 7.53e-05 -8.57e-22 -1.28e-04 2.42e-05 9.48e-05 -1.16e-04 1.04e-05 -7.82e-05 -3.66e-05 1.28e-04 5.70e-23 -7.45e-05 3.27e-04 -4.56e-05 4.06e-05 -6.66e-04 -1.07e-04 -2.42e-05 7.45e-05 -6.10e-22
    0.00e+00 -3.49e-04 1.21e-04 -3.44e-04 2.98e-05 4.26e-04 -9.85e-05 -1.38e-05 3.49e-04 -9.00e-22 -3.30e-04 -5.47e-05 -1.43e-04 5.55e-05 1.85e-04 3.79e-04 -1.21e-04 3.30e-04 -9.90e-22 3.21e-04 -8.88e-05 -2.30e-04 1.23e-04 -7.89e-05 3.44e-04 5.47e-05 -3.21e-04 -2.80e-22 -4.89e-04 1.77e-04 2.06e-04 5.66e-04 -2.98e-05 1.43e-04 8.88e-05 4.89e-04 -2.62e-21 5.59e-05 1.35e-05 -2.88e-04 -4.26e-04 -5.55e-05 2.30e-04 -1.77e-04 -5.59e-05 9.15e-22 -5.19e-05 -1.87e-04 9.85e-05 -1.85e-04 -1.23e-04 -2.06e-04 -1.35e-05 5.19e-05 -4.62e-22 1.48e-04 1.38e-05 -3.79e-04 7.89e-05 -5.66e-04 2.88e-04 1.87e-04 -1.48e-04 -7.04e-22
    0.00e+00 1.20e-04 1.53e-04 4.51e-05 -2.68e-04 7.58e-06 -1.25e-05 -6.60e-05 -1.20e-04 -1.17e-21 3.11e-04 6.03e-04 -1.11e-04 2.01e-04 1.43e-05 2.06e-05 -1.53e-04 -3.11e-04 3.58e-22 1.72e-05 -1.10e-04 -6.36e-05 -1.71e-04 3.35e-05 -4.51e-05 -6.03e-04 -1.72e-05 -5.35e-22 -1.47e-04 -3.54e-04 -5.98e-05 1.34e-05 2.68e-04 1.11e-04 1.10e-04 1.47e-04 -1.37e-21 4.75e-05 2.78e-04 3.27e-05 -7.58e-06 -2.01e-04 6.36e-05 3.54e-04 -4.75e-05 3.63e-21 -1.23e-04 4.50e-05 1.25e-05 -1.43e-05 1.71e-04 5.98e-05 -2.78e-04 1.23e-04 -2.41e-22 -1.94e-04 6.60e-05 -2.06e-05 -3.35e-05 -1.34e-05 -3.27e-05 -4.50e-05 1.94e-04 5.77e-22
    0.00e+00 -5.83e-05 -3.06e-04 -1.17e-04 -6.66e-05 1.10e-04 1.65e-04 -3.08e-04 5.83e-05 -1.98e-22 -3.82e-04 -2.64e-04 -7.12e-05 2.13e-04 -1.05e-06 -3.71e-04 3.06e-04 3.82e-04 1.27e-21 2.93e-04 1.15e-04 6.91e-05 -8.66e-05 -1.11e-04 1.17e-04 2.64e-04 -2.93e-04 -3.76e-21 -1.61e-04 -1.11e-05 -1.34e-04 -1.32e-04 6.66e-05 7.12e-05 -1.15e-04 1.61e-04 -3.06e-21 2.59e-05 1.86e-05 -2.36e-04 -1.10e-04 -2.13e-04 -6.91e-05 1.11e-05 -2.59e-05 4.67e-23 1.59e-04 -3.79e-05 -1.65e-04 1.05e-06 8.66e-05 1.34e-04 -1.86e-05 -1.59e-04 -7.27e-22 3.63e-04 3.08e-04 3.71e-04 1.11e-04 1.32e-04 2.36e-04 3.79e-05 -3.63e-04 -9.79e-22
    0.00e+00 -3.88e-04 2.01e-04 -7.77e-05 -3.34e-04 6.90e-04 -7.51e-05 -3.40e-04 3.88e-04 2.53e-22 -4.56e-05 -2.41e-04 -5.95e-05 1.04e-04 1.36e-04 1.12e-04 -2.01e-04 4.56e-05 -1.12e-21 4.60e-04 -7.62e-05 2.47e-05 -8.97e-06 4.09e-05 7.77e-05 2.41e-04 -4.60e-04 -1.59e-22 -2.22e-04 -1.80e-04 1.50e-04 7.75e-04 3.34e-04 5.95e-05 7.62e-05 2.22e-04 1.80e-22 1.14e-05 1.69e-05 -2.64e-06 -6.90e-04 -1.04e-04 -2.47e-05 1.80e-04 -1.14e-05 1.02e-21 -1.11e-04 -6.13e-05 7.51e-05 -1.36e-04 8.97e-06 -1.50e-04 -1.69e-05 1.11e-04 7.01e-23 -2.49e-05 3.40e-04 -1.12e-04 -4.09e-05 -7.75e-04 2.64e-06 6.13e-05 2.49e-05 2.18e-21
    0.00e+00 1.01e-04 -9.99e-05 2.28e-05 -8.99e-06 -1.09e-04 3.78e-04 -5.68e-05 -1.01e-04 2.74e-21 -7.57e-05 6.18e-04 -3.31e-04 -1.93e-04 -4.48e-04 -1.54e-04 9.99e-05 7.57e-05 1.83e-21 -3.50e-05 6.77e-05 4.87e-05 -1.52e-04 1.68e-04 -2.28e-05 -6.18e-04 3.50e-05 -1.86e-21 -1.65e-04 -2.53e-04 -1.69e-04 -1.21e-04 8.99e-06 3.31e-04 -6.77e-05 1.65e-04 2.74e-22 3.59e-05 1.07e-04 1.07e-06 1.09e-04 1.93e-04 -4.87e-05 2.53e-04 -3.59e-05 4.93e-21 -3.43e-04 3.58e-05 -3.78e-04 4.48e-04 1.52e-04 1.69e-04 -1.07e-04 3.43e-04 -6.79e-23 -2.28e-04 5.68e-05 1.54e-04 -1.68e-04 1.21e-04 -1.07e-06 -3.58e-05 2.28e-04 -2.59e-22
    0.00e+00 1.69e-04 3.02e-04 2.86e-04 -4.69e-04 -1.19e-04 -3.50e-04 -1.85e-04 -1.69e-04 -1.53e-21 2.17e-04 1.93e-04 2.00e-04 -4.33e-04 2.67e-04 -4.97e-05 -3.02e-04 -2.17e-04 -3.57e-22 7.84e-05 -1.77e-04 -6.02e-04 6.70e-05 2.04e-04 -2.86e-04 -1.93e-04 -7.84e-05 -3.48e-21 -8.71e-05 -4.92e-04 1.28e-04 2.36e-04 4.69e-04 -2.00e-04 1.77e-04 8.71e-05 2.36e-21 -8.60e-05 -1.20e-04 5.69e-06 1.19e-04 4.33e-04 6.02e-04 4.92e-04 8.60e-05 2.12e-22 -3.09e-04 -2.15e-04 3.50e-04 -2.67e-04 -6.70e-05 -1.28e-04 1.20e-04 3.09e-04 -1.75e-21 6.36e-05 1.85e-04 4.97e-05 -2.04e-04 -2.36e-04 -5.69e-06 2.15e-04 -6.36e-05 2.11e-22
    0.00e+00 -8.63e-05 1.67e-04 -3.30e-04 -1.12e-04 5.58e-04 2.42e-05 -2.04e-04 8.63e-05 4.46e-21 -1.20e-04 -4.96e-04 8.87e-05 1.49e-04 1.83e-04 1.10e-04 -1.67e-04 1.20e-04 -3.21e-22 5.68e-05 2.81e-05 -3.45e-04 5.37e-05 1.88e-04 3.30e-04 4.96e-04 -5.68e-05 4.01e-21 -4.34e-04 -3.89e-05 1.40e-04 7.01e-04 1.12e-04 -8.87e-05 -2.81e-05 4.34e-04 5.01e-22 -9.10e-06 -1.25e-04 -2.64e-04 -5.58e-04 -1.49e-04 3.45e-04 3.89e-05 9.10e-06 -2.43e-21 -4.75e-05 -3.38e-04 -2.42e-05 -1.83e-04 -5.37e-05 -1.40e-04 1.25e-04 4.75e-05 -3.19e-22 -5.08e-05 2.04e-04 -1.10e-04 -1.88e-04 -7.01e-04 2.64e-04 3.38e-04 5.08e-05 6.58e-22
    0.00e+00 1.35e-04 1.23e-04 1.32e-04 -4.18e-04 1.30e-04 -6.28e-05 -1.92e-05 -1.35e-04 -1.36e-22 -1.27e-05 1.85e-05 1.92e-04 -1.13e-05 4.36e-04 1.54e-04 -1.23e-04 1.27e-05 -1.16e-21 5.07e-05 -8.75e-05 1.79e-04 1.83e-04 1.84e-04 -1.32e-04 -1.85e-05 -5.07e-05 -1.21e-21 -7.65e-05 1.83e-04 2.22e-04 2.68e-04 4.18e-04 -1.92e-04 8.75e-05 7.65e-05 2.16e-22 -3.00e-04 -1.33e-04 -1.64e-04 -1.30e-04 1.13e-05 -1.79e-04 -1.83e-04 3.00e-04 6.43e-22 -3.73e-04 -1.44e-04 6.28e-05 -4.36e-04 -1.83e-04 -2.22e-04 1.33e-04 3.73e-04 3.60e-22 7.47e-05 1.92e-05 -1.54e-04 -1.84e-04 -2.68e-04 1.64e-04 1.44e-04 -7.47e-05 -2.49e-21
    0.00e+00 2.73e-04 -2.62e-04 -1.91e-05 -8.07e-05 -3.71e-04 3.61e-05 5.48e-05 -2.73e-04 6.19e-22 1.00e-04 2.75e-04 -1.61e-04 1.21e-04 -5.56e-04 3.33e-04 2.62e-04 -1.00e-04 1.18e-21 -1.30e-04 -3.45e-04 2.44e-04 -3.57e-04 1.18e-04 1.91e-05 -2.75e-04 1.30e-04 -1.26e-21 1.28e-04 1.04e-04 3.30e-04 -2.68e-04 8.07e-05 1.61e-04 3.45e-04 -1.28e-04 -2.44e-21 3.47e-05 -3.57e-04 1.33e-05 3.71e-04 -1.21e-04 -2.44e-04 -1.04e-04 -3.47e-05 -1.48e-21 3.77e-04 -1.15e-04 -3.61e-05 5.56e-04 3.57e-04 -3.30e-04 3.57e-04 -3.77e-04 7.16e-22 -4.30e-04 -5.48e-05 -3.33e-04 -1.18e-04 2.68e-04 -1.33e-05 1.15e-04 4.30e-04 5.05e-22
    0.00e+00 1.45e-04 1.45e-04 1.11e-04 -2.31e-04 2.35e-04 -8.39e-05 -3.07e-04 -1.45e-04 -2.76e-21 1.02e-04 3.13e-04 1.58e-04 6.23e-04 7.33e-04 5.32e-05 -1.45e-04 -1.02e-04 -2.14e-21 -1.49e-04 2.28e-04 -1.67e-04 1.42e-04 4.47e-04 -1.11e-04 -3.13e-04 1.49e-04 1.90e-22 -7.27e-05 2.98e-04 7.25e-05 6.83e-05 2.31e-04 -1.58e-04 -2.28e-04 7.27e-05 1.19e-21 -1.71e-04 5.58e-06 1.54e-05 -2.35e-04 -6.23e-04 1.67e-04 -2.98e-04 1.71e-04 1.56e-21 -2.55e-04 4.22e-04 8.39e-05 -7.33e-04 -1.42e-04 -7.25e-05 -5.58e-06 2.55e-04 4.13e-22 2.01e-04 3.07e-04 -5.32e-05 -4.47e-04 -6.83e-05 -1.54e-05 -4.22e-04 -2.01e-04 1.58e-21
    0.00e+00 -1.05e-04 -1.42e-04 -6.48e-06 1.47e-04 1.04e-06 1.72e-04 1.27e-04 1.05e-04 -1.04e-21 -1.36e-04 -1.59e-04 -1.50e-04 -8.70e-06 4.60e-05 1.17e-05 1.42e-04 1.36e-04 -8.35e-22 1.03e-05 1.69e-04 1.96e-05 1.23e-05 -4.93e-04 6.48e-06 1.59e-04 -1.03e-05 -4.45e-21 1.86e-04 2.28e-04 -1.24e-04 -5.26e-04 -1.47e-04 1.50e-04 -1.69e-04 -1.86e-04 -6.10e-22 5.23e-04 3.99e-05 -3.92e-04 -1.04e-06 8.70e-06 -1.96e-05 -2.28e-04 -5.23e-04 -4.63e-22 9.93e-05 4.83e-04 -1.72e-04 -4.60e-05 -1.23e-05 1.24e-04 -3.99e-05 -9.93e-05 -7.79e-21 5.61e-04 -1.27e-04 -1.17e-05 4.93e-04 5.26e-04 3.92e-04 -4.83e-04 -5.61e-04 1.42e-22
    0.00e+00 1.83e-04 -2.22e-04 5.27e-05 -1.15e-04 1.12e-04 -5.76e-05 -1.97e-04 -1.83e-04 4.22e-23 1.10e-04 2.64e-04 1.54e-04 7.00e-04 7.81e-04 -4.59e-05 2.22e-04 -1.10e-04 -2.06e-21 -1.67e-04 8.80e-05 -3.01e-04 1.63e-04 6.17e-04 -5.27e-05 -2.64e-04 1.67e-04 -4.23e-22 -1.62e-05 2.08e-04 1.13e-04 3.10e-06 1.15e-04 -1.54e-04 -8.80e-05 1.62e-05 -8.09e-23 -5.83e-05 3.62e-05 3.12e-05 -1.12e-04 -7.00e-04 3.01e-04 -2.08e-04 5.83e-05 1.73e-21 -1.22e-04 1.97e-04 5.76e-05 -7.81e-04 -1.63e-04 -1.13e-04 -3.62e-05 1.22e-04 1.77e-21 1.33e-04 1.97e-04 4.59e-05 -6.17e-04 -3.10e-06 -3.12e-05 -1.97e-04 -1.33e-04 1.31e-21
    0.00e+00 2.10e-04 -2.43e-04 9.81e-05 -9.11e-05 -8.96e-06 2.57e-04 6.54e-05 -2.10e-04 -2.12e-22 1.76e-04 -2.27e-04 1.39e-04 5.61e-04 6.51e-04 -6.36e-05 2.43e-04 -1.76e-04 3.08e-22 -2.55e-04 1.15e-04 -2.63e-04 5.65e-05 5.08e-04 -9.81e-05 2.27e-04 2.55e-04 -5.76e-22 -1.66e-05 3.05e-04 3.22e-04 -2.34e-05 9.11e-05 -1.39e-04 -1.15e-04 1.66e-05 7.73e-22 -1.49e-04 -1.12e-04 -2.34e-05 8.96e-06 -5.61e-04 2.63e-04 -3.05e-04 1.49e-04 -1.45e-21 -7.94e-05 1.59e-04 -2.57e-04 -6.51e-04 -5.65e-05 -3.22e-04 1.12e-04 7.94e-05 3.35e-21 -1.43e-04 -6.54e-05 6.36e-05 -5.08e-04 2.34e-05 2.34e-05 -1.59e-04 1.43e-04 -6.61e-22
    0.00e+00 -3.22e-04 4.24e-04 -2.50e-04 -1.98e-04 6.34e-04 -1.33e-05 -3.14e-04 3.22e-04 -1.55e-21 -1.28e-04 -3.37e-04 -7.93e-05 3.48e-05 1.86e-04 2.10e-04 -4.24e-04 1.28e-04 2.15e-21 3.75e-04 -8.08e-06 -7.61e-05 1.54e-05 -2.83e-05 2.50e-04 3.37e-04 -3.75e-04 2.08e-21 -4.02e-04 -5.58e-05 2.33e-04 8.48e-04 1.98e-04 7.93e-05 8.08e-06 4.02e-04 -6.43e-22 6.09e-05 -2.56e-05 -7.99e-05 -6.34e-04 -3.48e-05 7.61e-05 5.58e-05 -6.09e-05 3.02e-21 -6.93e-05 -1.30e-04 1.33e-05 -1.86e-04 -1.54e-05 -2.33e-04 2.56e-05 6.93e-05 -7.14e-22 -3.60e-05 3.14e-04 -2.10e-04 2.83e-05 -8.48e-04 7.99e-05 1.30e-04 3.60e-05 -9.53e-22
    0.00e+00 -3.51e-05 4.29e-04 2.62e-04 -2.24e-04 3.87e-04 7.07e-05 2.52e-04 3.51e-05 1.12e-22 1.72e-04 8.78e-05 -2.11e-04 -2.35e-04 1.50e-04 -9.68e-05 -4.29e-04 -1.72e-04 -6.29e-21 1.46e-04 -1.47e-04 -1.44e-04 1.04e-04 -7.18e-05 -2.62e-04 -8.78e-05 -1.46e-04 2.35e-21 8.60e-05 -6.17e-05 3.46e-05 7.94e-05 2.24e-04 2.11e-04 1.47e-04 -8.60e-05 1.53e-22 4.17e-04 -2.42e-04 3.78e-04 -3.87e-04 2.35e-04 1.44e-04 6.17e-05 -4.17e-04 -1.83e-21 -7.62e-05 1.94e-04 -7.07e-05 -1.50e-04 -1.04e-04 -3.46e-05 2.42e-04 7.62e-05 -4.09e-21 9.93e-05 -2.52e-04 9.68e-05 7.18e-05 -7.94e-05 -3.78e-04 -1.94e-04 -9.93e-05 1.71e-21
    0.00e+00 -8.12e-05 -2.75e-04 -3.62e-04 5.09e-05 3.74e-04 -9.99e-07 3.88e-04 8.12e-05 -3.30e-22 -9.74e-05 -4.27e-04 6.61e-06 2.52e-04 1.62e-04 3.12e-04 2.75e-04 9.74e-05 2.47e-21 1.15e-04 -8.99e-05 -3.67e-06 1.30e-04 3.11e-04 3.62e-04 4.27e-04 -1.15e-04 1.22e-21 -3.75e-04 7.84e-05 2.04e-04 4.58e-04 -5.09e-05 -6.61e-06 8.99e-05 3.75e-04 -1.72e-21 -1.33e-04 -4.93e-05 -4.72e-04 -3.74e-04 -2.52e-04 3.67e-06 -7.84e-05 1.33e-04 -7.87e-22 -1.38e-04 -2.57e-04 9.99e-07 -1.62e-04 -1.30e-04 -2.04e-04 4.93e-05 1.38e-04 -7.60e-21 2.65e-04 -3.88e-04 -3.12e-04 -3.11e-04 -4.58e-04 4.72e-04 2.57e-04 -2.65e-04 -6.84e-22
    0.00e+00 2.47e-04 -2.01e-04 -3.08e-04 -5.60e-05 2.26e-04 9.07e-05 4.71e-04 -2.47e-04 3.92e-22 -9.46e-05 -1.96e-04 2.33e-04 2.74e-04 8.64e-05 2.22e-05 2.01e-04 9.46e-05 1.79e-21 -9.79e-05 -8.21e-05 1.19e-04 3.21e-04 8.29e-04 3.08e-04 1.96e-04 9.79e-05 -1.92e-21 -2.73e-04 -1.44e-04 -9.19e-05 1.17e-04 5.60e-05 -2.33e-04 8.21e-05 2.73e-04 -4.68e-21 -5.19e-05 2.93e-06 -3.93e-04 -2.26e-04 -2.74e-04 -1.19e-04 1.44e-04 5.19e-05 2.27e-21 -4.07e-05 -3.07e-04 -9.07e-05 -8.64e-05 -3.21e-04 9.19e-05 -2.93e-06 4.07e-05 3.91e-21 1.68e-04 -4.71e-04 -2.22e-05 -8.29e-04 -1.17e-04 3.93e-04 3.07e-04 -1.68e-04 -1.11e-21
    0.00e+00 2.74e-04 -3.35e-04 -2.01e-04 -3.68e-05 1.35e-04 6.48e-05 3.26e-04 -2.74e-04 -3.18e-21 -1.13e-04 1.63e-04 2.05e-04 4.66e-04 3.48e-04 -2.49e-04 3.35e-04 1.13e-04 7.66e-22 -1.26e-04 -5.55e-05 -7.84e-05 3.58e-04 9.23e-04 2.01e-04 -1.63e-04 1.26e-04 -6.72e-22 -1.55e-04 -2.18e-04 -1.14e-04 1.02e-04 3.68e-05 -2.05e-04 5.55e-05 1.55e-04 -2.40e-22 3.31e-05 5.64e-05 -2.24e-04 -1.35e-04 -4.66e-04 7.84e-05 2.18e-04 -3.31e-05 2.95e-23 -2.74e-06 -3.52e-04 -6.48e-05 -3.48e-04 -3.58e-04 1.14e-04 -5.64e-05 2.74e-06 -8.62e-22 7.11e-05 -3.26e-04 2.49e-04 -9.23e-04 -1.02e-04 2.24e-04 3.52e-04 -7.11e-05 2.75e-22
    0.00e+00 -3.63e-04 1.30e-04 1.65e-04 -4.08e-04 7.24e-04 -4.03e-05 -2.41e-04 3.63e-04 -2.45e-21 9.02e-05 -2.53e-04 -8.25e-05 7.93e-05 1.18e-04 1.54e-05 -1.30e-04 -9.02e-05 7.83e-22 4.09e-04 -8.10e-05 1.51e-04 -1.28e-05 -3.06e-06 -1.65e-04 2.53e-04 -4.09e-04 8.13e-22 -1.35e-04 -2.05e-04 1.01e-04 6.76e-04 4.08e-04 8.25e-05 8.10e-05 1.35e-04 -7.39e-22 -8.68e-05 2.64e-05 6.52e-05 -7.24e-04 -7.93e-05 -1.51e-04 2.05e-04 8.68e-05 -6.62e-21 -1.36e-04 4.15e-05 4.03e-05 -1.18e-04 1.28e-05 -1.01e-04 -2.64e-05 1.36e-04 -4.51e-22 -9.68e-05 2.41e-04 -1.54e-05 3.06e-06 -6.76e-04 -6.52e-05 -4.15e-05 9.68e-05 1.05e-21
    0.00e+00 1.75e-04 1.95e-04 1.84e-04 -3.12e-04 2.83e-04 -4.30e-05 -1.91e-05 -1.75e-04 -1.06e-21 -4.78e-05 -6.59e-06 2.01e-04 3.78e-04 5.29e-04 2.30e-04 -1.95e-04 4.78e-05 2.99e-21 -5.86e-05 9.12e-05 2.82e-04 1.78e-04 2.64e-04 -1.84e-04 6.59e-06 5.86e-05 -2.52e-22 -7.03e-05 2.76e-04 2.38e-04 2.66e-04 3.12e-04 -2.01e-04 -9.12e-05 7.03e-05 -1.09e-21 -3.44e-04 -1.67e-04 -1.56e-04 -2.83e-04 -3.78e-04 -2.82e-04 -2.76e-04 3.44e-04 -2.67e-21 -3.98e-04 3.46e-05 4.30e-05 -5.29e-04 -1.78e-04 -2.38e-04 1.67e-04 3.98e-04 -8.22e-22 1.33e-04 1.91e-05 -2.30e-04 -2.64e-04 -2.66e-04 1.56e-04 -3.46e-05 -1.33e-04 -1.64e-21
    0.00e+00 -1.55e-04 -1.92e-04 -2.94e-04 -7.75e-05 1.42e-04 1.24e-04 -2.43e-04 1.55e-04 4.76e-22 -2.23e-04 -4.25e-04 -1.09e-04 1.60e-04 1.38e-04 -5.43e-04 1.92e-04 2.23e-04 1.74e-21 -2.52e-05 1.88e-04 9.19e-05 9.54e-05 -1.45e-04 2.94e-04 4.25e-04 2.52e-05 1.05e-21 -4.72e-05 -9.02e-05 -1.34e-04 -1.72e-04 7.75e-05 1.09e-04 -1.88e-04 4.72e-05 1.48e-21 1.16e-04 2.27e-05 -4.74e-04 -1.42e-04 -1.60e-04 -9.19e-05 9.02e-05 -1.16e-04 6.64e-22 1.55e-05 7.02e-05 -1.24e-04 -1.38e-04 -9.54e-05 1.34e-04 -2.27e-05 -1.55e-05 -1.47e-21 2.05e-04 2.43e-04 5.43e-04 1.45e-04 1.72e-04 4.74e-04 -7.02e-05 -2.05e-04 2.99e-21
    0.00e+00 -2.32e-06 2.51e-04 -9.88e-05 3.33e-04 3.29e-04 3.72e-04 3.47e-04 2.32e-06 -7.15e-22 -3.00e-05 -2.96e-04 -1.67e-04 -1.54e-04 8.60e-05 9.02e-05 -2.51e-04 3.00e-05 6.17e-21 2.65e-05 3.98e-05 2.94e-05 -2.89e-05 -1.83e-04 9.88e-05 2.96e-04 -2.65e-05 3.03e-21 3.25e-04 3.56e-04 5.19e-04 4.13e-04 -3.33e-04 1.67e-04 -3.98e-05 -3.25e-04 6.32e-22 3.53e-04 3.47e-04 2.94e-04 -3.29e-04 1.54e-04 -2.94e-05 -3.56e-04 -3.53e-04 -8.85e-22 8.91e-05 5.54e-05 -3.72e-04 -8.60e-05 2.89e-05 -5.19e-04 -3.47e-04 -8.91e-05 5.12e-22 1.04e-04 -3.47e-04 -9.02e-05 1.83e-04 -4.13e-04 -2.94e-04 -5.54e-05 -1.04e-04 3.91e-22
    0.00e+00 1.38e-04 9.04e-05 1.53e-04 -3.81e-04 3.16e-04 -1.07e-05 -1.13e-05 -1.38e-04 -3.41e-22 -5.04e-05 7.58e-05 1.84e-04 2.50e-04 3.91e-04 1.78e-04 -9.04e-05 5.04e-05 -2.57e-21 4.37e-05 -3.52e-05 3.00e-04 1.86e-04 2.19e-04 -1.53e-04 -7.58e-05 -4.37e-05 -7.25e-22 -3.63e-05 3.30e-04 1.65e-04 2.36e-04 3.81e-04 -1.84e-04 3.52e-05 3.63e-05 1.29e-21 -3.92e-04 -1.38e-04 -1.18e-04 -3.16e-04 -2.50e-04 -3.00e-04 -3.30e-04 3.92e-04 -2.55e-22 -3.09e-04 -2.94e-05 1.07e-05 -3.91e-04 -1.86e-04 -1.65e-04 1.38e-04 3.09e-04 -1.73e-21 1.47e-04 1.13e-05 -1.78e-04 -2.19e-04 -2.36e-04 1.18e-04 2.94e-05 -1.47e-04 -1.30e-21
    0.00e+00 6.99e-06 2.50e-04 1.80e-04 -4.12e-04 -1.37e-05 -3.60e-04 -8.72e-05 -6.99e-06 -1.99e-21 2.22e-04 2.24e-04 1.24e-04 -5.13e-04 3.84e-04 -1.49e-05 -2.50e-04 -2.22e-04 -8.74e-22 1.13e-06 -8.97e-05 -3.04e-04 1.02e-04 1.65e-04 -1.80e-04 -2.24e-04 -1.13e-06 -7.90e-22 -1.57e-04 -4.23e-04 2.19e-04 2.21e-04 4.12e-04 -1.24e-04 8.97e-05 1.57e-04 -6.10e-22 -5.74e-05 -5.03e-05 -8.59e-05 1.37e-05 5.13e-04 3.04e-04 4.23e-04 5.74e-05 3.22e-21 -3.86e-04 -2.78e-04 3.60e-04 -3.84e-04 -1.02e-04 -2.19e-04 5.03e-05 3.86e-04 -1.71e-22 -8.04e-06 8.72e-05 1.49e-05 -1.65e-04 -2.21e-04 8.59e-05 2.78e-04 8.04e-06 3.37e-21
    0.00e+00 -1.04e-04 -2.54e-04 -2.20e-04 -8.73e-05 1.45e-04 1.88e-04 -3.17e-04 1.04e-04 1.59e-22 -3.31e-04 -3.25e-04 -9.21e-05 2.18e-04 7.11e-05 -4.46e-04 2.54e-04 3.31e-04 -1.57e-21 2.68e-04 1.12e-04 6.86e-05 1.02e-06 -2.30e-04 2.20e-04 3.25e-04 -2.68e-04 -2.85e-21 -1.53e-04 -2.11e-05 -1.42e-04 -1.75e-04 8.73e-05 9.21e-05 -1.12e-04 1.53e-04 -4.41e-22 4.98e-05 -2.39e-05 -3.56e-04 -1.45e-04 -2.18e-04 -6.86e-05 2.11e-05 -4.98e-05 2.67e-22 1.05e-04 4.74e-05 -1.88e-04 -7.11e-05 -1.02e-06 1.42e-04 2.39e-05 -1.05e-04 -1.61e-22 3.88e-04 3.17e-04 4.46e-04 2.30e-04 1.75e-04 3.56e-04 -4.74e-05 -3.88e-04 -2.24e-23
    0.00e+00 -1.59e-04 -1.63e-04 -1.56e-04 -8.02e-05 1.72e-04 9.18e-05 -6.76e-05 1.59e-04 9.07e-22 -1.96e-04 -3.69e-04 -2.15e-04 1.17e-04 2.07e-04 -5.19e-04 1.63e-04 1.96e-04 2.22e-21 -4.59e-05 7.75e-05 5.91e-05 9.19e-05 -3.27e-04 1.56e-04 3.69e-04 4.59e-05 -4.95e-22 -6.44e-05 -1.09e-04 7.47e-05 -3.68e-04 8.02e-05 2.15e-04 -7.75e-05 6.44e-05 -8.97e-22 7.88e-05 4.53e-05 -5.36e-04 -1.72e-04 -1.17e-04 -5.91e-05 1.09e-04 -7.88e-05 9.93e-24 -2.89e-05 2.51e-04 -9.18e-05 -2.07e-04 -9.19e-05 -7.47e-05 -4.53e-05 2.89e-05 1.52e-21 1.73e-04 6.76e-05 5.19e-04 3.27e-04 3.68e-04 5.36e-04 -2.51e-04 -1.73e-04 -7.99e-22
    0.00e+00 -1.28e-04 -1.73e-04 4.18e-05 2.04e-05 1.43e-04 9.42e-06 9.90e-05 1.28e-04 9.16e-22 -1.37e-04 -4.27e-04 -1.43e-04 6.52e-05 1.52e-04 -4.69e-04 1.73e-04 1.37e-04 -3.83e-22 -3.09e-04 3.49e-05 -8.20e-05 2.02e-04 -3.70e-04 -4.18e-05 4.27e-04 3.09e-04 -1.67e-22 2.49e-05 -9.64e-05 6.71e-05 -3.48e-04 -2.04e-05 1.43e-04 -3.49e-05 -2.49e-05 1.23e-21 1.17e-04 8.14e-05 -6.35e-04 -1.43e-04 -6.52e-05 8.20e-05 9.64e-05 -1.17e-04 2.89e-22 -6.22e-05 2.38e-04 -9.42e-06 -1.52e-04 -2.02e-04 -6.71e-05 -8.14e-05 6.22e-05 -1.32e-21 -1.37e-05 -9.90e-05 4.69e-04 3.70e-04 3.48e-04 6.35e-04 -2.38e-04 1.37e-05 6.56e-22
    0.00e+00 7.19e-05 -1.93e-04 -3.18e-04 2.72e-06 4.56e-04 -4.58e-05 1.04e-04 -7.19e-05 7.80e-22 -1.16e-04 -5.12e-04 1.72e-04 2.41e-04 1.99e-04 1.60e-04 1.93e-04 1.16e-04 3.07e-21 4.13e-05 9.07e-06 -3.05e-04 1.65e-04 3.43e-04 3.18e-04 5.12e-04 -4.13e-05 7.29e-22 -4.66e-04 -5.13e-05 -2.89e-06 4.87e-04 -2.72e-06 -1.72e-04 -9.07e-06 4.66e-04 2.19e-22 -4.79e-05 -8.87e-05 -3.91e-04 -4.56e-04 -2.41e-04 3.05e-04 5.13e-05 4.79e-05 -8.56e-22 -1.03e-04 -5.22e-04 4.58e-05 -1.99e-04 -1.65e-04 2.89e-06 8.87e-05 1.03e-04 1.85e-22 6.38e-05 -1.04e-04 -1.60e-04 -3.43e-04 -4.87e-04 3.91e-04 5.22e-04 -6.38e-05 -9.47e-22
    0.00e+00 1.07e-05 2.06e-04 5.58e-05 -4.54e-04 7.42e-04 1.30e-04 7.62e-05 -1.07e-05 -2.37e-21 1.57e-04 2.00e-04 4.24e-05 1.76e-04 2.15e-04 -1.77e-04 -2.06e-04 -1.57e-04 -6.27e-21 1.60e-04 -4.17e-05 2.96e-04 1.25e-04 6.43e-05 -5.58e-05 -2.00e-04 -1.60e-04 -3.09e-21 -8.14e-05 1.40e-04 -1.49e-04 3.22e-04 4.54e-04 -4.24e-05 4.17e-05 8.14e-05 -2.21e-21 -3.13e-04 -5.31e-05 -6.63e-05 -7.42e-04 -1.76e-04 -2.96e-04 -1.40e-04 3.13e-04 -4.07e-22 -1.99e-04 2.94e-04 -1.30e-04 -2.15e-04 -1.25e-04 1.49e-04 5.31e-05 1.99e-04 2.65e-22 4.22e-05 -7.62e-05 1.77e-04 -6.43e-05 -3.22e-04 6.63e-05 -2.94e-04 -4.22e-05 1.23e-21
    0.00e+00 1.01e-04 -2.45e-04 -1.27e-04 -8.07e-05 -4.46e-05 1.94e-04 -2.15e-04 -1.01e-04 -1.04e-22 -3.31e-04 6.03e-05 -5.64e-05 1.50e-04 -4.87e-04 5.99e-05 2.45e-04 3.31e-04 -6.11e-22 -2.62e-05 -1.67e-04 1.43e-04 1.08e-04 8.06e-05 1.27e-04 -6.03e-05 2.62e-05 7.16e-22 -3.55e-05 -1.97e-04 2.60e-04 -1.50e-04 8.07e-05 5.64e-05 1.67e-04 3.55e-05 -2.12e-21 -8.08e-05 2.26e-04 -7.14e-05 4.46e-05 -1.50e-04 -1.43e-04 1.97e-04 8.08e-05 1.68e-21 1.91e-04 -2.03e-04 -1.94e-04 4.87e-04 -1.08e-04 -2.60e-04 -2.26e-04 -1.91e-04 -4.64e-21 -1.02e-04 2.15e-04 -5.99e-05 -8.06e-05 1.50e-04 7.14e-05 2.03e-04 1.02e-04 1.33e-21
    0.00e+00 -2.00e-04 -1.72e-04 -2.59e-04 9.32e-06 2.98e-04 -8.27e-06 3.75e-04 2.00e-04 1.60e-21 -1.44e-04 -2.23e-04 -1.83e-04 1.89e-04 1.64e-04 2.47e-04 1.72e-04 1.44e-04 8.63e-22 1.57e-04 -1.38e-04 6.38e-05 8.99e-05 1.04e-04 2.59e-04 2.23e-04 -1.57e-04 2.69e-21 -3.80e-04 1.13e-04 2.76e-04 2.62e-04 -9.32e-06 1.83e-04 1.38e-04 3.80e-04 -1.07e-21 -9.47e-05 -4.69e-07 -4.47e-04 -2.98e-04 -1.89e-04 -6.38e-05 -1.13e-04 9.47e-05 -1.70e-21 -9.58e-05 -1.11e-06 8.27e-06 -1.64e-04 -8.99e-05 -2.76e-04 4.69e-07 9.58e-05 1.10e-21 2.68e-04 -3.75e-04 -2.47e-04 -1.04e-04 -2.62e-04 4.47e-04 1.11e-06 -2.68e-04 -1.24e-21
    0.00e+00 -1.48e-04 -1.30e-04 -2.26e-04 -1.00e-04 2.15e-04 1.28e-04 1.54e-04 1.48e-04 1.90e-21 -1.06e-04 -2.55e-04 -1.75e-04 1.66e-04 2.22e-04 -7.78e-05 1.30e-04 1.06e-04 1.77e-21 1.87e-04 -1.25e-04 6.44e-05 1.04e-04 -1.70e-04 2.26e-04 2.55e-04 -1.87e-04 -5.00e-22 -2.63e-04 1.06e-05 1.51e-04 -2.24e-04 1.00e-04 1.75e-04 1.25e-04 2.63e-04 5.67e-22 -1.62e-05 -6.11e-06 -4.29e-04 -2.15e-04 -1.66e-04 -6.44e-05 -1.06e-05 1.62e-05 1.39e-21 -4.98e-05 2.31e-04 -1.28e-04 -2.22e-04 -1.04e-04 -1.51e-04 6.11e-06 4.98e-05 -2.48e-21 2.82e-04 -1.54e-04 7.78e-05 1.70e-04 2.24e-04 4.29e-04 -2.31e-04 -2.82e-04 -1.27e-21
    0.00e+00 -1.39e-04 3.75e-04 -2.69e-04 -1.02e-04 5.75e-04 9.34e-05 -2.68e-04 1.39e-04 -4.69e-22 -5.64e-05 -4.67e-04 5.65e-05 9.44e-05 1.78e-04 1.97e-05 -3.75e-04 5.64e-05 1.80e-21 1.30e-04 8.14e-05 -2.37e-04 3.78e-05 6.68e-05 2.69e-04 4.67e-04 -1.30e-04 -1.39e-21 -3.43e-04 -9.59e-05 1.79e-04 7.72e-04 1.02e-04 -5.65e-05 -8.14e-05 3.43e-04 -2.50e-21 -2.27e-05 -1.20e-04 -1.72e-04 -5.75e-04 -9.44e-05 2.37e-04 9.59e-05 2.27e-05 1.08e-21 -2.02e-06 -1.17e-04 -9.34e-05 -1.78e-04 -3.78e-05 -1.79e-04 1.20e-04 2.02e-06 3.01e-21 -1.09e-04 2.68e-04 -1.97e-05 -6.68e-05 -7.72e-04 1.72e-04 1.17e-04 1.09e-04 2.69e-21
    0.00e+00 3.71e-04 5.28e-04 5.52e-05 -4.03e-04 -3.92e-04 -7.80e-05 -1.19e-04 -3.71e-04 -9.85e-22 1.33e-04 1.88e-04 1.08e-04 -3.08e-04 2.08e-04 -1.16e-04 -5.28e-04 -1.33e-04 -9.48e-22 2.19e-04 -1.63e-04 -5.82e-04 1.50e-04 6.87e-05 -5.52e-05 -1.88e-04 -2.19e-04 6.90e-22 -3.49e-05 -5.67e-06 -2.00e-05 1.68e-04 4.03e-04 -1.08e-04 1.63e-04 3.49e-05 -9.69e-22 -1.97e-04 -1.45e-04 1.45e-04 3.92e-04 3.08e-04 5.82e-04 5.67e-06 1.97e-04 -1.82e-21 -6.78e-05 -1.80e-04 7.80e-05 -2.08e-04 -1.50e-04 2.00e-05 1.45e-04 6.78e-05 1.20e-21 6.18e-05 1.19e-04 1.16e-04 -6.87e-05 -1.68e-04 -1.45e-04 1.80e-04 -6.18e-05 2.94e-21
    0.00e+00 1.06e-04 1.92e-04 4.12e-05 -1.57e-04 1.86e-04 1.29e-04 4.02e-04 -1.06e-04 1.76e-21 6.98e-05 6.59e-05 -3.81e-05 -1.74e-04 1.40e-04 3.16e-05 -1.92e-04 -6.98e-05 -1.09e-22 2.55e-05 -1.34e-04 2.16e-05 -5.32e-05 -5.35e-05 -4.12e-05 -6.59e-05 -2.55e-05 6.90e-21 1.35e-04 2.48e-04 2.04e-04 3.01e-04 1.57e-04 3.81e-05 1.34e-04 -1.35e-04 1.65e-22 2.79e-04 3.91e-04 6.36e-04 -1.86e-04 1.74e-04 -2.16e-05 -2.48e-04 -2.79e-04 -4.24e-22 5.31e-05 1.70e-04 -1.29e-04 -1.40e-04 5.32e-05 -2.04e-04 -3.91e-04 -5.31e-05 -5.21e-23 7.82e-05 -4.02e-04 -3.16e-05 5.35e-05 -3.01e-04 -6.36e-04 -1.70e-04 -7.82e-05 -3.97e-23
    0.00e+00 5.61e-05 1.69e-04 1.59e-05 -4.38e-04 6.33e-04 7.62e-05 8.00e-05 -5.61e-05 2.93e-22 1.02e-04 3.25e-04 8.22e-05 1.87e-04 2.48e-04 -1.95e-05 -1.69e-04 -1.02e-04 1.12e-21 1.36e-04 -6.26e-05 3.30e-04 2.02e-04 1.95e-04 -1.59e-05 -3.25e-04 -1.36e-04 1.08e-21 -2.06e-05 3.75e-04 -8.87e-05 3.03e-04 4.38e-04 -8.22e-05 6.26e-05 2.06e-05 -3.20e-22 -4.26e-04 -6.45e-05 -8.85e-05 -6.33e-04 -1.87e-04 -3.30e-04 -3.75e-04 4.26e-04 -4.11e-22 -2.59e-04 2.40e-04 -7.62e-05 -2.48e-04 -2.02e-04 8.87e-05 6.45e-05 2.59e-04 -1.81e-21 1.05e-04 -8.00e-05 1.95e-05 -1.95e-04 -3.03e-04 8.85e-05 -2.40e-04 -1.05e-04 -1.52e-21
    0.00e+00 7.11e-05 -3.17e-04 -3.49e-04 2.30e-05 4.10e-04 -3.90e-05 3.41e-04 -7.11e-05 1.76e-21 -1.00e-04 -5.62e-04 1.40e-04 2.90e-04 1.54e-04 2.13e-04 3.17e-04 1.00e-04 5.68e-21 5.54e-05 -6.37e-05 -1.01e-04 1.42e-04 4.33e-04 3.49e-04 5.62e-04 -5.54e-05 -1.12e-21 -4.08e-04 -5.82e-05 2.04e-05 4.68e-04 -2.30e-05 -1.40e-04 6.37e-05 4.08e-04 -5.20e-22 -1.25e-04 -8.45e-05 -4.30e-04 -4.10e-04 -2.90e-04 1.01e-04 5.82e-05 1.25e-04 -5.76e-22 -1.56e-04 -4.71e-04 3.90e-05 -1.54e-04 -1.42e-04 -2.04e-05 8.45e-05 1.56e-04 -2.74e-22 1.90e-04 -3.41e-04 -2.13e-04 -4.33e-04 -4.68e-04 4.30e-04 4.71e-04 -1.90e-04 1.26e-21
    0.00e+00 2.64e-04 -2.37e-04 1.06e-04 5.36e-06 -3.38e-04 6.05e-05 1.56e-04 -2.64e-04 1.32e-21 1.47e-04 2.64e-04 -1.22e-04 6.66e-05 -5.30e-04 2.59e-04 2.37e-04 -1.47e-04 6.05e-22 -7.47e-05 -1.64e-04 3.60e-04 -3.37e-04 7.54e-05 -1.06e-04 -2.64e-04 7.47e-05 -4.08e-22 -3.34e-05 1.85e-05 3.47e-04 -2.11e-04 -5.36e-06 1.22e-04 1.64e-04 3.34e-05 3.24e-21 5.45e-05 -3.05e-04 1.40e-04 3.38e-04 -6.66e-05 -3.60e-04 -1.85e-05 -5.45e-05 2.43e-22 2.51e-04 6.92e-05 -6.05e-05 5.30e-04 3.37e-04 -3.47e-04 3.05e-04 -2.51e-04 4.63e-22 -5.21e-04 -1.56e-04 -2.59e-04 -7.54e-05 2.11e-04 -1.40e-04 -6.92e-05 5.21e-04 -3.75e-22
    0.00e+00 1.55e-04 1.31e-04 5.94e-05 -2.70e-05 -3.63e-04 -1.02e-04 2.65e-05 -1.55e-04 2.70e-20 7.26e-05 1.16e-04 8.59e-05 -3.61e-04 1.02e-04 1.25e-04 -1.31e-04 -7.26e-05 -5.01e-23 2.13e-05 -1.12e-04 -2.42e-05 -5.36e-05 4.12e-05 -5.94e-05 -1.16e-04 -2.13e-05 1.57e-22 -4.71e-07 1.12e-04 3.11e-05 4.41e-05 2.70e-05 -8.59e-05 1.12e-04 4.71e-07 -1.90e-22 -6.25e-04 8.85e-05 1.05e-04 3.63e-04 3.61e-04 2.42e-05 -1.12e-04 6.25e-04 2.84e-22 -1.70e-04 -1.62e-04 1.02e-04 -1.02e-04 5.36e-05 -3.11e-05 -8.85e-05 1.70e-04 4.36e-22 -6.54e-06 -2.65e-05 -1.25e-04 -4.12e-05 -4.41e-05 -1.05e-04 1.62e-04 6.54e-06 -7.54e-22
    0.00e+00 2.86e-04 -7.35e-05 9.06e-05 9.21e-05 -4.49e-04 4.49e-04 4.19e-05 -2.86e-04 1.33e-21 1.32e-04 5.76e-04 -3.97e-04 2.59e-04 -4.09e-05 1.08e-05 7.35e-05 -1.32e-04 2.49e-21 -4.02e-05 -1.10e-04 2.79e-04 -3.52e-04 6.67e-05 -9.06e-05 -5.76e-04 4.02e-05 -2.45e-22 -1.31e-04 -1.85e-04 2.03e-04 -3.01e-04 -9.21e-05 3.97e-04 1.10e-04 1.31e-04 -5.83e-22 -3.61e-05 2.55e-04 4.20e-05 4.49e-04 -2.59e-04 -2.79e-04 1.85e-04 3.61e-05 -2.75e-21 -2.52e-04 2.05e-04 -4.49e-04 4.09e-05 3.52e-04 -2.03e-04 -2.55e-04 2.52e-04 1.93e-22 -4.61e-04 -4.19e-05 -1.08e-05 -6.67e-05 3.01e-04 -4.20e-05 -2.05e-04 4.61e-04 2.66e-22
    0.00e+00 9.58e-05 -7.89e-05 -1.06e-05 -8.24e-07 -2.74e-04 -2.55e-05 1.54e-04 -9.58e-05 2.97e-21 -5.39e-05 2.54e-04 -2.49e-04 3.99e-04 -6.69e-05 -3.15e-05 7.89e-05 5.39e-05 -3.46e-21 1.25e-04 -3.34e-04 -8.35e-06 -3.72e-04 2.98e-04 1.06e-05 -2.54e-04 -1.25e-04 -1.67e-22 9.98e-05 8.74e-05 7.52e-05 -2.79e-04 8.24e-07 2.49e-04 3.34e-04 -9.98e-05 -1.18e-21 5.28e-05 1.56e-04 5.87e-05 2.74e-04 -3.99e-04 8.35e-06 -8.74e-05 -5.28e-05 -7.76e-23 -1.57e-04 -8.11e-05 2.55e-05 6.69e-05 3.72e-04 -7.52e-05 -1.56e-04 1.57e-04 -2.53e-22 5.86e-06 -1.54e-04 3.15e-05 -2.98e-04 2.79e-04 -5.87e-05 8.11e-05 -5.86e-06 2.91e-21
    0.00e+00 1.70e-04 3.36e-04 5.57e-05 -3.33e-04 4.25e-04 -4.77e-05 -7.08e-05 -1.70e-04 -4.26e-22 -9.53e-06 2.35e-04 1.87e-04 4.44e-04 4.74e-04 1.43e-04 -3.36e-04 9.53e-06 4.02e-22 -6.97e-05 1.45e-04 2.58e-04 1.79e-04 2.88e-04 -5.57e-05 -2.35e-04 6.97e-05 -3.41e-21 -5.16e-05 3.94e-04 1.01e-04 2.46e-04 3.33e-04 -1.87e-04 -1.45e-04 5.16e-05 -3.21e-22 -3.98e-04 -1.16e-04 -1.09e-04 -4.25e-04 -4.44e-04 -2.58e-04 -3.94e-04 3.98e-04 -1.63e-22 -3.64e-04 2.63e-04 4.77e-05 -4.74e-04 -1.79e-04 -1.01e-04 1.16e-04 3.64e-04 -3.81e-22 2.35e-04 7.08e-05 -1.43e-04 -2.88e-04 -2.46e-04 1.09e-04 -2.63e-04 -2.35e-04 2.86e-21
    0.00e+00 2.71e-04 -3.79e-04 -2.86e-04 -2.53e-05 2.26e-04 3.29e-05 3.32e-04 -2.71e-04 -7.71e-22 -1.05e-04 -9.36e-05 2.61e-04 3.66e-04 2.25e-04 -1.61e-04 3.79e-04 1.05e-04 -9.99e-22 -3.89e-05 -6.65e-05 -4.75e-05 3.37e-04 7.87e-04 2.86e-04 9.36e-05 3.89e-05 -2.61e-21 -2.62e-04 -2.78e-04 -1.74e-04 1.74e-04 2.53e-05 -2.61e-04 6.65e-05 2.62e-04 -9.28e-22 2.92e-06 3.95e-05 -3.21e-04 -2.26e-04 -3.66e-04 4.75e-05 2.78e-04 -2.92e-06 -1.83e-21 -8.43e-05 -4.61e-04 -3.29e-05 -2.25e-04 -3.37e-04 1.74e-04 -3.95e-05 8.43e-05 3.17e-21 1.06e-04 -3.32e-04 1.61e-04 -7.87e-04 -1.74e-04 3.21e-04 4.61e-04 -1.06e-04 -3.32e-21
    0.00e+00 -2.73e-04 3.88e-04 -1.08e-04 -2.96e-04 7.31e-04 6.91e-05 -2.61e-04 2.73e-04 -1.10e-21 3.24e-05 -4.00e-04 -5.15e-05 5.42e-05 1.50e-04 3.31e-05 -3.88e-04 -3.24e-05 2.30e-21 3.01e-04 -3.35e-05 3.92e-05 -1.79e-06 4.55e-06 1.08e-04 4.00e-04 -3.01e-04 -8.90e-23 -3.10e-04 -1.52e-04 1.78e-04 8.31e-04 2.96e-04 5.15e-05 3.35e-05 3.10e-04 -2.03e-21 -2.65e-05 -4.90e-05 -4.20e-05 -7.31e-04 -5.42e-05 -3.92e-05 1.52e-04 2.65e-05 8.87e-22 -8.02e-05 3.40e-05 -6.91e-05 -1.50e-04 1.79e-06 -1.78e-04 4.90e-05 8.02e-05 -9.59e-22 -1.35e-04 2.61e-04 -3.31e-05 -4.55e-06 -8.31e-04 4.20e-05 -3.40e-05 1.35e-04 6.07e-21
    0.00e+00 1.56e-04 -1.83e-04 -3.32e-04 -5.77e-05 2.85e-04 1.25e-04 4.66e-04 -1.56e-04 6.48e-22 -5.35e-05 -3.48e-04 1.56e-04 2.87e-04 5.69e-05 1.73e-04 1.83e-04 5.35e-05 -1.52e-21 -2.64e-05 -1.13e-04 1.27e-04 2.69e-04 7.02e-04 3.32e-04 3.48e-04 2.64e-05 -4.68e-21 -3.07e-04 -7.30e-05 -2.20e-05 1.53e-04 5.77e-05 -1.56e-04 1.13e-04 3.07e-04 -7.12e-22 -1.12e-04 -7.45e-06 -4.10e-04 -2.85e-04 -2.87e-04 -1.27e-04 7.30e-05 1.12e-04 1.25e-23 -5.03e-05 -2.26e-04 -1.25e-04 -5.69e-05 -2.69e-04 2.20e-05 7.45e-06 5.03e-05 -2.16e-21 2.29e-04 -4.66e-04 -1.73e-04 -7.02e-04 -1.53e-04 4.10e-04 2.26e-04 -2.29e-04 -2.71e-21
    0.00e+00 -7.85e-05 1.99e-04 1.07e-04 -5.02e-04 7.90e-04 1.16e-04 -2.11e-06 7.85e-05 -4.97e-21 1.92e-04 -5.09e-05 -3.39e-05 8.60e-05 2.02e-04 -2.71e-04 -1.99e-04 -1.92e-04 1.20e-21 2.24e-04 -6.73e-05 2.94e-04 7.04e-05 5.47e-05 -1.07e-04 5.09e-05 -2.24e-04 1.43e-21 -1.25e-04 -9.34e-06 -7.75e-05 5.13e-04 5.02e-04 3.39e-05 6.73e-05 1.25e-04 -4.85e-22 -2.46e-04 -2.52e-05 1.17e-05 -7.90e-04 -8.60e-05 -2.94e-04 9.34e-06 2.46e-04 1.41e-21 -1.58e-04 2.79e-04 -1.16e-04 -2.02e-04 -7.04e-05 7.75e-05 2.52e-05 1.58e-04 6.29e-22 -9.55e-05 2.11e-06 2.71e-04 -5.47e-05 -5.13e-04 -1.17e-05 -2.79e-04 9.55e-05 -2.85e-22
    0.00e+00 -1.64e-04 -2.09e-04 -1.93e-04 -1.13e-04 1.95e-04 1.35e-04 -6.68e-05 1.64e-04 -1.14e-22 -1.62e-04 -3.70e-04 -1.80e-04 1.70e-04 2.02e-04 -3.45e-04 2.09e-04 1.62e-04 7.81e-22 1.66e-04 -5.30e-05 4.04e-05 1.14e-04 -3.51e-04 1.93e-04 3.70e-04 -1.66e-04 1.34e-21 -1.93e-04 -6.34e-05 7.02e-05 -3.75e-04 1.13e-04 1.80e-04 5.30e-05 1.93e-04 -6.00e-22 2.29e-05 -2.12e-05 -4.32e-04 -1.95e-04 -1.70e-04 -4.04e-05 6.34e-05 -2.29e-05 -5.12e-22 -1.11e-05 2.62e-04 -1.35e-04 -2.02e-04 -1.14e-04 -7.02e-05 2.12e-05 1.11e-05 3.13e-22 2.98e-04 6.68e-05 3.45e-04 3.51e-04 3.75e-04 4.32e-04 -2.62e-04 -2.98e-04 6.18e-22
    0.00e+00 2.52e-04 -3.88e-04 -1.20e-04 -4.64e-05 9.73e-05 2.24e-05 1.83e-04 -2.52e-04 -1.68e-22 -4.35e-05 2.72e-04 1.74e-04 5.40e-04 4.69e-04 -2.87e-04 3.88e-04 4.35e-05 -4.20e-22 -1.16e-04 -5.64e-05 -2.07e-04 3.00e-04 8.50e-04 1.20e-04 -2.72e-04 1.16e-04 -2.08e-21 -1.14e-04 -1.43e-04 -1.03e-04 5.00e-05 4.64e-05 -1.74e-04 5.64e-05 1.14e-04 -1.18e-22 3.74e-05 6.46e-05 -1.43e-04 -9.73e-05 -5.40e-04 2.07e-04 1.43e-04 -3.74e-05 1.24e-21 -4.30e-05 -2.44e-04 -2.24e-05 -4.69e-04 -3.00e-04 1.03e-04 -6.46e-05 4.30e-05 4.12e-21 7.92e-05 -1.83e-04 2.87e-04 -8.50e-04 -5.00e-05 1.43e-04 2.44e-04 -7.92e-05 -1.01e-21
    0.00e+00 6.62e-05 3.04e-05 -2.88e-04 -9.78e-07 4.54e-04 -2.14e-05 -1.44e-04 -6.62e-05 -1.54e-21 -5.19e-05 -3.94e-04 2.06e-04 2.06e-04 2.26e-04 -3.01e-05 -3.04e-05 5.19e-05 2.05e-21 3.60e-06 1.21e-04 -3.82e-04 1.41e-04 3.20e-04 2.88e-04 3.94e-04 -3.60e-06 3.46e-21 -3.96e-04 -1.31e-04 -4.13e-05 5.15e-04 9.78e-07 -2.06e-04 -1.21e-04 3.96e-04 1.28e-21 -1.87e-05 -1.12e-04 -2.54e-04 -4.54e-04 -2.06e-04 3.82e-04 1.31e-04 1.87e-05 -5.48e-22 -4.00e-05 -3.66e-04 2.14e-05 -2.26e-04 -1.41e-04 4.13e-05 1.12e-04 4.00e-05 -3.03e-23 -6.68e-05 1.44e-04 3.01e-05 -3.20e-04 -5.15e-04 2.54e-04 3.66e-04 6.68e-05 1.98e-21
    0.00e+00 -4.28e-06 -5.60e-05 -3.26e-06 -5.57e-05 -8.97e-05 3.55e-04 -1.02e-04 4.28e-06 2.92e-22 -2.86e-05 7.06e-04 -3.14e-04 -2.26e-04 -3.77e-04 -1.33e-04 5.60e-05 2.86e-05 1.39e-21 -8.08e-05 8.23e-05 4.28e-05 -1.42e-04 1.50e-04 3.26e-06 -7.06e-04 8.08e-05 7.59e-22 -1.79e-04 -3.06e-04 -8.62e-05 -9.88e-05 5.57e-05 3.14e-04 -8.23e-05 1.79e-04 2.95e-22 5.94e-05 1.21e-04 -2.19e-05 8.97e-05 2.26e-04 -4.28e-05 3.06e-04 -5.94e-05 -5.08e-22 -3.36e-04 1.02e-05 -3.55e-04 3.77e-04 1.42e-04 8.62e-05 -1.21e-04 3.36e-04 -3.16e-21 -2.64e-04 1.02e-04 1.33e-04 -1.50e-04 9.88e-05 2.19e-05 -1.02e-05 2.64e-04 -1.33e-21
    0.00e+00 -2.36e-04 1.35e-04 1.70e-04 2.23e-04 -4.04e-04 8.63e-05 2.57e-04 2.36e-04 6.49e-22 6.47e-05 3.49e-05 1.40e-04 -3.91e-04 -2.20e-05 -8.01e-05 -1.35e-04 -6.47e-05 6.18e-22 3.54e-05 -7.06e-05 -5.05e-05 6.62e-05 -8.51e-05 -1.70e-04 -3.49e-05 -3.54e-05 1.68e-22 -8.03e-05 1.30e-04 9.11e-05 1.96e-06 -2.23e-04 -1.40e-04 7.06e-05 8.03e-05 -6.31e-23 -1.35e-04 3.47e-05 7.26e-04 4.04e-04 3.91e-04 5.05e-05 -1.30e-04 1.35e-04 9.21e-23 -4.91e-04 -3.97e-04 -8.63e-05 2.20e-05 -6.62e-05 -9.11e-05 -3.47e-05 4.91e-04 1.78e-22 4.91e-05 -2.57e-04 8.01e-05 8.51e-05 -1.96e-06 -7.26e-04 3.97e-04 -4.91e-05 -8.36e-22
    0.00e+00 -1.65e-04 1.95e-04 5.77e-06 4.28e-04 1.53e-04 1.14e-04 4.93e-04 1.65e-04 -2.23e-21 -1.01e-04 -3.40e-04 -2.96e-04 -2.06e-04 1.68e-04 -3.59e-05 -1.95e-04 1.01e-04 -1.99e-21 1.28e-04 2.98e-05 3.95e-05 -3.90e-05 -2.86e-04 -5.77e-06 3.40e-04 -1.28e-04 7.81e-22 -7.33e-05 2.61e-04 3.35e-04 -5.01e-05 -4.28e-04 2.96e-04 -2.98e-05 7.33e-05 -3.26e-22 5.11e-04 8.33e-05 -5.07e-05 -1.53e-04 2.06e-04 -3.95e-05 -2.61e-04 -5.11e-04 7.03e-22 6.91e-05 -2.26e-04 -1.14e-04 -1.68e-04 3.90e-05 -3.35e-04 -8.33e-05 -6.91e-05 -1.35e-21 2.48e-04 -4.93e-04 3.59e-05 2.86e-04 5.01e-05 5.07e-05 2.26e-04 -2.48e-04 -1.12e-21
    0.00e+00 1.86e-04 -3.82e-04 -3.39e-04 -6.91e-06 2.91e-04 5.99e-05 4.02e-04 -1.86e-04 -1.00e-22 -5.60e-05 -3.59e-04 2.62e-04 3.36e-04 1.29e-04 4.67e-05 3.82e-04 5.60e-05 -1.05e-21 1.93e-05 -9.19e-05 1.09e-05 3.39e-04 6.60e-04 3.39e-04 3.59e-04 -1.93e-05 1.34e-21 -3.43e-04 -2.01e-04 -1.27e-04 2.24e-04 6.91e-06 -2.62e-04 9.19e-05 3.43e-04 -4.21e-22 -7.27e-05 -2.22e-05 -4.13e-04 -2.91e-04 -3.36e-04 -1.09e-05 2.01e-04 7.27e-05 -4.11e-21 -1.01e-04 -4.84e-04 -5.99e-05 -1.29e-04 -3.39e-04 1.27e-04 2.22e-05 1.01e-04 -1.61e-21 1.84e-04 -4.02e-04 -4.67e-05 -6.60e-04 -2.24e-04 4.13e-04 4.84e-04 -1.84e-04 -9.56e-22
    0.00e+00 1.41e-04 4.22e-04 -3.45e-05 -4.40e-04 6.10e-04 5.89e-05 -1.11e-05 -1.41e-04 -9.68e-23 1.17e-04 3.62e-04 1.37e-04 2.72e-04 3.19e-04 -1.47e-04 -4.22e-04 -1.17e-04 3.45e-21 8.52e-06 8.90e-05 2.63e-04 1.69e-04 2.17e-04 3.45e-05 -3.62e-04 -8.52e-06 -1.03e-21 -1.04e-04 2.50e-04 -1.66e-04 3.17e-04 4.40e-04 -1.37e-04 -8.90e-05 1.04e-04 -2.03e-21 -3.53e-04 -5.88e-05 -8.39e-05 -6.10e-04 -2.72e-04 -2.63e-04 -2.50e-04 3.53e-04 -2.63e-21 -2.50e-04 4.52e-04 -5.89e-05 -3.19e-04 -1.69e-04 1.66e-04 5.88e-05 2.50e-04 8.66e-22 1.07e-04 1.11e-05 1.47e-04 -2.17e-04 -3.17e-04 8.39e-05 -4.52e-04 -1.07e-04 1.61e-21
    0.00e+00 -1.35e-04 -1.35e-04 -2.49e-04 -5.83e-05 2.65e-04 9.33e-05 3.29e-04 1.35e-04 6.97e-22 -8.18e-05 -2.85e-04 -1.48e-04 2.00e-04 1.92e-04 1.44e-04 1.35e-04 8.18e-05 -1.25e-21 1.68e-04 -1.62e-04 7.06e-05 1.40e-04 5.41e-05 2.49e-04 2.85e-04 -1.68e-04 8.81e-22 -3.38e-04 3.68e-05 1.84e-04 -1.85e-06 5.83e-05 1.48e-04 1.62e-04 3.38e-04 2.24e-21 -9.31e-05 -1.76e-05 -4.22e-04 -2.65e-04 -2.00e-04 -7.06e-05 -3.68e-05 9.31e-05 1.14e-21 -9.08e-05 1.41e-04 -9.33e-05 -1.92e-04 -1.40e-04 -1.84e-04 1.76e-05 9.08e-05 8.30e-22 3.28e-04 -3.29e-04 -1.44e-04 -5.41e-05 1.85e-06 4.22e-04 -1.41e-04 -3.28e-04 3.22e-22
    0.00e+00 1.05e-06 2.40e-04 7.74e-05 -3.69e-04 1.61e-04 -3.31e-04 -1.26e-04 -1.05e-06 -1.58e-21 1.77e-04 5.88e-04 -1.29e-05 -3.08e-04 2.23e-04 -6.83e-05 -2.40e-04 -1.77e-04 1.39e-21 -3.19e-05 -4.03e-05 -1.54e-04 -4.32e-05 1.19e-04 -7.74e-05 -5.88e-04 3.19e-05 -1.49e-21 -1.61e-04 -4.81e-04 6.51e-05 5.60e-05 3.69e-04 1.29e-05 4.03e-05 1.61e-04 9.95e-22 3.99e-05 1.29e-04 -5.71e-05 -1.61e-04 3.08e-04 1.54e-04 4.81e-04 -3.99e-05 2.35e-21 -3.01e-04 -1.00e-04 3.31e-04 -2.23e-04 4.32e-05 -6.51e-05 -1.29e-04 3.01e-04 -9.45e-22 -1.06e-04 1.26e-04 6.83e-05 -1.19e-04 -5.60e-05 5.71e-05 1.00e-04 1.06e-04 2.41e-21
    0.00e+00 -1.58e-04 -1.90e-04 -2.51e-04 -5.32e-05 1.33e-04 1.95e-04 -2.41e-04 1.58e-04 7.05e-22 -2.59e-04 -3.59e-04 -2.00e-04 1.74e-04 1.40e-04 -4.83e-04 1.90e-04 2.59e-04 -1.37e-22 1.12e-04 8.68e-05 9.14e-05 2.28e-05 -2.87e-04 2.51e-04 3.59e-04 -1.12e-04 -8.90e-22 -9.12e-05 -3.96e-05 -7.90e-05 -2.87e-04 5.32e-05 2.00e-04 -8.68e-05 9.12e-05 3.78e-22 9.88e-05 5.59e-07 -4.10e-04 -1.33e-04 -1.74e-04 -9.14e-05 3.96e-05 -9.88e-05 1.32e-21 1.47e-06 1.79e-04 -1.95e-04 -1.40e-04 -2.28e-05 7.90e-05 -5.59e-07 -1.47e-06 -3.34e-21 2.68e-04 2.41e-04 4.83e-04 2.87e-04 2.87e-04 4.10e-04 -1.79e-04 -2.68e-04 3.71e-21
    0.00e+00 1.17e-04 2.46e-04 3.66e-06 -3.98e-04 5.35e-04 3.63e-05 1.99e-05 -1.17e-04 -9.79e-22 3.34e-05 3.15e-04 1.67e-04 2.74e-04 3.49e-04 7.03e-05 -2.46e-04 -3.34e-05 -2.97e-22 4.88e-05 1.35e-05 2.92e-04 2.38e-04 2.49e-04 -3.66e-06 -3.15e-04 -4.88e-05 -2.02e-21 -5.15e-05 4.17e-04 -4.31e-05 2.84e-04 3.98e-04 -1.67e-04 -1.35e-05 5.15e-05 -2.67e-22 -4.37e-04 -1.10e-04 -1.19e-04 -5.35e-04 -2.74e-04 -2.92e-04 -4.17e-04 4.37e-04 -6.66e-22 -3.26e-04 2.20e-04 -3.63e-05 -3.49e-04 -2.38e-04 4.31e-05 1.10e-04 3.26e-04 9.46e-22 1.83e-04 -1.99e-05 -7.03e-05 -2.49e-04 -2.84e-04 1.19e-04 -2.20e-04 -1.83e-04 4.23e-22
    0.00e+00 1.68e-04 3.88e-04 3.98e-05 -3.16e-04 3.79e-04 -7.01e-05 -2.23e-04 -1.68e-04 -2.01e-21 4.69e-05 3.81e-04 1.69e-04 4.85e-04 5.94e-04 3.46e-05 -3.88e-04 -4.69e-05 1.91e-21 -9.98e-05 2.59e-04 6.25e-05 1.77e-04 3.33e-04 -3.98e-05 -3.81e-04 9.98e-05 2.58e-22 -8.99e-05 3.15e-04 -1.24e-05 1.70e-04 3.16e-04 -1.69e-04 -2.59e-04 8.99e-05 3.89e-22 -2.63e-04 -2.11e-05 -4.30e-05 -3.79e-04 -4.85e-04 -6.25e-05 -3.15e-04 2.63e-04 2.15e-21 -3.02e-04 5.01e-04 7.01e-05 -5.94e-04 -1.77e-04 1.24e-05 2.11e-05 3.02e-04 1.62e-22 2.07e-04 2.23e-04 -3.46e-05 -3.33e-04 -1.70e-04 4.30e-05 -5.01e-04 -2.07e-04 3.63e-21
    0.00e+00 1.97e-04 -3.10e-04 -4.61e-06 -8.75e-05 1.15e-04 -6.53e-05 -1.09e-04 -1.97e-04 2.32e-21 3.03e-05 3.20e-04 1.60e-04 6.52e-04 6.85e-04 -1.84e-04 3.10e-04 -3.03e-05 -3.92e-22 -1.38e-04 3.58e-05 -2.84e-04 2.24e-04 7.58e-04 4.61e-06 -3.20e-04 1.38e-04 -1.77e-22 -5.30e-05 6.31e-05 -7.37e-06 2.00e-05 8.75e-05 -1.60e-04 -3.58e-05 5.30e-05 4.14e-22 1.24e-06 8.52e-05 1.19e-05 -1.15e-04 -6.52e-04 2.84e-04 -6.31e-05 -1.24e-06 1.65e-21 -7.40e-05 4.43e-05 6.53e-05 -6.85e-04 -2.24e-04 7.37e-06 -8.52e-05 7.40e-05 -3.16e-21 1.02e-04 1.09e-04 1.84e-04 -7.58e-04 -2.00e-05 -1.19e-05 -4.43e-05 -1.02e-04 1.16e-22
    0.00e+00 -2.38e-04 1.98e-04 1.94e-04 -3.50e-04 7.15e-04 9.58e-05 -6.95e-05 2.38e-04 -1.13e-21 1.32e-04 -1.92e-04 -5.45e-05 1.11e-04 1.11e-04 -1.73e-04 -1.98e-04 -1.32e-04 1.34e-21 2.42e-04 -7.83e-05 1.41e-04 -2.07e-05 5.09e-05 -1.94e-04 1.92e-04 -2.42e-04 -1.37e-21 -9.30e-05 -1.46e-04 4.17e-05 4.92e-04 3.50e-04 5.45e-05 7.83e-05 9.30e-05 -1.70e-22 -1.33e-04 1.08e-05 -2.46e-05 -7.15e-04 -1.11e-04 -1.41e-04 1.46e-04 1.33e-04 1.27e-20 -6.66e-05 1.76e-04 -9.58e-05 -1.11e-04 2.07e-05 -4.17e-05 -1.08e-05 6.66e-05 -4.62e-22 -1.00e-04 6.95e-05 1.73e-04 -5.09e-05 -4.92e-04 2.46e-05 -1.76e-04 1.00e-04 -1.22e-21
    0.00e+00 3.35e-04 -1.41e-04 1.06e-04 1.32e-04 -4.56e-04 2.06e-04 1.58e-04 -3.35e-04 1.86e-23 1.59e-04 4.45e-04 -3.76e-04 2.69e-04 1.08e-05 5.82e-05 1.41e-04 -1.59e-04 2.12e-21 -5.36e-05 -1.52e-04 3.26e-04 -3.69e-04 1.05e-04 -1.06e-04 -4.45e-04 5.36e-05 1.11e-21 -8.73e-05 -1.03e-04 3.52e-04 -3.22e-04 -1.32e-04 3.76e-04 1.52e-04 8.73e-05 -1.09e-21 -1.64e-05 9.50e-05 6.78e-05 4.56e-04 -2.69e-04 -3.26e-04 1.03e-04 1.64e-05 -7.07e-23 -2.53e-04 2.33e-04 -2.06e-04 -1.08e-05 3.69e-04 -3.52e-04 -9.50e-05 2.53e-04 -9.10e-22 -4.19e-04 -1.58e-04 -5.82e-05 -1.05e-04 3.22e-04 -6.78e-05 -2.33e-04 4.19e-04 8.60e-22
    0.00e+00 1.66e-04 -9.90e-05 -2.05e-04 -8.70e-05 -3.52e-05 1.20e-04 4.42e-04 -1.66e-04 -2.58e-22 9.48e-06 2.24e-04 1.63e-04 4.48e-04 2.84e-04 -2.71e-04 9.90e-05 -9.48e-06 3.81e-21 -2.54e-04 -6.69e-05 -2.48e-04 1.84e-04 5.83e-04 2.05e-04 -2.24e-04 2.54e-04 -5.52e-22 -1.09e-04 3.15e-05 1.62e-04 -3.52e-05 8.70e-05 -1.63e-04 6.69e-05 1.09e-04 -3.23e-22 -2.22e-05 5.44e-06 -2.73e-04 3.52e-05 -4.48e-04 2.48e-04 -3.15e-05 2.22e-05 -1.74e-21 2.56e-04 -1.32e-04 -1.20e-04 -2.84e-04 -1.84e-04 -1.62e-04 -5.44e-06 -2.56e-04 -2.71e-21 -8.64e-05 -4.42e-04 2.71e-04 -5.83e-04 3.52e-05 2.73e-04 1.32e-04 8.64e-05 -2.69e-21
    0.00e+00 -2.58e-04 -5.21e-05 -1.88e-04 5.64e-05 3.09e-04 -8.82e-05 2.28e-04 2.58e-04 6.57e-22 -2.92e-04 -8.34e-05 -1.86e-04 9.94e-05 2.34e-04 2.76e-04 5.21e-05 2.92e-04 3.08e-22 1.64e-04 -8.48e-05 -1.93e-04 1.38e-04 -1.47e-04 1.88e-04 8.34e-05 -1.64e-04 -1.47e-21 -4.98e-04 1.97e-04 2.09e-04 3.02e-04 -5.64e-05 1.86e-04 8.48e-05 4.98e-04 2.39e-22 -1.45e-06 1.08e-05 -3.88e-04 -3.09e-04 -9.94e-05 1.93e-04 -1.97e-04 1.45e-06 3.28e-21 -6.79e-05 -6.83e-05 8.82e-05 -2.34e-04 -1.38e-04 -2.09e-04 -1.08e-05 6.79e-05 -1.08e-21 1.83e-04 -2.28e-04 -2.76e-04 1.47e-04 -3.02e-04 3.88e-04 6.83e-05 -1.83e-04 9.26e-22
    0.00e+00 2.54e-04 -2.95e-04 -1.22e-04 -1.42e-04 -3.48e-04 1.10e-04 -3.45e-05 -2.54e-04 1.34e-21 -5.78e-07 2.62e-04 -1.39e-04 1.75e-04 -5.15e-04 3.60e-04 2.95e-04 5.78e-07 4.03e-22 -1.31e-04 -4.26e-04 1.24e-04 -2.13e-04 1.38e-04 1.22e-04 -2.62e-04 1.31e-04 -9.10e-22 2.20e-04 1.25e-04 2.76e-04 -2.51e-04 1.42e-04 1.39e-04 4.26e-04 -2.20e-04 1.33e-21 2.53e-05 -1.78e-04 -8.34e-06 3.48e-04 -1.75e-04 -1.24e-04 -1.25e-04 -2.53e-05 1.65e-22 4.55e-04 -2.65e-04 -1.10e-04 5.15e-04 2.13e-04 -2.76e-04 1.78e-04 -4.55e-04 1.19e-21 -2.64e-04 3.45e-05 -3.60e-04 -1.38e-04 2.51e-04 8.34e-06 2.65e-04 2.64e-04 -1.33e-21
    0.00e+00 -1.09e-04 -2.41e-04 -2.12e-04 -1.23e-05 3.16e-04 5.17e-05 3.92e-04 1.09e-04 -4.55e-21 -9.51e-05 -4.03e-04 -9.50e-05 1.88e-04 2.02e-04 2.76e-04 2.41e-04 9.51e-05 -5.48e-21 1.33e-04 -9.02e-05 -1.08e-04 1.96e-04 5.16e-05 2.12e-04 4.03e-04 -1.33e-04 -1.82e-21 -4.90e-04 7.34e-05 8.82e-05 9.53e-05 1.23e-05 9.50e-05 9.02e-05 4.90e-04 -1.97e-21 -1.15e-04 -6.28e-05 -4.72e-04 -3.16e-04 -1.88e-04 1.08e-04 -7.34e-05 1.15e-04 1.81e-21 -1.36e-04 -1.21e-04 -5.17e-05 -2.02e-04 -1.96e-04 -8.82e-05 6.28e-05 1.36e-04 2.05e-21 2.98e-04 -3.92e-04 -2.76e-04 -5.16e-05 -9.53e-05 4.72e-04 1.21e-04 -2.98e-04 -2.37e-21
    0.00e+00 -1.41e-04 2.32e-04 8.57e-05 -5.00e-04 8.22e-04 1.03e-04 -9.04e-05 1.41e-04 -5.82e-21 2.09e-04 -2.32e-04 -4.62e-05 6.34e-05 1.79e-04 -2.32e-04 -2.32e-04 -2.09e-04 -2.86e-23 2.36e-04 -7.38e-05 2.25e-04 1.92e-05 3.85e-05 -8.57e-05 2.32e-04 -2.36e-04 1.23e-23 -1.73e-04 -9.83e-05 9.49e-06 6.19e-04 5.00e-04 4.62e-05 7.38e-05 1.73e-04 -3.15e-22 -1.47e-04 2.26e-06 2.49e-05 -8.22e-04 -6.34e-05 -2.25e-04 9.83e-05 1.47e-04 4.35e-22 -1.25e-04 2.38e-04 -1.03e-04 -1.79e-04 -1.92e-05 -9.49e-06 -2.26e-06 1.25e-04 -1.38e-22 -1.42e-04 9.04e-05 2.32e-04 -3.85e-05 -6.19e-04 -2.49e-05 -2.38e-04 1.42e-04 -8.50e-22
    0.00e+00 -6.23e-05 -1.67e-04 -3.60e-04 -6.11e-05 2.88e-04 7.92e-05 4.02e-04 6.23e-05 -2.41e-22 -4.11e-05 -3.16e-04 -2.66e-06 2.20e-04 1.53e-04 3.14e-04 1.67e-04 4.11e-05 -3.76e-21 1.53e-04 -1.91e-04 9.45e-05 9.32e-05 2.78e-04 3.60e-04 3.16e-04 -1.53e-04 3.84e-22 -2.90e-04 3.67e-05 1.90e-04 1.90e-04 6.11e-05 2.66e-06 1.91e-04 2.90e-04 -6.60e-22 -1.72e-04 -4.27e-05 -4.33e-04 -2.88e-04 -2.20e-04 -9.45e-05 -3.67e-05 1.72e-04 1.56e-21 -1.38e-04 1.58e-05 -7.92e-05 -1.53e-04 -9.32e-05 -1.90e-04 4.27e-05 1.38e-04 3.14e-22 3.30e-04 -4.02e-04 -3.14e-04 -2.78e-04 -1.90e-04 4.33e-04 -1.58e-05 -3.30e-04 2.48e-23
    0.00e+00 6.30e-05 -1.69e-05 -3.18e-04 1.74e-05 1.30e-04 -1.10e-05 3.24e-04 -6.30e-05 -3.17e-22 -2.28e-04 -2.05e-04 6.59e-05 1.14e-04 -2.90e-07 1.88e-04 1.69e-05 2.28e-04 -1.27e-21 -1.92e-04 -8.46e-05 1.95e-04 -1.39e-05 7.49e-04 3.18e-04 2.05e-04 1.92e-04 4.71e-23 -9.88e-05 3.44e-05 1.28e-04 2.98e-04 -1.74e-05 -6.59e-05 8.46e-05 9.88e-05 1.06e-22 -6.59e-05 3.72e-06 -2.27e-04 -1.30e-04 -1.14e-04 -1.95e-04 -3.44e-05 6.59e-05 -5.69e-22 -3.32e-05 -9.44e-05 1.10e-05 2.90e-07 1.39e-05 -1.28e-04 -3.72e-06 3.32e-05 -7.39e-22 1.31e-04 -3.24e-04 -1.88e-04 -7.49e-04 -2.98e-04 2.27e-04 9.44e-05 -1.31e-04 -4.24e-22
    0.00e+00 1.93e-04 3.10e-04 -2.43e-05 -5.28e-04 7.44e-05 -2.02e-04 -5.34e-05 -1.93e-04 -1.05e-21 -8.95e-06 2.74e-04 1.76e-04 -1.76e-04 4.87e-04 2.08e-05 -3.10e-04 8.95e-06 -8.83e-22 4.28e-05 -8.78e-05 9.94e-05 2.42e-04 1.80e-04 2.43e-05 -2.74e-04 -4.28e-05 1.80e-22 -1.33e-04 1.95e-04 1.44e-04 2.10e-04 5.28e-04 -1.76e-04 8.78e-05 1.33e-04 -2.57e-21 -2.81e-04 -1.10e-04 -6.44e-05 -7.44e-05 1.76e-04 -9.94e-05 -1.95e-04 2.81e-04 -3.62e-22 -3.50e-04 -6.39e-05 2.02e-04 -4.87e-04 -2.42e-04 -1.44e-04 1.10e-04 3.50e-04 3.97e-22 1.62e-04 5.34e-05 -2.08e-05 -1.80e-04 -2.10e-04 6.44e-05 6.39e-05 -1.62e-04 -5.37e-22
    0.00e+00 3.34e-04 4.09e-04 8.12e-05 -6.00e-04 -2.24e-04 -2.12e-04 -1.25e-04 -3.34e-04 1.20e-21 9.15e-05 2.93e-04 2.12e-04 -3.50e-04 3.08e-04 -7.15e-05 -4.09e-04 -9.15e-05 1.32e-21 1.57e-04 -2.09e-04 -3.92e-04 2.03e-04 1.10e-04 -8.12e-05 -2.93e-04 -1.57e-04 -6.02e-22 -9.48e-05 -8.92e-05 3.94e-05 2.02e-04 6.00e-04 -2.12e-04 2.09e-04 9.48e-05 1.94e-21 -1.90e-04 -1.15e-04 2.60e-05 2.24e-04 3.50e-04 3.92e-04 8.92e-05 1.90e-04 -2.02e-21 -2.07e-04 -2.71e-04 2.12e-04 -3.08e-04 -2.03e-04 -3.94e-05 1.15e-04 2.07e-04 1.02e-21 1.02e-04 1.25e-04 7.15e-05 -1.10e-04 -2.02e-04 -2.60e-05 2.71e-04 -1.02e-04 3.35e-22
    0.00e+00 5.94e-06 -3.49e-04 -1.26e-04 -1.45e-04 -1.81e-04 -6.70e-05 7.95e-05 -5.94e-06 2.21e-22 1.47e-05 2.73e-04 -3.73e-05 2.96e-04 -1.79e-04 -8.42e-05 3.49e-04 -1.47e-05 -3.09e-21 8.62e-05 -3.93e-04 -2.21e-04 -3.74e-04 -8.05e-05 1.26e-04 -2.73e-04 -8.62e-05 -5.45e-22 2.62e-05 6.89e-05 1.27e-04 -2.15e-04 1.45e-04 3.73e-05 3.93e-04 -2.62e-05 2.64e-22 -2.70e-05 1.67e-04 -5.77e-05 1.81e-04 -2.96e-04 2.21e-04 -6.89e-05 2.70e-05 3.68e-22 2.45e-04 -1.36e-04 6.70e-05 1.79e-04 3.74e-04 -1.27e-04 -1.67e-04 -2.45e-04 -1.96e-22 7.23e-05 -7.95e-05 8.42e-05 8.05e-05 2.15e-04 5.77e-05 1.36e-04 -7.23e-05 -4.07e-22
    0.00e+00 1.51e-04 2.59e-04 -4.39e-05 -5.06e-04 2.67e-04 -8.71e-05 -5.69e-05 -1.51e-04 -2.28e-22 -4.23e-05 2.42e-04 1.95e-04 7.39e-05 5.43e-04 7.33e-05 -2.59e-04 4.23e-05 5.69e-22 3.40e-05 -2.63e-05 2.55e-04 2.27e-04 2.73e-04 4.39e-05 -2.42e-04 -3.40e-05 -9.01e-22 -8.55e-05 3.79e-04 1.20e-04 2.91e-04 5.06e-04 -1.95e-04 2.63e-05 8.55e-05 -1.58e-21 -4.09e-04 -1.23e-04 -9.41e-05 -2.67e-04 -7.39e-05 -2.55e-04 -3.79e-04 4.09e-04 3.61e-21 -3.54e-04 7.72e-05 8.71e-05 -5.43e-04 -2.27e-04 -1.20e-04 1.23e-04 3.54e-04 -9.47e-23 1.68e-04 5.69e-05 -7.33e-05 -2.73e-04 -2.91e-04 9.41e-05 -7.72e-05 -1.68e-04 2.62e-21
    0.00e+00 -1.65e-04 -1.92e-04 1.73e-05 -1.61e-05 1.51e-04 3.80e-05 1.09e-04 1.65e-04 3.18e-23 -1.48e-04 -3.53e-04 -2.62e-04 8.75e-05 1.90e-04 -4.05e-04 1.92e-04 1.48e-04 8.95e-22 -1.93e-04 4.89e-05 -9.99e-05 1.84e-04 -3.56e-04 -1.73e-05 3.53e-04 1.93e-04 -8.89e-22 -1.76e-04 -6.31e-05 5.43e-05 -4.13e-04 1.61e-05 2.62e-04 -4.89e-05 1.76e-04 1.45e-21 7.85e-05 4.89e-05 -5.21e-04 -1.51e-04 -8.75e-05 9.99e-05 6.31e-05 -7.85e-05 1.84e-21 -6.05e-05 2.56e-04 -3.80e-05 -1.90e-04 -1.84e-04 -5.43e-05 -4.89e-05 6.05e-05 6.54e-21 1.26e-04 -1.09e-04 4.05e-04 3.56e-04 4.13e-04 5.21e-04 -2.56e-04 -1.26e-04 8.43e-22
    0.00e+00 -8.55e-06 2.27e-04 9.36e-05 -3.32e-04 2.92e-05 -1.28e-04 -8.93e-05 8.55e-06 -1.50e-21 -7.09e-05 2.32e-05 1.42e-04 1.36e-05 6.45e-04 1.58e-04 -2.27e-04 7.09e-05 7.97e-22 -1.48e-04 1.16e-04 1.81e-04 1.73e-04 2.36e-04 -9.36e-05 -2.32e-05 1.48e-04 1.46e-21 -1.47e-04 1.27e-04 3.02e-04 2.26e-04 3.32e-04 -1.42e-04 -1.16e-04 1.47e-04 7.70e-22 -2.21e-04 -8.83e-05 -1.52e-04 -2.92e-05 -1.36e-05 -1.81e-04 -1.27e-04 2.21e-04 -1.16e-21 -4.14e-04 -6.29e-05 1.28e-04 -6.45e-04 -1.73e-04 -3.02e-04 8.83e-05 4.14e-04 -2.78e-22 5.42e-05 8.93e-05 -1.58e-04 -2.36e-04 -2.26e-04 1.52e-04 6.29e-05 -5.42e-05 6.99e-22
    0.00e+00 2.52e-05 6.54e-04 8.34e-05 -2.35e-04 -1.72e-04 1.54e-04 -9.28e-05 -2.52e-05 -1.24e-21 1.74e-04 -9.87e-05 -2.17e-04 -2.43e-04 1.23e-04 -3.08e-05 -6.54e-04 -1.74e-04 -1.31e-21 1.96e-04 -2.58e-05 -3.53e-04 7.13e-05 -1.86e-04 -8.34e-05 9.87e-05 -1.96e-04 7.50e-22 2.35e-04 3.76e-05 2.50e-05 1.60e-04 2.35e-04 2.17e-04 2.58e-05 -2.35e-04 3.99e-22 1.51e-05 -1.25e-04 3.82e-04 1.72e-04 2.43e-04 3.53e-04 -3.76e-05 -1.51e-05 -2.02e-21 4.01e-05 1.41e-04 -1.54e-04 -1.23e-04 -7.13e-05 -2.50e-05 1.25e-04 -4.01e-05 3.43e-22 6.35e-07 9.28e-05 3.08e-05 1.86e-04 -1.60e-04 -3.82e-04 -1.41e-04 -6.35e-07 -6.34e-22
    0.00e+00 -8.96e-05 -2.86e-04 -1.82e-04 -1.37e-04 1.81e-04 2.00e-04 -1.67e-04 8.96e-05 5.22e-23 -1.27e-04 -3.24e-04 -8.29e-05 1.87e-04 8.16e-05 -2.33e-04 2.86e-04 1.27e-04 4.20e-22 4.14e-04 -2.16e-05 -2.78e-05 1.20e-04 -3.46e-04 1.82e-04 3.24e-04 -4.14e-04 1.06e-21 -2.88e-04 -1.33e-05 -1.00e-04 -3.13e-04 1.37e-04 8.29e-05 2.16e-05 2.88e-04 -9.74e-22 -5.01e-05 -9.39e-05 -2.26e-04 -1.81e-04 -1.87e-04 2.78e-05 1.33e-05 5.01e-05 -2.54e-22 6.64e-05 1.94e-04 -2.00e-04 -8.16e-05 -1.20e-04 1.00e-04 9.39e-05 -6.64e-05 -1.69e-21 4.50e-04 1.67e-04 2.33e-04 3.46e-04 3.13e-04 2.26e-04 -1.94e-04 -4.50e-04 -6.37e-22
    0.00e+00 8.36e-05 -9.00e-05 1.04e-04 -1.53e-04 1.97e-05 1.36e-04 -1.40e-04 -8.36e-05 -5.36e-22 1.01e-04 -1.28e-04 1.33e-04 3.37e-04 6.61e-04 1.03e-04 9.00e-05 -1.01e-04 1.35e-21 -2.39e-04 1.98e-04 -5.50e-05 8.48e-05 3.77e-04 -1.04e-04 1.28e-04 2.39e-04 -2.37e-22 -7.35e-05 3.18e-04 3.68e-04 1.67e-05 1.53e-04 -1.33e-04 -1.98e-04 7.35e-05 2.79e-22 -1.53e-04 -9.29e-05 -6.83e-06 -1.97e-05 -3.37e-04 5.50e-05 -3.18e-04 1.53e-04 -3.96e-21 -2.66e-04 1.60e-04 -1.36e-04 -6.61e-04 -8.48e-05 -3.68e-04 9.29e-05 2.66e-04 -7.72e-22 -9.25e-05 1.40e-04 -1.03e-04 -3.77e-04 -1.67e-05 6.83e-06 -1.60e-04 9.25e-05 3.87e-22
    0.00e+00 2.94e-06 -2.38e-04 -2.20e-04 -8.84e-06 3.39e-04 2.68e-05 3.38e-04 -2.94e-06 1.31e-21 -5.10e-05 -4.71e-04 5.06e-05 2.24e-04 1.82e-04 2.60e-04 2.38e-04 5.10e-05 -2.25e-22 1.16e-04 -4.08e-05 -2.44e-04 2.77e-04 1.20e-04 2.20e-04 4.71e-04 -1.16e-04 1.31e-21 -5.59e-04 -1.14e-05 -5.11e-05 1.67e-04 8.84e-06 -5.06e-05 4.08e-05 5.59e-04 -3.18e-22 -7.72e-05 -8.29e-05 -4.75e-04 -3.39e-04 -2.24e-04 2.44e-04 1.14e-05 7.72e-05 -3.52e-21 -9.02e-05 -3.47e-04 -2.68e-05 -1.82e-04 -2.77e-04 5.11e-05 8.29e-05 9.02e-05 -1.11e-21 2.24e-04 -3.38e-04 -2.60e-04 -1.20e-04 -1.67e-04 4.75e-04 3.47e-04 -2.24e-04 -1.61e-21
    0.00e+00 2.85e-05 3.37e-04 -2.36e-04 -1.17e-04 5.66e-04 3.48e-05 -2.49e-04 -2.85e-05 7.01e-22 1.79e-05 -3.48e-04 1.68e-04 1.61e-04 2.18e-04 -1.42e-04 -3.37e-04 -1.79e-05 1.89e-22 1.96e-05 1.51e-04 -2.77e-04 7.22e-05 1.86e-04 2.36e-04 3.48e-04 -1.96e-05 -6.40e-22 -3.47e-04 -1.30e-04 1.23e-05 6.31e-04 1.17e-04 -1.68e-04 -1.51e-04 3.47e-04 -1.63e-22 -1.32e-05 -1.27e-04 -1.73e-04 -5.66e-04 -1.61e-04 2.77e-04 1.30e-04 1.32e-05 -1.66e-21 -1.59e-05 -9.19e-05 -3.48e-05 -2.18e-04 -7.22e-05 -1.23e-05 1.27e-04 1.59e-05 -6.08e-23 -1.46e-04 2.49e-04 1.42e-04 -1.86e-04 -6.31e-04 1.73e-04 9.19e-05 1.46e-04 2.55e-21
    0.00e+00 -1.74e-04 -2.10e-04 -2.26e-04 6.09e-06 1.38e-04 2.40e-04 -2.53e-04 1.74e-04 1.94e-21 -3.47e-04 -3.34e-04 -1.91e-04 1.63e-04 4.28e-05 -4.38e-04 2.10e-04 3.47e-04 1.39e-22 1.14e-04 1.46e-04 8.58e-05 4.03e-05 -1.28e-04 2.26e-04 3.34e-04 -1.14e-04 1.82e-22 -8.86e-05 2.22e-06 -1.47e-04 -1.44e-04 -6.09e-06 1.91e-04 -1.46e-04 8.86e-05 6.26e-21 9.51e-05 1.04e-05 -3.38e-04 -1.38e-04 -1.63e-04 -8.58e-05 -2.22e-06 -9.51e-05 1.57e-21 3.21e-05 5.56e-05 -2.40e-04 -4.28e-05 -4.03e-05 1.47e-04 -1.04e-05 -3.21e-05 1.30e-21 2.52e-04 2.53e-04 4.38e-04 1.28e-04 1.44e-04 3.38e-04 -5.56e-05 -2.52e-04 4.09e-22
    0.00e+00 -2.49e-04 2.75e-04 4.40e-05 -4.14e-04 7.36e-04 4.95e-05 -2.03e-04 2.49e-04 1.04e-21 1.56e-04 -3.76e-04 -5.68e-05 2.51e-05 1.70e-04 -1.05e-04 -2.75e-04 -1.56e-04 2.51e-21 2.96e-04 -3.43e-05 1.96e-04 1.26e-05 -3.07e-05 -4.40e-05 3.76e-04 -2.96e-04 -2.24e-22 -2.26e-04 -1.72e-04 7.42e-05 7.31e-04 4.14e-04 5.68e-05 3.43e-05 2.26e-04 -3.16e-22 -7.69e-05 -4.42e-06 5.85e-05 -7.36e-04 -2.51e-05 -1.96e-04 1.72e-04 7.69e-05 -1.78e-21 -1.16e-04 1.34e-04 -4.95e-05 -1.70e-04 -1.26e-05 -7.42e-05 4.42e-06 1.16e-04 4.72e-22 -1.79e-04 2.03e-04 1.05e-04 3.07e-05 -7.31e-04 -5.85e-05 -1.34e-04 1.79e-04 2.88e-22
    0.00e+00 -8.61e-05 5.87e-04 5.16e-05 -2.48e-04 -1.72e-04 1.91e-04 -1.65e-04 8.61e-05 -9.42e-22 1.67e-04 -1.98e-04 -2.28e-04 -2.45e-04 7.67e-05 1.82e-05 -5.87e-04 -1.67e-04 2.47e-21 1.72e-04 -3.04e-05 -2.14e-04 7.80e-05 -2.66e-04 -5.16e-05 1.98e-04 -1.72e-04 -4.31e-21 2.59e-04 9.43e-05 7.39e-05 1.86e-04 2.48e-04 2.28e-04 3.04e-05 -2.59e-04 2.85e-21 -4.92e-05 -1.22e-04 3.92e-04 1.72e-04 2.45e-04 2.14e-04 -9.43e-05 4.92e-05 2.42e-21 3.21e-05 1.83e-04 -1.91e-04 -7.67e-05 -7.80e-05 -7.39e-05 1.22e-04 -3.21e-05 -6.22e-22 -2.19e-05 1.65e-04 -1.82e-05 2.66e-04 -1.86e-04 -3.92e-04 -1.83e-04 2.19e-05 9.52e-22
    0.00e+00 1.67e-04 -1.55e-05 3.06e-05 -1.68e-04 2.20e-04 -9.65e-05 -2.99e-04 -1.67e-04 5.31e-22 6.17e-05 3.56e-04 1.54e-04 6.08e-04 6.48e-04 -1.35e-04 1.55e-05 -6.17e-05 2.49e-21 -1.28e-04 1.62e-04 -2.92e-04 1.78e-04 5.44e-04 -3.06e-05 -3.56e-04 1.28e-04 -1.03e-21 -8.23e-05 1.19e-04 -6.20e-05 6.58e-05 1.68e-04 -1.54e-04 -1.62e-04 8.23e-05 -1.09e-21 -4.97e-05 6.63e-05 3.09e-05 -2.20e-04 -6.08e-04 2.92e-04 -1.19e-04 4.97e-05 4.32e-21 -1.52e-04 2.60e-04 9.65e-05 -6.48e-04 -1.78e-04 6.20e-05 -6.63e-05 1.52e-04 -1.26e-21 1.52e-04 2.99e-04 1.35e-04 -5.44e-04 -6.58e-05 -3.09e-05 -2.60e-04 -1.52e-04 7.97e-22
    0.00e+00 -1.05e-04 4.00e-04 -2.16e-05 1.62e-04 8.89e-06 1.49e-04 -3.25e-04 1.05e-04 2.03e-21 1.20e-04 -2.02e-04 -9.59e-05 -2.24e-04 1.14e-04 1.05e-05 -4.00e-04 -1.20e-04 -5.00e-21 1.49e-04 2.12e-04 -9.61e-05 3.12e-05 -2.17e-04 2.16e-05 2.02e-04 -1.49e-04 -2.64e-21 5.85e-05 -3.54e-05 8.51e-05 3.08e-04 -1.62e-04 9.59e-05 -2.12e-04 -5.85e-05 -1.19e-22 -1.94e-04 -1.17e-04 2.15e-04 -8.89e-06 2.24e-04 9.61e-05 3.54e-05 1.94e-04 1.17e-21 5.94e-05 1.18e-05 -1.49e-04 -1.14e-04 -3.12e-05 -8.51e-05 1.17e-04 -5.94e-05 -1.78e-22 -2.32e-05 3.25e-04 -1.05e-05 2.17e-04 -3.08e-04 -2.15e-04 -1.18e-05 2.32e-05 4.80e-22
    0.00e+00 9.95e-05 -2.92e-04 -3.90e-04 -4.10e-05 3.75e-04 4.89e-05 4.43e-04 -9.95e-05 1.53e-21 -1.94e-05 -4.53e-04 1.91e-04 3.40e-04 1.31e-04 2.21e-04 2.92e-04 1.94e-05 1.66e-21 3.42e-05 -1.40e-04 8.33e-05 2.05e-04 5.89e-04 3.90e-04 4.53e-04 -3.42e-05 1.37e-21 -3.28e-04 -8.66e-05 2.50e-05 3.17e-04 4.10e-05 -1.91e-04 1.40e-04 3.28e-04 1.38e-22 -1.61e-04 -3.96e-05 -4.56e-04 -3.75e-04 -3.40e-04 -8.33e-05 8.66e-05 1.61e-04 -1.79e-22 -1.52e-04 -2.97e-04 -4.89e-05 -1.31e-04 -2.05e-04 -2.50e-05 3.96e-05 1.52e-04 -1.72e-21 2.66e-04 -4.43e-04 -2.21e-04 -5.89e-04 -3.17e-04 4.56e-04 2.97e-04 -2.66e-04 -5.82e-22
    0.00e+00 -1.19e-04 4.87e-04 -1.70e-04 -2.63e-04 6.58e-04 1.24e-04 -2.42e-04 1.19e-04 5.92e-21 3.06e-05 -4.37e-04 2.90e-05 8.85e-05 2.04e-04 -1.18e-04 -4.87e-04 -3.06e-05 1.17e-21 1.51e-04 4.07e-05 -8.46e-05 7.87e-06 8.96e-06 1.70e-04 4.37e-04 -1.51e-04 3.52e-22 -3.10e-04 -1.66e-04 1.17e-04 7.37e-04 2.63e-04 -2.90e-05 -4.07e-05 3.10e-04 5.08e-21 -9.01e-05 -5.27e-05 -1.36e-04 -6.58e-04 -8.85e-05 8.46e-05 1.66e-04 9.01e-05 -7.78e-23 -2.26e-05 2.79e-05 -1.24e-04 -2.04e-04 -7.87e-06 -1.17e-04 5.27e-05 2.26e-05 6.10e-22 -1.76e-04 2.42e-04 1.18e-04 -8.96e-06 -7.37e-04 1.36e-04 -2.79e-05 1.76e-04 -1.18e-21
    0.00e+00 1.37e-05 5.21e-04 2.65e-05 -4.81e-05 -2.64e-04 1.39e-04 -1.64e-04 -1.37e-05 9.25e-22 1.35e-04 -5.51e-05 -2.36e-04 -2.33e-04 3.91e-05 -5.08e-05 -5.21e-04 -1.35e-04 1.52e-20 1.87e-04 7.34e-05 -3.74e-04 3.48e-05 -8.93e-05 -2.65e-05 5.51e-05 -1.87e-04 2.44e-21 2.16e-04 1.21e-04 -1.16e-05 1.80e-04 4.81e-05 2.36e-04 -7.34e-05 -2.16e-04 1.21e-21 -1.25e-04 -1.03e-04 3.19e-04 2.64e-04 2.33e-04 3.74e-04 -1.21e-04 1.25e-04 -6.62e-21 7.46e-06 8.31e-05 -1.39e-04 -3.91e-05 -3.48e-05 1.16e-05 1.03e-04 -7.46e-06 -1.28e-21 -1.75e-05 1.64e-04 5.08e-05 8.93e-05 -1.80e-04 -3.19e-04 -8.31e-05 1.75e-05 -3.04e-22
    0.00e+00 7.47e-05 4.43e-04 -6.08e-05 -4.35e-04 7.03e-04 1.17e-04 -1.06e-04 -7.47e-05 -3.26e-21 1.76e-04 3.75e-05 9.25e-05 1.84e-04 2.50e-04 -3.29e-04 -4.43e-04 -1.76e-04 1.23e-21 6.75e-05 1.06e-04 1.41e-04 9.67e-05 1.51e-04 6.08e-05 -3.75e-05 -6.75e-05 -2.39e-22 -2.11e-04 -3.41e-05 -1.49e-04 5.07e-04 4.35e-04 -9.25e-05 -1.06e-04 2.11e-04 1.23e-22 -2.09e-04 -4.34e-05 -6.42e-05 -7.03e-04 -1.84e-04 -1.41e-04 3.41e-05 2.09e-04 1.95e-21 -1.44e-04 3.99e-04 -1.17e-04 -2.50e-04 -9.67e-05 1.49e-04 4.34e-05 1.44e-04 -8.89e-22 -1.00e-04 1.06e-04 3.29e-04 -1.51e-04 -5.07e-04 6.42e-05 -3.99e-04 1.00e-04 -8.71e-22
    0.00e+00 -3.93e-05 -1.94e-04 -3.06e-04 -1.56e-05 2.89e-04 1.28e-04 4.59e-04 3.93e-05 -1.32e-22 -3.22e-05 -4.69e-04 -2.86e-05 2.09e-04 7.53e-05 2.87e-04 1.94e-04 3.22e-05 -4.69e-21 8.71e-05 -1.56e-04 1.39e-04 1.47e-04 4.27e-04 3.06e-04 4.69e-04 -8.71e-05 -1.77e-21 -3.06e-04 4.84e-05 1.04e-04 2.60e-04 1.56e-05 2.86e-05 1.56e-04 3.06e-04 1.41e-21 -1.44e-04 -1.41e-05 -4.36e-04 -2.89e-04 -2.09e-04 -1.39e-04 -4.84e-05 1.44e-04 1.94e-23 -1.38e-04 -1.30e-04 -1.28e-04 -7.53e-05 -1.47e-04 -1.04e-04 1.41e-05 1.38e-04 -5.19e-21 2.58e-04 -4.59e-04 -2.87e-04 -4.27e-04 -2.60e-04 4.36e-04 1.30e-04 -2.58e-04 3.11e-22
    0.00e+00 1.44e-04 -3.15e-04 -2.98e-04 -4.86e-05 3.54e-04 -1.37e-05 2.59e-04 -1.44e-04 6.71e-22 -1.38e-05 -4.34e-04 2.59e-04 2.72e-04 1.80e-04 7.44e-05 3.15e-04 1.38e-05 2.05e-21 8.25e-05 -2.52e-05 -2.33e-04 3.26e-04 4.03e-04 2.98e-04 4.34e-04 -8.25e-05 -1.01e-21 -4.42e-04 -1.69e-04 -1.54e-04 2.59e-04 4.86e-05 -2.59e-04 2.52e-05 4.42e-04 -2.56e-23 -7.92e-05 -2.33e-05 -3.87e-04 -3.54e-04 -2.72e-04 2.33e-04 1.69e-04 7.92e-05 -5.81e-22 -1.44e-04 -5.39e-04 1.37e-05 -1.80e-04 -3.26e-04 1.54e-04 2.33e-05 1.44e-04 -3.34e-22 1.14e-04 -2.59e-04 -7.44e-05 -4.03e-04 -2.59e-04 3.87e-04 5.39e-04 -1.14e-04 1.63e-21
    0.00e+00 2.58e-04 -4.03e-04 -2.07e-04 -4.80e-05 1.77e-04 -2.57e-05 1.31e-04 -2.58e-04 2.31e-21 -5.72e-05 9.10e-05 2.21e-04 4.76e-04 3.32e-04 -3.05e-04 4.03e-04 5.72e-05 -1.23e-21 -3.25e-05 -5.18e-05 -2.09e-04 2.88e-04 7.10e-04 2.07e-04 -9.10e-05 3.25e-05 -1.83e-22 -2.18e-04 -2.76e-04 -2.21e-04 1.44e-04 4.80e-05 -2.21e-04 5.18e-05 2.18e-04 1.31e-21 2.80e-05 8.02e-05 -1.98e-04 -1.77e-04 -4.76e-04 2.09e-04 2.76e-04 -2.80e-05 1.71e-21 -1.04e-04 -3.90e-04 2.57e-05 -3.32e-04 -2.88e-04 2.21e-04 -8.02e-05 1.04e-04 -1.17e-21 7.92e-05 -1.31e-04 3.05e-04 -7.10e-04 -1.44e-04 1.98e-04 3.90e-04 -7.92e-05 -1.05e-21
    0.00e+00 2.01e-04 -3.37e-04 -4.93e-05 -6.49e-05 1.21e-04 -6.26e-05 -8.33e-05 -2.01e-04 -3.08e-22 6.22e-06 3.24e-04 1.65e-04 5.97e-04 5.83e-04 -2.50e-04 3.37e-04 -6.22e-06 -1.24e-21 -8.27e-05 1.07e-05 -2.88e-04 2.68e-04 7.79e-04 4.93e-05 -3.24e-04 8.27e-05 -6.47e-21 -8.97e-05 -8.91e-05 -1.16e-04 5.04e-05 6.49e-05 -1.65e-04 -1.07e-05 8.97e-05 -2.63e-21 3.86e-05 9.38e-05 -5.94e-06 -1.21e-04 -5.97e-04 2.88e-04 8.91e-05 -3.86e-05 -6.78e-23 -4.88e-05 -7.37e-05 6.26e-05 -5.83e-04 -2.68e-04 1.16e-04 -9.38e-05 4.88e-05 1.42e-22 6.58e-05 8.33e-05 2.50e-04 -7.79e-04 -5.04e-05 5.94e-06 7.37e-05 -6.58e-05 1.00e-21
    0.00e+00 3.03e-05 -3.69e-04 -3.12e-05 -1.08e-04 -5.99e-05 1.34e-04 -1.02e-04 -3.03e-05 6.53e-22 -1.79e-04 -2.13e-04 -1.64e-04 1.80e-04 -4.20e-04 3.64e-05 3.69e-04 1.79e-04 2.19e-22 -1.11e-05 -2.68e-04 2.24e-04 5.22e-05 2.06e-04 3.12e-05 2.13e-04 1.11e-05 1.43e-21 -5.80e-05 -2.22e-04 3.22e-04 -3.16e-04 1.08e-04 1.64e-04 2.68e-04 5.80e-05 -7.28e-22 -1.63e-04 3.08e-04 -8.32e-05 5.99e-05 -1.80e-04 -2.24e-04 2.22e-04 1.63e-04 -1.17e-21 3.49e-04 -6.22e-05 -1.34e-04 4.20e-04 -5.22e-05 -3.22e-04 -3.08e-04 -3.49e-04 1.64e-21 -3.22e-04 1.02e-04 -3.64e-05 -2.06e-04 3.16e-04 8.32e-05 6.22e-05 3.22e-04 -1.87e-21
    0.00e+00 2.84e-04 3.36e-04 6.27e-06 -6.22e-04 -1.23e-04 -2.87e-04 -1.24e-04 -2.84e-04 4.27e-21 -3.75e-06 3.38e-04 1.94e-04 -3.47e-04 4.22e-04 -9.72e-05 -3.36e-04 3.75e-06 2.98e-23 1.42e-04 -1.48e-04 -1.43e-04 2.07e-04 1.11e-04 -6.27e-06 -3.38e-04 -1.42e-04 -1.32e-22 -1.31e-04 -1.54e-05 4.19e-05 1.72e-04 6.22e-04 -1.94e-04 1.48e-04 1.31e-04 1.33e-21 -2.14e-04 -5.79e-05 -1.03e-05 1.23e-04 3.47e-04 1.43e-04 1.54e-05 2.14e-04 -2.47e-21 -2.75e-04 -1.75e-04 2.87e-04 -4.22e-04 -2.07e-04 -4.19e-05 5.79e-05 2.75e-04 1.11e-21 9.66e-05 1.24e-04 9.72e-05 -1.11e-04 -1.72e-04 1.03e-05 1.75e-04 -9.66e-05 1.85e-22
    0.00e+00 -1.98e-04 -3.54e-05 -1.38e-04 3.11e-05 3.88e-04 -4.92e-05 1.74e-04 1.98e-04 -1.70e-21 -2.46e-04 -3.17e-04 -1.40e-04 8.88e-05 2.23e-04 3.16e-04 3.54e-05 2.46e-04 1.15e-21 3.22e-05 -5.90e-05 -3.76e-04 1.73e-04 -1.32e-04 1.38e-04 3.17e-04 -3.22e-05 2.29e-22 -4.94e-04 1.59e-04 1.48e-04 1.85e-04 -3.11e-05 1.40e-04 5.90e-05 4.94e-04 -6.60e-22 -4.55e-05 -1.05e-04 -4.63e-04 -3.88e-04 -8.88e-05 3.76e-04 -1.59e-04 4.55e-05 -2.49e-21 -8.22e-05 -1.99e-04 4.92e-05 -2.23e-04 -1.73e-04 -1.48e-04 1.05e-04 8.22e-05 -1.10e-22 1.87e-04 -1.74e-04 -3.16e-04 1.32e-04 -1.85e-04 4.63e-04 1.99e-04 -1.87e-04 9.13e-22
    0.00e+00 1.62e-04 4.74e-04 -7.24e-05 -4.31e-04 6.08e-04 5.29e-05 -9.31e-05 -1.62e-04 1.93e-22 1.44e-04 2.71e-04 1.41e-04 2.49e-04 3.38e-04 -2.98e-04 -4.74e-04 -1.44e-04 -8.35e-22 -1.20e-05 1.50e-04 1.39e-04 1.45e-04 2.14e-04 7.24e-05 -2.71e-04 1.20e-05 2.39e-21 -1.94e-04 1.06e-04 -2.28e-04 3.56e-04 4.31e-04 -1.41e-04 -1.50e-04 1.94e-04 -2.94e-22 -2.43e-04 -3.31e-05 -6.50e-05 -6.08e-04 -2.49e-04 -1.39e-04 -1.06e-04 2.43e-04 -1.58e-22 -1.92e-04 4.80e-04 -5.29e-05 -3.38e-04 -1.45e-04 2.28e-04 3.31e-05 1.92e-04 8.02e-22 2.61e-05 9.31e-05 2.98e-04 -2.14e-04 -3.56e-04 6.50e-05 -4.80e-04 -2.61e-05 1.19e-21
    0.00e+00 3.58e-04 -2.20e-04 1.24e-04 1.06e-04 -4.35e-04 1.43e-05 1.75e-04 -3.58e-04 -1.98e-22 1.80e-04 3.12e-04 -3.03e-04 1.77e-04 -2.77e-04 2.28e-04 2.20e-04 -1.80e-04 5.60e-22 -9.83e-05 -1.84e-04 3.85e-04 -3.53e-04 1.08e-04 -1.24e-04 -3.12e-04 9.83e-05 -1.86e-21 -3.32e-05 -5.71e-05 3.87e-04 -2.83e-04 -1.06e-04 3.03e-04 1.84e-04 3.32e-05 -4.37e-22 -2.78e-05 -2.20e-04 8.78e-05 4.35e-04 -1.77e-04 -3.85e-04 5.71e-05 2.78e-05 -8.20e-22 1.11e-05 1.56e-04 -1.43e-05 2.77e-04 3.53e-04 -3.87e-04 2.20e-04 -1.11e-05 -1.62e-22 -4.68e-04 -1.75e-04 -2.28e-04 -1.08e-04 2.83e-04 -8.78e-05 -1.56e-04 4.68e-04 -3.18e-22
    0.00e+00 1.62e-04 8.27e-05 1.95e-05 -3.05e-04 4.19e-04 7.16e-05 6.32e-05 -1.62e-04 -2.12e-22 8.13e-05 1.85e-04 1.65e-04 2.65e-04 2.21e-04 1.91e-04 -8.27e-05 -8.13e-05 1.06e-21 1.12e-04 -1.80e-04 2.40e-04 1.54e-04 1.93e-04 -1.95e-05 -1.85e-04 -1.12e-04 -9.00e-22 -3.87e-05 4.43e-04 1.28e-04 1.92e-04 3.05e-04 -1.65e-04 1.80e-04 3.87e-05 2.56e-23 -3.70e-04 -7.07e-05 -6.69e-05 -4.19e-04 -2.65e-04 -2.40e-04 -4.43e-04 3.70e-04 4.76e-22 -2.39e-04 -7.50e-07 -7.16e-05 -2.21e-04 -1.54e-04 -1.28e-04 7.07e-05 2.39e-04 -2.22e-21 1.55e-04 -6.32e-05 -1.91e-04 -1.93e-04 -1.92e-04 6.69e-05 7.50e-07 -1.55e-04 -1.34e-21
    0.00e+00 -7.92e-05 -8.47e-05 -7.37e-05 1.99e-04 -1.89e-07 1.15e-04 1.22e-04 7.92e-05 7.51e-22 -5.66e-05 -1.50e-04 -1.31e-04 -3.88e-05 8.03e-05 -3.06e-05 8.47e-05 5.66e-05 -6.12e-22 5.18e-05 1.47e-04 2.20e-05 -6.61e-06 -4.93e-04 7.37e-05 1.50e-04 -5.18e-05 -7.48e-22 1.73e-04 1.96e-04 -1.52e-05 -4.03e-04 -1.99e-04 1.31e-04 -1.47e-04 -1.73e-04 -5.24e-21 5.58e-04 5.95e-05 -1.35e-04 1.89e-07 3.88e-05 -2.20e-05 -1.96e-04 -5.58e-04 5.69e-22 5.05e-05 3.88e-04 -1.15e-04 -8.03e-05 6.61e-06 1.52e-05 -5.95e-05 -5.05e-05 4.18e-22 5.92e-04 -1.22e-04 3.06e-05 4.93e-04 4.03e-04 1.35e-04 -3.88e-04 -5.92e-04 1.14e-22
    0.00e+00 -1.16e-04 6.50e-05 7.21e-05 3.00e-04 1.77e-04 -3.24e-05 1.32e-04 1.16e-04 4.24e-23 -1.63e-04 -4.22e-04 4.42e-05 -8.90e-05 1.74e-04 -3.33e-04 -6.50e-05 1.63e-04 1.56e-22 -3.05e-04 -7.18e-05 -1.35e-04 1.39e-04 -3.61e-04 -7.21e-05 4.22e-04 3.05e-04 2.82e-22 2.86e-04 -1.49e-04 1.92e-04 -7.30e-05 -3.00e-04 -4.42e-05 7.18e-05 -2.86e-04 1.60e-21 1.75e-04 -8.07e-05 -5.34e-04 -1.77e-04 8.90e-05 1.35e-04 1.49e-04 -1.75e-04 1.12e-21 -6.19e-05 8.66e-05 3.24e-05 -1.74e-04 -1.39e-04 -1.92e-04 8.07e-05 6.19e-05 -1.00e-21 -9.38e-05 -1.32e-04 3.33e-04 3.61e-04 7.30e-05 5.34e-04 -8.66e-05 9.38e-05 8.01e-22
    0.00e+00 5.77e-05 -2.12e-04 -5.70e-05 -8.23e-05 -2.91e-04 -5.61e-05 6.50e-05 -5.77e-05 2.95e-21 1.44e-05 2.40e-04 -1.36e-04 2.84e-04 -1.65e-04 1.01e-04 2.12e-04 -1.44e-05 -1.65e-21 6.02e-05 -3.93e-04 -1.88e-05 -2.96e-04 1.68e-04 5.70e-05 -2.40e-04 -6.02e-05 6.35e-22 1.50e-04 1.28e-04 8.87e-05 -2.67e-04 8.23e-05 1.36e-04 3.93e-04 -1.50e-04 4.80e-22 1.98e-05 1.06e-04 4.93e-05 2.91e-04 -2.84e-04 1.88e-05 -1.28e-04 -1.98e-05 9.52e-22 8.55e-05 -1.66e-04 5.61e-05 1.65e-04 2.96e-04 -8.87e-05 -1.06e-04 -8.55e-05 -3.74e-21 6.90e-05 -6.50e-05 -1.01e-04 -1.68e-04 2.67e-04 -4.93e-05 1.66e-04 -6.90e-05 -7.87e-22
    0.00e+00 -3.10e-04 3.87e-04 -1.54e-04 -5.57e-05 4.78e-04 -1.54e-04 -2.19e-04 3.10e-04 -1.50e-21 -2.04e-04 -5.55e-05 -2.86e-05 -6.04e-05 1.99e-04 2.55e-04 -3.87e-04 2.04e-04 1.02e-21 3.72e-04 -3.68e-05 -6.68e-05 6.77e-05 -1.81e-04 1.54e-04 5.55e-05 -3.72e-04 5.95e-22 -3.34e-04 -6.97e-05 1.24e-04 4.40e-04 5.57e-05 2.86e-05 3.68e-05 3.34e-04 3.48e-22 1.16e-04 5.56e-06 -1.49e-04 -4.78e-04 6.04e-05 6.68e-05 6.97e-05 -1.16e-04 -2.55e-21 -4.23e-05 -1.98e-04 1.54e-04 -1.99e-04 -6.77e-05 -1.24e-04 -5.56e-06 4.23e-05 -7.73e-21 9.49e-05 2.19e-04 -2.55e-04 1.81e-04 -4.40e-04 1.49e-04 1.98e-04 -9.49e-05 -1.14e-21
    0.00e+00 1.52e-04 4.58e-04 -4.09e-05 -4.20e-04 5.23e-04 -3.18e-05 -8.98e-05 -1.52e-04 -2.59e-22 5.20e-05 3.96e-04 1.74e-04 3.47e-04 4.75e-04 -5.02e-05 -4.58e-04 -5.20e-05 -1.07e-21 -6.84e-05 1.42e-04 2.25e-04 1.93e-04 2.80e-04 4.09e-05 -3.96e-04 6.84e-05 -1.90e-23 -1.13e-04 3.38e-04 -9.91e-05 2.64e-04 4.20e-04 -1.74e-04 -1.42e-04 1.13e-04 -1.49e-21 -3.63e-04 -7.16e-05 -7.17e-05 -5.23e-04 -3.47e-04 -2.25e-04 -3.38e-04 3.63e-04 3.57e-22 -2.91e-04 4.55e-04 3.18e-05 -4.75e-04 -1.93e-04 9.91e-05 7.16e-05 2.91e-04 -1.96e-22 1.93e-04 8.98e-05 5.02e-05 -2.80e-04 -2.64e-04 7.17e-05 -4.55e-04 -1.93e-04 -5.64e-22
    0.00e+00 -1.77e-04 1.32e-04 2.53e-05 -2.91e-04 -3.69e-05 -1.35e-04 -1.54e-04 1.77e-04 -3.67e-22 4.98e-06 4.13e-04 -4.57e-05 -4.81e-04 2.25e-04 -7.84e-05 -1.32e-04 -4.98e-06 -2.53e-21 -1.22e-04 8.26e-05 -8.82e-06 5.56e-05 1.62e-04 -2.53e-05 -4.13e-04 1.22e-04 -1.48e-21 -2.15e-04 -3.43e-04 1.34e-04 1.14e-04 2.91e-04 4.57e-05 -8.26e-05 2.15e-04 -6.13e-22 -1.08e-05 7.58e-05 -1.00e-04 3.69e-05 4.81e-04 8.82e-06 3.43e-04 1.08e-05 -1.28e-21 -3.58e-04 -1.48e-04 1.35e-04 -2.25e-04 -5.56e-05 -1.34e-04 -7.58e-05 3.58e-04 1.16e-21 -1.81e-04 1.54e-04 7.84e-05 -1.62e-04 -1.14e-04 1.00e-04 1.48e-04 1.81e-04 -6.40e-22
    0.00e+00 -2.65e-05 -4.09e-04 -1.07e-04 -1.27e-04 9.71e-05 1.59e-04 -2.40e-04 2.65e-05 -5.00e-22 -8.23e-05 -2.92e-04 -4.94e-05 1.87e-04 -1.09e-04 -1.57e-04 4.09e-04 8.23e-05 1.13e-21 4.44e-04 -3.11e-05 -5.56e-05 1.78e-05 -2.55e-04 1.07e-04 2.92e-04 -4.44e-04 1.57e-22 -3.04e-04 -8.19e-05 -6.64e-05 -2.35e-04 1.27e-04 4.94e-05 3.11e-05 3.04e-04 7.41e-22 -8.64e-05 -8.11e-06 -6.10e-05 -9.71e-05 -1.87e-04 5.56e-05 8.19e-05 8.64e-05 9.26e-22 1.72e-04 8.31e-05 -1.59e-04 1.09e-04 -1.78e-05 6.64e-05 8.11e-06 -1.72e-04 -9.86e-22 3.22e-04 2.40e-04 1.57e-04 2.55e-04 2.35e-04 6.10e-05 -8.31e-05 -3.22e-04 -2.19e-21
    0.00e+00 1.78e-04 -1.97e-04 -2.48e-04 -8.03e-05 4.10e-04 -7.37e-05 2.51e-05 -1.78e-04 9.33e-22 -4.23e-05 -2.78e-04 2.88e-04 3.03e-04 2.34e-04 -1.20e-04 1.97e-04 4.23e-05 1.40e-21 -8.82e-06 5.18e-05 -3.17e-04 2.57e-04 5.02e-04 2.48e-04 2.78e-04 8.82e-06 -1.48e-21 -3.84e-04 -2.00e-04 -2.22e-04 3.51e-04 8.03e-05 -2.88e-04 -5.18e-05 3.84e-04 6.10e-22 2.42e-06 7.21e-06 -2.68e-04 -4.10e-04 -3.03e-04 3.17e-04 2.00e-04 -2.42e-06 -4.00e-22 -1.19e-04 -4.60e-04 7.37e-05 -2.34e-04 -2.57e-04 2.22e-04 -7.21e-06 1.19e-04 1.25e-22 1.91e-05 -2.51e-05 1.20e-04 -5.02e-04 -3.51e-04 2.68e-04 4.60e-04 -1.91e-05 2.29e-21
    0.00e+00 -5.74e-05 -2.71e-04 -1.49e-04 -8.34e-05 2.40e-04 1.82e-04 2.16e-04 5.74e-05 -4.34e-21 5.27e-09 -4.31e-04 -8.82e-05 1.78e-04 1.83e-04 7.25e-05 2.71e-04 -5.27e-09 1.54e-21 2.24e-04 -8.78e-05 -1.09e-04 2.71e-04 -1.53e-04 1.49e-04 4.31e-04 -2.24e-04 -1.33e-21 -4.75e-04 -4.50e-06 -4.77e-05 -2.08e-04 8.34e-05 8.82e-05 8.78e-05 4.75e-04 1.17e-21 -1.01e-04 -9.25e-05 -3.77e-04 -2.40e-04 -1.78e-04 1.09e-04 4.50e-06 1.01e-04 1.79e-21 -6.26e-05 8.75e-05 -1.82e-04 -1.83e-04 -2.71e-04 4.77e-05 9.25e-05 6.26e-05 1.93e-21 3.53e-04 -2.16e-04 -7.25e-05 1.53e-04 2.08e-04 3.77e-04 -8.75e-05 -3.53e-04 1.81e-21
    0.00e+00 5.32e-05 -8.09e-05 -3.38e-04 -1.20e-04 2.63e-04 1.38e-04 3.73e-04 -5.32e-05 -3.10e-21 4.04e-05 -3.02e-04 7.12e-05 2.60e-04 6.71e-05 2.57e-04 8.09e-05 -4.04e-05 1.52e-22 1.04e-04 -2.05e-04 9.24e-05 2.57e-04 3.77e-04 3.38e-04 3.02e-04 -1.04e-04 -1.77e-21 -2.54e-04 -1.94e-05 1.11e-04 -2.99e-08 1.20e-04 -7.12e-05 2.05e-04 2.54e-04 1.91e-22 -1.38e-04 -4.74e-05 -3.78e-04 -2.63e-04 -2.60e-04 -9.24e-05 1.94e-05 1.38e-04 1.17e-21 -3.49e-05 3.35e-05 -1.38e-04 -6.71e-05 -2.57e-04 -1.11e-04 4.74e-05 3.49e-05 -1.57e-21 2.86e-04 -3.73e-04 -2.57e-04 -3.77e-04 2.99e-08 3.78e-04 -3.35e-05 -2.86e-04 -5.69e-22
    0.00e+00 2.65e-05 3.21e-04 1.43e-05 -3.77e-04 1.54e-04 -1.95e-04 -2.15e-04 -2.65e-05 -1.47e-21 -9.50e-05 2.14e-04 1.82e-04 1.20e-04 7.22e-04 1.07e-04 -3.21e-04 9.50e-05 8.28e-22 -1.28e-04 2.00e-04 1.71e-04 2.07e-04 2.72e-04 -1.43e-05 -2.14e-04 1.28e-04 -6.44e-22 -1.43e-04 3.19e-04 2.18e-04 1.43e-04 3.77e-04 -1.82e-04 -2.00e-04 1.43e-04 -6.08e-22 -2.55e-04 -4.91e-05 -7.27e-05 -1.54e-04 -1.20e-04 -1.71e-04 -3.19e-04 2.55e-04 3.03e-22 -4.22e-04 2.36e-04 1.95e-04 -7.22e-04 -2.07e-04 -2.18e-04 4.91e-05 4.22e-04 -1.58e-21 1.50e-04 2.15e-04 -1.07e-04 -2.72e-04 -1.43e-04 7.27e-05 -2.36e-04 -1.50e-04 2.87e-21
    0.00e+00 -1.68e-04 4.39e-04 -8.15e-06 2.16e-04 2.76e-04 8.43e-05 -4.05e-04 1.68e-04 -6.18e-22 2.41e-05 -3.33e-04 3.45e-05 -2.51e-04 1.52e-04 6.36e-05 -4.39e-04 -2.41e-05 2.12e-21 1.53e-04 2.48e-04 2.79e-05 -5.92e-06 -3.04e-04 8.15e-06 3.33e-04 -1.53e-04 -1.07e-21 4.82e-05 -2.94e-04 1.56e-04 4.09e-04 -2.16e-04 -3.45e-05 -2.48e-04 -4.82e-05 2.04e-21 -3.68e-05 -1.58e-04 6.17e-05 -2.76e-04 2.51e-04 -2.79e-05 2.94e-04 3.68e-05 1.13e-21 2.65e-05 -2.26e-05 -8.43e-05 -1.52e-04 5.92e-06 -1.56e-04 1.58e-04 -2.65e-05 -4.14e-22 -2.37e-05 4.05e-04 -6.36e-05 3.04e-04 -4.09e-04 -6.17e-05 2.26e-05 2.37e-05 2.52e-21
    0.00e+00 -2.36e-05 3.10e-04 7.86e-05 -4.57e-04 4.03e-05 -3.53e-04 -1.39e-04 2.36e-05 2.96e-22 -1.68e-05 2.99e-04 9.93e-05 -4.10e-04 5.07e-04 6.30e-07 -3.10e-04 1.68e-05 -6.17e-22 -6.82e-05 3.53e-05 -5.45e-05 1.39e-04 1.82e-04 -7.86e-05 -2.99e-04 6.82e-05 5.58e-22 -1.73e-04 -2.31e-04 2.40e-04 1.78e-04 4.57e-04 -9.93e-05 -3.53e-05 1.73e-04 -4.54e-22 -8.44e-05 -5.41e-05 -1.11e-04 -4.03e-05 4.10e-04 5.45e-05 2.31e-04 8.44e-05 -7.73e-22 -4.19e-04 -2.31e-04 3.53e-04 -5.07e-04 -1.39e-04 -2.40e-04 5.41e-05 4.19e-04 -1.57e-21 1.16e-05 1.39e-04 -6.30e-07 -1.82e-04 -1.78e-04 1.11e-04 2.31e-04 -1.16e-05 -5.56e-22
    0.00e+00 9.61e-05 -4.71e-04 -8.67e-05 -1.64e-04 -1.09e-04 5.03e-05 -1.36e-04 -9.61e-05 -1.52e-21 -1.33e-04 -1.49e-04 -1.28e-04 2.25e-04 -4.47e-04 9.52e-06 4.71e-04 1.33e-04 -5.57e-22 1.70e-04 -1.68e-04 4.94e-05 9.54e-05 1.25e-04 8.67e-05 1.49e-04 -1.70e-04 -1.66e-21 -9.84e-05 -1.63e-04 3.28e-04 -2.70e-04 1.64e-04 1.28e-04 1.68e-04 9.84e-05 -1.17e-22 -1.43e-04 3.28e-04 -2.08e-05 1.09e-04 -2.25e-04 -4.94e-05 1.63e-04 1.43e-04 -4.18e-22 3.43e-04 -1.64e-04 -5.03e-05 4.47e-04 -9.54e-05 -3.28e-04 -3.28e-04 -3.43e-04 -6.56e-22 -1.19e-04 1.36e-04 -9.52e-06 -1.25e-04 2.70e-04 2.08e-05 1.64e-04 1.19e-04 3.31e-21
    0.00e+00 1.56e-04 5.68e-04 2.11e-06 -7.48e-05 -2.93e-04 8.26e-05 -1.14e-04 -1.56e-04 1.60e-22 1.70e-04 -7.28e-06 -1.58e-04 -2.30e-04 1.33e-04 -1.35e-04 -5.68e-04 -1.70e-04 7.06e-21 2.03e-04 3.67e-05 -4.62e-04 5.69e-05 -5.76e-05 -2.11e-06 7.28e-06 -2.03e-04 5.50e-22 1.37e-04 7.46e-05 -2.87e-05 2.15e-04 7.48e-05 1.58e-04 -3.67e-05 -1.37e-04 -7.13e-22 -1.79e-04 -7.99e-05 2.45e-04 2.93e-04 2.30e-04 4.62e-04 -7.46e-05 1.79e-04 -2.23e-21 3.18e-05 1.66e-05 -8.26e-05 -1.33e-04 -5.69e-05 2.87e-05 7.99e-05 -3.18e-05 -1.07e-21 -9.70e-07 1.14e-04 1.35e-04 5.76e-05 -2.15e-04 -2.45e-04 -1.66e-05 9.70e-07 8.08e-22
    0.00e+00 1.64e-04 3.80e-04 -9.42e-05 -4.85e-04 2.97e-04 -1.38e-04 -1.15e-04 -1.64e-04 2.93e-21 -1.11e-06 3.21e-04 2.10e-04 1.26e-04 5.76e-04 -8.06e-05 -3.80e-04 1.11e-06 -1.63e-21 -9.78e-05 1.31e-04 2.21e-04 2.60e-04 2.40e-04 9.42e-05 -3.21e-04 9.78e-05 -6.51e-22 -1.90e-04 3.91e-04 -5.78e-05 1.65e-04 4.85e-04 -2.10e-04 -1.31e-04 1.90e-04 -8.35e-22 -3.17e-04 -7.12e-05 -5.86e-05 -2.97e-04 -1.26e-04 -2.21e-04 -3.91e-04 3.17e-04 -1.85e-21 -3.02e-04 3.78e-04 1.38e-04 -5.76e-04 -2.60e-04 5.78e-05 7.12e-05 3.02e-04 1.31e-21 1.86e-04 1.15e-04 8.06e-05 -2.40e-04 -1.65e-04 5.86e-05 -3.78e-04 -1.86e-04 7.31e-22
    0.00e+00 2.03e-04 -4.25e-04 -3.23e-04 -3.21e-05 2.77e-04 1.46e-05 2.04e-04 -2.03e-04 -4.17e-22 -3.65e-05 -2.54e-04 2.74e-04 4.03e-04 2.09e-04 -1.60e-04 4.25e-04 3.65e-05 9.38e-22 1.15e-05 -5.50e-05 -1.79e-04 3.34e-04 6.08e-04 3.23e-04 2.54e-04 -1.15e-05 -6.76e-21 -3.12e-04 -2.91e-04 -1.86e-04 2.83e-04 3.21e-05 -2.74e-04 5.50e-05 3.12e-04 -1.85e-21 -2.49e-05 3.49e-05 -3.09e-04 -2.77e-04 -4.03e-04 1.79e-04 2.91e-04 2.49e-05 -4.63e-22 -1.04e-04 -5.51e-04 -1.46e-05 -2.09e-04 -3.34e-04 1.86e-04 -3.49e-05 1.04e-04 -2.67e-21 8.12e-05 -2.04e-04 1.60e-04 -6.08e-04 -2.83e-04 3.09e-04 5.51e-04 -8.12e-05 6.22e-22
    0.00e+00 2.22e-04 3.48e-04 -1.29e-04 -5.39e-04 3.64e-04 -2.65e-05 8.53e-07 -2.22e-04 1.00e-21 7.85e-05 2.71e-04 1.04e-04 -2.75e-05 3.81e-04 -2.34e-04 -3.48e-04 -7.85e-05 -1.88e-22 2.15e-05 1.48e-05 2.72e-04 2.30e-04 1.37e-04 1.29e-04 -2.71e-04 -2.15e-05 -1.08e-21 -1.47e-04 3.38e-04 -1.13e-04 2.60e-04 5.39e-04 -1.04e-04 -1.48e-05 1.47e-04 2.92e-22 -3.62e-04 -7.67e-05 -3.71e-05 -3.64e-04 2.75e-05 -2.72e-04 -3.38e-04 3.62e-04 4.40e-21 -2.39e-04 3.00e-04 2.65e-05 -3.81e-04 -2.30e-04 1.13e-04 7.67e-05 2.39e-04 -1.81e-21 8.45e-05 -8.53e-07 2.34e-04 -1.37e-04 -2.60e-04 3.71e-05 -3.00e-04 -8.45e-05 2.56e-21
    0.00e+00 1.44e-04 -4.32e-04 -3.49e-04 9.82e-07 2.97e-04 6.47e-05 2.90e-04 -1.44e-04 -1.10e-21 -6.35e-06 -4.26e-04 2.51e-04 3.63e-04 1.13e-04 5.77e-05 4.32e-04 6.35e-06 -8.73e-22 7.62e-05 -9.64e-05 -1.04e-04 3.41e-04 4.98e-04 3.49e-04 4.26e-04 -7.62e-05 1.48e-21 -3.65e-04 -2.70e-04 -8.76e-05 2.79e-04 -9.82e-07 -2.51e-04 9.64e-05 3.65e-04 -2.64e-21 -1.11e-04 -3.57e-07 -3.92e-04 -2.97e-04 -3.63e-04 1.04e-04 2.70e-04 1.11e-04 -2.98e-23 -7.71e-05 -5.43e-04 -6.47e-05 -1.13e-04 -3.41e-04 8.76e-05 3.57e-07 7.71e-05 -2.25e-22 1.52e-04 -2.90e-04 -5.77e-05 -4.98e-04 -2.79e-04 3.92e-04 5.43e-04 -1.52e-04 2.10e-21
    0.00e+00 1.44e-04 1.97e-04 -1.76e-04 -2.00e-04 5.20e-04 -2.48e-05 -2.33e-04 -1.44e-04 -2.67e-21 9.93e-06 -1.53e-04 2.25e-04 2.44e-04 2.74e-04 -2.58e-04 -1.97e-04 -9.93e-06 1.51e-22 -3.68e-05 1.47e-04 -3.16e-04 1.48e-04 3.67e-04 1.76e-04 1.53e-04 3.68e-05 -3.90e-22 -3.45e-04 -1.51e-04 -1.67e-04 4.44e-04 2.00e-04 -2.25e-04 -1.47e-04 3.45e-04 6.83e-22 -9.69e-06 -6.96e-06 -1.27e-04 -5.20e-04 -2.44e-04 3.16e-04 1.51e-04 9.69e-06 1.63e-21 -5.53e-05 -1.03e-04 2.48e-05 -2.74e-04 -1.48e-04 1.67e-04 6.96e-06 5.53e-05 6.17e-23 -8.42e-05 2.33e-04 2.58e-04 -3.67e-04 -4.44e-04 1.27e-04 1.03e-04 8.42e-05 -1.31e-21
    0.00e+00 1.78e-04 4.63e-04 -1.05e-04 -4.43e-04 4.35e-04 -8.34e-05 -1.19e-04 -1.78e-04 -5.92e-22 1.02e-04 3.04e-04 1.77e-04 2.44e-04 5.31e-04 -2.33e-04 -4.63e-04 -1.02e-04 -8.35e-22 -6.00e-05 1.65e-04 1.29e-04 2.05e-04 2.25e-04 1.05e-04 -3.04e-04 6.00e-05 9.75e-22 -2.09e-04 2.26e-04 -2.22e-04 2.76e-04 4.43e-04 -1.77e-04 -1.65e-04 2.09e-04 7.27e-22 -2.61e-04 -3.48e-05 -2.40e-06 -4.35e-04 -2.44e-04 -1.29e-04 -2.26e-04 2.61e-04 -5.08e-21 -2.16e-04 4.90e-04 8.34e-05 -5.31e-04 -2.05e-04 2.22e-04 3.48e-05 2.16e-04 1.56e-21 7.28e-05 1.19e-04 2.33e-04 -2.25e-04 -2.76e-04 2.40e-06 -4.90e-04 -7.28e-05 -2.24e-21
    0.00e+00 1.27e-04 3.39e-05 -6.08e-05 -4.52e-04 2.48e-04 9.21e-05 1.08e-05 -1.27e-04 1.51e-21 1.08e-04 2.27e-04 4.96e-05 -7.39e-05 1.84e-04 -2.02e-05 -3.39e-05 -1.08e-04 -8.50e-22 2.35e-04 -2.03e-04 1.09e-04 1.73e-04 1.10e-04 6.08e-05 -2.27e-04 -2.35e-04 -6.85e-22 -3.31e-05 3.70e-04 -1.51e-05 2.90e-04 4.52e-04 -4.96e-05 2.03e-04 3.31e-05 2.76e-22 -3.97e-04 -1.02e-04 -1.04e-05 -2.48e-04 7.39e-05 -1.09e-04 -3.70e-04 3.97e-04 9.57e-23 -2.02e-04 -3.51e-05 -9.21e-05 -1.84e-04 -1.73e-04 1.51e-05 1.02e-04 2.02e-04 -1.85e-21 9.06e-05 -1.08e-05 2.02e-05 -1.10e-04 -2.90e-04 1.04e-05 3.51e-05 -9.06e-05 -1.31e-22
    0.00e+00 2.08e-04 -6.46e-05 -2.47e-04 -1.04e-04 4.73e-05 1.33e-04 4.11e-04 -2.08e-04 1.92e-21 1.99e-05 8.17e-05 1.95e-04 3.22e-04 1.19e-04 -5.99e-05 6.46e-05 -1.99e-05 3.27e-22 -1.48e-04 -9.02e-05 -9.33e-05 3.76e-04 4.96e-04 2.47e-04 -8.17e-05 1.48e-04 -9.91e-22 -1.92e-04 -1.18e-04 5.10e-05 -2.85e-05 1.04e-04 -1.95e-04 9.02e-05 1.92e-04 -8.33e-22 -4.99e-05 -2.42e-05 -3.11e-04 -4.73e-05 -3.22e-04 9.33e-05 1.18e-04 4.99e-05 2.26e-21 2.03e-04 -1.19e-04 -1.33e-04 -1.19e-04 -3.76e-04 -5.10e-05 2.42e-05 -2.03e-04 -2.15e-21 8.21e-06 -4.11e-04 5.99e-05 -4.96e-04 2.85e-05 3.11e-04 1.19e-04 -8.21e-06 2.33e-21
    0.00e+00 -1.03e-04 5.42e-04 -2.38e-05 7.64e-05 3.72e-05 2.43e-04 -2.87e-04 1.03e-04 -2.05e-22 1.34e-04 -2.93e-04 -9.74e-05 -2.78e-04 1.24e-04 4.00e-05 -5.42e-04 -1.34e-04 2.77e-21 1.11e-04 1.88e-04 -6.94e-05 5.00e-05 -2.90e-04 2.38e-05 2.93e-04 -1.11e-04 -2.74e-21 1.63e-04 -3.03e-05 1.52e-04 2.86e-04 -7.64e-05 9.74e-05 -1.88e-04 -1.63e-04 7.45e-22 -1.59e-04 -1.71e-04 2.01e-04 -3.72e-05 2.78e-04 6.94e-05 3.03e-05 1.59e-04 -1.47e-21 7.52e-05 6.12e-05 -2.43e-04 -1.24e-04 -5.00e-05 -1.52e-04 1.71e-04 -7.52e-05 1.70e-21 -4.42e-05 2.87e-04 -4.00e-05 2.90e-04 -2.86e-04 -2.01e-04 -6.12e-05 4.42e-05 2.25e-21
    0.00e+00 1.61e-04 2.74e-04 2.58e-05 -2.79e-04 3.31e-04 -1.10e-04 -2.90e-04 -1.61e-04 -1.78e-23 2.63e-05 4.09e-04 1.67e-04 4.36e-04 5.85e-04 -1.39e-04 -2.74e-04 -2.63e-05 2.37e-22 -9.65e-05 2.42e-04 -1.59e-04 2.02e-04 4.09e-04 -2.58e-05 -4.09e-04 9.65e-05 5.01e-21 -1.18e-04 1.84e-04 -9.63e-05 8.42e-05 2.79e-04 -1.67e-04 -2.42e-04 1.18e-04 2.08e-22 -1.39e-04 3.69e-05 5.62e-06 -3.31e-04 -4.36e-04 1.59e-04 -1.84e-04 1.39e-04 3.43e-21 -2.68e-04 4.25e-04 1.10e-04 -5.85e-04 -2.02e-04 9.63e-05 -3.69e-05 2.68e-04 3.20e-21 1.79e-04 2.90e-04 1.39e-04 -4.09e-04 -8.42e-05 -5.62e-06 -4.25e-04 -1.79e-04 -6.21e-22
    0.00e+00 6.52e-05 -1.47e-04 2.07e-05 -1.30e-05 -1.34e-04 3.24e-04 3.79e-05 -6.52e-05 -2.72e-22 -8.18e-05 4.53e-04 -2.74e-04 -9.03e-05 -3.21e-04 -2.52e-04 1.47e-04 8.18e-05 2.27e-21 -2.12e-05 6.37e-05 -6.66e-05 -1.33e-04 2.63e-04 -2.07e-05 -4.53e-04 2.12e-05 -2.77e-22 -1.46e-04 -4.62e-05 -7.31e-05 -1.52e-04 1.30e-05 2.74e-04 -6.37e-05 1.46e-04 -1.13e-21 -4.89e-06 8.44e-05 -5.65e-07 1.34e-04 9.03e-05 6.66e-05 4.62e-05 4.89e-06 6.89e-22 -3.34e-04 9.38e-05 -3.24e-04 3.21e-04 1.33e-04 7.31e-05 -8.44e-05 3.34e-04 -1.99e-21 -2.84e-04 -3.79e-05 2.52e-04 -2.63e-04 1.52e-04 5.65e-07 -9.38e-05 2.84e-04 1.06e-21
    0.00e+00 1.36e-05 -2.43e-04 -2.57e-04 -6.90e-05 2.39e-04 2.06e-04 3.43e-04 -1.36e-05 -1.78e-21 7.59e-05 -4.49e-04 1.78e-05 2.59e-04 9.61e-05 2.11e-04 2.43e-04 -7.59e-05 8.80e-22 2.36e-04 -1.39e-04 3.44e-05 2.89e-04 1.28e-04 2.57e-04 4.49e-04 -2.36e-04 1.06e-21 -3.99e-04 -6.76e-05 -4.18e-05 -2.10e-05 6.90e-05 -1.78e-05 1.39e-04 3.99e-04 -6.69e-22 -1.40e-04 -5.86e-05 -3.91e-04 -2.39e-04 -2.59e-04 -3.44e-05 6.76e-05 1.40e-04 -1.84e-21 -5.78e-05 -9.03e-05 -2.06e-04 -9.61e-05 -2.89e-04 4.18e-05 5.86e-05 5.78e-05 1.15e-21 2.97e-04 -3.43e-04 -2.11e-04 -1.28e-04 2.10e-05 3.91e-04 9.03e-05 -2.97e-04 -2.35e-21
    0.00e+00 -9.07e-05 4.49e-04 -6.88e-05 -3.23e-04 7.19e-04 1.61e-04 -1.49e-04 9.07e-05 -1.63e-22 1.25e-04 -3.23e-04 3.48e-05 1.17e-04 1.86e-04 -2.18e-04 -4.49e-04 -1.25e-04 4.65e-21 1.61e-04 3.84e-05 6.85e-05 -1.92e-07 6.09e-05 6.88e-05 3.23e-04 -1.61e-04 -1.11e-21 -2.38e-04 -1.55e-04 3.67e-05 6.75e-04 3.23e-04 -3.48e-05 -3.84e-05 2.38e-04 1.13e-22 -1.21e-04 -7.98e-05 -4.50e-05 -7.19e-04 -1.17e-04 -6.85e-05 1.55e-04 1.21e-04 -1.64e-21 -5.09e-05 1.89e-04 -1.61e-04 -1.86e-04 1.92e-07 -3.67e-05 7.98e-05 5.09e-05 -2.22e-22 -1.83e-04 1.49e-04 2.18e-04 -6.09e-05 -6.75e-04 4.50e-05 -1.89e-04 1.83e-04 -2.27e-21
    0.00e+00 -6.92e-05 2.40e-04 8.23e-05 -3.61e-04 1.58e-06 -2.31e-04 -1.48e-04 6.92e-05 1.35e-22 -1.38e-04 1.72e-04 1.09e-04 -2.48e-04 5.45e-04 7.30e-05 -2.40e-04 1.38e-04 3.47e-21 -1.11e-04 9.32e-05 9.75e-05 1.46e-04 2.21e-04 -8.23e-05 -1.72e-04 1.11e-04 -6.83e-21 -1.59e-04 -1.42e-05 2.85e-04 1.48e-04 3.61e-04 -1.09e-04 -9.32e-05 1.59e-04 -1.28e-22 -1.61e-04 -4.54e-05 -1.08e-04 -1.58e-06 2.48e-04 -9.75e-05 1.42e-05 1.61e-04 4.04e-22 -3.94e-04 -1.12e-04 2.31e-04 -5.45e-04 -1.46e-04 -2.85e-04 4.54e-05 3.94e-04 7.42e-21 1.25e-05 1.48e-04 -7.30e-05 -2.21e-04 -1.48e-04 1.08e-04 1.12e-04 -1.25e-05 -1.69e-21
    0.00e+00 -8.39e-05 -2.61e-05 4.00e-05 -4.58e-04 5.16e-04 8.12e-05 3.09e-05 8.39e-05 -4.35e-21 2.11e-04 -1.04e-05 -9.86e-05 -5.06e-05 1.65e-04 -2.08e-04 2.61e-05 -2.11e-04 -8.18e-23 3.27e-04 -1.58e-04 2.65e-04 8.30e-05 7.70e-06 -4.00e-05 1.04e-05 -3.27e-04 7.14e-22 -4.52e-05 2.01e-04 -7.15e-05 4.57e-04 4.58e-04 9.86e-05 1.58e-04 4.52e-05 -1.23e-21 -3.45e-04 7.18e-06 8.88e-05 -5.16e-04 5.06e-05 -2.65e-04 -2.01e-04 3.45e-04 1.62e-22 -1.51e-04 1.59e-04 -8.12e-05 -1.65e-04 -8.30e-05 7.15e-05 -7.18e-06 1.51e-04 1.62e-21 -6.70e-05 -3.09e-05 2.08e-04 -7.70e-06 -4.57e-04 -8.88e-05 -1.59e-04 6.70e-05 1.24e-21
    0.00e+00 3.91e-04 1.93e-04 5.30e-05 -1.53e-04 -4.40e-04 2.05e-04 -1.19e-04 -3.91e-04 -1.67e-21 2.32e-04 5.26e-04 -1.12e-04 1.97e-04 1.48e-04 2.50e-04 -1.93e-04 -2.32e-04 5.71e-22 7.09e-05 -2.22e-04 1.51e-04 -1.35e-04 -7.72e-06 -5.30e-05 -5.26e-04 -7.09e-05 -1.35e-21 -1.68e-04 4.47e-06 6.04e-05 2.08e-05 1.53e-04 1.12e-04 2.22e-04 1.68e-04 -5.64e-22 -2.72e-04 4.77e-04 1.03e-04 4.40e-04 -1.97e-04 -1.51e-04 -4.47e-06 2.72e-04 -3.64e-21 1.70e-04 3.10e-05 -2.05e-04 -1.48e-04 1.35e-04 -6.04e-05 -4.77e-04 -1.70e-04 -1.73e-21 -1.97e-04 1.19e-04 -2.50e-04 7.72e-06 -2.08e-05 -1.03e-04 -3.10e-05 1.97e-04 8.54e-23
    0.00e+00 4.56e-05 -4.81e-04 -9.23e-05 -1.60e-04 4.22e-05 7.24e-05 -2.34e-04 -4.56e-05 -2.90e-22 -5.66e-05 -5.81e-05 -1.42e-05 1.89e-04 -2.29e-04 -3.06e-05 4.81e-04 5.66e-05 -1.65e-21 3.28e-04 -1.79e-04 -8.38e-05 -1.44e-04 -3.61e-04 9.23e-05 5.81e-05 -3.28e-04 -1.68e-21 -1.72e-04 -1.17e-04 1.87e-04 -2.22e-04 1.60e-04 1.42e-05 1.79e-04 1.72e-04 2.69e-21 -9.66e-05 2.05e-05 -6.01e-05 -4.22e-05 -1.89e-04 8.38e-05 1.17e-04 9.66e-05 -5.54e-22 2.44e-04 6.39e-05 -7.24e-05 2.29e-04 1.44e-04 -1.87e-04 -2.05e-05 -2.44e-04 -1.64e-22 2.13e-04 2.34e-04 3.06e-05 3.61e-04 2.22e-04 6.01e-05 -6.39e-05 -2.13e-04 1.40e-21
    0.00e+00 -1.64e-04 -2.15e-04 -3.02e-05 -4.67e-05 1.68e-04 7.23e-05 1.35e-04 1.64e-04 8.17e-22 -1.37e-04 -3.55e-04 -2.81e-04 1.24e-04 2.06e-04 -2.54e-04 2.15e-04 1.37e-04 -7.86e-22 -8.41e-05 3.08e-05 -1.40e-04 1.97e-04 -2.98e-04 3.02e-05 3.55e-04 8.41e-05 -1.23e-21 -3.18e-04 3.66e-07 1.57e-05 -3.78e-04 4.67e-05 2.81e-04 -3.08e-05 3.18e-04 7.95e-22 2.79e-05 -9.94e-06 -4.41e-04 -1.68e-04 -1.24e-04 1.40e-04 -3.66e-07 -2.79e-05 1.73e-21 -8.21e-05 2.12e-04 -7.23e-05 -2.06e-04 -1.97e-04 -1.57e-05 9.94e-06 8.21e-05 3.42e-22 2.10e-04 -1.35e-04 2.54e-04 2.98e-04 3.78e-04 4.41e-04 -2.12e-04 -2.10e-04 -8.51e-22
    0.00e+00 2.18e-04 -3.70e-04 -1.19e-04 -1.92e-04 -3.85e-04 9.71e-05 -2.18e-05 -2.18e-04 -2.52e-21 -1.93e-05 2.69e-04 -1.27e-04 2.17e-04 -3.81e-04 3.18e-04 3.70e-04 1.93e-05 -6.09e-21 -1.01e-04 -4.53e-04 3.99e-05 -9.66e-05 1.27e-04 1.19e-04 -2.69e-04 1.01e-04 2.65e-22 2.26e-04 1.36e-04 2.97e-04 -2.95e-04 1.92e-04 1.27e-04 4.53e-04 -2.26e-04 -1.68e-22 -1.64e-05 1.28e-05 4.98e-06 3.85e-04 -2.17e-04 -3.99e-05 -1.36e-04 1.64e-05 -5.32e-22 3.82e-04 -2.67e-04 -9.71e-05 3.81e-04 9.66e-05 -2.97e-04 -1.28e-05 -3.82e-04 -1.78e-22 -9.54e-05 2.18e-05 -3.18e-04 -1.27e-04 2.95e-04 -4.98e-06 2.67e-04 9.54e-05 1.65e-21
    0.00e+00 -1.26e-04 -1.78e-04 -6.82e-05 2.05e-04 2.06e-05 1.66e-04 1.66e-04 1.26e-04 3.26e-21 -1.26e-04 -3.20e-04 -2.08e-04 1.50e-05 9.17e-05 2.54e-06 1.78e-04 1.26e-04 -8.56e-22 -6.72e-05 -9.51e-05 -5.26e-05 1.11e-04 -1.83e-04 6.82e-05 3.20e-04 6.72e-05 3.94e-21 2.11e-05 2.35e-04 2.33e-05 2.31e-04 -2.05e-04 2.08e-04 9.51e-05 -2.11e-05 -1.34e-21 1.57e-04 -1.43e-04 1.40e-04 -2.06e-05 -1.50e-05 5.26e-05 -2.35e-04 -1.57e-04 -2.32e-23 2.46e-06 8.37e-06 -1.66e-04 -9.17e-05 -1.11e-04 -2.33e-05 1.43e-04 -2.46e-06 1.03e-21 9.00e-05 -1.66e-04 -2.54e-06 1.83e-04 -2.31e-04 -1.40e-04 -8.37e-06 -9.00e-05 1.65e-21
    0.00e+00 -7.53e-05 2.46e-04 -1.97e-04 2.60e-04 1.82e-04 2.69e-04 3.57e-04 7.53e-05 -1.91e-21 -6.66e-05 -3.43e-04 -1.22e-04 -1.49e-04 1.64e-04 1.22e-04 -2.46e-04 6.66e-05 -2.07e-21 7.68e-05 4.89e-05 4.37e-05 -6.67e-05 -2.20e-04 1.97e-04 3.43e-04 -7.68e-05 -3.99e-22 2.75e-04 2.47e-04 5.35e-04 3.64e-04 -2.60e-04 1.22e-04 -4.89e-05 -2.75e-04 2.42e-21 1.98e-04 1.90e-04 6.76e-05 -1.82e-04 1.49e-04 -4.37e-05 -2.47e-04 -1.98e-04 -1.82e-21 8.51e-05 4.45e-05 -2.69e-04 -1.64e-04 6.67e-05 -5.35e-04 -1.90e-04 -8.51e-05 -9.04e-22 1.58e-04 -3.57e-04 -1.22e-04 2.20e-04 -3.64e-04 -6.76e-05 -4.45e-05 -1.58e-04 6.23e-21
    0.00e+00 3.09e-04 -2.13e-04 1.05e-04 1.22e-04 -4.45e-04 4.26e-05 2.14e-04 -3.09e-04 3.70e-22 1.63e-04 3.91e-04 -3.19e-04 2.81e-04 -7.11e-05 1.16e-04 2.13e-04 -1.63e-04 -3.33e-21 -7.68e-05 -1.86e-04 3.23e-04 -3.92e-04 1.41e-04 -1.05e-04 -3.91e-04 7.68e-05 -6.96e-23 -2.09e-05 -5.04e-05 3.80e-04 -3.39e-04 -1.22e-04 3.19e-04 1.86e-04 2.09e-05 3.71e-22 2.76e-05 -7.43e-05 7.00e-05 4.45e-04 -2.81e-04 -3.23e-04 5.04e-05 -2.76e-05 -2.85e-22 -1.81e-04 1.87e-04 -4.26e-05 7.11e-05 3.92e-04 -3.80e-04 7.43e-05 1.81e-04 9.47e-22 -4.05e-04 -2.14e-04 -1.16e-04 -1.41e-04 3.39e-04 -7.00e-05 -1.87e-04 4.05e-04 1.37e-21
    0.00e+00 1.22e-04 8.88e-06 5.59e-05 8.74e-05 -3.61e-04 2.17e-06 2.86e-05 -1.22e-04 -9.33e-22 4.28e-05 -1.79e-05 -1.04e-05 -4.01e-04 -3.09e-04 -6.79e-06 -8.88e-06 -4.28e-05 -3.46e-22 -1.73e-05 -6.41e-05 1.19e-05 5.57e-05 -5.60e-05 -5.59e-05 1.79e-05 1.73e-05 -2.47e-22 -5.39e-06 -1.63e-04 3.92e-05 -4.53e-05 -8.74e-05 1.04e-05 6.41e-05 5.39e-06 5.33e-22 -4.56e-04 -7.95e-05 2.44e-05 3.61e-04 4.01e-04 -1.19e-05 1.63e-04 4.56e-04 9.01e-22 2.05e-04 2.72e-05 -2.17e-06 3.09e-04 -5.57e-05 -3.92e-05 7.95e-05 -2.05e-04 -3.75e-22 -4.80e-04 -2.86e-05 6.79e-06 5.60e-05 4.53e-05 -2.44e-05 -2.72e-05 4.80e-04 -2.07e-21
    0.00e+00 1.70e-04 4.47e-04 -1.45e-04 -3.68e-04 5.04e-04 1.84e-05 -2.10e-04 -1.70e-04 2.16e-21 1.27e-04 7.13e-05 1.68e-04 2.28e-04 3.78e-04 -3.55e-04 -4.47e-04 -1.27e-04 1.80e-21 -4.79e-05 1.93e-04 -7.28e-05 1.50e-04 2.20e-04 1.45e-04 -7.13e-05 4.79e-05 2.08e-22 -3.06e-04 -4.38e-05 -2.54e-04 4.38e-04 3.68e-04 -1.68e-04 -1.93e-04 3.06e-04 -1.24e-22 -1.26e-04 -2.93e-05 -2.78e-05 -5.04e-04 -2.28e-04 7.28e-05 4.38e-05 1.26e-04 -1.34e-21 -1.31e-04 3.41e-04 -1.84e-05 -3.78e-04 -1.50e-04 2.54e-04 2.93e-05 1.31e-04 -3.45e-22 -7.60e-05 2.10e-04 3.55e-04 -2.20e-04 -4.38e-04 2.78e-05 -3.41e-04 7.60e-05 -6.74e-22
    0.00e+00 1.04e-04 4.03e-04 -1.78e-04 -3.09e-04 6.40e-04 3.97e-05 -2.61e-04 -1.04e-04 -5.14e-21 7.07e-05 -1.74e-04 1.68e-04 2.49e-04 3.26e-04 -3.03e-04 -4.03e-04 -7.07e-05 1.37e-21 -4.26e-06 1.52e-04 -1.64e-04 1.16e-04 2.79e-04 1.78e-04 1.74e-04 4.26e-06 2.44e-21 -3.53e-04 -1.53e-04 -1.52e-04 5.74e-04 3.09e-04 -1.68e-04 -1.52e-04 3.53e-04 3.13e-22 -5.42e-05 -3.99e-05 -6.39e-05 -6.40e-04 -2.49e-04 1.64e-04 1.53e-04 5.42e-05 -1.57e-21 -5.95e-05 9.95e-05 -3.97e-05 -3.26e-04 -1.16e-04 1.52e-04 3.99e-05 5.95e-05 -5.43e-23 -1.28e-04 2.61e-04 3.03e-04 -2.79e-04 -5.74e-04 6.39e-05 -9.95e-05 1.28e-04 6.06e-22
    0.00e+00 2.61e-04 3.72e-04 1.79e-04 -4.39e-04 -3.70e-05 -3.75e-04 -8.38e-05 -2.61e-04 -7.35e-21 1.71e-04 2.86e-04 1.62e-04 -2.62e-04 2.73e-04 -8.52e-05 -3.72e-04 -1.71e-04 1.18e-23 1.01e-04 -1.17e-04 -2.46e-04 1.18e-04 1.87e-05 -1.79e-04 -2.86e-04 -1.01e-04 4.61e-22 -7.69e-05 -3.02e-04 5.18e-05 9.75e-05 4.39e-04 -1.62e-04 1.17e-04 7.69e-05 -2.45e-22 3.54e-05 -1.19e-04 9.64e-05 3.70e-05 2.62e-04 2.46e-04 3.02e-04 -3.54e-05 6.44e-21 -1.67e-04 -1.98e-04 3.75e-04 -2.73e-04 -1.18e-04 -5.18e-05 1.19e-04 1.67e-04 -1.45e-21 1.48e-04 8.38e-05 8.52e-05 -1.87e-05 -9.75e-05 -9.64e-05 1.98e-04 -1.48e-04 1.07e-22
    0.00e+00 -2.03e-04 -1.60e-04 -9.86e-05 1.81e-05 2.49e-04 -3.33e-05 3.16e-04 2.03e-04 -2.90e-24 -2.00e-04 -1.80e-04 -2.29e-04 1.09e-04 2.75e-04 9.57e-05 1.60e-04 2.00e-04 -2.68e-21 4.49e-05 -6.37e-05 -1.51e-04 1.76e-04 -2.05e-04 9.86e-05 1.80e-04 -4.49e-05 1.38e-21 -4.52e-04 1.20e-04 1.50e-04 -8.31e-05 -1.81e-05 2.29e-04 6.37e-05 4.52e-04 -5.10e-21 -2.83e-05 -1.28e-06 -4.55e-04 -2.49e-04 -1.09e-04 1.51e-04 -1.20e-04 2.83e-05 6.15e-22 -1.12e-04 1.03e-04 3.33e-05 -2.75e-04 -1.76e-04 -1.50e-04 1.28e-06 1.12e-04 2.34e-22 2.47e-04 -3.16e-04 -9.57e-05 2.05e-04 8.31e-05 4.55e-04 -1.03e-04 -2.47e-04 -8.58e-22
    0.00e+00 8.62e-05 -1.86e-04 4.12e-06 -1.92e-05 -1.15e-04 2.96e-04 8.37e-05 -8.62e-05 3.04e-21 -6.40e-06 2.58e-04 -1.48e-04 1.72e-04 -1.07e-04 -2.92e-04 1.86e-04 6.40e-06 1.21e-21 -4.72e-05 2.57e-05 -2.60e-04 -1.60e-04 3.97e-04 -4.12e-06 -2.58e-04 4.72e-05 -5.30e-22 -8.97e-05 1.79e-04 2.13e-05 -1.46e-04 1.92e-05 1.48e-04 -2.57e-05 8.97e-05 2.95e-22 -5.73e-05 5.24e-05 -5.73e-06 1.15e-04 -1.72e-04 2.60e-04 -1.79e-04 5.73e-05 1.10e-21 -2.12e-04 1.22e-04 -2.96e-04 1.07e-04 1.60e-04 -2.13e-05 -5.24e-05 2.12e-04 4.01e-22 -2.64e-04 -8.37e-05 2.92e-04 -3.97e-04 1.46e-04 5.73e-06 -1.22e-04 2.64e-04 4.15e-22
    0.00e+00 1.03e-04 -2.24e-04 4.91e-05 -4.69e-05 -1.10e-04 2.55e-04 2.79e-05 -1.03e-04 -2.01e-21 6.91e-05 1.82e-05 2.92e-05 2.55e-04 2.68e-04 -2.16e-04 2.24e-04 -6.91e-05 1.84e-22 -1.57e-04 8.54e-05 -2.72e-04 -3.30e-05 4.75e-04 -4.91e-05 -1.82e-05 1.57e-04 2.25e-22 -6.47e-05 2.36e-04 2.10e-04 -8.92e-05 4.69e-05 -2.92e-05 -8.54e-05 6.47e-05 -2.29e-22 -9.38e-05 -6.42e-05 -2.45e-05 1.10e-04 -2.55e-04 2.72e-04 -2.36e-04 9.38e-05 9.97e-22 -1.59e-04 1.38e-04 -2.55e-04 -2.68e-04 3.30e-05 -2.10e-04 6.42e-05 1.59e-04 -1.31e-21 -2.32e-04 -2.79e-05 2.16e-04 -4.75e-04 8.92e-05 2.45e-05 -1.38e-04 2.32e-04 1.61e-21
    0.00e+00 -2.39e-04 8.00e-05 1.33e-04 -3.43e-04 4.48e-04 -1.75e-05 -1.91e-04 2.39e-04 -2.42e-22 2.02e-04 -2.44e-04 -1.39e-04 -8.52e-05 1.57e-04 -1.23e-04 -8.00e-05 -2.02e-04 4.17e-22 3.57e-04 -2.90e-05 2.05e-04 1.08e-05 -9.66e-05 -1.33e-04 2.44e-04 -3.57e-04 -5.42e-22 -7.26e-05 -6.69e-05 -1.44e-05 5.54e-04 3.43e-04 1.39e-04 2.90e-05 7.26e-05 -1.92e-21 -1.80e-04 4.39e-05 2.06e-04 -4.48e-04 8.52e-05 -2.05e-04 6.69e-05 1.80e-04 -1.18e-21 -1.12e-04 1.39e-04 1.75e-05 -1.57e-04 -1.08e-05 1.44e-05 -4.39e-05 1.12e-04 6.55e-22 -1.20e-04 1.91e-04 1.23e-04 9.66e-05 -5.54e-04 -2.06e-04 -1.39e-04 1.20e-04 -6.68e-22
    0.00e+00 1.93e-04 -1.78e-04 -4.34e-05 -1.22e-04 2.30e-04 -1.04e-04 -2.25e-04 -1.93e-04 1.76e-21 -3.46e-06 3.10e-04 1.64e-04 5.49e-04 5.11e-04 -2.76e-04 1.78e-04 3.46e-06 8.90e-22 -9.52e-05 9.50e-05 -3.20e-04 2.13e-04 6.62e-04 4.34e-05 -3.10e-04 9.52e-05 1.48e-21 -1.27e-04 -7.31e-05 -1.88e-04 8.68e-05 1.22e-04 -1.64e-04 -9.50e-05 1.27e-04 9.03e-22 9.54e-06 1.05e-04 -9.96e-06 -2.30e-04 -5.49e-04 3.20e-04 7.31e-05 -9.54e-06 4.58e-22 -9.71e-05 8.12e-06 1.04e-04 -5.11e-04 -2.13e-04 1.88e-04 -1.05e-04 9.71e-05 -2.18e-21 9.15e-05 2.25e-04 2.76e-04 -6.62e-04 -8.68e-05 9.96e-06 -8.12e-06 -9.15e-05 -7.46e-22
    0.00e+00 2.18e-04 -2.77e-04 -1.23e-04 -7.25e-05 1.93e-04 -1.04e-04 -1.15e-04 -2.18e-04 -1.87e-22 -2.63e-05 1.86e-04 1.97e-04 4.80e-04 4.02e-04 -3.02e-04 2.77e-04 2.63e-05 1.28e-21 -2.69e-05 4.05e-05 -3.06e-04 2.56e-04 6.65e-04 1.23e-04 -1.86e-04 2.69e-05 8.33e-21 -1.83e-04 -2.28e-04 -2.53e-04 1.26e-04 7.25e-05 -1.97e-04 -4.05e-05 1.83e-04 -1.81e-21 6.68e-05 1.09e-04 -5.97e-05 -1.93e-04 -4.80e-04 3.06e-04 2.28e-04 -6.68e-05 -6.24e-22 -7.71e-05 -2.17e-04 1.04e-04 -4.02e-04 -2.56e-04 2.53e-04 -1.09e-04 7.71e-05 9.11e-22 3.52e-05 1.15e-04 3.02e-04 -6.65e-04 -1.26e-04 5.97e-05 2.17e-04 -3.52e-05 3.29e-21
    0.00e+00 1.36e-05 -2.40e-04 -1.63e-04 -1.28e-04 -1.20e-04 -3.40e-05 2.15e-04 -1.36e-05 2.34e-21 6.35e-05 2.49e-04 5.74e-05 3.59e-04 -1.81e-05 -2.47e-04 2.40e-04 -6.35e-05 -1.79e-22 -5.04e-06 -2.42e-04 -3.06e-04 -2.65e-04 4.58e-05 1.63e-04 -2.49e-04 5.04e-06 5.70e-22 -6.84e-05 5.75e-05 1.51e-04 -1.41e-04 1.28e-04 -5.74e-05 2.42e-04 6.84e-05 -3.91e-22 -2.36e-05 1.25e-04 -1.47e-04 1.20e-04 -3.59e-04 3.06e-04 -5.75e-05 2.36e-05 -6.42e-22 2.60e-04 -1.16e-04 3.40e-05 1.81e-05 2.65e-04 -1.51e-04 -1.25e-04 -2.60e-04 1.59e-21 -2.97e-05 -2.15e-04 2.47e-04 -4.58e-05 1.41e-04 1.47e-04 1.16e-04 2.97e-05 -3.49e-21
    0.00e+00 -4.82e-05 -4.27e-04 -1.85e-05 -2.56e-05 2.82e-05 7.02e-05 -1.62e-04 4.82e-05 -5.78e-22 -2.08e-04 -2.61e-04 -1.05e-04 1.69e-04 -2.04e-04 -9.11e-05 4.27e-04 2.08e-04 1.28e-21 1.99e-04 -6.72e-05 6.37e-05 -5.43e-05 1.81e-04 1.85e-05 2.61e-04 -1.99e-04 2.59e-21 -2.32e-04 -8.09e-05 1.72e-05 -1.78e-04 2.56e-05 1.05e-04 6.72e-05 2.32e-04 1.79e-21 -6.36e-05 9.71e-05 3.78e-05 -2.82e-05 -1.69e-04 -6.37e-05 8.09e-05 6.36e-05 -2.82e-21 1.54e-04 -7.30e-05 -7.02e-05 2.04e-04 5.43e-05 -1.72e-05 -9.71e-05 -1.54e-04 2.40e-21 7.43e-05 1.62e-04 9.11e-05 -1.81e-04 1.78e-04 -3.78e-05 7.30e-05 -7.43e-05 8.60e-22
    0.00e+00 7.94e-05 -1.92e-04 -1.96e-04 -1.05e-04 -8.81e-05 4.57e-05 3.37e-04 -7.94e-05 -2.30e-21 8.75e-05 2.41e-04 1.23e-04 4.21e-04 1.41e-04 -2.81e-04 1.92e-04 -8.75e-05 -6.14e-22 -8.55e-05 -1.51e-04 -3.20e-04 -5.08e-05 2.47e-04 1.96e-04 -2.41e-04 8.55e-05 2.04e-21 -1.03e-04 2.14e-05 1.79e-04 -1.10e-04 1.05e-04 -1.23e-04 1.51e-04 1.03e-04 1.71e-22 -2.91e-05 6.65e-05 -2.07e-04 8.81e-05 -4.21e-04 3.20e-04 -2.14e-05 2.91e-05 -8.60e-23 2.77e-04 -1.21e-04 -4.57e-05 -1.41e-04 5.08e-05 -1.79e-04 -6.65e-05 -2.77e-04 1.73e-21 -6.38e-05 -3.37e-04 2.81e-04 -2.47e-04 1.10e-04 2.07e-04 1.21e-04 6.38e-05 -3.73e-22
    0.00e+00 4.32e-05 3.11e-04 8.08e-05 -3.68e-04 9.70e-05 -4.32e-04 -1.01e-04 -4.32e-05 -3.33e-21 1.44e-04 4.44e-04 5.17e-05 -3.87e-04 2.97e-04 -7.36e-05 -3.11e-04 -1.44e-04 -9.23e-22 3.11e-05 -6.38e-05 -7.79e-05 2.74e-05 8.39e-05 -8.08e-05 -4.44e-04 -3.11e-05 -5.16e-22 -1.47e-04 -4.25e-04 9.24e-05 8.69e-05 3.68e-04 -5.17e-05 6.38e-05 1.47e-04 4.28e-21 -1.43e-05 8.93e-05 1.09e-05 -9.70e-05 3.87e-04 7.79e-05 4.25e-04 1.43e-05 -5.87e-21 -2.18e-04 -1.67e-04 4.32e-04 -2.97e-04 -2.74e-05 -9.24e-05 -8.93e-05 2.18e-04 -7.78e-23 4.71e-06 1.01e-04 7.36e-05 -8.39e-05 -8.69e-05 -1.09e-05 1.67e-04 -4.71e-06 -1.78e-21
    0.00e+00 1.49e-04 -4.48e-04 -1.61e-04 -2.25e-04 -2.85e-04 8.84e-05 -9.47e-05 -1.49e-04 1.43e-22 -3.61e-05 1.71e-04 -8.40e-05 2.59e-04 -4.10e-04 2.09e-04 4.48e-04 3.61e-05 -2.37e-21 5.65e-05 -4.16e-04 -5.66e-05 2.07e-05 -4.72e-06 1.61e-04 -1.71e-04 -5.65e-05 -3.94e-22 1.21e-04 -5.13e-06 3.29e-04 -3.05e-04 2.25e-04 8.40e-05 4.16e-04 -1.21e-04 -2.32e-22 -9.39e-05 2.39e-04 -6.57e-06 2.85e-04 -2.59e-04 5.66e-05 5.13e-06 9.39e-05 4.43e-22 3.91e-04 -2.29e-04 -8.84e-05 4.10e-04 -2.07e-05 -3.29e-04 -2.39e-04 -3.91e-04 2.92e-21 1.21e-06 9.47e-05 -2.09e-04 4.72e-06 3.05e-04 6.57e-06 2.29e-04 -1.21e-06 1.59e-21
    0.00e+00 7.36e-05 -4.65e-04 -1.59e-04 -1.48e-04 -6.94e-05 7.79e-06 -9.58e-05 -7.36e-05 2.08e-23 -2.17e-05 1.50e-04 9.81e-06 2.35e-04 -3.09e-04 3.59e-05 4.65e-04 2.17e-05 -1.85e-21 1.65e-04 -3.59e-04 -1.76e-04 -2.29e-04 -2.55e-04 1.59e-04 -1.50e-04 -1.65e-04 -1.53e-21 -5.64e-05 -7.80e-05 1.97e-04 -1.62e-04 1.48e-04 -9.81e-06 3.59e-04 5.64e-05 2.82e-22 -6.41e-05 1.36e-04 -4.30e-05 6.94e-05 -2.35e-04 1.76e-04 7.80e-05 6.41e-05 1.39e-21 2.69e-04 -9.43e-05 -7.79e-06 3.09e-04 2.29e-04 -1.97e-04 -1.36e-04 -2.69e-04 -5.76e-22 1.48e-04 9.58e-05 -3.59e-05 2.55e-04 1.62e-04 4.30e-05 9.43e-05 -1.48e-04 -1.17e-21
    0.00e+00 -1.24e-05 5.63e-04 -1.18e-04 -3.29e-04 6.37e-04 1.58e-04 -1.72e-04 1.24e-05 -1.01e-21 1.33e-04 -2.64e-04 4.75e-05 3.80e-05 2.24e-04 -3.02e-04 -5.63e-04 -1.33e-04 2.22e-22 8.52e-05 8.15e-05 2.48e-05 3.34e-05 -4.62e-05 1.18e-04 2.64e-04 -8.52e-05 -2.21e-21 -2.65e-04 -1.51e-04 -1.62e-05 5.63e-04 3.29e-04 -4.75e-05 -8.15e-05 2.65e-04 -8.05e-21 -1.50e-04 -4.74e-05 -8.24e-05 -6.37e-04 -3.80e-05 -2.48e-05 1.51e-04 1.50e-04 -1.14e-20 -5.50e-05 2.37e-04 -1.58e-04 -2.24e-04 -3.34e-05 1.62e-05 4.74e-05 5.50e-05 -1.91e-22 -1.47e-04 1.72e-04 3.02e-04 4.62e-05 -5.63e-04 8.24e-05 -2.37e-04 1.47e-04 -3.93e-21
    0.00e+00 1.90e-04 1.98e-04 -1.16e-04 -2.46e-04 3.81e-04 -9.32e-05 -2.36e-04 -1.90e-04 2.15e-22 3.99e-05 1.35e-04 2.23e-04 3.39e-04 4.02e-04 -3.30e-04 -1.98e-04 -3.99e-05 -3.86e-21 -6.23e-05 1.96e-04 -2.55e-04 1.93e-04 3.80e-04 1.16e-04 -1.35e-04 6.23e-05 6.43e-23 -2.67e-04 -7.84e-05 -3.02e-04 2.54e-04 2.46e-04 -2.23e-04 -1.96e-04 2.67e-04 2.08e-22 -3.67e-05 3.62e-05 -4.65e-05 -3.81e-04 -3.39e-04 2.55e-04 7.84e-05 3.67e-05 -1.58e-21 -1.34e-04 1.38e-04 9.32e-05 -4.02e-04 -1.93e-04 3.02e-04 -3.62e-05 1.34e-04 8.41e-22 -6.95e-06 2.36e-04 3.30e-04 -3.80e-04 -2.54e-04 4.65e-05 -1.38e-04 6.95e-06 1.82e-21
    0.00e+00 1.68e-04 3.13e-04 -4.32e-05 -2.89e-04 4.14e-04 -9.41e-05 -2.84e-04 -1.68e-04 5.38e-22 3.19e-05 3.35e-04 1.71e-04 3.82e-04 4.52e-04 -2.98e-04 -3.13e-04 -3.19e-05 -1.64e-22 -1.10e-04 2.27e-04 -1.78e-04 1.95e-04 4.06e-04 4.32e-05 -3.35e-04 1.10e-04 -3.32e-21 -1.80e-04 5.54e-05 -2.76e-04 1.50e-04 2.89e-04 -1.71e-04 -2.27e-04 1.80e-04 -1.84e-22 -1.07e-04 1.68e-05 2.26e-05 -4.14e-04 -3.82e-04 1.78e-04 -5.54e-05 1.07e-04 -1.16e-21 -1.33e-04 3.64e-04 9.41e-05 -4.52e-04 -1.95e-04 2.76e-04 -1.68e-05 1.33e-04 -1.04e-21 6.63e-05 2.84e-04 2.98e-04 -4.06e-04 -1.50e-04 -2.26e-05 -3.64e-04 -6.63e-05 -7.89e-22
    0.00e+00 -1.04e-04 4.56e-04 7.73e-05 -1.49e-04 -6.47e-05 1.94e-04 -2.30e-04 1.04e-04 1.14e-21 1.70e-04 -2.28e-04 -1.51e-04 -2.43e-04 1.28e-04 1.47e-05 -4.56e-04 -1.70e-04 2.72e-21 1.28e-04 1.85e-05 -3.66e-05 1.07e-04 -2.61e-04 -7.73e-05 2.28e-04 -1.28e-04 1.77e-21 1.85e-04 2.29e-05 1.22e-04 2.09e-04 1.49e-04 1.51e-04 -1.85e-05 -1.85e-04 -1.47e-21 -1.39e-04 -1.16e-04 2.67e-04 6.47e-05 2.43e-04 3.66e-05 -2.29e-05 1.39e-04 -1.19e-21 3.07e-05 1.02e-04 -1.94e-04 -1.28e-04 -1.07e-04 -1.22e-04 1.16e-04 -3.07e-05 -1.29e-21 -3.66e-05 2.30e-04 -1.47e-05 2.61e-04 -2.09e-04 -2.67e-04 -1.02e-04 3.66e-05 -4.86e-22
    0.00e+00 -2.48e-04 3.66e-04 -3.69e-05 1.51e-04 3.87e-04 -7.65e-05 -2.40e-04 2.48e-04 -1.63e-22 -1.94e-04 -1.84e-04 1.05e-04 -1.53e-04 1.36e-04 9.37e-05 -3.66e-04 1.94e-04 -2.07e-21 1.92e-04 1.01e-04 -4.30e-05 4.79e-05 -3.12e-04 3.69e-05 1.84e-04 -1.92e-04 -1.20e-21 -7.94e-05 -2.61e-04 1.32e-04 2.19e-04 -1.51e-04 -1.05e-04 -1.01e-04 7.94e-05 2.34e-22 1.55e-04 -7.71e-05 -1.97e-04 -3.87e-04 1.53e-04 4.30e-05 2.61e-04 -1.55e-04 -4.62e-22 -4.13e-05 -8.49e-05 7.65e-05 -1.36e-04 -4.79e-05 -1.32e-04 7.71e-05 4.13e-05 -2.25e-22 1.09e-05 2.40e-04 -9.37e-05 3.12e-04 -2.19e-04 1.97e-04 8.49e-05 -1.09e-05 3.43e-21
    0.00e+00 -9.25e-05 5.00e-04 -7.24e-05 3.11e-04 2.13e-04 2.17e-04 -2.92e-04 9.25e-05 5.20e-22 -6.01e-05 -3.35e-04 9.43e-05 -2.53e-04 1.24e-04 7.41e-05 -5.00e-04 6.01e-05 -1.27e-21 -2.08e-05 2.33e-04 -3.46e-05 1.85e-05 -3.17e-04 7.24e-05 3.35e-04 2.08e-05 1.24e-21 6.69e-05 -2.29e-04 2.55e-04 2.91e-04 -3.11e-04 -9.43e-05 -2.33e-04 -6.69e-05 2.42e-21 -1.54e-05 -1.66e-04 -3.89e-05 -2.13e-04 2.53e-04 3.46e-05 2.29e-04 1.54e-05 -2.40e-22 8.94e-05 -3.79e-05 -2.17e-04 -1.24e-04 -1.85e-05 -2.55e-04 1.66e-04 -8.94e-05 3.38e-22 -3.88e-05 2.92e-04 -7.41e-05 3.17e-04 -2.91e-04 3.89e-05 3.79e-05 3.88e-05 3.28e-21
    0.00e+00 -2.51e-05 -2.80e-04 -1.16e-04 -1.15e-04 1.41e-04 2.66e-04 -1.39e-05 2.51e-05 1.08e-22 2.63e-05 -3.94e-04 -5.39e-05 1.94e-04 6.28e-05 5.42e-06 2.80e-04 -2.63e-05 -2.22e-22 3.79e-04 -5.80e-05 -1.21e-04 2.37e-04 -2.43e-04 1.16e-04 3.94e-04 -3.79e-04 5.21e-22 -3.83e-04 -6.33e-05 -1.68e-04 -2.52e-04 1.15e-04 5.39e-05 5.80e-05 3.83e-04 5.87e-22 -7.67e-05 -9.11e-05 -1.92e-04 -1.41e-04 -1.94e-04 1.21e-04 6.33e-05 7.67e-05 -3.62e-22 6.24e-05 1.16e-04 -2.66e-04 -6.28e-05 -2.37e-04 1.68e-04 9.11e-05 -6.24e-05 -7.89e-22 2.85e-04 1.39e-05 -5.42e-06 2.43e-04 2.52e-04 1.92e-04 -1.16e-04 -2.85e-04 2.02e-21
    0.00e+00 -1.90e-04 -2.12e-04 -2.04e-05 6.37e-06 2.08e-04 4.90e-05 1.97e-04 1.90e-04 -4.02e-21 -1.62e-04 -3.24e-04 -2.86e-04 1.01e-04 1.85e-04 -9.44e-05 2.12e-04 1.62e-04 1.54e-21 -6.61e-05 1.66e-05 -1.72e-04 1.95e-04 -2.41e-04 2.04e-05 3.24e-04 6.61e-05 -9.72e-22 -3.82e-04 7.23e-05 2.80e-05 -2.60e-04 -6.37e-06 2.86e-04 -1.66e-05 3.82e-04 -9.60e-22 -9.39e-07 -3.49e-05 -4.30e-04 -2.08e-04 -1.01e-04 1.72e-04 -7.23e-05 9.39e-07 2.45e-22 -7.60e-05 1.17e-04 -4.90e-05 -1.85e-04 -1.95e-04 -2.80e-05 3.49e-05 7.60e-05 -7.11e-21 2.40e-04 -1.97e-04 9.44e-05 2.41e-04 2.60e-04 4.30e-04 -1.17e-04 -2.40e-04 1.58e-23
    0.00e+00 2.20e-04 -1.56e-04 -1.77e-04 -1.11e-04 3.04e-04 -1.21e-04 -1.10e-04 -2.20e-04 -1.63e-21 -1.42e-05 3.94e-06 2.55e-04 4.30e-04 3.40e-04 -2.84e-04 1.56e-04 1.42e-05 3.12e-21 -2.05e-05 8.97e-05 -3.29e-04 2.53e-04 5.48e-04 1.77e-04 -3.94e-06 2.05e-05 -9.60e-22 -2.81e-04 -2.17e-04 -3.06e-04 2.45e-04 1.11e-04 -2.55e-04 -8.97e-05 2.81e-04 -5.47e-22 2.90e-05 5.47e-05 -1.25e-04 -3.04e-04 -4.30e-04 3.29e-04 2.17e-04 -2.90e-05 -4.14e-22 -1.16e-04 -2.74e-04 1.21e-04 -3.40e-04 -2.53e-04 3.06e-04 -5.47e-05 1.16e-04 -5.62e-22 1.60e-05 1.10e-04 2.84e-04 -5.48e-04 -2.45e-04 1.25e-04 2.74e-04 -1.60e-05 -8.12e-23
    0.00e+00 -5.59e-05 1.52e-04 -6.76e-05 3.39e-04 1.67e-04 1.74e-04 4.54e-04 5.59e-05 3.31e-22 -4.67e-05 -2.99e-04 -2.34e-04 -1.18e-04 1.09e-04 8.73e-05 -1.52e-04 4.67e-05 6.35e-22 1.21e-04 6.92e-05 6.90e-05 -8.21e-05 -2.11e-04 6.76e-05 2.99e-04 -1.21e-04 -2.42e-21 2.26e-04 2.99e-04 4.23e-04 1.65e-04 -3.39e-04 2.34e-04 -6.92e-05 -2.26e-04 -2.52e-21 3.11e-04 2.46e-04 6.04e-05 -1.67e-04 1.18e-04 -6.90e-05 -2.99e-04 -3.11e-04 -1.46e-21 9.21e-05 -6.44e-05 -1.74e-04 -1.09e-04 8.21e-05 -4.23e-04 -2.46e-04 -9.21e-05 1.49e-22 1.94e-04 -4.54e-04 -8.73e-05 2.11e-04 -1.65e-04 -6.04e-05 6.44e-05 -1.94e-04 -6.02e-21
    0.00e+00 -1.38e-04 2.36e-04 -3.13e-05 4.61e-04 1.89e-04 9.84e-06 1.68e-04 1.38e-04 4.91e-23 -1.62e-04 -3.55e-04 6.53e-05 -2.01e-04 2.01e-04 -1.18e-04 -2.36e-04 1.62e-04 -7.11e-22 -8.79e-05 2.05e-06 -5.59e-05 3.38e-05 -2.92e-04 3.13e-05 3.55e-04 8.79e-05 6.80e-23 3.11e-06 -1.78e-04 2.48e-04 -7.46e-05 -4.61e-04 -6.53e-05 -2.05e-06 -3.11e-06 -5.65e-22 2.29e-04 -2.01e-04 -3.57e-04 -1.89e-04 2.01e-04 5.59e-05 1.78e-04 -2.29e-04 2.73e-22 7.81e-07 -9.96e-05 -9.84e-06 -2.01e-04 -3.38e-05 -2.48e-04 2.01e-04 -7.81e-07 5.00e-22 6.92e-06 -1.68e-04 1.18e-04 2.92e-04 7.46e-05 3.57e-04 9.96e-05 -6.92e-06 -1.28e-21
    0.00e+00 1.88e-04 2.51e-05 -4.11e-05 -1.68e-04 2.79e-04 -1.18e-04 -2.88e-04 -1.88e-04 8.71e-24 1.95e-05 3.22e-04 1.77e-04 5.22e-04 4.75e-04 -2.72e-04 -2.51e-05 -1.95e-05 -1.34e-22 -8.28e-05 1.63e-04 -3.22e-04 2.14e-04 5.26e-04 4.11e-05 -3.22e-04 8.28e-05 6.03e-22 -1.50e-04 -4.95e-05 -2.61e-04 1.25e-04 1.68e-04 -1.77e-04 -1.63e-04 1.50e-04 -1.93e-22 -1.73e-05 5.51e-05 7.86e-06 -2.79e-04 -5.22e-04 3.22e-04 4.95e-05 1.73e-05 1.21e-21 -1.22e-04 1.36e-04 1.18e-04 -4.75e-04 -2.14e-04 2.61e-04 -5.51e-05 1.22e-04 1.38e-22 9.49e-05 2.88e-04 2.72e-04 -5.26e-04 -1.25e-04 -7.86e-06 -1.36e-04 -9.49e-05 3.30e-22
    0.00e+00 3.40e-04 2.16e-04 9.96e-05 -1.72e-04 -4.32e-04 1.09e-04 -6.25e-05 -3.40e-04 -5.34e-22 1.67e-04 2.76e-04 -6.36e-06 -1.92e-04 4.46e-05 2.01e-04 -2.16e-04 -1.67e-04 1.12e-22 2.90e-05 -2.31e-04 1.60e-04 -7.14e-05 -7.85e-06 -9.96e-05 -2.76e-04 -2.90e-05 -1.70e-22 -1.08e-04 1.31e-04 1.24e-04 5.76e-05 1.72e-04 6.36e-06 2.31e-04 1.08e-04 1.01e-20 -4.14e-04 3.64e-04 1.32e-04 4.32e-04 1.92e-04 -1.60e-04 -1.31e-04 4.14e-04 4.03e-22 1.11e-04 -6.52e-05 -1.09e-04 -4.46e-05 7.14e-05 -1.24e-04 -3.64e-04 -1.11e-04 1.26e-21 -1.42e-04 6.25e-05 -2.01e-04 7.85e-06 -5.76e-05 -1.32e-04 6.52e-05 1.42e-04 -1.51e-21
    0.00e+00 2.19e-04 -2.29e-04 -1.42e-04 -5.78e-05 2.77e-04 -1.31e-04 -1.13e-04 -2.19e-04 8.58e-22 -2.60e-05 9.93e-05 2.29e-04 5.05e-04 3.44e-04 -3.00e-04 2.29e-04 2.60e-05 -1.32e-22 -8.04e-06 4.18e-05 -3.31e-04 2.51e-04 6.49e-04 1.42e-04 -9.93e-05 8.04e-06 -9.55e-22 -2.20e-04 -2.78e-04 -3.16e-04 2.07e-04 5.78e-05 -2.29e-04 -4.18e-05 2.20e-04 2.17e-22 8.48e-07 4.51e-05 -9.70e-05 -2.77e-04 -5.05e-04 3.31e-04 2.78e-04 -8.48e-07 1.26e-22 -1.02e-04 -2.58e-04 1.31e-04 -3.44e-04 -2.51e-04 3.16e-04 -4.51e-05 1.02e-04 -2.01e-21 5.49e-05 1.13e-04 3.00e-04 -6.49e-04 -2.07e-04 9.70e-05 2.58e-04 -5.49e-05 4.14e-25
    0.00e+00 -2.62e-04 3.50e-04 -1.84e-04 3.83e-04 1.13e-04 1.18e-04 2.95e-04 2.62e-04 -7.09e-22 -1.37e-04 -3.99e-04 -1.17e-04 -2.55e-04 2.22e-04 6.74e-05 -3.50e-04 1.37e-04 1.10e-21 8.00e-05 -8.57e-07 4.27e-05 -5.78e-05 -3.02e-04 1.84e-04 3.99e-04 -8.00e-05 -9.53e-22 -9.93e-05 -9.14e-05 3.77e-04 -2.54e-05 -3.83e-04 1.17e-04 8.57e-07 9.93e-05 1.57e-21 2.25e-04 -1.58e-04 -3.22e-04 -1.13e-04 2.55e-04 -4.27e-05 9.14e-05 -2.25e-04 1.30e-21 4.70e-05 -1.70e-04 -1.18e-04 -2.22e-04 5.78e-05 -3.77e-04 1.58e-04 -4.70e-05 -4.60e-22 1.74e-04 -2.95e-04 -6.74e-05 3.02e-04 2.54e-05 3.22e-04 1.70e-04 -1.74e-04 -3.38e-22
    0.00e+00 8.83e-05 1.69e-04 1.28e-04 -1.65e-04 4.89e-04 -2.94e-07 7.73e-05 -8.83e-05 -3.64e-21 1.22e-04 1.28e-04 1.22e-04 3.30e-04 1.68e-04 9.96e-05 -1.69e-04 -1.22e-04 -1.72e-21 6.84e-05 -5.97e-05 2.45e-04 4.74e-05 5.39e-05 -1.28e-04 -1.28e-04 -6.84e-05 -1.82e-21 6.97e-05 1.92e-04 1.08e-04 1.07e-04 1.65e-04 -1.22e-04 5.97e-05 -6.97e-05 -1.72e-21 -1.86e-04 -7.48e-05 1.11e-05 -4.89e-04 -3.30e-04 -2.45e-04 -1.92e-04 1.86e-04 -5.81e-22 -1.80e-04 8.97e-05 2.94e-07 -1.68e-04 -4.74e-05 -1.08e-04 7.48e-05 1.80e-04 -3.02e-21 1.61e-04 -7.73e-05 -9.96e-05 -5.39e-05 -1.07e-04 -1.11e-05 -8.97e-05 -1.61e-04 1.97e-21
    0.00e+00 -1.41e-04 3.32e-04 -7.39e-05 3.39e-04 2.48e-04 3.19e-05 -4.13e-05 1.41e-04 8.68e-23 -1.60e-04 -3.04e-04 1.32e-04 -2.13e-04 2.03e-04 -8.29e-08 -3.32e-04 1.60e-04 -1.19e-21 -3.72e-05 7.67e-05 -6.03e-05 5.12e-05 -3.08e-04 7.39e-05 3.04e-04 3.72e-05 -5.15e-22 1.84e-05 -2.44e-04 2.70e-04 6.09e-05 -3.39e-04 -1.32e-04 -7.67e-05 -1.84e-05 1.71e-22 1.88e-04 -1.85e-04 -3.19e-04 -2.48e-04 2.13e-04 6.03e-05 2.44e-04 -1.88e-04 -6.69e-22 2.02e-05 -6.22e-05 -3.19e-05 -2.03e-04 -5.12e-05 -2.70e-04 1.85e-04 -2.02e-05 -1.38e-21 -1.75e-05 4.13e-05 8.29e-08 3.08e-04 -6.09e-05 3.19e-04 6.22e-05 1.75e-05 -8.56e-22
    0.00e+00 -1.24e-04 2.85e-04 -2.44e-04 3.11e-04 1.33e-04 2.15e-04 3.47e-04 1.24e-04 1.50e-21 -7.03e-05 -3.77e-04 -7.98e-05 -1.96e-04 1.94e-04 1.12e-04 -2.85e-04 7.03e-05 -1.05e-21 9.42e-05 5.99e-05 5.52e-05 -6.20e-05 -2.69e-04 2.44e-04 3.77e-04 -9.42e-05 1.81e-22 1.46e-04 1.15e-04 5.28e-04 2.62e-04 -3.11e-04 7.98e-05 -5.99e-05 -1.46e-04 -4.79e-21 1.64e-04 8.95e-06 -1.63e-04 -1.33e-04 1.96e-04 -5.52e-05 -1.15e-04 -1.64e-04 -2.05e-21 9.11e-05 -1.71e-05 -2.15e-04 -1.94e-04 6.20e-05 -5.28e-04 -8.95e-06 -9.11e-05 1.05e-21 1.77e-04 -3.47e-04 -1.12e-04 2.69e-04 -2.62e-04 1.63e-04 1.71e-05 -1.77e-04 -7.06e-22
    0.00e+00 1.25e-04 1.16e-04 -1.12e-04 -4.87e-04 3.33e-04 9.86e-05 6.66e-05 -1.25e-04 -1.10e-21 1.97e-04 1.31e-04 -1.21e-05 -8.79e-05 2.29e-04 -2.39e-04 -1.16e-04 -1.97e-04 -1.26e-21 2.25e-04 -1.40e-04 2.39e-04 1.47e-04 2.59e-05 1.12e-04 -1.31e-04 -2.25e-04 -1.87e-21 -8.10e-05 2.82e-04 -1.27e-04 3.64e-04 4.87e-04 1.21e-05 1.40e-04 8.10e-05 5.43e-22 -3.62e-04 -5.95e-05 4.45e-05 -3.33e-04 8.79e-05 -2.39e-04 -2.82e-04 3.62e-04 -7.19e-22 -1.63e-04 1.52e-04 -9.86e-05 -2.29e-04 -1.47e-04 1.27e-04 5.95e-05 1.63e-04 9.27e-22 -2.05e-05 -6.66e-05 2.39e-04 -2.59e-05 -3.64e-04 -4.45e-05 -1.52e-04 2.05e-05 -2.60e-21
    0.00e+00 1.30e-04 -3.00e-04 -2.78e-04 -1.08e-04 1.72e-04 1.48e-04 4.00e-04 -1.30e-04 6.99e-22 1.24e-04 -3.27e-04 1.50e-04 2.66e-04 9.44e-05 9.54e-05 3.00e-04 -1.24e-04 -1.49e-21 2.67e-04 -1.19e-04 -8.29e-05 3.77e-04 1.52e-04 2.78e-04 3.27e-04 -2.67e-04 -9.70e-22 -3.94e-04 -2.19e-04 -1.22e-04 -1.58e-05 1.08e-04 -1.50e-04 1.19e-04 3.94e-04 -9.04e-22 -1.20e-04 -2.57e-05 -3.55e-04 -1.72e-04 -2.66e-04 8.29e-05 2.19e-04 1.20e-04 -4.82e-23 -3.46e-05 -2.16e-04 -1.48e-04 -9.44e-05 -3.77e-04 1.22e-04 2.57e-05 3.46e-05 -6.96e-23 2.39e-04 -4.00e-04 -9.54e-05 -1.52e-04 1.58e-05 3.55e-04 2.16e-04 -2.39e-04 9.85e-22
    0.00e+00 -4.44e-05 5.32e-05 7.96e-05 -2.28e-04 4.07e-05 -5.00e-05 -2.26e-04 4.44e-05 3.47e-22 -3.72e-05 1.09e-04 8.60e-05 1.13e-04 6.00e-04 1.13e-04 -5.32e-05 3.72e-05 -2.72e-21 -1.77e-04 2.16e-04 4.31e-05 1.30e-04 3.35e-04 -7.96e-05 -1.09e-04 1.77e-04 1.48e-21 -1.07e-04 2.57e-04 3.12e-04 2.89e-05 2.28e-04 -8.60e-05 -2.16e-04 1.07e-04 -4.36e-22 -1.48e-04 -9.04e-06 -1.16e-05 -4.07e-05 -1.13e-04 -4.31e-05 -2.57e-04 1.48e-04 4.15e-22 -3.97e-04 1.81e-04 5.00e-05 -6.00e-04 -1.30e-04 -3.12e-04 9.04e-06 3.97e-04 1.46e-21 -9.49e-06 2.26e-04 -1.13e-04 -3.35e-04 -2.89e-05 1.16e-05 -1.81e-04 9.49e-06 -2.34e-21
    0.00e+00 1.86e-04 -5.85e-05 -9.71e-05 -1.13e-04 2.17e-04 -1.38e-04 -2.59e-04 -1.86e-04 -4.93e-22 1.28e-05 2.07e-04 1.89e-04 4.45e-04 3.93e-04 -3.30e-04 5.85e-05 -1.28e-05 4.35e-22 -3.50e-05 1.26e-04 -3.47e-04 2.08e-04 4.92e-04 9.71e-05 -2.07e-04 3.50e-05 2.18e-22 -1.84e-04 -1.81e-04 -3.02e-04 1.35e-04 1.13e-04 -1.89e-04 -1.26e-04 1.84e-04 -4.96e-22 3.43e-05 5.29e-05 1.61e-05 -2.17e-04 -4.45e-04 3.47e-04 1.81e-04 -3.43e-05 -4.12e-22 -6.78e-05 -3.45e-05 1.38e-04 -3.93e-04 -2.08e-04 3.02e-04 -5.29e-05 6.78e-05 -9.07e-22 1.46e-05 2.59e-04 3.30e-04 -4.92e-04 -1.35e-04 -1.61e-05 3.45e-05 -1.46e-05 -1.90e-22
    0.00e+00 1.38e-04 2.28e-04 1.11e-04 -3.66e-04 2.12e-04 -2.41e-04 -4.48e-06 -1.38e-04 -1.05e-22 2.46e-04 4.36e-04 7.92e-05 1.28e-04 2.83e-04 4.58e-05 -2.28e-04 -2.46e-04 -3.11e-23 6.94e-05 -1.17e-04 -5.99e-05 -7.99e-05 3.00e-05 -1.11e-04 -4.36e-04 -6.94e-05 -1.90e-24 -5.67e-05 -2.54e-04 -3.37e-07 8.53e-05 3.66e-04 -7.92e-05 1.17e-04 5.67e-05 3.77e-21 7.42e-05 1.26e-04 6.88e-05 -2.12e-04 -1.28e-04 5.99e-05 2.54e-04 -7.42e-05 2.92e-22 -5.04e-05 4.43e-05 2.41e-04 -2.83e-04 7.99e-05 3.37e-07 -1.26e-04 5.04e-05 1.09e-21 -1.49e-05 4.48e-06 -4.58e-05 -3.00e-05 -8.53e-05 -6.88e-05 -4.43e-05 1.49e-05 3.17e-22 ];

Qmat  = zeros(8,8,314);
for n = 1:314
    Qmat(:,:,n) = reshape(complex(locQmatr(n,:),locQmati(n,:)),[8,8]);
end

end

function Qmat = Load_default_SARglobal_Qmat


gloQmatr = [2.23e-04 -4.26e-05 1.17e-05 -1.64e-05 -1.25e-06 7.21e-06 -6.51e-06 -1.76e-05 -4.26e-05 2.52e-04 -1.41e-05 1.38e-05 9.67e-06 1.01e-05 4.44e-06 -2.20e-05 1.17e-05 -1.41e-05 2.40e-04 -3.35e-05 2.57e-06 -2.82e-05 -5.63e-06 -3.55e-06 -1.64e-05 1.38e-05 -3.35e-05 2.16e-04 6.73e-06 -2.53e-05 -1.65e-05 3.56e-06 -1.25e-06 9.67e-06 2.57e-06 6.73e-06 2.44e-04 1.09e-05 1.82e-05 -5.62e-06 7.21e-06 1.01e-05 -2.82e-05 -2.53e-05 1.09e-05 2.53e-04 -1.16e-05 -6.35e-06 -6.51e-06 4.44e-06 -5.63e-06 -1.65e-05 1.82e-05 -1.16e-05 2.81e-04 -1.35e-07 -1.76e-05 -2.20e-05 -3.55e-06 3.56e-06 -5.62e-06 -6.35e-06 -1.35e-07 2.26e-04 ];

gloQmati = [0.00e+00 1.61e-05 9.51e-06 -2.05e-08 -1.03e-05 -9.75e-06 1.76e-05 5.88e-06 -1.61e-05 0.00e+00 2.57e-06 7.14e-07 -8.91e-06 -4.60e-06 1.83e-06 -1.10e-05 -9.51e-06 -2.57e-06 0.00e+00 4.59e-06 -1.37e-05 -1.39e-05 -1.31e-06 -5.70e-06 2.05e-08 -7.14e-07 -4.59e-06 0.00e+00 -1.68e-05 -2.96e-06 1.87e-05 -6.04e-06 1.03e-05 8.91e-06 1.37e-05 1.68e-05 0.00e+00 -1.46e-06 4.21e-06 3.73e-06 9.75e-06 4.60e-06 1.39e-05 2.96e-06 1.46e-06 0.00e+00 -6.44e-06 -3.38e-06 -1.76e-05 -1.83e-06 1.31e-06 -1.87e-05 -4.21e-06 6.44e-06 0.00e+00 -1.28e-06 -5.88e-06 1.10e-05 5.70e-06 6.04e-06 -3.73e-06 3.38e-06 1.28e-06 0.00e+00 ];

Qmat = reshape(complex(gloQmatr,gloQmati),[8,8]);
end