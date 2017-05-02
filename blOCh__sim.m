function varargout = blOCh__sim(Fun,spc,khr,opt,varargin)
% function [sim,Msg] = blOCh__sim(spc,khr,opt,varargin)
%
%   This script simulates an output from blOCh
%
%
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
%

if isempty(Fun)
    warning off
    Status = 0;
    varargout{1} = -1;
    
    if spc.pTx > 1
        valid_sim_Txs = {'all','sep'};
        def_sim_Tx = valid_sim_Txs{1};
        
    else
        valid_sim_Txs = {'all'};
        def_sim_Tx = valid_sim_Txs{1};
    end
    
    valid_sim_ks = {'last','all','firstlast'};
    def_sim_k = valid_sim_ks{1};
    
    valid_sim_ns = {'last','all'};
    def_sim_n = valid_sim_ns{1};
    
    valid_opt_ns = {'opt1','opt2'};
    def_opt_n = 'opt1';
    valid_sim_Rlxs = [1,0];
    def_sim_Rlx = valid_sim_Rlxs(2);
    
    valid_Shows = [1,0];
    def_Show = valid_Shows(1);
    
    
    
    valid_SaveFigs = [1,0];
    def_SaveFig = valid_SaveFigs(1);
    
    valid_Nuc      = {'1H','13C','19F'};
    
    def_Nuc = '1H';
    def_B1inhom_lim = [1,1];
    def_B1inhom_N = 1;
    def_B0inhom_lim = [0;0];
    def_B0inhom_N = 1;
    def_dt  = [];
    def_u = [];
    def_v = [];
    def_g = [];
    
    p = inputParser;
    
    try p.addRequired('Fun', @(x)isstruct(x)||isempty(x));
        try p.addRequired('spc', @(x)isstruct(x));
            try p.addRequired('khr', @(x)isstruct(x)||isempty(x));
                try p.addRequired('opt', @(x)isstruct(x)||isempty(x));
                    try p.addParamValue('opt_n',def_opt_n,@(x)any(valid_opt_ns));
                        try p.addParamValue('sim_Rlx',def_sim_Rlx,@(x)any(valid_sim_Rlxs));
                            try p.addParamValue('sim_Tx',def_sim_Tx,@(x)any(strcmpi(x,valid_sim_Txs)));
                                try p.addParamValue('sim_k',def_sim_k,@(x)Validate_ks(x));
                                    try p.addParamValue('sim_n',def_sim_n,@(x)any(strcmpi(x,valid_sim_ns))); % ||validateattributes(x,{'numeric'},{'size',[1,1],'real','integer','>=',1,'<=',opt.N})
                                        try p.addParamValue('Show',def_Show,@(x)any(valid_Shows));
                                            
                                            try p.addParamValue('B1inhom_lim',def_B1inhom_lim,@(x)validateattributes(x,{'numeric'},{'size',[1,2],'real','>=',0,'<=',2}));
                                                try p.addParamValue('B1inhom_N',def_B1inhom_N,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>=',1,'real','finite'}));
                                                    try p.addParamValue('B0inhom_lim',def_B0inhom_lim,@(x)validateattributes(x,{'numeric'},{'size',[1,2],'real','finite'}));
                                                        try p.addParamValue('B0inhom_N',def_B0inhom_N,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>=',1,'real','finite'}));
                                                            
                                                            try p.addParamValue('dt',def_dt,@(x)validateattributes(x,{'numeric'},{'size',[NaN,NaN],'real','finite','>',0}));
                                                                try p.addParamValue('u',def_u,@(x)validateattributes(x,{'numeric'},{'size',[NaN,NaN,NaN],'real','finite'}));
                                                                    try p.addParamValue('v',def_v,@(x)validateattributes(x,{'numeric'},{'size',[NaN,NaN,NaN],'real','finite'}));
                                                                        try p.addParamValue('g',def_g,@(x)validateattributes(x,{'numeric'},{'size',[NaN,NaN,NaN],'real','finite'}));
                                                                            try p.addParamValue('Nuc',def_Nuc,@(x)any(strcmpi(x,valid_Nuc)));
                                                                                try p.addParamValue('SaveFig',def_SaveFig,@(x)any(valid_SaveFigs));
                                                                                    Status = 1;
                                                                                    
                                                                                catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                                                                            catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                                                                        catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                                                                    catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                                                                catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                                                            catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                                                            
                                                        catch  me;Display_Message(['blOCh__4_0__sim: ',me.message,' B0inhom_N: Resolution to simulate B0 inhomogeneity in'],2);end
                                                    catch  me;Display_Message(['blOCh__4_0__sim: ',me.message,' B0inhom_lim: should be a [a,b] array of units Hz like [-200,200] to simulate RF inhomogeneity from -200 Hz to 200 Hz'],2);end
                                                catch  me;Display_Message(['blOCh__4_0__sim: ',me.message,' B1inhom_N: Resolution to simulate RF inhomogeneity in'],2);end
                                            catch  me;Display_Message(['blOCh__4_0__sim: ',me.message,' B1inhom_lim: should be a [a,b] array like [0.9,1.1] to simulate RF inhomogeneity of 90% to 110%'],2);end
                                        catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                                    catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                                catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
                            catch  me;Display_Message(['blOCh__sim: ',me.message,' sim_n: the time frame(s) the simulation stores; ''last'', ''all'', or # (number of dt between frames)'],2); end
                        catch  me;Display_Message(['blOCh__sim: ',me.message,' sim_k: the iteration(s) the simulation does; ''last'', ''all'', ''firstlast'', or # (a specific iteration)'],2); end
                    catch  me;Display_Message(['blOCh__sim: ',me.message,' sim_Tx:'],2);end
                catch  me;Display_Message(['blOCh__sim: ',me.message,' sim_Rlx:'],2);end
            catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
        catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
    catch  me;Display_Message(['blOCh__sim: ',me.message],2);end
    
    if Status
        try
            p.parse(Fun,spc,khr,opt,varargin{:})
            
            sim.sim_n = p.Results.sim_n;
            sim.opt_n = p.Results.opt_n;
            sim.sim_Tx = p.Results.sim_Tx;
            sim.sim_k = p.Results.sim_k;
            sim.sim_Rlx = p.Results.sim_Rlx;
            sim.Show = p.Results.Show;
            sim.SaveFig = p.Results.SaveFig;
            
            sim.B1inhom_lim = p.Results.B1inhom_lim;
            sim.B1inhom_N = p.Results.B1inhom_N;
            sim.B0inhom_lim = p.Results.B0inhom_lim;
            sim.B0inhom_N = p.Results.B0inhom_N;
            
            
            sim.temp.u = p.Results.u;
            sim.temp.v = p.Results.v;
            sim.temp.g = p.Results.g;
            sim.temp.dt = p.Results.dt;
            sim.temp.Nuc = p.Results.Nuc;
            switch sim.sim_Tx
                case 'all'
                    sim.S = 1;
                    sim.s = 1;
                case 'sep'
                    sim.S = spc.pTx+1;
                    sim.s = [1:sim.S];
            end
            
        catch me
            Display_Message(['blOCh__sim: ',me.message],2)
            
        end
        
        if ~isempty(opt) && ~isempty(khr)
            
            if isdef(varargin,'g')
                Display_Message(sprintf('blOCh__sim: both khr and g were specified. g will be ignored'),1);
            end
            if isdef(varargin,'u') || isdef(varargin,'v')
                Display_Message(sprintf('blOCh__sim: both opt and u or v were specified. u,v will be ignored'),1);
            end
            if isdef(varargin,'dt')
                Display_Message(sprintf('blOCh__sim: both opt and dt were specified. dt will be ignored'),1);
            end
            if isdef(varargin,'Nuc')
                Display_Message(sprintf('blOCh__sim: both khr and Nuc were specified. Nuc will be ignored'),1);
            end
            
            if strcmp(sim.opt_n,'opt1')
                sim.uo = opt.opt1.uo;
                sim.vo = opt.opt1.vo;
                Kact = opt.opt1.ksafe;
                sim.MaxIter = opt.opt1.MaxIter;
                sim.mon = opt.opt1.mon;
                sim.g = opt.opt1.g;
                
                sim.OptNum = 1;
            elseif strcmp(sim.opt_n,'opt2')
                sim.uo = opt.opt2.uo;
                sim.vo = opt.opt2.vo;
                Kact = opt.opt2.ksafe;
                sim.MaxIter = opt.opt2.MaxIter;
                sim.mon = opt.opt2.mon;
                sim.g = opt.opt2.g;
                
                sim.OptNum = 2;
            end
            sim.Kact = Kact;
            sim.dt = opt.dt;
            sim.gamma = khr.gamma;
            Nact = opt.N;
            
        elseif isempty(opt) && ~isempty(khr)
            
            if isdef(varargin,'g')
                Display_Message(sprintf('blOCh__sim: both khr and g were specified. g will be ignored'),1);
            end
            if isdef(varargin,'u') && isdef(varargin,'v')
                Display_Message(sprintf('blOCh__sim: both u and v were specified. '),1);
            elseif isdef(varargin,'u') && ~isdef(varargin,'v')
                Display_Message(sprintf('blOCh__sim: only u was specified. Aborting '),2);
                return
            elseif ~isdef(varargin,'u') && isdef(varargin,'v')
                Display_Message(sprintf('blOCh__sim: only v was specified. Aborting '),2);
                return
            end
            if isdef(varargin,'dt')
                Display_Message(sprintf('blOCh__sim: dt was specified.'),1);
            end
            if isdef(varargin,'Nuc')
                Display_Message(sprintf('blOCh__sim: both khr and Nuc were specified. Nuc will be ignored'),1);
            end
            
            
            %             sim.dt = sim.temp.dt; % bug 170309
            sim.dt = khr.dt;
            sim.uo = sim.temp.u;
            sim.vo = sim.temp.v;
            sim.Nuc = khr.Nuc;
            sim.gamma = khr.gamma;
            Nact = size(sim.uo,2);
            Kact = size(sim.uo,3);
            sim.Kact = Kact;
            sim.MaxIter = Kact-1;
            sim.N = Nact;
            sim.OptNum = 1;
            sim.Mask = 0;
            sim = blOCh__opt('Prepare_khr_4_opt',[],[],[],[],sim,khr);
            
            
        elseif isempty(opt) && isempty(khr)
            
            if isdef(varargin,'g')
                Display_Message(sprintf('blOCh__sim: g was specified.'),1);
            else
                Display_Message(sprintf('blOCh__sim: g needs to be specified. Aborting '),2);
                return
            end
            if isdef(varargin,'u') && isdef(varargin,'v')
                Display_Message(sprintf('blOCh__sim: both u and v were specified. '),1);
            elseif isdef(varargin,'u') && ~isdef(varargin,'v')
                Display_Message(sprintf('blOCh__sim: only u was specified. Aborting '),2);
                return
            elseif ~isdef(varargin,'u') && isdef(varargin,'v')
                Display_Message(sprintf('blOCh__sim: only v was specified. Aborting '),2);
                return
            end
            if isdef(varargin,'dt')
                Display_Message(sprintf('blOCh__sim: dt was specified.'),1);
            else
                Display_Message(sprintf('blOCh__sim: dt needs to be specified. Aborting '),2);
                return
            end
            if isdef(varargin,'Nuc')
                Display_Message(sprintf('blOCh__sim: Nuc was specified. '),1);
            else
                Display_Message(sprintf('blOCh__sim: Nuc wasn''t specified. Assuming it''s 1H'),1);
                
            end
            
            
            
            
            
            sim.dt = sim.temp.dt;
            sim.g = sim.temp.g;
            sim.uo = sim.temp.u;
            sim.vo = sim.temp.v;
            sim.Nuc = sim.temp.Nuc;
            Nact = size(sim.uo,2);
            Kact = size(sim.uo,3);
            sim.MaxIter = 1;
            sim.OptNum = 1;
            sim.mon = true(1,Nact);
            %             khr.Nuc = sim.Nuc;
            switch sim.Nuc
                case '1H'
                    sim.gamma = 2.6751289e8;
                case '13C'
                    sim.gamma = 67.262e6;
                case '19F'
                    sim.gamma = 251.662e6;
            end
            khr.gamma = sim.gamma;
            khr.GmHW = 1;
            khr.Gmpct = 1;
            khr.SmHW = 1;
            khr.Smpct = 1;
            
        end
        
        switch sim.sim_n
            
            case 'last'
                sim.M = 1;
                sim.N = Nact;
            case 'all'
                sim.M = Nact+1;
                sim.N = Nact;
                
        end
        
        
        switch sim.sim_k
            
            case 'last'
                sim.K = 1;
                sim.k = Kact;
            case 'all'
                sim.K = Kact;
                sim.k = [1:Kact];
            case 'firstlast'
                sim.K =  2;
                sim.k = [1,Kact];
            otherwise
                sim.K = 1;
                sim.k = sim.sim_k;
        end
        
        
        if sim.sim_Rlx == 1
            sim.Rlx = 1;
        else
            sim.Rlx = 0;
        end
        
        sim.M_t_big = zeros(3*spc.P,sim.M,sim.K,sim.S,sim.B1inhom_N,sim.B0inhom_N);
        
        if exist('blOCh__Get_Relaxors','file') && sim.sim_Rlx == 1
            
            [sim.FRlxT1,sim.FRlxT2,sim.CRlxT1] = blOCh__Get_Relaxors(spc.T1map,spc.T2map,sim.dt,spc.P);
        end
        
        for b0i = 1:sim.B0inhom_N
            
            sim.B0inhom_offsets = linspace(sim.B0inhom_lim(1),sim.B0inhom_lim(2),sim.B0inhom_N);
            sim.B0inhom_offset = sim.B0inhom_offsets(b0i);
            
            for rfi = 1:sim.B1inhom_N
                
                sim.B1inhom_scales = linspace(sim.B1inhom_lim(1),sim.B1inhom_lim(2),sim.B1inhom_N);
                sim.B1inhom_scale = sim.B1inhom_scales(rfi);
                
                
                
                k_counter = 1;
                for k = sim.k
                    if size(sim.g,3)>1
                        g = sim.g(:,:,k);
                    else
                        g = sim.g;
                    end
                    for s = sim.s
                        
                        switch sim.sim_Tx
                            
                            case 'all'
                                sim.u = sim.uo(:,:,k);
                                sim.v = sim.vo(:,:,k);
                            case 'sep'
                                if s == sim.S
                                    sim.u = sim.uo(:,:,k);
                                    sim.v = sim.vo(:,:,k);
                                else
                                    
                                    sim.u = zeros(size(sim.uo,1),size(sim.uo,2));
                                    sim.v = zeros(size(sim.vo,1),size(sim.vo,2));
                                    sim.u(s,:) = sim.uo(s,:,k);
                                    sim.v(s,:) = sim.vo(s,:,k);
                                end
                        end
                        
                        
                        sim = blOCh__opt('Allocate_Variables',[],[],[],[],spc,khr,sim);
                        sim = rmfield(sim,'L_t');
                        sim = rmfield(sim,'Durations');
                        sim = rmfield(sim,'Fun');
                        sim = rmfield(sim,'Eff');
                        sim = rmfield(sim,'Pen');
                        sim = rmfield(sim,'dFun');
                        %                         sim = rmfield(sim,'MaxIter');
                        
                        
                        
                        if exist('blOCh__Rotate_Relax','file') && sim.sim_Rlx == 1
                            for n = 1:sim.N
                                [R11f,R12f,R13f,R21f,R22f,R23f,R31f,R32f,R33f] = blOCh__opt('Get_Rotator',[],[],[],[],sim.u.*sim.B1inhom_scale,sim.v.*sim.B1inhom_scale,g,sim.w0+2*pi*sim.B0inhom_offset,sim.yx,sim.yy,sim.yz,spc.pTx,sim.sr,sim.si,sim.dt,n);
                                sim.M_t(:,n+1) = blOCh__Rotate_Relax(sim.M_t(:,n),R11f,R12f,R13f,R21f,R22f,R23f,R31f,R32f,R33f,sim.FRlxT1,sim.FRlxT2,sim.CRlxT1,'Forward');
                            end
                            
                            
                        else
                            for n = 1:sim.N
                                [R11f,R12f,R13f,R21f,R22f,R23f,R31f,R32f,R33f] = blOCh__opt('Get_Rotator',[],[],[],[],sim.u.*sim.B1inhom_scale,sim.v.*sim.B1inhom_scale,g,sim.w0+2*pi*sim.B0inhom_offset,sim.yx,sim.yy,sim.yz,spc.pTx,sim.sr,sim.si,sim.dt,n);
                                sim.M_t(:,n+1) = blOCh__opt('Rotate',[],[],[],[],sim.M_t(:,n),R11f,R12f,R13f,R21f,R22f,R23f,R31f,R32f,R33f,'Forward');
                                
                            end
                            
                        end
                        
                        
                        sim.idxtot = spc.idxtot;
                        switch sim.sim_n
                            
                            case 'last'
                                sim.M_t_big(:,1,k_counter,s,rfi,b0i) = sim.M_t(:,end);
                            case 'all'
                                sim.M_t_big(:,:,k_counter,s,rfi,b0i) = sim.M_t;
                        end
                        
                    end
                    k_counter = k_counter + 1;
                end
            end
        end
        sim.M_t = sim.M_t_big;
        sim = rmfield(sim,'M_t_big');
        %         size(sim.M_t_big)
        sim.Dim = spc.Dim;
        sim.pTx = spc.pTx;
        sim.R = spc.R;
        sim.Rv = spc.Rv;
        sim.D = spc.D;
        sim.Dv = spc.Dv;
        sim.Md = spc.Md;
        sim.idxtot=spc.idxtot;
        
        if sim.Show > 0
            sim.fig = Show_Sim(sim);
        end
        
        
        
        
        
        varargout{1} = sim;
    end
else
    Nargout = nargout;
    fun = str2func(Fun);
    test = version('-release');
    if strcmp(test,'2015a')
        v = [];
    else
        v = {};
    end
    
    switch Nargout
        case 0
            fun(varargin{:});
        case 1
            v{1} = fun(varargin{:});
        case 2
            [v{1},v{2}] = fun(varargin{:});
        case 3
            [v{1},v{2},v{3}] = fun(varargin{:});
        case 4
            [v{1},v{2},v{3},v{4}] = fun(varargin{:});
        case 5
            [v{1},v{2},v{3},v{4},v{5}] = fun(varargin{:});
        case 6
            [v{1},v{2},v{3},v{4},v{5},v{6}] = fun(varargin{:});
        case 7
            [v{1},v{2},v{3},v{4},v{5},v{6},v{7}] = fun(varargin{:});
        case 8
            [v{1},v{2},v{3},v{4},v{5},v{6},v{7},v{8}] = fun(varargin{:});
        case 9
            [v{1},v{2},v{3},v{4},v{5},v{6},v{7},v{8},v{9}] = fun(varargin{:});
        case 10
            [v{1},v{2},v{3},v{4},v{5},v{6},v{7},v{8},v{9},v{10}] = fun(varargin{:});
        otherwise
            Display_Message(['spc: Fun: This switch needs to be expanded'],2)
            
    end
    
    varargout = v;
    
end
end

function test = Validate_rf(x,spc,khr,RFinputtype)

test = false;

if isempty(x) && RFinputtype == 1
    
    test = true;
    
elseif ~isempty(x) && RFinputtype == 2
    
    
    
    
    
end




end

function Out = isdef(In,Case)
Out = false;
for n = 1:length(In)
    if ischar(In{n})
        if strcmp(In{n},Case)
            Out = true;
        end
    end
end


end

function test = Validate_ks(x)

test = false;

if ischar(x)
    
    switch x
        
        case {'last','all','firstlast'}
            
            test = true;
        otherwise
            
    end
    
elseif isnumeric(x)
    
    
    
    test = true;
    
end

end


function Display_Message(Msg,Type)
if nargin == 1
    Type = 0;
end

if iscell(Msg)
    
    for n = 1:length(Msg)
        Msg_temp = Msg{n};
        Msg_temp = Repair_Msg(Msg_temp);
        
        if Type == 1
            
            fprintf(1,[Msg_temp,'\n']);
            
        elseif Type == 2
            
            fprintf(2,[Msg_temp,'\n']);
        elseif Type == 3
            
            fprintf(1,[Msg_temp]);
            
        else
            
            fprintf(Type,[Msg,'\n']);
        end
    end
else
    Msg_temp = Msg;
    Msg_temp = Repair_Msg(Msg_temp);
    
    if Type == 1
        
        fprintf(1,[Msg_temp,'\n']);
        
    elseif Type == 2
        
        fprintf(2,[Msg_temp,'\n']);
    elseif Type == 0
        
    elseif Type == 3
        
        fprintf(1,[Msg]);
    else
        
        fprintf(1,[Msg,'\n']);
    end
end


end

function Msg = Repair_Msg(Msg)
% function Msg = Repair_Msg(Msg)

Msg = regexprep(Msg,'\\','msv8');

Msg = regexprep(Msg,'msv8','\\\');
end

function fig = Show_Sim(varargin)
global hsim

hsim.sim = varargin{1};


hsim.fig = figure('Visible','on');


mp = get(0, 'MonitorPositions');
if size(mp,1) > 1
    mp = mp(1,:);
end

set(hsim.fig,'Units','pixels')
set(hsim.fig,'Position',[1,1,mp(3).*0.9,mp(4).*0.9])
WHratio = get(hsim.fig,'Position');
WHratio = WHratio(3)/WHratio(4);
set(hsim.fig,'Units','normalized')

%%



clear varargin
%%

hsim.axes_Magn = axes('Parent',hsim.fig,'Visible','on');
hsim.axes_Colorbar = axes('Parent',hsim.fig,'Visible','on');
hsim.axes_Pulse = axes('Parent',hsim.fig,'Visible','on');


hsim.listbox1= uicontrol('Style','listbox','Callback',@listbox1_Callback);

hsim.slider_n= uicontrol('Style','slider','Callback',@slider_n_Callback);
hsim.slider_Sl= uicontrol('Style','slider','Callback',@slider_Sl_Callback);
hsim.slider_Tx= uicontrol('Style','slider','Callback',@slider_Tx_Callback);
hsim.slider_freq= uicontrol('Style','slider','Callback',@slider_freq_Callback);
hsim.slider_k= uicontrol('Style','slider','Callback',@slider_k_Callback);
hsim.slider_B1= uicontrol('Style','slider','Callback',@slider_B1_Callback);
hsim.slider_B0= uicontrol('Style','slider','Callback',@slider_B0_Callback);
hsim.popupmenu_Magn= uicontrol('Style','popupmenu','Callback',@popupmenu_Magn_Callback);
hsim.popupmenu_Pulse= uicontrol('Style','popupmenu','Callback',@popupmenu_Pulse_Callback);

hsim.pushbutton_playx1= uicontrol('Style','pushbutton','Callback',@pushbutton_playx1_Callback);
hsim.pushbutton_playx10= uicontrol('Style','pushbutton','Callback',@pushbutton_playx10_Callback);
hsim.pushbutton_playx100= uicontrol('Style','pushbutton','Callback',@pushbutton_playx100_Callback);

hsim.text1= uicontrol('Style','text');
hsim.text2= uicontrol('Style','text');
hsim.text_n= uicontrol('Style','text');
hsim.text_Sl= uicontrol('Style','text');
hsim.text_Tx= uicontrol('Style','text');
hsim.text_freq= uicontrol('Style','text');
hsim.text_k= uicontrol('Style','text');
hsim.text_B1= uicontrol('Style','text');
hsim.text_B0= uicontrol('Style','text');


Sliderheight = 0.02;
Textheight = 0.02;
Textwidth = 0.04;
Spacer  = 0.001;

set(hsim.slider_n,'Units','normalized')
set(hsim.slider_Sl,'Units','normalized')
set(hsim.slider_Tx,'Units','normalized')
set(hsim.slider_freq,'Units','normalized')
set(hsim.slider_k,'Units','normalized')
set(hsim.slider_B1,'Units','normalized')
set(hsim.slider_B0,'Units','normalized')

set(hsim.pushbutton_playx1,'Units','normalized')
set(hsim.pushbutton_playx10,'Units','normalized')
set(hsim.pushbutton_playx100,'Units','normalized')


set(hsim.text1,'Units','normalized')
set(hsim.text2,'Units','normalized')
set(hsim.text_n,'Units','normalized')
set(hsim.text_Sl,'Units','normalized')
set(hsim.text_Tx,'Units','normalized')
set(hsim.text_freq,'Units','normalized')
set(hsim.text_k,'Units','normalized')
set(hsim.text_B1,'Units','normalized')
set(hsim.text_B0,'Units','normalized')

set(hsim.axes_Magn,'Units','normalized')
set(hsim.axes_Colorbar,'Units','normalized')
set(hsim.axes_Pulse,'Units','normalized')
set(hsim.listbox1,'Units','normalized')
% set(hsim.listboxPulse,'Units','normalized')

set(hsim.popupmenu_Magn,'Units','normalized')

set(hsim.popupmenu_Pulse,'Units','normalized')


%%

hsim.r.x_slider_n = 0.69;
hsim.r.x_slider_Sl = 0.69;
hsim.r.x_slider_Tx = 0.69;
hsim.r.x_slider_freq = 0.69;
hsim.r.x_slider_k = 0.69;
hsim.r.x_slider_B1 = 0.69;
hsim.r.x_slider_B0 = 0.69;



hsim.r.y_slider_B0 = 0.6;
hsim.r.y_slider_B1 = hsim.r.y_slider_B0+Sliderheight+0.001;
hsim.r.y_slider_k = hsim.r.y_slider_B1+Sliderheight+0.001;
hsim.r.y_slider_freq = hsim.r.y_slider_k+Sliderheight+0.001;
hsim.r.y_slider_Tx = hsim.r.y_slider_freq+Sliderheight+0.001;
hsim.r.y_slider_Sl = hsim.r.y_slider_Tx+Sliderheight+0.001;
hsim.r.y_slider_n = hsim.r.y_slider_Sl+Sliderheight+0.001;

hsim.r.y_text_B0 = hsim.r.y_slider_B0;
hsim.r.y_text_B1 = hsim.r.y_slider_B1;
hsim.r.y_text_k = hsim.r.y_slider_k;
hsim.r.y_text_freq = hsim.r.y_slider_freq;
hsim.r.y_text_Tx = hsim.r.y_slider_Tx;
hsim.r.y_text_Sl = hsim.r.y_slider_Sl;
hsim.r.y_text_n = hsim.r.y_slider_n;


hsim.r.y_text1 = 0.91;%0.095+Sliderheight+Spacer;
hsim.r.y_text2 = 0.2;% 0.05+Sliderheight+Spacer;

hsim.r.w_slider_n = 0.3;
hsim.r.w_slider_Sl = 0.3;
hsim.r.w_slider_Tx = 0.3;
hsim.r.w_slider_freq = 0.3;
hsim.r.w_slider_k = 0.3;
hsim.r.w_slider_B1 = 0.3;
hsim.r.w_slider_B0 = 0.3;

hsim.r.h_slider_n = Sliderheight;
hsim.r.h_slider_Sl = Sliderheight;
hsim.r.h_slider_Tx = Sliderheight;
hsim.r.h_slider_freq = Sliderheight;
hsim.r.h_slider_k = Sliderheight;
hsim.r.h_slider_B1 = Sliderheight;
hsim.r.h_slider_B0 = Sliderheight;

hsim.r.x_axes_Magn = 0.05;
hsim.r.xo_axes_Magn = 0.01;
hsim.r.yo_axes_Magn = 0.25;
hsim.r.y_axes_Magn = 0.3;
hsim.r.ho_axes_Magn = 0.7;
hsim.r.h_axes_Magn = 0.65;
hsim.r.wo_axes_Magn = hsim.r.ho_axes_Magn/WHratio;
hsim.r.w_axes_Magn = hsim.r.h_axes_Magn/WHratio;
% WHratio
hsim.r.x_axes_Colorbar = 0.625;
hsim.r.y_axes_Colorbar = 0.25;
hsim.r.h_axes_Colorbar = 0.65;
hsim.r.w_axes_Colorbar = 0.02;

hsim.r.x_axes_Pulse = 0.05;
hsim.r.xo_axes_Pulse = 0.01;
hsim.r.yo_axes_Pulse = 0.01;
hsim.r.y_axes_Pulse = 0.05;
hsim.r.ho_axes_Pulse = 0.20;
hsim.r.h_axes_Pulse = 0.19;
hsim.r.wo_axes_Pulse = 0.95;
hsim.r.w_axes_Pulse = 0.9;






hsim.r.w_text1 =0.05;% 0.98;
hsim.r.w_text2 =0.05;% 0.48;
hsim.r.w_text_n =0.19;% 0.48;
hsim.r.w_text_Sl =0.19;% 0.48;
hsim.r.w_text_Tx =0.19;% 0.48;
hsim.r.w_text_freq =0.19;% 0.2;
hsim.r.w_text_k =0.19;% 0.2;
hsim.r.w_text_B1 =0.19;% 0.2;
hsim.r.w_text_B0 =0.19;% 0.2;

hsim.r.x_text1 = 0.62;
hsim.r.x_text2 = 0.62;
hsim.r.x_text_n = 0.49;
hsim.r.x_text_Sl = 0.49;
hsim.r.x_text_Tx = 0.49;
hsim.r.x_text_freq = 0.49;
hsim.r.x_text_k = 0.49;
hsim.r.x_text_B1 = 0.49;
hsim.r.x_text_B0 = 0.49;

hsim.r.h_text1 = Textheight;
hsim.r.h_text2 = Textheight;
hsim.r.h_text_n = Textheight;
hsim.r.h_text_Sl = Textheight;
hsim.r.h_text_Tx = Textheight;
hsim.r.h_text_freq = Textheight;
hsim.r.h_text_k = Textheight;
hsim.r.h_text_B1 = Textheight;
hsim.r.h_text_B0 = Textheight;

hsim.r.x_popupmenu_Magn = 0.005;
hsim.r.y_popupmenu_Magn = 0.94;
hsim.r.w_popupmenu_Magn = 0.1;
hsim.r.h_popupmenu_Magn = 0.05;

hsim.r.x_popupmenu_Pulse = 0.005+0.1+0.01;
hsim.r.y_popupmenu_Pulse = 0.94;
hsim.r.w_popupmenu_Pulse = 0.1;
hsim.r.h_popupmenu_Pulse = 0.05;


hsim.r.x_listbox1 = 0.69;
hsim.r.y_listbox1 = 0.8;%0.2+0.266667*2;
hsim.r.w_listbox1 = 0.3;
hsim.r.h_listbox1 = 0.18;% 0.26;



hsim.r.x_pushbutton_playx1 = 0.69;
hsim.r.y_pushbutton_playx1 = 0.77;
hsim.r.w_pushbutton_playx1 = 0.1;
hsim.r.h_pushbutton_playx1 = 0.02;

hsim.r.x_pushbutton_playx10 = 0.69+0.1;
hsim.r.y_pushbutton_playx10 = 0.77;
hsim.r.w_pushbutton_playx10 = 0.1;
hsim.r.h_pushbutton_playx10 = 0.02;

hsim.r.x_pushbutton_playx100 = 0.69+0.2;
hsim.r.y_pushbutton_playx100 = 0.77;
hsim.r.w_pushbutton_playx100 = 0.1;
hsim.r.h_pushbutton_playx100 = 0.02;

set(hsim.slider_n,'Position',[hsim.r.x_slider_n,hsim.r.y_slider_n,hsim.r.w_slider_n,hsim.r.h_slider_n])
set(hsim.slider_Sl,'Position',[hsim.r.x_slider_Sl,hsim.r.y_slider_Sl,hsim.r.w_slider_Sl,hsim.r.h_slider_Sl])
set(hsim.slider_Tx,'Position',[hsim.r.x_slider_Tx,hsim.r.y_slider_Tx,hsim.r.w_slider_Tx,hsim.r.h_slider_Tx])
set(hsim.slider_freq,'Position',[hsim.r.x_slider_freq,hsim.r.y_slider_freq,hsim.r.w_slider_freq,hsim.r.h_slider_freq])
set(hsim.slider_k,'Position',[hsim.r.x_slider_k,hsim.r.y_slider_k,hsim.r.w_slider_k,hsim.r.h_slider_k])
set(hsim.slider_B1,'Position',[hsim.r.x_slider_B1,hsim.r.y_slider_B1,hsim.r.w_slider_B1,hsim.r.h_slider_B1])
set(hsim.slider_B0,'Position',[hsim.r.x_slider_B0,hsim.r.y_slider_B0,hsim.r.w_slider_B0,hsim.r.h_slider_B0])

set(hsim.text1,'Position',[hsim.r.x_text1,hsim.r.y_text1,hsim.r.w_text1,hsim.r.h_text1])
set(hsim.text2,'Position',[hsim.r.x_text2,hsim.r.y_text2,hsim.r.w_text2,hsim.r.h_text2])
set(hsim.text_n,'Position',[hsim.r.x_text_n,hsim.r.y_text_n,hsim.r.w_text_n,hsim.r.h_text_n],'HorizontalAlignment','right')
set(hsim.text_Sl,'Position',[hsim.r.x_text_Sl,hsim.r.y_text_Sl,hsim.r.w_text_Sl,hsim.r.h_text_Sl],'HorizontalAlignment','right')
set(hsim.text_Tx,'Position',[hsim.r.x_text_Tx,hsim.r.y_text_Tx,hsim.r.w_text_Tx,hsim.r.h_text_Tx],'HorizontalAlignment','right')
set(hsim.text_freq,'Position',[hsim.r.x_text_freq,hsim.r.y_text_freq,hsim.r.w_text_freq,hsim.r.h_text_freq],'HorizontalAlignment','right')
set(hsim.text_k,'Position',[hsim.r.x_text_k,hsim.r.y_text_k,hsim.r.w_text_k,hsim.r.h_text_k],'HorizontalAlignment','right')
set(hsim.text_B1,'Position',[hsim.r.x_text_B1,hsim.r.y_text_B1,hsim.r.w_text_B1,hsim.r.h_text_B1],'HorizontalAlignment','right')
set(hsim.text_B0,'Position',[hsim.r.x_text_B0,hsim.r.y_text_B0,hsim.r.w_text_B0,hsim.r.h_text_B0],'HorizontalAlignment','right')

set(hsim.axes_Magn,'OuterPosition',[hsim.r.xo_axes_Magn,hsim.r.yo_axes_Magn,hsim.r.wo_axes_Magn,hsim.r.ho_axes_Magn])
set(hsim.axes_Magn,'Position',[hsim.r.x_axes_Magn,hsim.r.y_axes_Magn,hsim.r.w_axes_Magn,hsim.r.h_axes_Magn])

set(hsim.axes_Pulse,'OuterPosition',[hsim.r.xo_axes_Pulse,hsim.r.yo_axes_Pulse,hsim.r.wo_axes_Pulse,hsim.r.ho_axes_Pulse])
set(hsim.axes_Pulse,'Position',[hsim.r.x_axes_Pulse,hsim.r.y_axes_Pulse,hsim.r.w_axes_Pulse,hsim.r.h_axes_Pulse])

set(hsim.axes_Colorbar,'Position',[hsim.r.x_axes_Colorbar,hsim.r.y_axes_Colorbar,hsim.r.w_axes_Colorbar,hsim.r.h_axes_Colorbar])
set(hsim.listbox1,'Position',[hsim.r.x_listbox1,hsim.r.y_listbox1,hsim.r.w_listbox1,hsim.r.h_listbox1])

set(hsim.popupmenu_Magn,'Position',[hsim.r.x_popupmenu_Magn,hsim.r.y_popupmenu_Magn,hsim.r.w_popupmenu_Magn,hsim.r.h_popupmenu_Magn])
set(hsim.popupmenu_Pulse,'Position',[hsim.r.x_popupmenu_Pulse,hsim.r.y_popupmenu_Pulse,hsim.r.w_popupmenu_Pulse,hsim.r.h_popupmenu_Pulse])

set(hsim.pushbutton_playx1,'Position',[hsim.r.x_pushbutton_playx1,hsim.r.y_pushbutton_playx1,hsim.r.w_pushbutton_playx1,hsim.r.h_pushbutton_playx1],'String','Play x1')
set(hsim.pushbutton_playx10,'Position',[hsim.r.x_pushbutton_playx10,hsim.r.y_pushbutton_playx10,hsim.r.w_pushbutton_playx10,hsim.r.h_pushbutton_playx10],'String','Play x10')
set(hsim.pushbutton_playx100,'Position',[hsim.r.x_pushbutton_playx100,hsim.r.y_pushbutton_playx100,hsim.r.w_pushbutton_playx100,hsim.r.h_pushbutton_playx100],'String','Play x100')

if str2double(hsim.sim.Dim(1)) > 1
    set(hsim.axes_Colorbar,'Visible','on')
    set(hsim.text1,'Visible','on')
    set(hsim.text2,'Visible','on')
else
    set(hsim.axes_Colorbar,'Visible','off')
    set(hsim.text1,'Visible','off')
    set(hsim.text2,'Visible','off')
end
set(hsim.text_freq,'Visible','on')
set(hsim.text_k,'Visible','on')
%%
if strcmpi(hsim.sim.sim_n,'last')
    set(hsim.pushbutton_playx1,'Visible','off')
    set(hsim.pushbutton_playx10,'Visible','off')
    set(hsim.pushbutton_playx100,'Visible','off')
    set(hsim.slider_n,'Visible','off','Value',1)
    set(hsim.text_n,'Visible','off');
    hsim.nNo = 1;
    set(hsim.pushbutton_playx1,'Visible','off')
elseif strcmpi(hsim.sim.sim_n,'all')
    set(hsim.pushbutton_playx1,'Visible','on')
    set(hsim.pushbutton_playx10,'Visible','on')
    set(hsim.pushbutton_playx100,'Visible','on')
    set(hsim.slider_n,'Min',1)
    set(hsim.slider_n,'Max',hsim.sim.N)
    TimeSliderStep = [1, 1] / (hsim.sim.N - 1);
    set(hsim.slider_n,'SliderStep',TimeSliderStep)
    hsim.nNo = hsim.sim.N;
    set(hsim.slider_n,'Visible','on','Value',hsim.nNo)
    set(hsim.text_n,'Visible','on','String',sprintf('Timeframe %i dt (%1.1e s) of %i dt (%1.1e s)',hsim.nNo,hsim.nNo*hsim.sim.dt,hsim.sim.N,hsim.sim.N.*hsim.sim.dt))
    
end

%%
if hsim.sim.Dim(1) == '3'
    set(hsim.slider_Sl,'Min',1)
    set(hsim.slider_Sl,'Max',hsim.sim.R(3))
    SliceSliderStep = [1, 1] / (hsim.sim.R(3) - 1);
    set(hsim.slider_Sl,'SliderStep',SliceSliderStep)
    hsim.SlNo = round(hsim.sim.R(3)/2);
    set(hsim.slider_Sl,'Visible','on','Value',hsim.SlNo)
    set(hsim.text_Sl,'Visible','on','String',sprintf('Slice %i of %i',hsim.SlNo,hsim.sim.R(3)))
    
else
    set(hsim.slider_Sl,'Visible','off','Value',1)
    set(hsim.text_Sl,'Visible','off');
    hsim.SlNo = 1;
end
%%


if hsim.sim.pTx > 1 && strcmpi(hsim.sim.sim_Tx,'sep')
    set(hsim.slider_Tx,'Min',1)
    set(hsim.slider_Tx,'Max',hsim.sim.pTx+1)
    ChannelSliderStep = [1, 1] / (hsim.sim.pTx+1 - 1);
    set(hsim.slider_Tx,'SliderStep',ChannelSliderStep)
    hsim.pTxNo = hsim.sim.pTx+1;
    set(hsim.slider_Tx,'Visible','on','Value',hsim.pTxNo)
    if hsim.pTxNo == hsim.sim.pTx+1
        set(hsim.text_Tx,'Visible','on','String',sprintf('All %i Tx',hsim.sim.pTx))
    else
        set(hsim.text_Tx,'Visible','on','String',sprintf('Tx %i of %i',hsim.pTxNo,hsim.sim.pTx))
    end
else
    set(hsim.slider_Tx,'Visible','off','Value',1)
    set(hsim.text_Tx,'Visible','off');
    hsim.pTxNo = 1;
end

%%

if hsim.sim.Dim(2) == '+' && hsim.sim.Dim(1) ~= '0'
    hsim.freq = linspace(hsim.sim.Dv(1),hsim.sim.Dv(2),hsim.sim.Rv);
    
    set(hsim.slider_freq,'Min',1)
    set(hsim.slider_freq,'Max',hsim.sim.Rv)
    FreqSliderStep = [1, 1] / (hsim.sim.Rv - 1);
    set(hsim.slider_freq,'SliderStep',FreqSliderStep)
    hsim.freqNo = round(hsim.sim.Rv/2);
    
    if hsim.sim.Dim(1) == '1'
        
        set(hsim.slider_freq,'Visible','off')
        set(hsim.text_freq,'Visible','off')
    else
        set(hsim.slider_freq,'Visible','on','Value',hsim.freqNo)
        set(hsim.text_freq,'Visible','on','String',sprintf('Frequency: %.2f [Hz]',hsim.freq(hsim.freqNo)))
        
    end
else
    set(hsim.slider_freq,'Visible','off','Value',1)
    set(hsim.text_freq,'Visible','off');
    hsim.freqNo = 1;
end
%%
if hsim.sim.K > 1
    set(hsim.slider_k,'Min',1)
    set(hsim.slider_k,'Max',hsim.sim.K)
    IterationSliderStep = [1, 1] / (hsim.sim.K - 1);
    set(hsim.slider_k,'SliderStep',IterationSliderStep)
    hsim.kmNo = hsim.sim.K;
    hsim.kpNo = hsim.sim.Kact;
    set(hsim.slider_k,'Visible','on','Value',hsim.kmNo)
    
    set(hsim.text_k,'Visible','on','String',sprintf('Iteration %i of %i',hsim.kmNo-1,hsim.sim.Kact-1))
else
    set(hsim.slider_k,'Visible','off','Value',1)
    set(hsim.text_k,'Visible','off');
    hsim.kmNo = 1;
    hsim.kpNo = hsim.sim.Kact;
end
%%

%%
if hsim.sim.B1inhom_N==1
    set(hsim.slider_B1,'Visible','off','Value',1)
    set(hsim.text_B1,'Visible','off');
    hsim.nB1inh = 1;
else
    set(hsim.slider_B1,'Min',1)
    set(hsim.slider_B1,'Max',hsim.sim.B1inhom_N)
    B1SliderStep = [1, 1] / (hsim.sim.B1inhom_N - 1);
    set(hsim.slider_B1,'SliderStep',B1SliderStep)
    hsim.nB1inh = round(hsim.sim.B1inhom_N/2);
    set(hsim.slider_B1,'Visible','on','Value',hsim.nB1inh)
    set(hsim.text_B1,'Visible','on','String',sprintf('B1 scale %i of %i (%1.1e %%)',hsim.nB1inh,hsim.sim.B1inhom_N,hsim.sim.B1inhom_scales(hsim.nB1inh)))
    
end
%%
if hsim.sim.B0inhom_N==1
    set(hsim.slider_B0,'Visible','off','Value',1)
    set(hsim.text_B0,'Visible','off');
    hsim.nB0inh = 1;
else
    set(hsim.slider_B0,'Min',1)
    set(hsim.slider_B0,'Max',hsim.sim.B1inhom_N)
    B0SliderStep = [1, 1] / (hsim.sim.B0inhom_N - 1);
    set(hsim.slider_B0,'SliderStep',B0SliderStep)
    hsim.nB0inh = round(hsim.sim.B0inhom_N/2);
    set(hsim.slider_B0,'Visible','on','Value',hsim.nB0inh)
    set(hsim.text_B0,'Visible','on','String',sprintf('B0 offset %i of %i (%1.1e Hz)',hsim.nB0inh,hsim.sim.B0inhom_N,hsim.sim.B0inhom_offsets(hsim.nB0inh)))
    
end
%%
List1 = Populate_Listbox(hsim.sim);

set(hsim.listbox1,'String',List1);


%%

String = cell(1,9);
String{1} = 'Magnetization, |Mxy| [M0]';
String{2} = 'Magnetization, arg(Mxy) [rad]';

String{3} = 'Magnetization, Mx [M0]';
String{4} = 'Magnetization, My [M0]';
String{5} = 'Magnetization, Mz [M0]';

String{6} = 'Desired Magnetization, |Mxy| [M0]';
String{7} = 'Desired Magnetization, Mx [M0]';
String{8} = 'Desired Magnetization, My [M0]';
String{9} = 'Desired Magnetization, Mz [M0]';

set(hsim.popupmenu_Magn, 'String', String);
set(hsim.popupmenu_Magn, 'Value',1);

%%

String = cell(1,14);
String{1} = 'RF Pulse Amplitude [rad/s]';
String{2} = 'RF Pulse Phase [rad]';
String{3} = 'RF Pulse Real [rad/s]';
String{4} = 'RF Pulse Imaginary [rad/s]';
String{5} = 'Gradient x [T/m]';
String{6} = 'Gradient y [T/m]';
String{7} = 'Gradient z [T/m]';
String{8} = 'Gradient x, y, and z [T/m]';
String{9} = 'RF Pulse Amplitude & Gradients x,y, and z';
String{10} = 'kx [1/m]';
String{11} = 'ky [1/m]';
String{12} = 'kz [1/m]';
String{13} = 'kx, ky, and kz [1/m]';
String{14} = 'RF Pulse Amplitude & kx, ky, and kz ';

set(hsim.popupmenu_Pulse, 'String', String);
set(hsim.popupmenu_Pulse, 'Value',1);
%%
[hsim.colmap.Jet,hsim.colmap.Gray] = GetColormaps;

hsim.Colmap = hsim.colmap.Jet;
colormap jet

hsim = Plotting(hsim);

fig = hsim.fig;
end

function popupmenu_Magn_Callback(hOb, ed)
global hsim
hsim.popupmenu_Magnselect = get(hsim.popupmenu_Magn,'Value');
hsim = Plotting(hsim);
end

function popupmenu_Pulse_Callback(hOb, ed)
global hsim

hsim.popupmenu_Pulseselect = get(hsim.popupmenu_Pulse,'Value');
hsim = Plotting(hsim);
end

function listbox1_Callback(hOb, ed)
end


function slider_n_Callback(hOb, ed)
global hsim
hsim.nNo = round(get(hOb, 'Value'));
set(hOb, 'Value',hsim.nNo);

set(hsim.slider_n,'Visible','on','Value',hsim.nNo)
set(hsim.text_n,'Visible','on','String',sprintf('Timeframe %i dt (%1.1e s) of %i dt (%1.1e s)',hsim.nNo,hsim.nNo*hsim.sim.dt,hsim.sim.N,hsim.sim.N.*hsim.sim.dt))
hsim = Plotting(hsim);
end

function slider_Sl_Callback(hOb, ed)
global hsim
hsim.SlNo = round(get(hOb, 'Value'));
% hsim.SlNo
set(hsim.text_Sl,'String',sprintf('Slice %i of %i',hsim.SlNo,hsim.sim.R(3)))
hsim = Plotting(hsim);
end

function slider_Tx_Callback(hOb, ed)
global hsim
hsim.pTxNo = round(get(hOb, 'Value'));
hsim.pTxNo
if hsim.pTxNo == hsim.sim.pTx+1
    set(hsim.text_Tx,'Visible','on','String',sprintf('All %i Tx',hsim.sim.pTx))
else
    set(hsim.text_Tx,'Visible','on','String',sprintf('Tx %i of %i',hsim.pTxNo,hsim.sim.pTx))
end

hsim = Plotting(hsim);
end

function slider_freq_Callback(hOb, ed)
global hsim
hsim.freqNo = round(get(hOb, 'Value'));
set(hOb, 'Value',hsim.freqNo);
set(hsim.text_freq,'Visible','on','String',sprintf('Frequency: %.2f [Hz]',hsim.freq(hsim.freqNo)))
hsim = Plotting(hsim);
end

function slider_k_Callback(hOb, ed)
global hsim

switch hsim.sim.sim_k
    
    case 'last'
        hsim.kmNo = hsim.sim.k;
        hsim.kpNo = hsim.sim.k;
        set(hOb, 'Value',hsim.kmNo);
    case 'all'
        hsim.kmNo = round(get(hOb, 'Value'));
        hsim.kpNo = hsim.kmNo;
        set(hOb, 'Value',hsim.kmNo);
        set(hsim.text_k,'Visible','on','String',sprintf('Iteration %i of %i',hsim.kmNo-1,hsim.sim.Kact-1))
    case 'firstlast'
        if round(get(hOb, 'Value')) == 1
            hsim.kmNo = 1;
            hsim.kpNo = 1;
        elseif round(get(hOb, 'Value')) == 2
            hsim.kmNo = 2;
            hsim.kpNo = hsim.sim.Kact;
            
        end
        
        set(hOb, 'Value',hsim.kmNo);
        set(hsim.text_k,'Visible','on','String',sprintf('Iteration %i of %i',hsim.kpNo-1,hsim.sim.Kact-1))
        
end




hsim = Plotting(hsim);
end

function slider_B1_Callback(hOb, ed)
global hsim
hsim.nB1inh = round(get(hOb, 'Value'));
set(hOb, 'Value',hsim.nB1inh);
set(hsim.text_B1,'Visible','on','String',sprintf('B1 scale %i of %i (%1.1e %%)',hsim.nB1inh,hsim.sim.B1inhom_N,hsim.sim.B1inhom_scales(hsim.nB1inh)))

hsim = Plotting(hsim);
end

function slider_B0_Callback(hOb, ed)
global hsim
hsim.nB0inh = round(get(hOb, 'Value'));
set(hOb, 'Value',hsim.nB0inh);
set(hsim.text_B0,'Visible','on','String',sprintf('B0 offset %i of %i (%1.1e Hz)',hsim.nB0inh,hsim.sim.B0inhom_N,hsim.sim.B0inhom_offsets(hsim.nB0inh)))

hsim = Plotting(hsim);
end

function List = Populate_Listbox(Struct)
global hsim
names = fieldnames(Struct);
vals = cell(length(names),1);
for n = 1:length(names)
    
    if ischar(getfield(Struct,names{n}))
        vals{n} = getfield(Struct,names{n});
    elseif isnumeric(getfield(Struct,names{n}))
        [A,B,C,D,E] = size(getfield(Struct,names{n}));
        % 		[A,B,C,D,E]
        
        if E == 1
            
            if D == 1 && C == 1 && B == 1 && A == 1
                vals{n} = num2str(getfield(Struct,names{n}));
            elseif D == 1 && C == 1 && B ~= 1 && A == 1
                if B < 6
                    temp = getfield(Struct,names{n});
                    
                    Str = '[';
                    for b = 1:B
                        
                        Str = [Str,num2str(temp(b))];
                        if b<B
                            Str = [Str,sprintf('\t\t,\t\t')];
                        end
                    end
                    Str = [Str,']'];
                    vals{n} = Str;
                else
                    
                    vals{n} = sprintf('<%ix%i>',A,B);
                end
            elseif D == 1 && C == 1 && B ~= 1 && A ~= 1
                if B < 4 && A < 3
                    temp = getfield(Struct,names{n});
                    
                    Str = '[';
                    for a = 1:A
                        for b = 1:B
                            
                            Str = [Str,num2str(temp(a,b))];
                            if b<B
                                Str = [Str,sprintf(',')];
                            end
                        end
                        if a ~=  A
                            Str = [Str,sprintf(';')];
                            
                        end
                    end
                    Str = [Str,']'];
                    vals{n} = Str;
                else
                    
                    vals{n} = sprintf('<%ix%i>',A,B);
                end
            elseif D == 1 && C == 1 && B == 1 && A ~= 1
                if A < 3
                    temp = getfield(Struct,names{n});
                    
                    Str = '[';
                    for a = 1:A
                        for b = 1:B
                            
                            Str = [Str,num2str(temp(a,b))];
                            if b<B
                                Str = [Str,sprintf(',')];
                            end
                        end
                        if a ~=  A
                            Str = [Str,sprintf(';')];
                            
                        end
                    end
                    Str = [Str,']'];
                    vals{n} = Str;
                else
                    
                    vals{n} = sprintf('<%ix%i>',A,B);
                end
            else
                
                vals{n} = sprintf('<%ix%ix%ix%i>',A,B,C,D);
            end
            
        else
            vals{n} = sprintf('<%ix%ix%ix%ix%i>',A,B,C,D,E);
        end
        
        
    elseif iscell(getfield(Struct,names{n}))
        [A,B,C,D,E] = size(getfield(Struct,names{n}));
        if E == 1
            
            if D == 1 && C == 1 && B == 1 && A == 1
                vals{n} = num2str(getfield(Struct,names{n}));
            elseif D == 1 && C == 1 && B ~= 1 && A == 1
                
                vals{n} = sprintf('<%ix%i cell>',A,B);
                
            elseif D == 1 && C == 1 && B ~= 1 && A ~= 1
                
                vals{n} = sprintf('<%ix%i cell>',A,B);
                
            elseif D == 1 && C == 1 && B == 1 && A ~= 1
                
                vals{n} = sprintf('<%ix%i cell>',A,B);
                
            else
                
                vals{n} = sprintf('<%ix%ix%ix%i cell>',A,B,C,D);
            end
            
        else
            vals{n} = sprintf('<%ix%ix%ix%ix%i cell>',A,B,C,D,E);
        end
        
    elseif isstruct(getfield(Struct,names{n}))
        [A,B] = size(getfield(Struct,names{n}));
        vals{n} = sprintf('<%ix%i struct>',A,B);
    end
end

List = cell(length(names),1);

for n = 1:length(names)
    List{n} = sprintf('%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s',names{n},vals{n});
    
end
end

function [Jet,Gray] = GetColormaps

jett = colormap(jet(512));

Jet = zeros(512,10,3);
Jet(:,:,1) = repmat(jett(:,1),[1,10,1]);
Jet(:,:,2) = repmat(jett(:,2),[1,10,1]);
Jet(:,:,3) = repmat(jett(:,3),[1,10,1]);


grayy = colormap(gray(512));

Gray = zeros(512,10,3);
Gray(:,:,1) = repmat(grayy(:,1),[1,10,1]);
Gray(:,:,2) = repmat(grayy(:,2),[1,10,1]);
Gray(:,:,3) = repmat(grayy(:,3),[1,10,1]);
end

function pushbutton_playx1_Callback(hOb, ed)
global hsim
for n = 1:hsim.sim.N
    hsim.nNo = n;
    set(hsim.slider_n,'Value',hsim.nNo)
    set(hsim.text_n,'String',sprintf('Timeframe %i dt (%1.1e s) of %i dt (%1.1e s)',hsim.nNo,hsim.nNo*hsim.sim.dt,hsim.sim.N,hsim.sim.N.*hsim.sim.dt))
    hsim = Plotting(hsim);
    pause(0.05)
    drawnow
end
end

function pushbutton_playx10_Callback(hOb, ed)
global hsim

n_ = round(linspace(1,hsim.sim.N,hsim.sim.N/10));

for n = n_;
    hsim.nNo = n;
    set(hsim.slider_n,'Value',hsim.nNo)
    set(hsim.text_n,'String',sprintf('Timeframe %i dt (%1.1e s) of %i dt (%1.1e s)',hsim.nNo,hsim.nNo*hsim.sim.dt,hsim.sim.N,hsim.sim.N.*hsim.sim.dt))
    hsim = Plotting(hsim);
    pause(0.05)
    drawnow
end
end

function pushbutton_playx100_Callback(hOb, ed)
global hsim

n_ = round(linspace(1,hsim.sim.N,hsim.sim.N/100));

for n = n_;
    hsim.nNo = n;
    set(hsim.slider_n,'Value',hsim.nNo)
    set(hsim.text_n,'String',sprintf('Timeframe %i dt (%1.1e s) of %i dt (%1.1e s)',hsim.nNo,hsim.nNo*hsim.sim.dt,hsim.sim.N,hsim.sim.N.*hsim.sim.dt))
    hsim = Plotting(hsim);
    pause(0.05)
    drawnow
end
end

function hsim = what2plot(hsim)
% tic
M_t = blOCh__spc('list2grid',[],[],[],[],hsim.sim.M_t(:,hsim.nNo,hsim.kmNo,hsim.pTxNo,hsim.nB1inh,hsim.nB0inh),[hsim.sim.R,hsim.sim.Rv],hsim.sim.idxtot,3);
Md = blOCh__spc('list2grid',[],[],[],[],hsim.sim.Md,[hsim.sim.R,hsim.sim.Rv],hsim.sim.idxtot,3);
Mx_max = max(max(hsim.sim.M_t(1:3:end,:)));
Mx_min = min(min(hsim.sim.M_t(1:3:end,:)));
My_max = max(max(hsim.sim.M_t(2:3:end,:)));
My_min = min(min(hsim.sim.M_t(2:3:end,:)));
Mz_max = max(max(hsim.sim.M_t(3:3:end,:)));
Mz_min = min(min(hsim.sim.M_t(3:3:end,:)));
Mxy_max = max(max(abs(complex(hsim.sim.M_t(1:3:end,:),hsim.sim.M_t(2:3:end,:)))));
Mxy_min = 0;

Mdx_max = max(max(hsim.sim.Md(1:3:end,:)));
Mdx_min = min(min(hsim.sim.Md(1:3:end,:)));
Mdy_max = max(max(hsim.sim.Md(2:3:end,:)));
Mdy_min = min(min(hsim.sim.Md(2:3:end,:)));
Mdz_max = max(max(hsim.sim.Md(3:3:end,:)));
Mdz_min = min(min(hsim.sim.Md(3:3:end,:)));
Mdxy_max = max(max(abs(complex(hsim.sim.Md(1:3:end,:),hsim.sim.Md(2:3:end,:)))));
Mdxy_min = 0;

% toc
switch get(hsim.popupmenu_Magn,'Value')
    
    case 1
        
        
        
        ordinate = abs(complex(M_t(:,:,:,:,1),M_t(:,:,:,:,2)));
        
        
        hsim.axm.ordinate_mn = Mxy_min;
        hsim.axm.ordinate_mx = Mxy_max;
        %       'Desired Magnetization, |Mxy| [M0]';
    case 2
        ordinate = angle(complex(M_t(:,:,:,:,1),M_t(:,:,:,:,2)));
        
        
        hsim.axm.ordinate_mn = -pi;
        hsim.axm.ordinate_mx = pi;
        %       'Desired Magnetization, arg(Mxy) [M0]';
    case 3
        %       'Desired Magnetization, Mx [M0]';
        ordinate = M_t(:,:,:,:,1);
        
        hsim.axm.ordinate_mn = Mx_min;
        hsim.axm.ordinate_mx = Mx_max;
    case 4
        %       'Desired Magnetization, My [M0]';
        ordinate = M_t(:,:,:,:,2);
        hsim.axm.ordinate_mn = My_min;
        hsim.axm.ordinate_mx = My_max;
    case 5
        %       'Desired Magnetization, Mz [M0]';
        ordinate = M_t(:,:,:,:,3);
        hsim.axm.ordinate_mn = Mz_min;
        hsim.axm.ordinate_mx = Mz_max;
    case 6
        
        
        ordinate = abs(complex(Md(:,:,:,:,1),Md(:,:,:,:,2)));
        
        
        hsim.axm.ordinate_mn = Mdxy_min;
        hsim.axm.ordinate_mx = Mdxy_max;
        %       'Desired Magnetization, |Mxy| [M0]';
    case 7
        %       'Desired Magnetization, Mx [M0]';
        ordinate = Md(:,:,:,:,1);
        
        hsim.axm.ordinate_mn = Mdx_min;
        hsim.axm.ordinate_mx = Mdx_max;
    case 8
        %       'Desired Magnetization, My [M0]';
        ordinate = Md(:,:,:,:,2);
        hsim.axm.ordinate_mn = Mdy_min;
        hsim.axm.ordinate_mx = Mdy_max;
    case 9
        %       'Desired Magnetization, Mz [M0]';
        ordinate = Md(:,:,:,:,3);
        hsim.axm.ordinate_mn = Mdz_min;
        hsim.axm.ordinate_mx = Mdz_max;
end
if hsim.axm.ordinate_mx == hsim.axm.ordinate_mn
    hsim.axm.ordinate_mx = hsim.axm.ordinate_mx + eps;
end
% hsim.axm.ordinate = ordinate;

switch hsim.sim.Dim
    
    case '1DSI'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,3),hsim.sim.D(2,3),5);
        
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(3),5);
        hsim.axm.Dim1label = ' z [m]';
        hsim.axm.ordinate = squeeze(ordinate);
        
        
        
        
    case '1DAP'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,2),hsim.sim.D(2,2),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(2),5);
        hsim.axm.Dim1label = ' y [m]';
        hsim.axm.ordinate = squeeze(ordinate);
    case '1DRL'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,1),hsim.sim.D(2,1),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(1),5);
        hsim.axm.Dim1label = ' x [m]';
        hsim.axm.ordinate = squeeze(ordinate);
    case '1+1DSI'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,3),hsim.sim.D(2,3),5);
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.Dv(1),hsim.sim.Dv(2),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(3),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.Rv,5);
        hsim.axm.Dim1label = ' z [m]';
        hsim.axm.Dim2label = ' v [Hz]';
        hsim.axm.ordinate = squeeze(ordinate);
    case '1+1DAP'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,2),hsim.sim.D(2,2),5);
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.Dv(1),hsim.sim.Dv(2),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(2),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.Rv,5);
        hsim.axm.Dim1label = ' y [m]';
        hsim.axm.Dim2label = ' v [Hz]';
        hsim.axm.ordinate = squeeze(ordinate);
    case '1+1DRL'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,1),hsim.sim.D(2,1),5);
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.Dv(1),hsim.sim.Dv(2),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(1),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.Rv,5);
        hsim.axm.Dim1label = ' x [m]';
        hsim.axm.Dim2label = ' v [Hz]';
        hsim.axm.ordinate = squeeze(ordinate);
    case '2DAx'
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.D(1,1),hsim.sim.D(2,1),5);
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,2),hsim.sim.D(2,2),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.R(1),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(2),5);
        hsim.axm.Dim1label = ' y [m]';
        hsim.axm.Dim2label = ' x [m]';
        hsim.axm.ordinate = squeeze(ordinate);
    case '2DCo'
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.D(1,1),hsim.sim.D(2,1),5);
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,3),hsim.sim.D(2,3),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.R(1),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(3),5);
        
        hsim.axm.Dim2label = ' x [m]';
        hsim.axm.Dim1label = ' z [m]';
        hsim.axm.ordinate = squeeze(ordinate);
    case '2DSa'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,2),hsim.sim.D(2,2),5);
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.D(1,3),hsim.sim.D(2,3),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(2),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.R(3),5);
        hsim.axm.Dim2label = ' y [m]';
        hsim.axm.Dim1label = ' z [m]';
        hsim.axm.ordinate = squeeze(ordinate);
    case '2+1DAx'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,1),hsim.sim.D(2,1),5);
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.D(1,2),hsim.sim.D(2,2),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(1),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.R(2),5);
        hsim.axm.Dim2label = ' x [m]';
        hsim.axm.Dim1label = ' y [m]';
        
        hsim.axm.ordinate = squeeze(ordinate(:,:,:,hsim.freqNo));
    case '2+1DCo'
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.D(1,1),hsim.sim.D(2,1),5);
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,3),hsim.sim.D(2,3),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.R(1),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(3),5);
        
        hsim.axm.Dim2label = ' x [m]';
        hsim.axm.Dim1label = ' z [m]';
        hsim.axm.ordinate = squeeze(ordinate(:,:,:,hsim.freqNo));
    case '2+1DSa'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,2),hsim.sim.D(2,2),5);
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.D(1,3),hsim.sim.D(2,3),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(2),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.R(3),5);
        hsim.axm.Dim2label = ' y [m]';
        hsim.axm.Dim1label = ' z [m]';
        hsim.axm.ordinate = squeeze(ordinate(:,:,:,hsim.freqNo));
    case '3D'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,1),hsim.sim.D(2,1),5);
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.D(1,2),hsim.sim.D(2,2),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(1),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.R(2),5);
        hsim.axm.Dim2label = ' x [m]';
        hsim.axm.Dim1label = ' y [m]';
        hsim.axm.ordinate = squeeze(ordinate(:,:,hsim.SlNo,1));
    case '3+1D'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.D(1,1),hsim.sim.D(2,1),5);
        hsim.axm.Dim2ticklabel = linspace(hsim.sim.D(1,2),hsim.sim.D(2,2),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.R(1),5);
        hsim.axm.Dim2axis = linspace(1,hsim.sim.R(2),5);
        hsim.axm.Dim2label = ' x [m]';
        hsim.axm.Dim1label = ' y [m]';
        hsim.axm.ordinate = squeeze(ordinate(:,:,hsim.SlNo,hsim.freqNo));
    case '0+1D'
        hsim.axm.Dim1ticklabel = linspace(hsim.sim.Dv(1),hsim.sim.Dv(2),5);
        hsim.axm.Dim1axis = linspace(1,hsim.sim.Rv,5);
        hsim.axm.Dim1label = ' v [Hz]';
        hsim.axm.ordinate = squeeze(ordinate);
end
switch get(hsim.popupmenu_Pulse,'Value')
    
    case 1
        if hsim.pTxNo == hsim.sim.pTx+1
            temp = winter(hsim.sim.pTx);
            hsim.axp.colorpulse = temp;
        else
            temp = winter(hsim.sim.pTx);
            hsim.axp.colorpulse = temp(hsim.pTxNo,:);
        end
        if hsim.pTxNo == hsim.sim.pTx+1
            hsim.axp.ordinate = abs(complex(hsim.sim.uo(:,:,hsim.kpNo),hsim.sim.vo(:,:,hsim.kpNo))).';
            
        else
            hsim.axp.ordinate = abs(complex(hsim.sim.uo(hsim.pTxNo,:,hsim.kpNo),hsim.sim.vo(hsim.pTxNo,:,hsim.kpNo))).';
            
        end
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(0,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(1,Mx,3);
        hsim.axp.Dim2label = ' |\omega| [rad/s]';
        hsim.axp.Dim1label = ' Time frame [#]';
        
    case 2
        if hsim.pTxNo == hsim.sim.pTx+1
            temp = winter(hsim.sim.pTx);
            hsim.axp.colorpulse = temp;
        else
            temp = winter(hsim.sim.pTx);
            hsim.axp.colorpulse = temp(hsim.pTxNo,:);
        end
        if hsim.pTxNo == hsim.sim.pTx+1
            hsim.axp.ordinate = angle(complex(hsim.sim.uo(:,:,hsim.kpNo),hsim.sim.vo(:,:,hsim.kpNo))).';
        else
            hsim.axp.ordinate = angle(complex(hsim.sim.uo(hsim.pTxNo,:,hsim.kpNo),hsim.sim.vo(hsim.pTxNo,:,hsim.kpNo))).';
        end
        
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = {'-pi','0','pi'};
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(-pi,pi,3);
        hsim.axp.Dim2label = ' \angle\omega [rad]';
        hsim.axp.Dim1label = ' Time frame [#]';
        
        
    case 3
        if hsim.pTxNo == hsim.sim.pTx+1
            temp = winter(hsim.sim.pTx);
            hsim.axp.colorpulse = temp;
        else
            temp = winter(hsim.sim.pTx);
            hsim.axp.colorpulse = temp(hsim.pTxNo,:);
        end
        if hsim.pTxNo == hsim.sim.pTx+1
            hsim.axp.ordinate = real(complex(hsim.sim.uo(:,:,hsim.kpNo),hsim.sim.vo(:,:,hsim.kpNo))).';
        else
            hsim.axp.ordinate = real(complex(hsim.sim.uo(hsim.pTxNo,:,hsim.kpNo),hsim.sim.vo(hsim.pTxNo,:,hsim.kpNo))).';
        end
        Mn = min(hsim.axp.ordinate(:))-eps;
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(Mn,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(Mn,Mx,3);
        hsim.axp.Dim2label = ' Re(\omega) [rad/s]';
        hsim.axp.Dim1label = ' Time frame [#]';
        
        
    case 4
        if hsim.pTxNo == hsim.sim.pTx+1
            temp = winter(hsim.sim.pTx);
            hsim.axp.colorpulse = temp;
        else
            temp = winter(hsim.sim.pTx);
            hsim.axp.colorpulse = temp(hsim.pTxNo,:);
        end
        if hsim.pTxNo == hsim.sim.pTx+1
            hsim.axp.ordinate = imag(complex(hsim.sim.uo(:,:,hsim.kpNo),hsim.sim.vo(:,:,hsim.kpNo))).';
        else
            hsim.axp.ordinate = imag(complex(hsim.sim.uo(hsim.pTxNo,:,hsim.pmNo),hsim.sim.vo(hsim.pTxNo,:,hsim.kpNo))).';
        end
        Mn = min(hsim.axp.ordinate(:))-eps;
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(Mn,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(Mn,Mx,3);
        hsim.axp.Dim2label = ' Im(\omega) [rad/s]';
        hsim.axp.Dim1label = ' Time frame [#]';
    case 5
        temp = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        hsim.axp.colorpulse = temp(1,:);
        hsim.axp.ordinate = 1e3.*hsim.sim.g(1,:).';
        Mn = min(hsim.axp.ordinate(:))-eps;
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(Mn,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(Mn,Mx,3);
        hsim.axp.Dim2label = ' G_x [mT/m]';
        hsim.axp.Dim1label = ' Time frame [#]';
    case 6
        temp = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        hsim.axp.colorpulse = temp(2,:);
        hsim.axp.ordinate = 1e3.*hsim.sim.g(2,:).';
        Mn = min(hsim.axp.ordinate(:))-eps;
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(Mn,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(Mn,Mx,3);
        hsim.axp.Dim2label = ' G_y [mT/m]';
        hsim.axp.Dim1label = ' Time frame [#]';
    case 7
        temp = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        hsim.axp.colorpulse = temp(3,:);
        hsim.axp.ordinate = 1e3.*hsim.sim.g(3,:).';
        Mn = min(hsim.axp.ordinate(:))-eps;
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(Mn,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(Mn,Mx,3);
        hsim.axp.Dim2label = ' G_z [mT/m]';
        hsim.axp.Dim1label = ' Time frame [#]';
    case 8
        tempc = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        
        hsim.axp.colorpulse = tempc;
        temp = [hsim.sim.g(1,:).',hsim.sim.g(2,:).',hsim.sim.g(3,:).'];
        temp = temp./max(abs(temp(:)));
        temp(:,2) = temp(:,2)+2.2;
        temp(:,1) = temp(:,1)+4.2;
        hsim.axp.ordinate = temp;
        %         Mn = min(hsim.axp.ordinate(:))-eps;
        %         Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = {'G_z','G_y','G_x'};
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = [0,2.2,4.2];
        hsim.axp.Dim2label = ' [arb.unit.]';
        hsim.axp.Dim1label = ' Time frame [#]';
    case 9
        if hsim.pTxNo == hsim.sim.pTx+1
            temp = winter(hsim.sim.pTx);
        else
            temp = winter(hsim.sim.pTx);
            temp=  temp(hsim.pTxNo,:);
        end
        temp2 = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        hsim.axp.colorpulse = [temp2;temp];
        
        if hsim.pTxNo == hsim.sim.pTx+1
            temprf = [abs(complex(hsim.sim.uo(:,:,hsim.kpNo),hsim.sim.vo(:,:,hsim.kpNo))).'./max(max(abs(complex(hsim.sim.uo(:,:,hsim.kpNo),hsim.sim.vo(:,:,hsim.kpNo)))))];
            
        else
            temprf = [abs(complex(hsim.sim.uo(hsim.pTxNo,:,hsim.kpNo),hsim.sim.vo(hsim.pTxNo,:,hsim.kpNo))).'./max(max(abs(complex(hsim.sim.uo(hsim.pTxNo,:,hsim.kpNo),hsim.sim.vo(hsim.pTxNo,:,hsim.kpNo)))))];
            
        end
        
        tempg = [hsim.sim.g(1,:).',hsim.sim.g(2,:).',hsim.sim.g(3,:).']./max(abs(hsim.sim.g(:)));
        tempg(:,2) = tempg(:,2)+2.2;
        tempg(:,1) = tempg(:,1)+4.2;
        
        
        hsim.axp.ordinate = [tempg,temprf+5.2];
        
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = {'G_z','G_y','G_x','|\omega|'};
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = [0,2.2,4.2,5.2];
        hsim.axp.Dim2label = ' [arb.unit.]';
        hsim.axp.Dim1label = ' Time frame [#]';
        
        
    case 10
        temp = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        hsim.axp.colorpulse = temp(1,:);
        hsim.axp.ordinate = cumsum(hsim.sim.g(1,:).').*hsim.sim.dt*hsim.sim.gamma/2/pi - sum(hsim.sim.g(1,:).').*hsim.sim.dt*hsim.sim.gamma/2/pi;
        
        Mn = min(hsim.axp.ordinate(:))-eps;
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(Mn,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(Mn,Mx,3);
        hsim.axp.Dim2label = ' k_x [1/m]';
        hsim.axp.Dim1label = ' Time frame [#]';
        
        
    case 11
        temp = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        hsim.axp.colorpulse = temp(2,:);
        hsim.axp.ordinate = cumsum(hsim.sim.g(2,:).').*hsim.sim.dt*hsim.sim.gamma/2/pi - sum(hsim.sim.g(2,:).').*hsim.sim.dt*hsim.sim.gamma/2/pi;
        
        Mn = min(hsim.axp.ordinate(:))-eps;
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(Mn,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(Mn,Mx,3);
        hsim.axp.Dim2label = ' k_y [1/m]';
        hsim.axp.Dim1label = ' Time frame [#]';
    case 12
        temp = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        hsim.axp.colorpulse = temp(3,:);
        hsim.axp.ordinate = cumsum(hsim.sim.g(3,:).').*hsim.sim.dt*hsim.sim.gamma/2/pi - sum(hsim.sim.g(3,:).').*hsim.sim.dt*hsim.sim.gamma/2/pi;
        
        Mn = min(hsim.axp.ordinate(:))-eps;
        Mx = max(hsim.axp.ordinate(:))+eps;
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = linspace(Mn,Mx,3);
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = linspace(Mn,Mx,3);
        hsim.axp.Dim2label = ' k_z [1/m]';
        hsim.axp.Dim1label = ' Time frame [#]';
        
    case 13
        temp = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        
        hsim.axp.colorpulse = temp;
        %         temp = [hsim.sim.g(1,:).',hsim.sim.g(2,:).',hsim.sim.g(3,:).']./max(abs(hsim.sim.g(:)));
        temp = cumsum([hsim.sim.g(1,:).',hsim.sim.g(2,:).',hsim.sim.g(3,:).']);
        temp = temp./max(abs(temp(:)));
        temp(:,2) = temp(:,2)+2.2;
        temp(:,1) = temp(:,1)+4.2;
        hsim.axp.ordinate = temp;
        
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = {'k_z','k_y','k_x'};
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = [0,2.2,4.2];
        hsim.axp.Dim2label = ' [arb.unit.]';
        hsim.axp.Dim1label = ' Time frame [#]';
        
    case 14
        if hsim.pTxNo == hsim.sim.pTx+1
            temp = winter(hsim.sim.pTx);
        else
            temp = winter(hsim.sim.pTx);
            temp=  temp(hsim.pTxNo,:);
        end
        temp2 = [0,0,0;0.2,0.2,0.2;0.4,0.4,0.4];
        hsim.axp.colorpulse = [temp2;temp];
        
        if hsim.pTxNo == hsim.sim.pTx+1
            temprf = [abs(complex(hsim.sim.uo(:,:,hsim.kpNo),hsim.sim.vo(:,:,hsim.kpNo))).'./max(max(abs(complex(hsim.sim.uo(:,:,hsim.kpNo),hsim.sim.vo(:,:,hsim.kpNo)))))];
            
        else
            temprf = [abs(complex(hsim.sim.uo(hsim.pTxNo,:,hsim.kpNo),hsim.sim.vo(hsim.pTxNo,:,hsim.kpNo))).'./max(max(abs(complex(hsim.sim.uo(hsim.pTxNo,:,hsim.kpNo),hsim.sim.vo(hsim.pTxNo,:,hsim.kpNo)))))];
            
        end
        
        tempg = cumsum([hsim.sim.g(1,:).',hsim.sim.g(2,:).',hsim.sim.g(3,:).'])-repmat(sum(hsim.sim.g,2).',[hsim.sim.N,1]);
        tempg = tempg./max(abs(tempg(:)));
        tempg(:,2) = tempg(:,2)+2.2;
        tempg(:,1) = tempg(:,1)+4.2;
        
        
        hsim.axp.ordinate = [tempg,temprf+5.2];
        
        hsim.axp.Dim1ticklabel = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2ticklabel = {'k_z','k_y','k_x','|\omega|'};
        hsim.axp.Dim1axis = linspace(1,hsim.sim.N,5);
        hsim.axp.Dim2axis = [0,2.2,4.2,5.2];
        hsim.axp.Dim2label = ' [arb.unit.]';
        hsim.axp.Dim1label = ' Time frame [#]';
end


end

function hsim=Plotting(hsim)

hsim = what2plot(hsim);



switch hsim.sim.Dim
    case '1DSI'
        axes(hsim.axes_Magn)
        
        plot(hsim.axm.ordinate,'linewidth',2)
        axis([1,hsim.sim.R(3),hsim.axm.ordinate_mn,hsim.axm.ordinate_mx])
        set(hsim.axes_Magn,'YTick',linspace(hsim.axm.ordinate_mn,hsim.axm.ordinate_mx,5))
        set(hsim.axes_Magn,'XTickLabel',hsim.axm.Dim1ticklabel)
        xlabel(hsim.axm.Dim1label)
    case '1DAP'
        axes(hsim.axes_Magn)
        
        plot(hsim.axm.ordinate,'linewidth',2)
        axis([1,hsim.sim.R(2),hsim.axm.ordinate_mn,hsim.axm.ordinate_mx])
        set(hsim.axes_Magn,'YTick',linspace(hsim.axm.ordinate_mn,hsim.axm.ordinate_mx,5))
        set(hsim.axes_Magn,'XTickLabel',hsim.axm.Dim1ticklabel)
        xlabel(hsim.axm.Dim1label)
        
    case '1DRL'
        axes(hsim.axes_Magn)
        
        plot(hsim.axm.ordinate,'linewidth',2)
        axis([1,hsim.sim.R(1),hsim.axm.ordinate_mn,hsim.axm.ordinate_mx])
        set(hsim.axes_Magn,'YTick',linspace(hsim.axm.ordinate_mn,hsim.axm.ordinate_mx,5))
        set(hsim.axes_Magn,'XTickLabel',hsim.axm.Dim1ticklabel)
        xlabel(hsim.axm.Dim1label)
        
        
    case {'1+1DSI','1+1DRL','1+1DAP','2DAx','2DCo','2DSa','2+1DAx','2+1DCo','2+1DSa','3D','3+1D'}
        axes(hsim.axes_Magn)
        imagesc(hsim.axm.ordinate,[hsim.axm.ordinate_mn,hsim.axm.ordinate_mx])
        set(hsim.axes_Magn,'XTick',hsim.axm.Dim2axis,'YTick',hsim.axm.Dim1axis)
        set(hsim.axes_Magn,'XTickLabel',hsim.axm.Dim2ticklabel,'YTickLabel',hsim.axm.Dim1ticklabel)
        set(hsim.axes_Magn,'Ydir','normal')
        xlabel(hsim.axm.Dim2label)
        ylabel(hsim.axm.Dim1label)
        axes(hsim.axes_Colorbar)
        imagesc(permute(hsim.Colmap,[1,2,3])), axis off
        set(hsim.text1,'Visible','on','string',sprintf('%1.10f',hsim.axm.ordinate_mn))
        set(hsim.text2,'Visible','on','string',sprintf('%1.10f',hsim.axm.ordinate_mx))
end

axes(hsim.axes_Pulse)
for n = 1:size(hsim.axp.ordinate,2)
    stairs(hsim.axp.ordinate(:,n),'color',hsim.axp.colorpulse(n,:),'linewidth',2)
    hold on
end
hold on
% val = max(abs(,abs(max(hsim.axp.ordinate(:))))

switch get(hsim.popupmenu_Pulse,'Value')
    case {8,13}
        line([1,hsim.nNo],[0,0],'linewidth',1,'color','black')
        line([1,hsim.nNo],[2.2,2.2],'linewidth',1,'color','black')
        line([1,hsim.nNo],[4.2,4.2],'linewidth',1,'color','black')
    case {9,14}
        line([1,hsim.nNo],[0,0],'linewidth',1,'color','black')
        line([1,hsim.nNo],[2.2,2.2],'linewidth',1,'color','black')
        line([1,hsim.nNo],[4.2,4.2],'linewidth',1,'color','black')
        line([1,hsim.nNo],[5.2,5.2],'linewidth',1,'color','black')
end

line([hsim.nNo,hsim.nNo],[min(hsim.axp.ordinate(:))-0.1*abs(min(hsim.axp.ordinate(:))),max(hsim.axp.ordinate(:))+0.1*abs(max(hsim.axp.ordinate(:)))],'linewidth',1,'color','red')
hold off
axis tight

set(hsim.axes_Pulse,'XTick',(hsim.axp.Dim1axis),'YTick',hsim.axp.Dim2axis)
set(hsim.axes_Pulse,'XTickLabel',round(hsim.axp.Dim1ticklabel),'YTickLabel',hsim.axp.Dim2ticklabel)
xlabel(hsim.axp.Dim1label)
ylabel(hsim.axp.Dim2label)


end



