function varargout = blOCh__opt(Fun,spc,khr,Method,Par,varargin)
% function varargout = blOCh__opt(Fun,spc,khr,Method,Par,varargin)
%
%   This blOCh script conducts the optimization
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



if isempty(Fun)
    DEBUG = 0;
    
    if DEBUG
        ccc
        opt = [];
        Msg = [];
        
        Par.dt = 1;
        Par.lambda = 1e-6;
        Par.epsilon = 2e-2;
        Par.Grad = '1st';
        Par.Grad = 'Schirmer';
        
        Par = {Par,Par};
        spc = blOCh__spc([],'2DAx','','','','pTx',1,'Show',1,'L',[0.192,0.192],'R',[5,5],'v0fac',1,'Mdflip',[10,0]);
        
        Park.T = 1;
        tic
        khr = blOCh__khr([],spc,'2Dvdspii',Park,'Show',1,'Sys','Sie7T');
        
        Method = {'GRAPE_Khaneja','GRAPE_Khaneja'};
        varargin = {'dt',10e-6,'Init','Random','Handover',1,'MaxIter',[10,10],};
        E = mfilename;
        def_Name = mfilename;
        def_Save2 = [mfilename,filesep];
        Caller = mfilename;
        Fun = [];
    else
        E = dbstack('-completenames');
        def_Name = E(2,1).name;
        def_Save2 = [E(2,1).file(1:end-2),filesep];
        Caller = E(2,1).file;
    end
    Now = now;
    %% Defining Defaults and Valids
    ValidMethods = {'GRAPE_Khaneja','MC_MadayTurinici','QN_LBFGS'};
    
    def_ksaveintermediate = 10;
    
    def_Show = 1;
    
    def_Save = '000000';
    
    def_Init = {'Random',NaN};
    
    def_RFipct = 10;
    
    def_RFmpct = 100;
    
    val_Handovers = [1,0];
    
    def_Handover = val_Handovers(2);
    
    val_HandoverWhat = {'last','best'};
    
    def_HandoverWhat = val_HandoverWhat{1};
    
    val_deriv_checks = {'off','on'};
    
    def_deriv_check = val_deriv_checks{1};
    
    val_Masks      = [0,1];
    
    def_Mask = val_Masks(2);
    
    def_par_Ncores = 1;
    
    def_par_Type = 1;
    
    valid_Shows = [1,0];
    def_Show = valid_Shows(1);
    
    
    
    
    %% Starting to read input
    p = inputParser;
    test = false;
    try p.addRequired('Fun', @(x)isempty(x));
        try p.addRequired('spc', @(x)isstruct(x));
            try p.addRequired('khr', @(x)isstruct(x)); %#ok<*ALIGN>
                try def_dt = khr.dt; % Here I should error check the attributes of dt
                    try def_NOC = khr.N; % Here I should error check the attributes of N
                        try p.addRequired('Method',@(x)Validate_Method(x,ValidMethods));
                            if iscell(Method)
                                if numel(Method) == 2
                                    opt.NumMethods = 2;
                                    
                                    def_TolCon = ones(1,2).*1e-6;
                                    def_TolX = ones(1,2).*1e-6;
                                    def_ObjLim = ones(1,2).*1e-10;
                                    def_TolFun = ones(1,2).*1e-6;
                                    def_MaxIter = ones(1,2).*5000;
                                    def_MaxFunEvals = ones(1,2).*Inf;
                                else
                                    opt.NumMethods = 1;
                                    
                                    def_TolCon = 1e-6;
                                    def_TolX = 1e-6;
                                    def_ObjLim = 1e-10;
                                    def_TolFun = 1e-6;
                                    def_MaxIter = 5000;
                                    def_MaxFunEvals = Inf;
                                end
                            else
                                opt.NumMethods = 1;
                                
                                def_TolCon = 1e-6;
                                def_TolX = 1e-6;
                                def_ObjLim = 1e-10;
                                def_TolFun = 1e-6;
                                def_MaxIter = 5000;
                                def_MaxFunEvals = Inf;
                            end
                            
                            try p.addRequired('Par',@(x)Validate_Par(x,opt.NumMethods));
                                
                                Display_Message('opt: There is currently no error checking of the content(s) of the Par input(s)',1);
                                
                                try p.addParamValue('dt',def_dt,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','finite','nonnegative'})); %#ok<*NVREPL>
                                    try p.addParamValue('NOC',def_NOC,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','integer','nonnegative'}));
                                        try p.addParamValue('Init',def_Init,@(x)Validate_Init(x));
                                            try p.addParamValue('Handover',def_Handover,@(x)any(val_Handovers));
                                                try p.addParamValue('HandoverWhat',def_HandoverWhat,@(x)any(strcmpi(x,val_HandoverWhats)));
                                                    try p.addParamValue('Show',def_Show,@(x)any(valid_Shows));
                                                        try p.addParamValue('RFipct',def_RFipct,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','finite','nonnegative'}));
                                                            try p.addParamValue('RFmpct',def_RFmpct,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','finite','nonnegative'}));
                                                                try p.addParamValue('Save2',def_Save2,@(x)validateattributes(x,{'char'},{'row'}));
                                                                    try p.addParamValue('Save',def_Save,@(x)Validate_Save(x));
                                                                        try p.addParamValue('Name',def_Name,@(x)validateattributes(x,{'char'},{'row'}));
                                                                            try p.addParamValue('ksaveintermediate',def_ksaveintermediate,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','integer','nonnegative','>',1}));
                                                                                try p.addParamValue('Mask',def_Mask,@(x)any(val_Masks));
                                                                                    try p.addParamValue('deriv_check',def_deriv_check,@(x)any(strcmpi(x,val_deriv_checks)));
                                                                                        try p.addParamValue('TolCon',def_TolCon,@(x)validateattributes(x,{'numeric'},{'size',[1,opt.NumMethods],'real','finite','nonnegative'}));
                                                                                            try p.addParamValue('TolFun',def_TolFun,@(x)validateattributes(x,{'numeric'},{'size',[1,opt.NumMethods],'real','finite','nonnegative'}));
                                                                                                try p.addParamValue('ObjLim',def_ObjLim,@(x)validateattributes(x,{'numeric'},{'size',[1,opt.NumMethods],'real','finite','nonnegative'}));
                                                                                                    try p.addParamValue('TolX',def_TolX,@(x)validateattributes(x,{'numeric'},{'size',[1,opt.NumMethods],'real','finite','nonnegative'}));
                                                                                                        try p.addParamValue('MaxFunEvals',def_MaxFunEvals,@(x)validateattributes(x,{'numeric'},{'size',[1,opt.NumMethods],'real','finite','nonnegative'}));
                                                                                                            try p.addParamValue('MaxIter',def_MaxIter,@(x)validateattributes(x,{'numeric'},{'size',[1,opt.NumMethods],'real','finite','nonnegative'}));
                                                                                                                try p.addParamValue('par_Ncores',def_par_Ncores,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','integer','finite'}));
                                                                                                                    try p.addParamValue('par_Type',def_par_Type,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','finite'}));
                                                                                                                        
                                                                                                                        
                                                                                                                        test= true;
                                                                                                                        
                                                                                                                        
                                                                                                                        
                                                                                                                    catch  me;Display_Message(['blOCh__opt: par_Type: ',me.message],2); end
                                                                                                                catch  me;Display_Message(['blOCh__opt: par_Ncores: ',me.message],2); end
                                                                                                                
                                                                                                            catch  me;Display_Message(['blOCh__opt: MaxIter: ',me.message],2); end
                                                                                                        catch  me;Display_Message(['blOCh__opt: MaxFunEvals: ',me.message],2); end
                                                                                                    catch  me;Display_Message(['blOCh__opt: TolX: ',me.message],2); end
                                                                                                catch  me;Display_Message(['blOCh__opt: ObjLim: ',me.message],2); end
                                                                                            catch  me;Display_Message(['blOCh__opt: TolFun: ',me.message],2); end
                                                                                        catch  me;Display_Message(['blOCh__opt: TolCon: ',me.message],2); end
                                                                                    catch  me;Display_Message(['blOCh__opt: deriv_check: ',me.message],2); end
                                                                                catch  me;Display_Message(['blOCh__opt: Mask: ',me.message],2); end
                                                                            catch  me;Display_Message(['blOCh__opt: ksaveintermediate: ',me.message],2); end
                                                                        catch  me;Display_Message(['blOCh__opt: Name: ',me.message],2); end
                                                                    catch  me;Display_Message(['blOCh__opt: Save: ',me.message],2); end
                                                                catch  me;Display_Message(['blOCh__opt: Save2:',me.message],2); end
                                                            catch  me;Display_Message(['blOCh__opt: RFmpct: ',me.message],2);end
                                                        catch  me;Display_Message(['blOCh__opt: RFipct: ',me.message],2);end
                                                    catch  me;Display_Message(['blOCh__opt: Show: ',me.message],2); end
                                                catch  me;Display_Message(['blOCh__opt: HandoverWhat: ',me.message],2); end
                                            catch  me;Display_Message(['blOCh__opt: Handover: ',me.message],2);end
                                            
                                        catch  me;Display_Message(['blOCh__opt: Init: ',me.message],2); end
                                    catch  me;Display_Message(['blOCh__opt: NOC: ',me.message],2); end
                                catch  me;Display_Message(['blOCh__opt: dt:',me.message],2); end
                            catch  me;Display_Message(['blOCh__opt: Par: ',me.message],2);end
                        catch  me;Display_Message(['blOCh__opt: Method: ',me.message],2);end
                    catch  me;Display_Message(['blOCh__opt: def_NOC: ',me.message],2); end
                catch  me;Display_Message(['blOCh__opt: def_dt: ',me.message],2); end
            catch  me;Display_Message(['blOCh__opt: khr: ',me.message],2);end
        catch  me;Display_Message(['blOCh__opt: spc: ',me.message],2);end
    catch  me;Display_Message(['blOCh__opt: Fun: ',me.message],2);end
    if test
        % try p.parse(spc,khr,Method,Par);
        test = false;
        try p.parse(Fun,spc,khr,Method,Par,varargin{:});
            test = true;
        catch  me; Display_Message(['blOCh__opt: ',me.message],2); end
    end
    
    
    
    %% Compiling input
    if test
    opt.TimeStamp               = [datestr(Now,'yymmdd'),'_at_',datestr(Now,'HHMMSS')];
    opt.Method                  = p.Results.Method;
    opt.Par                  = p.Results.Par;
    opt.dt                      = p.Results.dt;
    opt.N                       = p.Results.NOC;
    opt.ArgOrder                = varargin;
    [opt.dt,opt.N]            = Set_NOC_dt(opt.N,opt.dt,opt.ArgOrder,khr);
    
    
    opt.Handover                = p.Results.Handover;
    opt.HandoverWhat            = p.Results.HandoverWhat;
    
    if opt.Handover == 1
        if iscell(  opt.Method)
            if numel(opt.Method) == 1
                opt.Handover = 0;
            end
        else
            opt.Handover = 0;
        end
    end
    
    opt.Save                    = p.Results.Save;
    opt.Save2                   = p.Results.Save2;
    opt.Name                    = p.Results.Name;
    
    opt.par_Ncores              = p.Results.par_Ncores;
    opt.par_Type                = p.Results.par_Type;
    
    opt.MaxIter                 = p.Results.MaxIter;
    opt.TolX                    = p.Results.TolX;
    opt.ObjLim                   = p.Results.ObjLim;
    opt.TolFun                  = p.Results.TolFun;
    opt.TolCon                  = p.Results.TolCon;
    opt.MaxFunEvals             = p.Results.MaxFunEvals;
    
    opt.ksaveintermediate       = p.Results.ksaveintermediate;
    
    
    
    
    opt.Show                    = p.Results.Show;
    opt.Caller                  = Caller;
    opt.Mask                    = p.Results.Mask;
    opt.deriv_check             = p.Results.deriv_check;
    opt.Init                    = p.Results.Init;
    opt.RFipct                  = p.Results.RFipct;
    opt.RFmpct                  = p.Results.RFmpct;
    
    
    if ispc
        home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    else
        home = getenv('HOME');
        
        if ~isempty(strfind(opt.Save2,'~'))
            opt.Save2           = strrep(opt.Save2,'~',home);
        end
        
    end
    
    
    opt = Run_this(spc,khr,opt);
    
    varargout{1} = opt;
    else
        varargout{1} = -1;
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
            Display_Message(['blOCh__opt: Fun: This switch needs to be expanded'],2)
            
    end
    
    varargout = v;
    
end
end

%% Validation Functions
function test = Validate_Method(x,ValidMethods)
DEBUG = 0;
if DEBUG
    x = {'GRAPE_Khaneja'};
    x = {'GRAPE_Khaneja','GRAPE_Khaneja'};

    ValidMethods = {'GRAPE_Khaneja','MC_MadayTurinici'};
    
end


test = false;

if iscell(x)
    Nmethods = numel(x);
    if Nmethods > 2
        Display_Message('Validate_Method: Maximum number of Methods is 2',2);
        return
    elseif Nmethods == 2
        Display_Message('Validate_Method: The number of Methods is 2',1);
    elseif Nmethods == 1
        Display_Message('Validate_Method: 1 Method intended',1);
    else
        Display_Message('Validate_Method: Something wrong with the Method',2);
        return
    end
    Tests = zeros(1,Nmethods);
    
    for n = 1:Nmethods
        Tests(n) = Validate_Single_Method(x{n},ValidMethods);
    end
    if sum(Tests) == Nmethods
        test = true;
    end
else
    Display_Message('Validate_Method: 1 Method intended',1);
    test = Validate_Single_Method(x,ValidMethods);
end
end

function test = Validate_Single_Method(x,ValidMethods)
if ischar(x)
    if ismember(x,ValidMethods)
        test = true;
        Display_Message(sprintf('Validate_Method: Method "%s" is OK',x),1);
    else
        
        if exist(x,'file')
            test = true;
            Display_Message(sprintf('Validate_Single_Method: File for Method "%s" exists',x),1);
        else
            test = false;
            Display_Message(sprintf('Validate_Single_Method: Somthing is wrong with Method "%s"',x),2);
            
        end
    end
else
    Display_Message(sprintf('Validate_Single_Method: The Method must be a literal string'),2);
    test = false;
end

end

function test = Validate_Par(x,NumMethods)
DEBUG = 0;
if DEBUG
    x.dt= 1;

    NumMethods = 2;
end


test = false;

if iscell(x)
    Npars = numel(x);
    if Npars > 2
        Display_Message('Validate_Par: Maximum number of Par inputs is 2',2);
        
    elseif Npars == 2 && Npars == NumMethods
        Display_Message('Validate_Par: The number of Par inputs is 2',1);
    elseif Npars ~= NumMethods
        Display_Message('Validate_Par: The number of Par inputs should match the number Methods',2);
        return
        
    elseif Npars == 1
        Display_Message('Validate_Par: 1 Par input intended',1);
    else
        Display_Message('Validate_Par: Something wrong with the Par',2);
        return
    end
    Tests = zeros(1,Npars);
    
    for n = 1:Npars
        Tests(n) = Validate_Single_Par(x{n});
    end
    if sum(Tests) == Npars
        test = true;
    end
else
    Npars = 1;
    if Npars ~= NumMethods
        Display_Message('Validate_Par: The number of Par inputs should match the number Methods',2);
        return
    else
        Display_Message('Validate_Par: 1 Par input intended',1);
        test = Validate_Single_Par(x);
    end
end
end

function test = Validate_Single_Par(x)
if isstruct(x)
    test = true;
else
    Display_Message(sprintf('Validate_Single_Par: The Par input must be a struct'),2);
    test = false;
end

end

function test = Validate_Init(x)

test = false;
if iscell(x)
    
    if numel(x) == 2;
        x1 = x{1};
        x2 = x{2};
        
        if ischar(x1)
            
            switch x1
                
                case {'Random','0s'}
                    
                    test = true;
                    % Doesn't care what x2 might be...
                otherwise
                    if exist(x1,'file')
                        if isnumeric(x2)
                            if isnan(x2)
                                test = true;
                            elseif x2 > 0
                                if rem(x2,2) == 0
                                    test = true;
                                end
                            end
                        elseif ischar(x2)
                            switch x2
                                case 'opt1best'
                                    test = true;
                                case 'opt2best'
                                    test = true;
                            end
                        else
                            return
                        end
                    end
            end
            
        end
        
    elseif numel(x) == 1
        
        x = x{1};
        
        switch x
            
            case {'Random','0s'}
                
                test = true;
                
            otherwise
                
                if exist(x,'file')
                    test = true;
                end
                
        end
        
    end
elseif ischar(x)
    
    switch x
        
        case {'Random','0s'}
            
            test = true;
            
        otherwise
            
            if exist(x,'file')
                test = true;
            end
            
    end
    
    
    
end


if test
    Display_Message('Validate_Init: Passed',1);
    
end
end


function test = Validate_Save(x)

test = false;

nn = numel(x);
nvalid = 6;
if nn ~= nvalid
    error('Save must be %i characters',nvalid)
end

t1 = {'0','1','2'};
t2 = {'0','1','2'};
t3 = {'0','1','2','3','4','5'};
t4 = {'0','1'};
t5 = {'0','1'};
t6 = {'0','1','2'};


OK = zeros(1,nvalid);

for n = 1:nvalid
    
    t = eval(['t',num2str(n)]);
    fun = @(xx)any(strcmpi(xx,t));
    OK(n) = fun(x(n));
end
OK = sum(OK);

if OK == nvalid
    test = true;
end

end

%% Allocation Functions
function opt = Get_Initial_RF(spc,khr,opt)


if opt.OptNum == 1
    Init = opt.Init;
    if iscell(Init)
        
        
        if numel(Init) == 2
            
            I1 = Init{1};
            I2 = Init{2};
            
            switch I1
                case 'Random'
                    [opt.u,opt.v,opt.uo,opt.vo] = Get_Random_RF(spc.pTx,opt.Non,opt.mon,opt.N,opt.RFipct/100*khr.RFmHW,opt.MaxIter);
                    
                case '0s'
                    [opt.u,opt.v,opt.uo,opt.vo] = Get_0s_RF(spc.pTx,opt.Non,opt.mon,opt.N,opt.RFipct/100*khr.RFmHW,opt.MaxIter);
                    
                otherwise
                    
                    [opt.u,opt.v,opt.uo,opt.vo,Msg,Nfo] = Load_RF(I1,I2,opt,spc.pTx);
                    
            end
            
        elseif numel(Init) == 1
            I1 = Init{1};
            I2 = NaN;
            switch I1
                case 'Random'
                    [opt.u,opt.v,opt.uo,opt.vo] = Get_Random_RF(spc.pTx,opt.Non,opt.mon,opt.N,opt.RFipct/100*khr.RFmHW,opt.MaxIter);
                    
                case '0s'
                    [opt.u,opt.v,opt.uo,opt.vo] = Get_0s_RF(spc.pTx,opt.Non,opt.mon,opt.N,opt.RFipct/100*khr.RFmHW,opt.MaxIter);
                    
                otherwise
                    
                    [opt.u,opt.v,opt.uo,opt.vo,Msg,Nfo] = Load_RF(I1,I2,opt,spc.pTx);
                    
            end
        end
    elseif ischar(Init)
        
        switch Init
            
            case 'Random'
                
                
                
                [opt.u,opt.v,opt.uo,opt.vo] = Get_Random_RF(spc.pTx,opt.Non,opt.mon,opt.N,opt.RFipct/100*khr.RFmHW,opt.MaxIter);
                
                
            case '0s' % pTx,Non,mon,N,K
                [opt.u,opt.v,opt.uo,opt.vo] = Get_0s_RF(spc.pTx,opt.Non,opt.mon,opt.N,opt.MaxIter);
                
                
                
            otherwise
                [opt.u,opt.v,opt.uo,opt.vo,Msg,Nfo] = Load_RF(Init,NaN,opt,spc.pTx);
                
                
        end
    end
    
    
    
else
    N = opt.N;
    
    K = opt.MaxIter;
    opt.u = opt.HandoverControls1{1};
    opt.v = opt.HandoverControls1{2};
    opt.uo = zeros(spc.pTx,N,K+1);
    opt.vo = zeros(spc.pTx,N,K+1);
    opt.uo(:,:,1) = opt.u;
    opt.vo(:,:,1) = opt.v;
end
end

function [wx,wy,wxo,wyo] = Get_Random_RF(pTx,Non,mon,N,RFi,K)
rand('seed',sum(100*clock()));
wx_ = (rand(pTx,Non)-0.5)*2.*RFi;
rand('seed',sum(200*clock()));
wy_ = (rand(pTx,Non)-0.5)*2.*RFi;

wx = zeros(pTx,N);
wy = zeros(pTx,N);

wx(:,mon) = wx_;
wy(:,mon) = wy_;


wxo = zeros(pTx,N,K+1);
wyo = zeros(pTx,N,K+1);
wxo(:,:,1) = wx;
wyo(:,:,1) = wy;

end

function [wx,wy,wxo,wyo] = Get_0s_RF(pTx,Non,mon,N,K)

wx_ = zeros(pTx,Non);
wy_ = zeros(pTx,Non);

wx = zeros(pTx,N);
wy = zeros(pTx,N);

wx(:,mon) = wx_;
wy(:,mon) = wy_;


wxo = zeros(pTx,N,K+1);
wyo = zeros(pTx,N,K+1);
wxo(:,:,1) = wx;
wyo(:,:,1) = wy;

end

function [wx,wy,wxo,wyo,Msg,Nfo] = Load_RF(I1,I2,opt,pTx)
N = opt.N;
K = opt.MaxIter;
[pathstr, name, extension] = fileparts(I1);

Msg = [];

switch extension
    
    case {'bnd','.bnd'}
        
        old = load(I1,'-mat');
        
        if isfield(old,'opt')
            
            names = fieldnames(old.opt);
            
            if numel(names) == 1
                names = names{1};
                if isfield(eval(sprintf('old.opt.%s',names)),'uo') && isfield(eval(sprintf('old.opt.%s',names)),'vo')
                    
                    if ischar(I2)
                        
                        switch I2
                            case 'best'
                                [Mx,idx] = max(eval(sprintf('old.opt.%s.Fun(:)',names)));
                                
                                wx_ = eval(sprintf('old.opt.%s.uo(:,:,%i)',names,idx));
                                wy_ = eval(sprintf('old.opt.%s.vo(:,:,%i)',names,idx));
                        end
                        
                    elseif isnumeric(I2)
                        
                        if isnan(I2)
                            
                            [Mx,idx] = max(eval(sprintf('old.opt.%s.Fun(:)',names)));
                            
                            wx_ = eval(sprintf('old.opt.%s.uo(:,:,%i)',names,idx));
                            wy_ = eval(sprintf('old.opt.%s.vo(:,:,%i)',names,idx));
                            
                        else
                            [A,B,C] = size(eval(sprintf('old.opt.%s.uo',names)));
                            
                            if I2 > C,I2 = C;
                            end
                            wx_ = eval(sprintf('old.opt.%s.uo(:,:,%i)',names,I2));
                            wy_ = eval(sprintf('old.opt.%s.vo(:,:,%i)',names,I2));
                        end
                        
                    end
                    
                else
                    Msg = sprintf('Init bundle opt.%s doesn''t contain uo and vo',names);
                end
            else
                Msg = sprintf('Init opt struct contains multiple structs');
            end
            
        else
            Msg = sprintf('Init bundle misses an opt struct');
        end
        
    case {'txt','.txt'}
        
        old = load(Init,'-ascii');
        
        [A,B] = size(old);
        
        if B/2 ~= spc.pTx
            Msg = sprintf('For the chosen pTx, Init:\n\t%s\ndoesn''t contain enough columns (u(s=1),...,u(s=pTx),v(s=1),...,v(s=pTx))',Init);
        else
            temp1 = old(:,1:B/2);
            temp2 = old(:,B/2+1:end);
            
            wx_ = permute(temp1,[2,1]);
            wy_ = permute(temp2,[2,1]);
        end
        
        
        
    otherwise
        
end

wx = zeros(pTx,N);
wy = zeros(pTx,N);


if size(wx_,2) ~= N
    for s = 1:pTx
        wx(s,:) = interp1(linspace(1,N,size(wx_,2)),wx_(s,:),[1:N]);
    end
else
    for s = 1:pTx
        wx(s,:) = wx_(s,:);
    end
end
if size(wy_,2) ~= N
    for s = 1:pTx
        wy(s,:) = interp1(linspace(1,N,size(wy_,2)),wy_(s,:),[1:N]);
    end
else
    for s = 1:pTx
        wy(s,:) = wy_(s,:);
    end
end

wx(isnan(wx)) = 0;
wy(isnan(wx)) = 0;
wx(isinf(abs(wx))) = 0;
wy(isinf(abs(wx))) = 0;

wxo = zeros(pTx,N,K+1);
wyo = zeros(pTx,N,K+1);
wxo(:,:,1) = wx;
wyo(:,:,1) = wy;

Nfo = '';




end

function [dt,NOC] = Set_NOC_dt(NOC,dt,ArgOrder,khr)

idx_dt = find(strcmpi(ArgOrder,'dt'));
idx_no = find(strcmpi(ArgOrder,'NOC'));
if ~isempty(idx_dt) || ~isempty(idx_no)
    
    % Check NOC vs dt
    test1 = NOC ~= khr.N;
    test2 = dt ~= khr.dt;
    temp1 = 0;
    temp2 = 0;
    if numel(test1) > 1 || test1 == 1
        
        dt = khr.T./NOC;
        Display_Message('Set_NOC_dt: NOC was specified different from khr.N => dt adjusted accordingly',1);
        
        temp1 = 1;
    end
    if numel(test2) > 1 || test2 == 1
        
        NOC = round(khr.T./dt);
        Display_Message('Set_NOC_dt: dt was specified different from khr.dt => NOC adjusted accordingly',1);
        temp2 = 1;
    end
    
    if temp1 && temp2
        Display_Message('Set_NOC_dt: Both NOC and dt were specified different from khr. N and khr.dt => Both are adjusted accordingly and NOC overrules dt',1);
        dt = khr.T./NOC;
    end
end
end

function opt = Prepare_khr_4_opt(opt,khr)

temp = linspace(1,khr.N,khr.N);
temp2 = linspace(1,khr.N,opt.N);

if opt.OptNum == 1
    
    opt.g = zeros(3,opt.N);
    
    
    
    for n = 1:3
        opt.g(n,:) = interp1(temp,khr.g(n,:),temp2,'linear');
    end
    % opt
    opt.g(isnan(opt.g)) = 0;
    opt.g(isinf(opt.g)) = 0;
    
    
    
elseif opt.OptNum == 2
    
    for n = 1:3
        opt.g(n,:) = opt.HandoverControls1{n+2};
    end
    
    
end
if opt.Mask
    opt.m = logical(interp1(temp,double(khr.m),temp2,'nearest'));
else
    opt.m = true(1,opt.N);
end
opt.Non = sum(opt.m);
opt.Noff = sum(~opt.m);
opt.mon = find(opt.m);
opt.moff = find(~opt.m);


end

function opt = Allocate_Variables(spc,khr,opt)

if opt.OptNum >= 1
    
    if ~isreal(spc.B1map)
        opt.sr = real(spc.B1map);
        opt.si = imag(spc.B1map);
    else
        opt.sr = spc.B1map;
        opt.si = zeros(size(spc.B1map));
    end
    opt.Durations = zeros(opt.MaxIter+1,1);
    opt.Fun = zeros(opt.MaxIter+1,1);
    opt.Eff = zeros(opt.MaxIter+1,1);
    opt.Pen = zeros(opt.MaxIter+1,1);
    opt.dFun = nan(opt.MaxIter+1,1);
    
    
    opt.w0 = spc.w0map+spc.V*2*pi;
    opt.yx = spc.X.*khr.gamma;
    opt.yy = spc.Y.*khr.gamma;
    opt.yz = spc.Z.*khr.gamma;
    
    
    opt.Gm = khr.GmHW*khr.Gmpct/100;
    opt.Sm = khr.SmHW*khr.Smpct/100;
    
    opt.M_t = zeros(3*spc.P,opt.N+1);
    
    
    
    
    opt.M_t(:,1:opt.mon(1)) = repmat(spc.M0,1,opt.mon(1));
    opt.L_t = zeros(3*spc.P,opt.N+1);
    opt.L_t(:,end-opt.mon(end)+1:end) = repmat(spc.Md,1,opt.mon(end));
    
    
    
    
end
end


%% Organization Functions
function opt = Run_this(spc,khr,opt)

opt.OptNum = 1;
opt2 = [];
%% First optimization
opt1 = opt;

opt1.MaxIter                 = opt.MaxIter(1);
opt1.TolX                    = opt.TolX(1);
opt1.TolFun                  = opt.TolFun(1);
opt1.TolCon                  = opt.TolCon(1);
opt1.ObjLim                  = opt.ObjLim(1);
opt1.MaxFunEvals             = opt.MaxFunEvals(1);
if iscell(opt.Method)
    opt1.Method = opt.Method{1};
else
    opt1.Method = opt.Method;
end
if iscell(opt.Par)
    opt1.Par = opt.Par{1};
else
    opt1.Par = opt.Par;
end



try
    opt1 = Prepare_khr_4_opt(opt1,khr);
    try
        opt1 = Get_Initial_RF(spc,khr,opt1);
        try
            switch opt1.Method
                
                case 'MC_MadayTurinici'
                    opt1 = MC_MadayTurinici(spc,khr,opt1);
                case 'GRAPE_Khaneja'
                    opt1 = GRAPE_Khaneja(spc,khr,opt1);
                    
                case 'QN_LBFGS'
                    opt1= QN_LBFGS(spc,khr,opt1);
                otherwise
                    
                    fun = str2func(opt1.Method);
                    opt1 = fun(spc,khr,opt1);
                    
            end
            try
                opt1 = Save_Job(spc,khr,opt1);
                %
                
                %
            catch me;Display_Message(['Run_this: Save_Job: ',me.message],2);end
        catch me;Display_Message(sprintf('Run_this: %s: %s',opt1.Method,me.message),2);end
    catch me;Display_Message(['Run_this: Get_Initial_RF: ',me.message],2);end
catch me;Display_Message(['Run_this: Prepare_khr_4_opt: ',me.message],2);end



%% Second Optimization
if opt.Handover
    
    opt = Handover(opt,opt1);
    opt.OptNum = 2;
    
    opt2 = opt;
    
    opt2.MaxIter                 = opt.MaxIter(2);
    opt2.TolX                    = opt.TolX(2);
    opt2.TolFun                  = opt.TolFun(2);
    opt2.TolCon                  = opt.TolCon(2);
    opt2.ObjLim                  = opt.ObjLim (2);
    opt2.MaxFunEvals             = opt.MaxFunEvals(2);
    if iscell(opt.Method)
        opt2.Method = opt.Method{2};
    end
    if iscell(opt.Par)
        opt2.Par = opt.Par{2};
        
    end
    
    try
        opt2 = Prepare_khr_4_opt(opt2,khr);
        try
            opt2 = Get_Initial_RF(spc,khr,opt2);
            
            try
                switch opt2.Method
                    
                    case 'MC_MadayTurinici'
                        opt2 = MC_MadayTurinici(spc,khr,opt2);
                    case 'GRAPE_Khaneja'
                        opt2 = GRAPE_Khaneja(spc,khr,opt2);
                        
                    case 'QN_LBFGS'
                        opt2 = QN_LBFGS(spc,khr,opt2);
                    otherwise
                        
                        fun = str2func(opt2.Method);
                        opt2 = fun(spc,khr,opt2);
                        %                         eval(sprintf('opt2 = %s(spc,khr,opt2);',opt2.Method))
                        
                end
                 try 
                    opt2 = Save_Job(spc,khr,opt2);

                 catch me;Display_Message(['Run_this: Save_Job: ',me.message],2);end
             catch me;Display_Message(sprintf('Run_this: %s: %s',opt2.Method,me.message),2);end
        catch me;Display_Message(['Run_this: Get_Initial_RF: ',me.message],2);end
     catch me;Display_Message(['Run_this: Prepare_khr_4_opt: ',me.message],2);end
    
    
    opt.opt1 = opt1;
    opt.opt2 = opt2;
else
    opt.opt1 = opt1;
end
try
opt.Go = 0;
% just to be sure
opt = Save_Job(spc,khr,opt);
%

%
catch me;Display_Message(['Run_this: Save_Job: ',me.message],2);end

end

function opt = Handover(opt,opt1)

if opt.OptNum == 1 && opt.Handover
    
    switch opt.HandoverWhat
        
        case 'best'
            
            
            [opt.HandoverFun,opt.HandoverIdx] = max(opt1.Fun);
            if isfield(opt1,'gxo')
                opt.HandoverControls1 = {opt1.uo(:,:,opt.HandoverIdx);opt1.vo(:,:,opt.HandoverIdx);opt1.gxo(1,:,opt.HandoverIdx);opt1.gyo(1,:,opt.HandoverIdx);opt1.gzo(1,:,opt.HandoverIdx)};
            else
                opt.HandoverControls1 = {opt1.uo(:,:,opt.HandoverIdx);opt1.vo(:,:,opt.HandoverIdx);opt1.g(1,:);opt1.g(2,:);opt1.g(3,:)};
                
            end
        case 'last'
            opt.HandoverIdx = numel(opt1.Fun);
            opt.HandoverFun = opt1.Fun(end);
            if isfield(opt1,'gxo')
                opt.HandoverControls1 = {opt1.uo(:,:,end);opt1.vo(:,:,end);opt1.gxo(1,:,end);opt1.gyo(1,:,end);opt1.gzo(1,:,end)};
            else
                opt.HandoverControls1 = {opt1.uo(:,:,end);opt1.vo(:,:,end);opt1.g(1,:);opt1.g(2,:);opt1.g(3,:)};
                
            end
    end
    
    
    
    
    
end

end

%% Rotation stuff

function [R11,R12,R13,R21,R22,R23,R31,R32,R33] = Get_Rotator(u,v,g,w0,yx,yy,yz,pTx,sr,si,dt,n)

%BUG:
wx = zeros(size(w0));
wy = zeros(size(w0));
for s = 1:pTx
    wx = wx+(sr(:,s).*u(s,n)-si(:,s).*v(s,n)).*dt;
    wy = wy+(si(:,s).*u(s,n)+sr(:,s).*v(s,n)).*dt;
end
wz = (yx*g(1,n) + yy*g(2,n)+yz*g(3,n)+w0).*dt;





wx = wx+eps;
wy = wy+eps;
wz = wz+eps;

wx = -wx;
wy = -wy;
wz = -wz;


theta = sqrt(wx.^2+wy.^2+wz.^2);

x = wx./theta;
y = wy./theta;
z = wz./theta;



c = cos(theta);
s = sin(theta);

C =1-c;

zs = z.*s;
xs = x.*s;
ys = y.*s;

xyC = x.*y.*C;
yzC = y.*z.*C;
xzC = x.*z.*C;

R11 = x.*x.*C+c;
R12 = xyC-zs;
R13 = xzC+ys;

R21 = xyC+zs;
R22 = y.*y.*C+c;
R23 = yzC-xs;

R31 = xzC-ys;
R32 = yzC+xs;
R33 = z.*z.*C+c;

end

function ML_out = Rotate(ML_in,R11,R12,R13,R21,R22,R23,R31,R32,R33,Type)



[P,N] = size(R11);

% become row
ML = ML_in(:)';

% become
% [Mx1 Mx2     MxP]
% [My1 My2 ... MyP]
% [Mz1 Mz2     MzP];

ML = reshape(ML,3,P);

% become
% [Mx1 My1 Mz1;
%  Mx2 My2 Mz2;
%    .
%    .
%    .
%  MxP MyP MzP];

ML = permute(ML,[2 1]);

MLx = ML(:,1);
MLy = ML(:,2);
MLz = ML(:,3);


if strcmp(Type,'Forward')
    
    if N > 1
        
        for n = 1:N
            
            MLx_ = R11(:,n).*MLx + R12(:,n).*MLy + R13(:,n).*MLz;
            MLy_ = R21(:,n).*MLx + R22(:,n).*MLy + R23(:,n).*MLz;
            MLz_ = R31(:,n).*MLx + R32(:,n).*MLy + R33(:,n).*MLz;
            
            MLx = MLx_;
            MLy = MLy_;
            MLz = MLz_;
            
            
            
        end
        
        
    else
        
        MLx_ = R11.*MLx + R12.*MLy + R13.*MLz;
        MLy_ = R21.*MLx + R22.*MLy + R23.*MLz;
        MLz_ = R31.*MLx + R32.*MLy + R33.*MLz;
        
        MLx = MLx_;
        MLy = MLy_;
        MLz = MLz_;
        
    end
    
    
    
elseif strcmp(Type,'Backward')
    
    if N > 1
        
        for n = N:-1:1
            
            MLx_ = R11(:,n).*MLx + R21(:,n).*MLy + R31(:,n).*MLz;
            MLy_ = R12(:,n).*MLx + R22(:,n).*MLy + R32(:,n).*MLz;
            MLz_ = R13(:,n).*MLx + R23(:,n).*MLy + R33(:,n).*MLz;
            
            MLx = MLx_;
            MLy = MLy_;
            MLz = MLz_;
            
            
            
        end
        
        
    else
        
        MLx_ = R11.*MLx + R21.*MLy + R31.*MLz;
        MLy_ = R12.*MLx + R22.*MLy + R32.*MLz;
        MLz_ = R13.*MLx + R23.*MLy + R33.*MLz;
        
        MLx = MLx_;
        MLy = MLy_;
        MLz = MLz_;
        
        
    end
    
    
    
end

ML(:,1) = MLx;
ML(:,2) = MLy;
ML(:,3) = MLz;

ML = permute(ML,[2 1]);

ML = reshape(ML,3*P,1);

ML_out = ML(:);
end

%% Control Update stuff

function [Gu,Gv] = RF_Update_Terms(spc,opt,n)


if ~isempty(n)
    if ismember(n,opt.mon)
        
        switch opt.Par.Grad
            case 'Schirmer'
                
                [Gu,Gv] =    Get_pTx_u_v_Schirmer_Gradients(opt.u,opt.v,opt.M_t,opt.L_t,opt.g,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.w0,opt.dt,n);
                
            case '1st'
                [Gu,Gv] = Get_pTx_u_v_1st_Gradients(opt.M_t,opt.L_t,opt.sr,opt.si,spc.pTx,n);
                
            case 'type1'
                [Gu,Gv] = Get_pTx_u_v_type1_Gradients(opt.M_t,opt.L_t,opt.sr,opt.si,spc.pTx,n);
                
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
                    [Gu(:,n),Gv(:,n)] =    Get_pTx_u_v_Schirmer_Gradients(opt.u,opt.v,opt.M_t,opt.L_t,opt.g,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.w0,opt.dt,n,opt.par_Ncores,opt.par_Type);
                    end
                end
            elseif opt.par_Ncores > 1 && opt.par_Type == 2
                % Type 2 parallelizes Gu,Gv per pTx channel
                for n = 1:opt.N
                    if ismember(n,opt.mon)
                    [Gu(:,n),Gv(:,n)] =    Get_pTx_u_v_Schirmer_Gradients(opt.u,opt.v,opt.M_t,opt.L_t,opt.g,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.w0,opt.dt,n,opt.par_Ncores,opt.par_Type);
                end
            end
            else
                
                for n = 1:opt.N
                    if ismember(n,opt.mon)
                        [Gu(:,n),Gv(:,n)] =    Get_pTx_u_v_Schirmer_Gradients(opt.u,opt.v,opt.M_t,opt.L_t,opt.g,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.w0,opt.dt,n,opt.par_Ncores,opt.par_Type);
                    end
                end
            end
            
        case '1st'
            
            if opt.par_Ncores > 1
                parfor (n = 1:opt.N,opt.par_Ncores)
                    if ismember(n,opt.mon)
                        [Gu(:,n),Gv(:,n)] = Get_pTx_u_v_1st_Gradients(opt.M_t,opt.L_t,opt.sr,opt.si,spc.pTx,n);
                    end
                end
            else
                for n = 1:opt.N
                    if ismember(n,opt.mon)
                        [Gu(:,n),Gv(:,n)] = Get_pTx_u_v_1st_Gradients(opt.M_t,opt.L_t,opt.sr,opt.si,spc.pTx,n);
                    end
                end
            end
            
            
            
        case 'type1'
            if opt.par_Ncores > 1
                parfor (n = 1:opt.N,opt.par_Ncores)
                    if ismember(n,mon)
                        [Gu(:,n),Gv(:,n)] = Get_pTx_u_v_type1_Gradients(opt.M_t,opt.L_t,opt.sr,opt.si,spc.pTx,n);
                    end
                    
                end
            else
                for n = 1:opt.N
                    if ismember(n,opt.mon)
                        [Gu(:,n),Gv(:,n)] = Get_pTx_u_v_type1_Gradients(opt.M_t,opt.L_t,opt.sr,opt.si,spc.pTx,n);
                    end
                end
                
            end
            
    end
    
end
Gu(isnan(Gu)) = 0;
Gv(isnan(Gv)) = 0;




end

function [Gu,Gv] = Get_pTx_u_v_Schirmer_Gradients(u,v,M_t,L_t,g,yx,yy,yz,pTx,sr,si,w0,dt,n,par_Ncores,par_Type)

P = size(w0,1);
z = (yx*g(1,n) + yy*g(2,n)+yz*g(3,n)+w0);
x = sum((sr.*repmat(u(:,n).',[P,1])-si.*repmat(v(:,n).',[P,1])),2);
y = sum((si.*repmat(u(:,n).',[P,1])+sr.*repmat(v(:,n).',[P,1])),2);


N_block = P;


idx1 = repmat([1:3],[1,N_block]);
idx2 = repmat([0:6:(N_block-1)*6],[3,1]);
idx2 = idx2(:);
idx = idx1+idx2.';

vector1 = M_t(:,n);
vector1 = reshape(vector1(:).',[3,P]);
vector1 = repmat(vector1,[2,1]);
vector1 = vector1(:);

vector2 = M_t(:,n+1);
vector2 = reshape(vector2(:).',[3,P]);
vector2 = repmat(vector2,[2,1]);
vector2 = vector2(:);

vector3 = L_t(:,n+1);

G = zeros(pTx,2);
for g = 1:2
    
    switch g
        
        case 1
            
            a = si; b = -sr;
        case 2
            a = sr; b = si;
            
    end
    
    
    if par_Ncores > 1 && par_Type == 2
        parfor (s = 1:pTx, par_Ncores)
            Dm2 = zeros(P,6);
            k = 1;
            Dm2(:,k) = y;
            Dm2(:,k+3) = y;
            
            Dm1 = zeros(P,6);
            k = 1;
            Dm1(:,k) = -z;
            Dm1(:,k+1) = -x;
            Dm1(:,k+3) = -z;
            Dm1(:,k+4) = -x;
            Dp1 = zeros(P,6);
            k = 1;
            Dp1(:,k) = z;
            Dp1(:,k+1) = x;
            Dp1(:,k+2) = a(:,s);
            Dp1(:,k+3) = z;
            Dp1(:,k+4) = x;
            Dp2 = zeros(P,6);
            k = 1;
            Dp2(:,k) = -y;
            
            Dp2(:,k+2) = b(:,s);
            Dp2(:,k+3) = -y;
            
            Dp4 = zeros(P,6);
            k = 1;
            Dp4(:,k+1) = -b(:,s);
            
            Dp5 = zeros(P,6);
            k = 1;
            Dp5(:,k) = -a(:,s);
            
            
            Dm2_ = Dm2.';Dm2_ = Dm2_(:);
            Dm1_ = Dm1.';Dm1_ = Dm1_(:);
            Dp1_ = Dp1.';Dp1_ = Dp1_(:);
            Dp2_ = Dp2.';Dp2_ = Dp2_(:);
            Dp4_ = Dp4.';Dp4_ = Dp4_(:);
            Dp5_ = Dp5.';Dp5_ = Dp5_(:);
            
            Dm2_ = [Dm2_;0;0;0;0;0]; %#ok<*AGROW>
            Dm1_ = [Dm1_;0;0;0;0;0];
            Dp1_ = [0;Dp1_;0;0;0;0];
            Dp2_ = [0;0;Dp2_;0;0;0];
            Dp4_ = [0;0;0;0;Dp4_;0];
            Dp5_ = [0;0;0;0;0;Dp5_];
            
            matrix = spdiags([Dm2_,Dm1_,Dp1_,Dp2_,Dp4_,Dp5_],[-2,-1,1,2,4,5],N_block*6,N_block*6);
            
            [vector_prop1, err, hump] = expv3(dt,matrix,vector1,1e-4,10);
            vector_prop = vector_prop1 - vector2;
            
            vector_prop2 = vector_prop(idx);
            vector_prop2 = vector_prop2(:).';
            vector_prop2 = reshape(vector_prop2,[3*P,1]);
            
            
            
            G(s,g) = vector3'*vector_prop2;
            
            
            
        end
    else
        for s = 1:pTx
            
            Dm2 = zeros(P,6);
            k = 1;
            Dm2(:,k) = y;
            Dm2(:,k+3) = y;
            
            Dm1 = zeros(P,6);
            k = 1;
            Dm1(:,k) = -z;
            Dm1(:,k+1) = -x;
            Dm1(:,k+3) = -z;
            Dm1(:,k+4) = -x;
            
            Dp1 = zeros(P,6);
            k = 1;
            Dp1(:,k) = z;
            Dp1(:,k+1) = x;
            Dp1(:,k+2) = a(:,s);
            Dp1(:,k+3) = z;
            Dp1(:,k+4) = x;
            Dp2 = zeros(P,6);
            k = 1;
            Dp2(:,k) = -y;
            
            Dp2(:,k+2) = b(:,s);
            Dp2(:,k+3) = -y;
            
            Dp4 = zeros(P,6);
            k = 1;
            Dp4(:,k+1) = -b(:,s);
            
            Dp5 = zeros(P,6);
            k = 1;
            Dp5(:,k) = -a(:,s);
            
            
            Dm2_ = Dm2.';Dm2_ = Dm2_(:);
            Dm1_ = Dm1.';Dm1_ = Dm1_(:);
            Dp1_ = Dp1.';Dp1_ = Dp1_(:);
            Dp2_ = Dp2.';Dp2_ = Dp2_(:);
            Dp4_ = Dp4.';Dp4_ = Dp4_(:);
            Dp5_ = Dp5.';Dp5_ = Dp5_(:);
            
            Dm2_ = [Dm2_;0;0;0;0;0]; %#ok<*AGROW>
            Dm1_ = [Dm1_;0;0;0;0;0];
            Dp1_ = [0;Dp1_;0;0;0;0];
            Dp2_ = [0;0;Dp2_;0;0;0];
            Dp4_ = [0;0;0;0;Dp4_;0];
            Dp5_ = [0;0;0;0;0;Dp5_];
            
            matrix = spdiags([Dm2_,Dm1_,Dp1_,Dp2_,Dp4_,Dp5_],[-2,-1,1,2,4,5],N_block*6,N_block*6);
            [vector_prop1, err, hump] = expv3(dt,matrix,vector1,1e-4,10);
            vector_prop = vector_prop1 - vector2;
            vector_prop2 = vector_prop(idx);
            vector_prop2 = vector_prop2(:).';
            vector_prop2 = reshape(vector_prop2,[3*P,1]);
            
            
            
            G(s,g) = vector3'*vector_prop2;
            
            
            
        end
    end
    
end
Gu = G(:,1);
Gv = G(:,2);
end

function [Gu,Gv] = Get_pTx_u_v_1st_Gradients(M,L,sr,si,pTx,n)
% function [Gu,Gv] = blOCh__4_0__Get_pTx_u_v_1st_Gradients(M,L,sr,si,pTx,n)
%
%   This script gets the 1st order u,v pTx control gradients

Gu = zeros(pTx,1);
Gv = zeros(pTx,1);
for s = 1:pTx
    Gu(s) = sum(si(:,s).*L(3:3:end,n+1).*M(1:3:end,n))-sum(sr(:,s).*L(3:3:end,n+1).*M(2:3:end,n))-sum(si(:,s).*L(1:3:end,n+1).*M(3:3:end,n))+sum(sr(:,s).*L(2:3:end,n+1).*M(3:3:end,n));
    Gv(s) = sum(si(:,s).*L(3:3:end,n+1).*M(2:3:end,n))+sum(sr(:,s).*L(3:3:end,n+1).*M(1:3:end,n))-sum(si(:,s).*L(2:3:end,n+1).*M(3:3:end,n))-sum(sr(:,s).*L(1:3:end,n+1).*M(3:3:end,n));
end

end

function [Gu,Gv] = Get_pTx_u_v_type1_Gradients(M,L,sr,si,pTx,n)

Gu_left = zeros(pTx,1);
Gv_left = zeros(pTx,1);
Gu_right = zeros(pTx,1);
Gv_right = zeros(pTx,1);
for s = 1:pTx
    Gu_left(s)  = sum(si(:,s).*L(3:3:end,n).*M(1:3:end,n))-sum(sr(:,s).*L(3:3:end,n).*M(2:3:end,n))-sum(si(:,s).*L(1:3:end,n).*M(3:3:end,n))+sum(sr(:,s).*L(2:3:end,n).*M(3:3:end,n));
    Gv_left(s)  = sum(si(:,s).*L(3:3:end,n).*M(2:3:end,n))+sum(sr(:,s).*L(3:3:end,n).*M(1:3:end,n))-sum(si(:,s).*L(2:3:end,n).*M(3:3:end,n))-sum(sr(:,s).*L(1:3:end,n).*M(3:3:end,n));
    Gu_right(s) = sum(si(:,s).*L(3:3:end,n+1).*M(1:3:end,n+1))-sum(sr(:,s).*L(3:3:end,n+1).*M(2:3:end,n+1))-sum(si(:,s).*L(1:3:end,n+1).*M(3:3:end,n+1))+sum(sr(:,s).*L(2:3:end,n+1).*M(3:3:end,n+1));
    Gv_right(s) = sum(si(:,s).*L(3:3:end,n+1).*M(2:3:end,n+1))+sum(sr(:,s).*L(3:3:end,n+1).*M(1:3:end,n+1))-sum(si(:,s).*L(2:3:end,n+1).*M(3:3:end,n+1))-sum(sr(:,s).*L(1:3:end,n+1).*M(3:3:end,n+1));
end
Gu = (Gu_left+Gu_right)./2;
Gv = (Gv_left+Gv_right)./2;
end

function [w, err, hump] = expv3( t, A, v, tol, m )
% function [w, err, hump] = expv3( t, A, v, tol, m )
%  [w, err, hump] = expv( t, A, v, tol, m )
%  EXPV computes an approximation of w = exp(t*A)*v for a
%  general matrix A using Krylov subspace  projection techniques.
%  It does not compute the matrix exponential in isolation but instead,
%  it computes directly the action of the exponential operator on the
%  operand vector. This way of doing so allows for addressing large
%  sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%
%  w = expv( t, A, v )
%  computes w = exp(t*A)*v using a default tol = 1.0e-7 and m = 30.
%
%  [w, err] = expv( t, A, v )
%  renders an estimate of the error on the approximation.
%
%  [w, err] = expv( t, A, v, tol )
%  overrides default tolerance.
%
%  [w, err, hump] = expv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace,
%  and renders an approximation of the `hump'.
%
%  The hump is defined as:
%          hump = max||exp(sA)||, s in [0,t]  (or s in [t,0] if t < 0).
%  It is used as a measure of the conditioning of the matrix exponential
%  problem. The matrix exponential is well-conditioned if hump = 1,
%  whereas it is poorly-conditioned if hump >> 1. However the solution
%  can still be relatively fairly accurate even when the hump is large
%  (the hump is an upper bound), especially when the hump and
%  ||w(t)||/||v|| are of the same order of magnitude (further details in
%  reference below).
%
%  Example 1:
%  ----------
%    n = 100;
%    A = rand(n);
%    v = eye(n,1);
%    w = expv(1,A,v);
%
%  Example 2:
%  ----------
%    % generate a random sparse matrix
%    n = 100;
%    A = rand(n);
%    for j = 1:n
%        for i = 1:n
%            if rand < 0.5, A(i,j) = 0; end;
%        end;
%    end;
%    v = eye(n,1);
%    A = sparse(A); % invaluable for a large and sparse matrix.
%
%    tic
%    [w,err] = expv(1,A,v);
%    toc
%
%    disp('w(1:10) ='); disp(w(1:10));
%    disp('err =');     disp(err);
%
%    tic
%    w_matlab = expm(full(A))*v;
%    toc
%
%    disp('w_matlab(1:10) ='); disp(w_matlab(1:10));
%    gap = norm(w-w_matlab)/norm(w_matlab);
%    disp('||w-w_matlab|| / ||w_matlab|| ='); disp(gap);
%
%  In the above example, n could have been set to a larger value,
%  but the computation of w_matlab will be too long (feel free to
%  discard this computation).
%
%  See also MEXPV, EXPOKIT.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials.
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

[n,n] = size(A);
if nargin == 3,
    tol = 1.0e-7;
    m = min(n,30);
end;
if nargin == 4,
    m = min(n,30);
end;

anorm = norm(A,'inf');
mxrej = 10;  btol  = 1.0e-7;
gamma = 0.9; delta = 1.2;
mb    = m; t_out   = abs(t);
nstep = 0; t_new   = 0;
t_now = 0; s_error = 0;
rndoff= anorm*eps;

k1 = 2; xm = 1/m; normv = norm(v); beta = normv;
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s;
sgn = sign(t); nstep = 0;

w = v;
hump = normv;
while t_now < t_out
    nstep = nstep + 1;
    t_step = min( t_out-t_now,t_new );
    V = zeros(n,m+1);
    H = zeros(m+2,m+2);
    
    V(:,1) = (1/beta)*w;
    for j = 1:m
        p = A*V(:,j);
        for i = 1:j
            H(i,j) = V(:,i)'*p;
            p = p-H(i,j)*V(:,i);
        end;
        s = norm(p);
        if s < btol,
            k1 = 0;
            mb = j;
            t_step = t_out-t_now;
            break;
        end;
        H(j+1,j) = s;
        V(:,j+1) = (1/s)*p;
    end;
    if k1 ~= 0,
        H(m+2,m+1) = 1;
        avnorm = norm(A*V(:,m+1));
    end;
    ireject = 0;
    while ireject <= mxrej,
        mx = mb + k1;
        F = expm(sgn*t_step*H(1:mx,1:mx));
        if k1 == 0,
            err_loc = btol;
            break;
        else
            phi1 = abs( beta*F(m+1,1) );
            phi2 = abs( beta*F(m+2,1) * avnorm );
            if phi1 > 10*phi2,
                err_loc = phi2;
                xm = 1/m;
            elseif phi1 > phi2,
                err_loc = (phi1*phi2)/(phi1-phi2);
                xm = 1/m;
            else
                err_loc = phi1;
                xm = 1/(m-1);
            end;
        end;
        if err_loc <= delta * t_step*tol,
            break;
        else
            t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
            s = 10^(floor(log10(t_step))-1);
            t_step = ceil(t_step/s) * s;
            if ireject == mxrej,
                error('The requested tolerance is too high.');
            end;
            ireject = ireject + 1;
        end;
    end;
    mx = mb + max( 0,k1-1 );
    w = V(:,1:mx)*(beta*F(1:mx,1));
    beta = norm( w );
    hump = max(hump,beta);
    
    t_now = t_now + t_step;
    t_new = gamma * t_step * (t_step*tol/err_loc)^xm;
    s = 10^(floor(log10(t_new))-1);
    t_new = ceil(t_new/s) * s;
    
    err_loc = max(err_loc,rndoff);
    s_error = s_error + err_loc;
end;
err = s_error;
hump = hump / normv;



end

%% Utilities

function Display_Message(Msg,Type)



if nargin == 1
    Type = 1;
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

function opt = Save_Job(spc,khr,opt,varargin)
% function opt = Save_Job(spc,khr,opt,varargin)
%
%   Stores the different structures, and prints output to a datafile
%   All dependent on the choices described below.
%
%   If:
%       opt.Save =   0;     Stores nothing
%       opt.Save =   ABCDEF;
%       A = 1: matlab mat-file with reduced spc,khr,opt and sim-structs
%       A = 2: matlab mat-file with complete spc,khr,opt and sim-structs
%       B = 1: ascii txt-file with important data values stored in the end of an optimization instance
%       B = 2: ascii txt-file with important data values stored every ksaveintermediate iteration
%       C = 1: ascii txt-file of the last controls stored in the end of an optimization instance
%       C = 2: ascii txt-file of the controls stored every ksaveintermediate iteration
%       C = 3: ascii controls txt-file of the controls stored every ksaveintermediate iteration and the last controls stored in the end of an optimization instance
%       C = 4: ascii controls txt-file of the controls stored every ksaveintermediate iteration and all the controls stored in the end of an optimization instance
%       C = 5: ascii txt-file of all the controls stored in the end of an optimization instance
%       D = 1: png images of displayed figures
%       E = 1: system specific control files ... well not supported yet.
%       F = 1: Saves the main script from which the optimization is run.
%       F = 2: As with F=1, but also a zip of the folder with the blOCh
%       scripts, but not package scripts yet.

if sum(str2num(opt.Save)) > 0
    
    SaveThis = SaveThis(opt.Save);
    
    Type = My_Type(opt);
    
    if length(varargin) >= 1
        if isstruct(varargin{1})
            sim = varargin{1};
            AppendName = '_sim';
        else
            sim = [];
            AppendName = '';
        end
    else
        sim = [];
        AppendName = '';
    end
    
    % About bundles
    
    if numel(SaveThis) >=1
        if SaveThis(1) > 0
            
            opt.Folder_bnd = [opt.Save2,opt.TimeStamp,filesep,'bundles',filesep];
            opt.File_bnd = sprintf('%s%s_on%i.bnd',Type,AppendName,opt.OptNum);
        end
    end
    
    % abouyt ascii data
    
    if numel(SaveThis) >=2
        
        if opt.Go
            if SaveThis(2) == 2
                opt.Folder_dat = [opt.Save2,opt.TimeStamp,filesep,'data',filesep];
                opt.File_dat_intermediate = sprintf('%s%s_ik%04.f.txt',Type,AppendName,opt.k);
                opt.Hdr_dat_intermediate = sprintf('%s%s_ik%04.f_hdr.txt',Type,AppendName,opt.k);
            end
        else
            if SaveThis(2) >= 1
                opt.Folder_dat = [opt.Save2,opt.TimeStamp,filesep,'data',filesep];
                opt.File_dat = sprintf('%s%s_on%i.txt',Type,AppendName,opt.OptNum);
                opt.Hdr_dat = sprintf('%s%s_on%i_hdr.txt',Type,AppendName,opt.OptNum);
            end
        end
        
        
    end
    
    % about ascii controls
    
    if numel(SaveThis) >=3
        
        if opt.Go
            
            
            if SaveThis(3) == 2 || SaveThis(3) == 3 || SaveThis(3) == 4
                
                opt.Folder_ctrl = [opt.Save2,opt.TimeStamp,filesep,'controls',filesep];
                
                opt.File_ctrl_RF_intermediate  = sprintf('%s%s_on%i_ik%04.f.rf.txt',Type,AppendName,opt.OptNum,opt.k);
                
            end
            
            
        else
            
            if SaveThis(3) == 1 || SaveThis(3) == 3
                
                opt.Folder_ctrl = [opt.Save2,opt.TimeStamp,filesep,'controls',filesep];
                
                opt.File_ctrl_RF_last = sprintf('%s%s_on%i_last.rf.txt',Type,AppendName,opt.OptNum);
                
            end
            
            if SaveThis(3) == 4 || SaveThis(3) == 5
                
                rfx = opt.uo;
                rfy = opt.vo;
                
                [A,B,C] = size(rfx);
                
                opt.Folder_ctrl = [opt.Save2,opt.TimeStamp,filesep,'controls',filesep];
                opt.File_ctrl_RF_all = cell(C,1);
                for c = 1:C
                    opt.File_ctrl_RF_all{c} = sprintf('%s%s_%i_k%04.f.rf.txt',Type,AppendName,opt.OptNum,c);
                end
            end
        end
        
        
        
        
    end
    
    % about pngs
    
    opt.Folder_pics = [opt.Save2,opt.TimeStamp,filesep,'pics',filesep];
    
    if numel(SaveThis) >=4
        
        if SaveThis(4)
            
            opt.File_pic_khr = sprintf('khr_%s%s_on%i.png',Type,AppendName,opt.OptNum);
            opt.File_pic_spc = sprintf('spc_%s%s_on%i.png',Type,AppendName,opt.OptNum);
            opt.File_pic_opt = sprintf('opt_%s%s_on%i.png',Type,AppendName,opt.OptNum);
            opt.File_pic_sim = sprintf('sim_%s%s_on%i.png',Type,AppendName,opt.OptNum);
        end
    end
    
    % about system related files
    % not supported yet
%     if numel(SaveThis) >= 5
%         if SaveThis(5)
%             opt.Folder_sys = [opt.Save2,opt.TimeStamp,filesep,'sys',filesep];
%             opt.File_system_RF  = sprintf('%s%s_on%i_k%04.f',Type,AppendName,opt.OptNum,opt.k);
%         end
%         
%     end
    
    % abouts main and scripts files
    
    if numel(SaveThis) >=6
        if SaveThis(6)
            opt.Folder_scripts = [opt.Save2,opt.TimeStamp,filesep,'scripts',filesep];
        end
        
    end
    %% Figures
    % NB: This must appear before Save_Bundle because it erases the
    % fig-handles from spc and khr.
    try
        Save_Figures(spc,khr,opt,sim,SaveThis);
    catch me;Msg = ['Save_Job ',me.message];
        Display_Message(Msg,2);
    end
    %% Bundle
    try
        Save_Bundle(spc,khr,opt,sim,SaveThis);
    catch me;Msg = ['Save_Job ',me.message];
        Display_Message(Msg,2);
    end
    
    
    %% Data
    try
        Save_Data(spc,khr,opt,sim,SaveThis);
    catch me;Msg = ['Save_Job ',me.message];
        Display_Message(Msg,2);
    end
    %% Controls
    try
        Save_Controls(spc,khr,opt,sim,SaveThis);
    catch me;Msg = ['Save_Job ',me.message];
        Display_Message(Msg,2);
    end
    
    %% System specific
%     try
%         Save_System(spc,khr,opt,sim,SaveThis);
%     catch me;Msg = ['Save_Job ',me.message];
%         Display_Message(Msg,2);
%     end
    %% Scripts
    try
        Save_Scripts(spc,khr,opt,sim,SaveThis);
    catch me;Msg = ['Save_Job ',me.message];
        Display_Message(Msg,2);
    end
    %    %% Remove what cannot exist for next optimization
    
    
end
end

function Type = My_Type(opt)
if opt.OptNum == 1
    if iscell(opt.Method)
        Type = opt.Method{1};
    else
        Type = opt.Method;
    end
else
    if iscell(opt.Method)
        Type = opt.Method{2};
    else
        Type = opt.Method;
    end
end
end

function SaveThis = SaveThis(Save)

N = numel(Save);
% N = numel(num2str(fix(abs(Save))));

% B = num2str(fix(abs(Save)));

SaveThis = zeros(1,N);

for n = 1:N
    
    SaveThis(n) = str2double(Save(n));
    
    
end


end

function [spc,khr,opt,sim] = Remove_unimportant_stuff(spc,khr,opt,sim)

if isfield(opt,'opt1')
    temp = opt.opt1;
    
    temp = Remove_unimportant_from_opt(temp);
    opt.opt1 = temp;
end
if isfield(opt,'opt2')
    temp = opt.opt2;
    
    temp = Remove_unimportant_from_opt(temp);
    opt.opt2 = temp;
end

spc = Remove_unimportant_from_spc(spc);
sim = Remove_unimportant_from_sim(sim);
khr = Remove_unimportant_from_khr(khr);

end

function spc = Remove_unimportant_from_spc(spc)
if isfield(spc,'fig')
    spc = rmfield(spc,'fig');
    
end
if isfield(spc,'sar')
    spc = rmfield(spc,'sar');
end
end

function khr = Remove_unimportant_from_khr(khr)
if isfield(khr,'fig')
    khr = rmfield(khr,'fig');
end
end

function sim = Remove_unimportant_from_sim(sim)
if isfield(sim,'fig')
    sim = rmfield(sim,'fig');
end
end

function opt = Remove_unimportant_from_opt(opt)

if isfield(opt,'M_t')
    opt = rmfield(opt,'M_t');
end
if isfield(opt,'fig')
    opt = rmfield(opt,'fig');
end
if isfield(opt,'M_T')
    opt = rmfield(opt,'M_T');
end
if isfield(opt,'L_t')
    opt = rmfield(opt,'L_t');
end
if isfield(opt,'L_T')
    opt = rmfield(opt,'L_T');
end
if isfield(opt,'w0')
    opt = rmfield(opt,'w0');
end
if isfield(opt,'sr')
    opt = rmfield(opt,'sr');
end
if isfield(opt,'si')
    opt = rmfield(opt,'si');
end
if isfield(opt,'yx')
    opt = rmfield(opt,'yx');
end
if isfield(opt,'yy')
    opt = rmfield(opt,'yy');
end
if isfield(opt,'yz')
    opt = rmfield(opt,'yz');
end
if isfield(opt,'u')
    opt = rmfield(opt,'u');
end
if isfield(opt,'v')
    opt = rmfield(opt,'v');
end



if isfield(opt,'gt')
    opt = rmfield(opt,'gt');
end

if isfield(opt,'T1')
    opt = rmfield(opt,'T1');
end
if isfield(opt,'T2')
    opt = rmfield(opt,'T2');
end
if isfield(opt,'ut')
    opt = rmfield(opt,'ut');
end
if isfield(opt,'vt')
    opt = rmfield(opt,'vt');
end

if isfield(opt,'LocUCor')
    opt = rmfield(opt,'LocUCor');
end
if isfield(opt,'h')
    opt = rmfield(opt,'h');
end




end

function [spc,khr,opt,sim] = Save_Bundle(spc,khr,opt,sim,SaveThis)


opt = orderfields(opt); % order fields
tempdum = opt;					% safe store opt

if SaveThis(1) == 1				% remove unimportant for save
    [spc,khr,opt,sim] = Remove_unimportant_stuff(spc,khr,opt,sim);
end


%now save what shall be saved
if SaveThis(1) > 0
    % first create folder

    if ~exist(tempdum.Folder_bnd','dir')
        mkdir(tempdum.Folder_bnd)
    end
    % filename
    temp2 = [tempdum.Folder_bnd,tempdum.File_bnd];

    save(temp2,'spc','khr','opt','sim','-mat','-v7.3')
    % message
    Display_Message(sprintf('Saved: %s',temp2),1);
    
end

% restore opt
opt = tempdum;
% also pass opt type specific


end

function [spc,khr,opt,sim] = Save_Data(spc,khr,opt,sim,SaveThis)

if numel(SaveThis >=2)
    
    if opt.Go
        if SaveThis(2) == 2
            
            K = size(opt.uo,3);
            
            Fun = opt.Fun;
            Eff = opt.Eff;
            Pen = opt.Pen;
            
            
            
            if ~exist(opt.Folder_dat,'dir')
                mkdir(opt.Folder_dat)
            end
            
            k = opt.k-1; % subtract 1 because where this is called (Eval*_Prog*) Go = 1 and k has added 1 already
            
            stat = 1;
            
            fh = fopen([opt.Folder_dat,opt.Hdr_dat_intermediate],'w');
            fprintf(fh,'k\nFun\nEff\nPen');
            fclose(fh);
            
            f = fopen([opt.Folder_dat,opt.File_dat_intermediate],'w');
            fprintf(f,'%i %e %e %e\n',k,Fun(k),Eff(k),Pen(k));
            stat = fclose(f);
            
            
            if stat == 0
                Display_Message(sprintf('Saved: %s',[opt.Folder_dat,opt.File_dat_intermediate]),1);
            end
        end
    else
        
        
        
        
        if SaveThis(2) >= 1
            
            if ~exist(opt.Folder_dat,'dir')
                mkdir(opt.Folder_dat)
            end
            
            if isfield(opt,'opt1')
                K = size(opt.opt1.uo,3);
                
                Fun = opt.opt1.Fun;
                Eff = opt.opt1.Eff;
                Pen = opt.opt1.Pen;
            elseif isfield(opt,'opt2')
                K = size(opt.opt2.uo,3);
                
                Fun = opt.opt2.Fun;
                Eff = opt.opt2.Eff;
                Pen = opt.opt2.Pen;
            else
                K = size(opt.uo,3);
                
                Fun = opt.Fun;
                Eff = opt.Eff;
                Pen = opt.Pen;
            end

            fh = fopen([opt.Folder_dat,opt.Hdr_dat],'w');
            fprintf(fh,'k\nFun\nEff\nPen');
            fclose(fh);
            
            
            f = fopen([opt.Folder_dat,opt.File_dat],'w');
            for k = 1:K
                fprintf(f,'%i %e %e %e\n',k,Fun(k),Eff(k),Pen(k));
            end
            fclose(f);
            
            Display_Message(sprintf('Saved: %s',[opt.Folder_dat,opt.File_dat]),1);
            
        end
    end
    
    
end

end

function [spc,khr,opt,sim] = Save_Controls(spc,khr,opt,sim,SaveThis)

if numel(SaveThis >=3)
    
    if opt.Go
        if SaveThis(3) == 2 || SaveThis(3) == 3 || SaveThis(3) == 4
            rf = opt.xintermed;
            
            if ~exist(opt.Folder_ctrl,'dir')
                mkdir(opt.Folder_ctrl)
            end
            
            fid = fopen([opt.Folder_ctrl,opt.File_ctrl_RF_intermediate],'w');
            for nnn = 1:length(rf)
                fprintf(fid,'%f ',rf);
                
            end
            fclose(fid);
            
            if isunix
                if ismac
                    zip([opt.Folder_ctrl,opt.File_ctrl_RF_intermediate,'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_intermediate])
                    delete([opt.Folder_ctrl,opt.File_ctrl_RF_intermediate])
                else
                    zip([opt.Folder_ctrl,opt.File_ctrl_RF_intermediate,'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_intermediate])
                    delete([opt.Folder_ctrl,opt.File_ctrl_RF_intermediate])
                end
            else
                zip([opt.Folder_ctrl,opt.File_ctrl_RF_intermediate,'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_intermediate])
                delete([opt.Folder_ctrl,opt.File_ctrl_RF_intermediate])
            end
            
            Display_Message(sprintf('Saved: %s',[opt.Folder_ctrl,opt.File_ctrl_RF_intermediate,'.zip']),1);
            
        end
    else
        
        if SaveThis(3) == 1 || SaveThis(3) == 3
            
            if isfield(opt,'opt1')
                                rfx = opt.opt1.uo;
                rfy = opt.opt1.vo;
            elseif isfield(opt,'opt2')
                                rfx = opt.opt2.uo;
                rfy = opt.opt2.vo;
            else
                rfx = opt.uo;
                rfy = opt.vo;
            end

            [A,B,C] = size(rfx);
            
            stat = 1;
            if ~exist(opt.Folder_ctrl,'dir')
                mkdir(opt.Folder_ctrl)
            end
            fid = fopen([opt.Folder_ctrl,opt.File_ctrl_RF_last],'w');
            for b = 1:B
                
                for a = 1:A
                    fprintf(fid,'%f ',rfx(a,b));
                end
                for a = 1:A
                    fprintf(fid,'%f ',rfy(a,b));
                end
                fprintf(fid,'\n ');
            end
            stat = fclose(fid);
            
            if isunix
                if ismac
                    zip([opt.Folder_ctrl,opt.File_ctrl_RF_last,'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_last])
                    delete([opt.Folder_ctrl,opt.File_ctrl_RF_last])
                else
                    zip([opt.Folder_ctrl,opt.File_ctrl_RF_last,'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_last])
                    delete([opt.Folder_ctrl,opt.File_ctrl_RF_last])
                end
            else
                zip([opt.Folder_ctrl,opt.File_ctrl_RF_last,'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_last])
                delete([opt.Folder_ctrl,opt.File_ctrl_RF_last])
            end
            
            
            
            if ~stat
                Display_Message(sprintf('Saved: %s',[opt.Folder_ctrl,opt.File_ctrl_RF_last,'.zip']),1);
            end
            
            
        end
        
        
        
        
        
        
        if SaveThis(3) == 4 || SaveThis(3) == 5
            if opt.OptNum == 1
                
                rfx = opt.opt1.uo(:,:,end);
                rfy = opt.opt1.vo(:,:,end);
                
            elseif OptNum == 2
                
                rfx = opt.opt2.uo(:,:,end);
                rfy = opt.opt2.vo(:,:,end);
                
            end
            
            [A,B,C] = size(rfx);
            
            
            if ~exist(opt.Folder_ctrl,'dir')
                mkdir(opt.Folder_ctrl)
            end
            
            for c = 1:C
                stat = 1;
                fid = fopen([opt.Folder_ctrl,opt.File_ctrl_RF_all{c}],'w');
                for b = 1:B
                    
                    for a = 1:A
                        fprintf(fid,'%f ',rfx(a,b,c));
                    end
                    for a = 1:A
                        fprintf(fid,'%f ',rfy(a,b,c));
                    end
                    fprintf(fid,'\n ');
                end
                stat = fclose(fid);
                
                
                if isunix
                    if ismac
                        zip([opt.Folder_ctrl,opt.File_ctrl_RF_all{c},'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_all{c}])
                        delete([opt.Folder_ctrl,opt.File_ctrl_RF_all{c}])
                    else
                        zip([opt.Folder_ctrl,opt.File_ctrl_RF_all,'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_all{c}])
                        delete([opt.Folder_ctrl,opt.File_ctrl_RF_all{c}])
                    end
                else
                    zip([opt.Folder_ctrl,opt.File_ctrl_RF_all{c},'.zip'],[opt.Folder_ctrl,opt.File_ctrl_RF_all{c}])
                    delete([opt.Folder_ctrl,opt.File_ctrl_RF_all{c}])
                end
                
                
                if ~stat
                    Display_Message(sprintf('Saved: %s',[opt.Folder_ctrl,opt.File_ctrl_RF_all{c},'.zip']),1);
                end
            end
            
            
        end
        
        
        
    end
    
    
    
    
    
    
end

end

function [spc,khr,opt,sim] = Save_Figures(spc,khr,opt,sim,SaveThis)


if numel(SaveThis >=4)
    if ~exist(opt.Folder_pics,'dir')
        mkdir(opt.Folder_pics)
    end
    if SaveThis(4)
        if ~opt.Go
        if opt.OptNum == 1
            
            if isfield(khr,'fig')
                print(khr.fig,[opt.Folder_pics,opt.File_pic_khr],'-dpng')
                Display_Message(sprintf('Saved: %s',[opt.Folder_pics,opt.File_pic_khr]),1);
            end
            if isfield(spc,'fig')
                print(spc.fig,[opt.Folder_pics,opt.File_pic_spc],'-dpng')
                Display_Message(sprintf('Saved: %s',[opt.Folder_pics,opt.File_pic_spc]),1);
            end
            
            if isfield(opt,'fig')
                print(opt.fig,[opt.Folder_pics,opt.File_pic_opt],'-dpng')
                Display_Message(sprintf('Saved: %s',[opt.Folder_pics,opt.File_pic_opt]),1);
            end
            
        else
            
            if isfield(opt,'fig')
                print(opt.fig,[opt.Folder_pics,opt.File_pic_opt],'-dpng')
                Display_Message(sprintf('Saved: %s',[opt.Folder_pics,opt.File_pic_opt]),1);
            end
            
        end
        if ~isempty(sim)
            
            if isfield(sim,'fig')
                print(sim.fig,[opt.Folder_pics,opt.File_pic_sim],'-dpng')
                Display_Message(sprintf('Saved: %s',[opt.Folder_pics,opt.File_pic_sim]),1);
            end
        end
        end
        
    end
end

end


function [spc,khr,opt,sim] = Save_Scripts(spc,khr,opt,sim,SaveThis)


if numel(SaveThis) >=6
    if SaveThis(6) > 0
        % BUG: Perhaps Caller should be more that dbstack to 2nd level...
        if ~exist(opt.Folder_scripts,'dir')
            mkdir(opt.Folder_scripts)
        end
        if opt.OptNum == 1
            if SaveThis(6) == 1
                
                temp = strsplit(opt.Caller,'/');
                
                if strcmp(temp{end},'parallel_function.m')
                    
                else
                    copyfile(opt.Caller,[opt.Folder_scripts,'/',temp{end},'.run'])
                    Display_Message(sprintf('Saved main file to %s',opt.Folder_scripts),1);
                    
                end
                
                
                
                
                
                
                
            elseif SaveThis(6) == 2
                
                
                
                Path = strrep(mfilename('fullpath'),mfilename,'');
                
                % 			  copyfile(Path,opt.Folder_scripts)
                
                temp = strsplit(opt.Caller,'/');
                
                
                
                try
                    if strcmp(temp{end},'parallel_function.m')
                        Display_Message(sprintf('BUG: when running parallel job, blOCh is not able to save the main file...yet'),1);
                    else
                        copyfile(opt.Caller,[opt.Folder_scripts,'/',temp{end},'.run'])
                        Display_Message(sprintf('Saved main file to %s',opt.Folder_scripts),1);
                        
                    end
                catch me;
                    Display_Message(sprintf(['Was not abve to save main file ',me.message]),1);
                end
                
                try
                    if isunix
                        if ismac
                            zip([opt.Folder_scripts,'/scripts'],Path)
                        else
                            zip([opt.Folder_scripts,'/scripts'],Path)
                        end
                    else
                        zip([opt.Folder_scripts,'\scripts'],Path)
                    end
                    Display_Message(sprintf('Saved all blOCh-scripts to %s',opt.Folder_scripts),1);
                    
                catch me;
                    Display_Message(sprintf(['Was not able to save and zip scripts: ',me.message]),1);
                    
                end
                
                
                
                
                
                
                
            end
        end
    end
end
end
%% GRAPE_Khaneja

function opt = GRAPE_Khaneja(spc,khr,opt)

tic_Opt_Wall_Clock_Time = tic;

opt = Allocate_Variables(spc,khr,opt);
tic_One_Iteration = tic;
opt.k = 1;
% if opt.OptNum == 1
opt = Iteration_GRAPE_Khaneja(spc,opt,0);


% else
%     opt = Iteration_GRAPE_Khaneja(spc,opt,1);
% end
opt.M_T = opt.M_t(:,end);
opt.M_T_store(:,opt.k) = opt.M_T;






Eff = Efficiency_GRAPE_Khaneja(spc,opt);
opt.uo(:,:,opt.k) = opt.u;
opt.vo(:,:,opt.k) = opt.v;

Pen = opt.Par.lambda*sum((opt.u(:)).^2+(opt.v(:)).^2)*opt.dt;
Fun  = Eff-Pen;

opt.dFun(opt.k) = Fun;
opt.Fun(opt.k) = Fun;
opt.Eff(opt.k) = Eff;
opt.Pen(opt.k) = Pen;

Print_GRAPE_Khaneja(opt.k,Fun,Eff,Pen);

opt = Progress_GRAPE_Khaneja(spc,khr,opt);
opt.k = 2;

opt.Durations(opt.k-1) = toc(tic_One_Iteration);

while opt.Go
    tic_One_Iteration = tic;
    opt = Iteration_GRAPE_Khaneja(spc,opt,1);
    
    opt.M_T = opt.M_t(:,end);
    opt.M_T_store(:,opt.k) = opt.M_T;
    opt.uo(:,:,opt.k) = opt.u;
    opt.vo(:,:,opt.k) = opt.v;
    
    Eff = Efficiency_GRAPE_Khaneja(spc,opt);
    Pen = opt.Par.lambda*sum((opt.u(:)).^2+(opt.v(:)).^2)*opt.dt;
    
    Fun  = Eff-Pen;
    
    
    
    opt.dFun(opt.k) = Fun-opt.Fun(opt.k-1);
    opt.Fun(opt.k) = Fun;
    opt.Eff(opt.k) = Eff;
    opt.Pen(opt.k) = Pen;
    
    Print_GRAPE_Khaneja(opt.k,Fun,Eff,Pen);
    
    opt = Progress_GRAPE_Khaneja(spc,khr,opt);
    opt.Durations(opt.k-1) = toc(tic_One_Iteration);
end



opt.uo(:,:,opt.ksafe+1:end) = [];
opt.vo(:,:,opt.ksafe+1:end) = [];


opt.Opt_Wall_Clock_Time = toc(tic_Opt_Wall_Clock_Time);

end

function Eff = Efficiency_GRAPE_Khaneja(spc,opt,M_T)

if nargin ~= 3
    M_T = opt.M_T;
end

A_x = spc.Md(1:3:end);
A_y = spc.Md(2:3:end);
A_z = spc.Md(3:3:end);

B_x = M_T(1:3:end);
B_y = M_T(2:3:end);
B_z = M_T(3:3:end);
Phi = sum(A_x.*B_x+A_y.*B_y+A_z.*B_z);


Eff = Phi/spc.P;
end

function opt = Iteration_GRAPE_Khaneja(spc,opt,Option)

if nargin == 2
    Option = 1; % This propagates state vectors and updates controls
end
R11 = cell(1,opt.N);R12 = cell(1,opt.N);R13 = cell(1,opt.N);
R21 = cell(1,opt.N);R22 = cell(1,opt.N);R23 = cell(1,opt.N);
R31 = cell(1,opt.N);R32 = cell(1,opt.N);R33 = cell(1,opt.N);

for n = 1:opt.N
    [R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n}] = Get_Rotator(opt.u,opt.v,opt.g,opt.w0,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.dt,n);
end

for n = 1:opt.N
    opt.M_t(:,n+1) = Rotate(opt.M_t(:,n),R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n},'Forward');
end

for n = opt.N:-1:1
    opt.L_t(:,n) = Rotate(opt.L_t(:,n+1),R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n},'Backward');
end


if Option
    
    [Gu,Gv] = RF_Update_Terms(spc,opt,[]);
    
    %     [Gu(13),Gv(13)]
    for n = 1:opt.N
        
        opt = Update_GRAPE_Khaneja(spc,opt,Gu,Gv,n);
    end
    
end
end

function Print_GRAPE_Khaneja(k,Fun,Eff,Pen)

if k == 1
    
    Display_Message(sprintf('\nIter Fun Eff Pen'),1);
    
    Display_Message(sprintf('----------------------'),1);
    
    
    
    Display_Message(sprintf('%i %6.7g% 6.7g% 6.7g%',k-1,Fun,Eff,Pen),1);
    
    
    
    
else
    Display_Message(sprintf('%i %6.7g% 6.7g% 6.7g%',k-1,Fun,Eff,Pen),1);
    
end
end

function opt = Progress_GRAPE_Khaneja(spc,khr,opt)


opt.Conv = 0;

% Convergence, by Fun change


if ~isempty(opt.dFun)
    
    if abs(opt.dFun(opt.k)) > opt.TolFun
        opt.Conv = 0;
    elseif abs(opt.dFun(opt.k)) <= opt.TolFun
        
        if opt.k > 1
            
            
            
            
            
            
            if opt.Fun(opt.k) < opt.Fun(1)
                opt.Conv = 0;
                Display_Message(sprintf('Progress_GRAPE_Khaneja: Convergence was obtained (dFun < TolFun), but Fun(k) < Fun(1), so we continue...'),1);
            else
                opt.Conv = 1;
                Display_Message(sprintf('Progress_GRAPE_Khaneja: Convergence is obtained: dFun < TolFun'),1);
            end
            
        else
            opt.Conv = 1;
            Display_Message(sprintf('Progress_GRAPE_Khaneja: Convergence is obtained: dFun < TolFun'),1);
            
            
        end
        
        
    end
    
else
    opt.Conv = 0;
end






% Convergence, by Fun value

if opt.Fun(opt.k) >= 1-opt.ObjLim
    opt.Bingo = 1;
    
    Display_Message(sprintf('Progress_GRAPE_Khaneja: Reached desired performance: Fun >= 1-ObjLim'),1);
else
    opt.Bingo = 0;
end





if opt.k == opt.MaxIter+1
    opt.Done = 1;
    Display_Message(sprintf('Progress_GRAPE_Khaneja: Reached Iteration MaxIter'),1);
else
    opt.Done = 0;
end

if opt.Bingo == 0
    if opt.Done == 0
        if opt.Conv == 0
            opt.Go = true;
            opt.ksafe = opt.k;
            opt.k = opt.k+1;
            
            
            
        else
            opt.Go = false;
            opt.ksafe = opt.k;
        end
        
    else
        opt.Go = false;
        opt.ksafe = opt.k;
    end
    
else
    opt.Go = false;
    opt.ksafe = opt.k;
end


% save intermediate
if mod(opt.k,opt.ksaveintermediate) == 0
    
    
    tempsave = opt.Save; % take backup
    %
    opt.Save = '000000'; % reset
    opt.Save(2) = tempsave(2); % get initial intention back for data saving
    opt.Save(3) = tempsave(3); % get initial intention back for controls saving
    % thus, if intermediate saving of data and controls were intended it
    % will be saved. But other (bundle, pics, scripts etc.) will not, at
    % least not before the end if intended.
    try opt = Save_Job(spc,khr,opt);
    catch me; Display_Message(['Progress_GRAPE_Khaneja: Save_Job ',me.message],2); end
    
    opt.Save = tempsave;
    
end
end

function opt = Update_GRAPE_Khaneja(spc,opt,Gu,Gv,n)


switch opt.Par.Grad
    case {'1st','type1'}
        opt.u(:,n) = opt.u(:,n)+opt.Par.epsilon.*Gu(:,n)./(2*opt.Par.lambda*spc.P);
        opt.v(:,n) = opt.v(:,n)+opt.Par.epsilon.*Gv(:,n)./(2*opt.Par.lambda*spc.P);
    case {'Schirmer'}
        opt.u(:,n) = opt.u(:,n)+opt.Par.epsilon.*Gu(:,n)./(2*opt.Par.lambda);
        opt.v(:,n) = opt.v(:,n)+opt.Par.epsilon.*Gv(:,n)./(2*opt.Par.lambda);
        
end
end

%% MC_MadayTurinici


function opt = MC_MadayTurinici(spc,khr,opt)
% function opt = Do_MonCon(spc,khr,opt)
%
%   This script does Monotonic Convergent optimisation of Maday and
%   Turinici.


tic_Opt_Wall_Clock_Time = tic;

tic_One_Iteration = tic;
opt = Allocate_Variables(spc,khr,opt);
opt.ut = opt.u;
opt.vt = opt.v;
opt.k = 1;
opt = Iteration_GRAPE_Khaneja(spc,opt,0);

opt.M_T = opt.M_t(:,end);

Eff = Efficiency_GRAPE_Khaneja(spc,opt);


Pen = opt.Par.lambda*sum((opt.u(:)).^2+(opt.v(:)).^2)*opt.dt;
Fun  = Eff-Pen;


opt.dFun(opt.k) = Fun;


opt.Fun(opt.k) = Fun;
opt.Eff(opt.k) = Eff;
opt.Pen(opt.k) = Pen;


Print_GRAPE_Khaneja(opt.k,Fun,Eff,Pen);

opt = Progress_MC_MadayTurinici(spc,khr,opt);


opt.Durations(opt.k-1) = toc(tic_One_Iteration);


tic_Opt_Wall_Clock_Time = tic;

while opt.Go
    tic_One_Iteration = tic;
    
    opt = Iteration_MC_MadayTurinici(spc,khr,opt);
    
    opt.M_T = opt.M_t(:,end);
    
    opt.uo(:,:,opt.k) = opt.u;
    opt.vo(:,:,opt.k) = opt.v;
    
    
    Eff = Efficiency_GRAPE_Khaneja(spc,opt);
    
    Pen = opt.Par.lambda*sum((opt.u(:)).^2+(opt.v(:)).^2)*opt.dt;
    
    Fun  = Eff-Pen;
    
    Print_GRAPE_Khaneja(opt.k,Fun,Eff,Pen);
    
    opt.Fun(opt.k) = Fun;
    
    opt.Eff(opt.k) = Eff;
    
    opt.dFun(opt.k) = Fun-opt.Fun(opt.k-1);
    
    opt.Pen(opt.k) = Pen;
    
    opt = Progress_MC_MadayTurinici(spc,khr,opt);
    
    opt.Durations(opt.k-1) = toc(tic_One_Iteration);
end



if opt.ksafe < 1
    opt.ksafe = 1;
end
opt.uo(:,:,opt.ksafe+1:end) = [];
opt.vo(:,:,opt.ksafe+1:end) = [];
opt.Fun(opt.ksafe+1:end) = [];
opt.Eff(opt.ksafe+1:end) = [];
opt.Pen(opt.ksafe+1:end) = [];
opt.dFun(opt.ksafe+1:end) = [];

opt.Durations(opt.ksafe+1:end) = [];


opt.Opt_Wall_Clock_Time = toc(tic_Opt_Wall_Clock_Time);
end

function opt = Iteration_MC_MadayTurinici(spc,khr,opt)
% function [u,v,ut,vt,M_t,L_t] = One_MonCon_Iteration(spc,khr,opt)
%
%   This script does one MonCon-based iteration, i.e., a full propation forward
%   and a full propagation backward as described by Maday and Turinici.
%

N = opt.N;

[Gu,Gv] = RF_Update_Terms(spc,opt,1);

opt = Update_MC_MadayTurinici(spc,opt,Gu,Gv,1,'Forward');

[R11f,R12f,R13f,R21f,R22f,R23f,R31f,R32f,R33f] = Get_Rotator(opt.u,opt.v,opt.g,opt.w0,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.dt,1);

for n = 1:N-1
    
    opt.M_t(:,n+1) = Rotate(opt.M_t(:,n),R11f,R12f,R13f,R21f,R22f,R23f,R31f,R32f,R33f,'Forward');
    
    [Gu,Gv] = RF_Update_Terms(spc,opt,n+1);
    
    opt = Update_MC_MadayTurinici(spc,opt,Gu,Gv,n+1,'Forward');
    
    [R11f,R12f,R13f,R21f,R22f,R23f,R31f,R32f,R33f] = Get_Rotator(opt.u,opt.v,opt.g,opt.w0,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.dt,n+1);
    
end

opt.M_t(:,N+1) = Rotate(opt.M_t(:,N),R11f,R12f,R13f,R21f,R22f,R23f,R31f,R32f,R33f,'Forward');

[Gu,Gv] = RF_Update_Terms(spc,opt,N);

opt = Update_MC_MadayTurinici(spc,opt,Gu,Gv,N,'Backward');

[R11b,R12b,R13b,R21b,R22b,R23b,R31b,R32b,R33b] = Get_Rotator(opt.ut,opt.vt,opt.g,opt.w0,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.dt,N);


for n = N:-1:2
    
    opt.L_t(:,n) = Rotate(opt.L_t(:,n+1),R11b,R12b,R13b,R21b,R22b,R23b,R31b,R32b,R33b,'Backward');
    
    [Gu,Gv] = RF_Update_Terms(spc,opt,n-1);
    
    opt = Update_MC_MadayTurinici(spc,opt,Gu,Gv,n-1,'Backward');
    
    [R11b,R12b,R13b,R21b,R22b,R23b,R31b,R32b,R33b] = Get_Rotator(opt.ut,opt.vt,opt.g,opt.w0,opt.yx,opt.yy,opt.yz,spc.pTx,opt.sr,opt.si,opt.dt,n-1);
    
    
end

opt.L_t(:,1) = Rotate(opt.L_t(:,2),R11b,R12b,R13b,R21b,R22b,R23b,R31b,R32b,R33b,'Backward');



end

function opt = Progress_MC_MadayTurinici(spc,khr,opt)
% function opt = Evaluate_Progress(spc,khr,opt)
%
%   This script evaluates the overall progress of the optimization based on threshholds,
%   iteration limits etc.
%   For the MonCon-based algorithm the monotonicity is
%   is also a criterium.





% Monotonicity
idx = find(~isnan(opt.dFun));
if opt.dFun(max(idx)) > 0
    opt.Mono = 1;
else
    opt.Mono = 0;
    
    Display_Message(sprintf('Progress_MC_MadayTurinici: Monotonicity was lost: dFun < 0'),1);
end


% Convergence

opt.Conv = 0;

% Convergence, by Fun change


if ~isempty(opt.dFun)
    
    if abs(opt.dFun(opt.k)) > opt.TolFun
        opt.Conv = 0;
    elseif abs(opt.dFun(opt.k)) <= opt.TolFun
        
        if opt.k > 1
            
            
            if opt.Fun(opt.k) < opt.Fun(1)
                opt.Conv = 0;
                Display_Message(sprintf('Progress_MC_MadayTurinici: Convergence was obtained (dFun < TolFun), but Fun(k) < Fun(1), so we continue...'),1);
            else
                opt.Conv = 1;
                Display_Message(sprintf('Progress_MC_MadayTurinici: Convergence is obtained: dFun < TolFun'),1);
            end
            
        else
            opt.Conv = 1;
            Display_Message(sprintf('Progress_MC_MadayTurinici: Convergence is obtained: dFun < TolFun'),1);
            
            
        end
        
        
    end
    
else
    opt.Conv = 0;
end






% Convergence, by Fun value

if opt.Fun(opt.k) >= 1-opt.ObjLim
    opt.Bingo = 1;
    
    Display_Message(sprintf('Progress_MC_MadayTurinici: Reached desired performance: Fun >= 1-ObjLim'),1);
else
    opt.Bingo = 0;
end





if opt.k == opt.MaxIter+1
    opt.Done = 1;
    Display_Message(sprintf('Progress_MC_MadayTurinici: Reached Iteration limit MaxIter'),1);
else
    opt.Done = 0;
end





if abs(opt.Mono) == 1
    if opt.Bingo == 0
        if opt.Done == 0
            if opt.Conv == 0
                opt.Go = true;
                opt.ksafe = opt.k;
                opt.k = opt.k+1;
                
            else
                opt.Go = false;
                opt.ksafe = opt.k;
            end
            
        else
            opt.Go = false;
            opt.ksafe = opt.k;
        end
        
    else
        opt.Go = false;
        opt.ksafe = opt.k;
    end
else
    opt.Go = false;
    %     opt.ksafe = opt.ksafe - 1;
end

% save intermediate
if mod(opt.k,opt.ksaveintermediate) == 0
    
    
    tempsave = opt.Save; % take backup
    %
    opt.Save = '000000'; % reset
    opt.Save(2) = tempsave(2); % get initial intention back for data saving
    opt.Save(3) = tempsave(3); % get initial intention back for controls saving
    % thus, if intermediate saving of data and controls were intended it
    % will be saved. But other (bundle, pics, scripts etc.) will not, at
    % least not before the end if intended.
    try opt = Save_Job(spc,khr,opt);
    catch me; Display_Message(['Progress_MC_MadayTurinici: Save_Job ',me.message],2); end
    
    opt.Save = tempsave;
    
end


end

function opt = Update_MC_MadayTurinici(spc,opt,Gu,Gv,n,Dir)
% function opt = Update_Controls_4_MonCon(spc,khr,opt,Gu,Gv,n,Dir)
%
%   This updates the RF controls for MC

switch opt.Par.Grad
    
    
    case {'1st','type1'}
        
        Gu = Gu./(2*opt.Par.lambda*spc.P);
        Gv = Gv./(2*opt.Par.lambda*spc.P);
    otherwise
        Display_Message(sprintf('Update_MC_MadayTurinici: Currently only support 1st and type1 gradients'),1);
        Gu = 0;
        Gv = 0;
        
end


switch Dir
    
    case 'Forward'
        opt.u(:,n) = (1-opt.Par.delta)*opt.ut(:,n)+opt.Par.delta*Gu;
        opt.v(:,n) = (1-opt.Par.delta)*opt.vt(:,n)+opt.Par.delta*Gv;
        
        
        
    case 'Backward'
        opt.ut(:,n) = (1-opt.Par.eta)*opt.u(:,n)+opt.Par.eta*Gu;
        opt.vt(:,n) = (1-opt.Par.eta)*opt.v(:,n)+opt.Par.eta*Gv;
        
end

end

%% QN_LBFGS


function opt = QN_LBFGS(spc,khr,opt)
% function opt = Do_QNewton(spc,khr,opt)
%
%   This script does QNewton optimisation

global history opt2
tic_Opt_Wall_Clock_Time = tic;

Par.Grad = 'Schirmer';
Par.Constr.Cope = 3; % Coping

Par.Constr.RFpeak.Type = 'c';
Par.Constr.RFpeak.Lim = 1000;
Par.Constr.RFpeak.Unit = 'Hz';
Par.Constr.RFpeak.Power = 1;



opt.Par = blOCh__khr('Get_NewPar',[],[],[],opt.Par,Par);


opt.k = 1;

opt = Arrange_Constraints_4_QN_LBFGS(spc,khr,opt);
opt = Allocate_Variables(spc,khr,opt);
if opt.par_Ncores > 1
    delete(gcp)
    parpool('local', opt.par_Ncores)
end


x0 = real(Rearrange_controls(opt.u,opt.v));
x0(isinf(x0)) = 0;
x0(isnan(x0)) = 0;
xlim = ones(1,spc.pTx*opt.N*2)*opt.w1m;
temp = zeros(1,opt.N);
temp(opt.mon) = 1;

temp = repmat(temp,[1,spc.pTx*2]);
xlim = xlim.*temp;
xlim(isnan(xlim)) = eps;

opt.maxtest = max(abs(xlim));
opt2 = opt;
xlim = xlim./opt.maxtest;
x0 = x0./opt.maxtest;

options=optimset('TolX',opt.TolX,'TolFun',opt.TolFun,'Display','off','MaxIter',opt.MaxIter,'SubproblemAlgorithm','cg',...
    'MaxFunEvals',opt.MaxFunEvals,'Algorithm','interior-point','Hessian',{'lbfgs',10},'ScaleProblem','none',...
    'GradObj','on','DerivativeCheck',opt.deriv_check,'FinDiffType','central','LargeScale','off','GradConstr','on','OutputFcn',@Output_QN_LBFGS);% ,'OutputFcn',@DisplayIterFun
%,'OutputFcn',@Output_FMINs_Fcns
if opt.Par.Constr.Cope > 2
    if strcmp(Par.Constr.RFpeak.Type,'c')
        [x,fval,opt.exitflag,opt.output] =fmincon(@Iteration_QN_LBFGS,x0,[],[],[],[],-xlim,xlim,@Constr2_QN_LBFGS,options,spc,khr,opt);
    else
        
        if strcmp(Par.Constr.RFpeak.Type,'0')
            [x,fval,opt.exitflag,opt.output] =fmincon(@Iteration_QN_LBFGS,x0,[],[],[],[],-1*ones(size(x0))*Inf,ones(size(x0))*Inf,@Constr2_QN_LBFGS,options,spc,khr,opt);
        else
            [x,fval,opt.exitflag,opt.output] =fmincon(@Iteration_QN_LBFGS,x0,[],[],[],[],[],[],@Constr2_QN_LBFGS,options,spc,khr,opt);
        end
    end
    
    
    
else
    [x,fval,opt.exitflag,opt.output] =fmincon(@Iteration_QN_LBFGS,x0,[],[],[],[],[],[],[],options,spc,khr,opt);
    
end



Display_Message(opt.output.message,1);
opt.Go = true;
opt.k = 1;
opt.ksafe = history.safeiter;
opt.Par.Constr.RFpeak.Viol = history.Constr.RFpeak.Viol;
if opt.Par.Constr.Cope == 2
    if opt.Par.Constr.RFpeak.Viol > 0
        opt.ksafe = opt.ksafe-1;
        if opt.ksafe < 1
            opt.ksafe = 1;
        end
    end
    
end

opt.Fun = history.Fun(1:opt.ksafe);

opt.con.Peak   = zeros(1,opt.ksafe);

opt.con.Peak(1)   = history.con(1).Peak;
if isfield(history,'Durations')
    opt.Durations = history.Durations;
else
    opt.Durations = -1;
end
for p = 2:opt.ksafe
    
    opt.k = p;
    if isfield(opt,'maxtest')
        [opt.u,opt.v] = Rearrange_controls(history.x(p,:).*opt.maxtest,opt.N,spc.pTx);
    else
        [opt.u,opt.v] = Rearrange_controls(history.x(p,:),opt.N,spc.pTx);
        
    end
    opt.uo(:,:,p) = opt.u;
    opt.vo(:,:,p) = opt.v;
    
    
    opt.con.Peak(p)   = history.con(p).Peak;
    
    
    
end
opt.uo(:,:,opt.ksafe+1:end) = [];
opt.vo(:,:,opt.ksafe+1:end) = [];
opt.Durations(opt.ksafe+1:end) = [];
clearvars -global history
opt.Opt_Wall_Clock_Time = toc(tic_Opt_Wall_Clock_Time);
opt.Go = opt2.Go;
end

function [J,Grad] = Iteration_QN_LBFGS(x,spc,khr,opt)

global M_T Duration
tic_One_Iteration = tic;
x_prog = x.*opt.maxtest;



[opt.u,opt.v] = Rearrange_controls(x_prog,opt.N,spc.pTx);

opt = Iteration_GRAPE_Khaneja(spc,opt,0);

[Gu,Gv] = RF_Update_Terms(spc,opt,[]);


Gu = -opt.maxtest.*Gu./(spc.P);
Gv = -opt.maxtest.*Gv./(spc.P);
% [Gu(13),Gv(13)]

M_T = opt.M_t(:,end);

opt.Fun = Efficiency_QN_LBFGS(spc,opt,M_T);

J = -opt.Fun;

if nargout == 2
    Grad = Rearrange_controls(Gu,Gv);
end
Duration = toc(tic_One_Iteration);
end

function Print_QN_LBFGS(opt,Type,Fun,con)
% function Print_Evaluations(spc,khr,opt,Type,Fun,Eff,Pen,con,pen)
%
%   This script prints the relevant output.


if nargin == 2
    Type = 'all';
end


switch Type
    case 'all'
        
        A = opt.k-1;
        B = Fun;
        
        J = con.Peak;
        
        if opt.k == 1
            
            Msg = sprintf('\nIter  | Fun          | Constraint  |');
            Display_Message(Msg,1);
            Msg = sprintf(  '------+--------------+-------------+');
            
            Display_Message(Msg,1);
            Msg = sprintf(  '      |              |  Peak       |');
            Display_Message(Msg,1);
            
            
            
            
            Msg = sprintf('% 5i | %6.10g | %6.10g |',A,B,J);
            Display_Message(Msg,1);
            
            
            
        else
            Msg = sprintf('% 5i | %6.10g | %6.10g |',A,B,J);
            Display_Message(Msg,1);
            
        end
        
        
    case 'con'
        
        F = con.Peak;
        Msg = sprintf('      |              | %6.10g |',F);
        Display_Message(Msg,1);
end

end

function [con,Grad] = Constr1_QN_LBFGS(spc,khr,opt)

if opt.Par.Constr.Cope > 0
    
    x = Rearrange_controls(opt.u,opt.v);
    
    [c,ceq,dc,dceq,con]=Constr2_QN_LBFGS(x,spc,khr,opt);
    
    
    if nargout == 2
        Grad = dc;
    end
else
    
    con.Peak = 0;
    Grad = zeros(1,spc.pTx*opt.N*2);
end
end

function [c,ceq,dc,dceq,con]=Constr2_QN_LBFGS(x,spc,khr,opt)

ceq = [];
dceq   = [];
x = x.*opt.maxtest;

[rfx_,rfy_] = Rearrange_controls(x,opt.N,spc.pTx);

rf = complex(rfx_,rfy_);
rf_long = rf(:);


%% Power

% peak
if strcmp(opt.Par.Constr.RFpeak.Type,'1') || strcmp(opt.Par.Constr.RFpeak.Type,'c')
    
    if opt.Par.Constr.RFpeak.Power == 2
        c3 = abs(complex(x(1:opt.N*spc.pTx),x(opt.N*spc.pTx+1:end))).^2 - opt.Par.Constr.RFpeak.LimCor;
        
        dc3 = [2*diag(x(1:opt.N*spc.pTx));2*diag(x(opt.N*spc.pTx+1:end))].*opt.maxtest;
    else
        
        c3 = abs(complex(x(1:opt.N*spc.pTx),x(opt.N*spc.pTx+1:end))).^2 - opt.Par.Constr.RFpeak.LimCor.^2;
        dc3 = [2*diag(x(1:opt.N*spc.pTx));2*diag(x(opt.N*spc.pTx+1:end))].*opt.maxtest;
        
    end
    c3 = c3(:);
    
    
    
else
    c3 = [];
    dc3 = [];
end


%% Compile
if strcmp(opt.Par.Constr.RFpeak.Type,'c')
    c= [];
    dc = [];
else
    c= c3;
    dc = sparse(dc3);
end

%% Status
if spc.Print || narargout == 5
    
    if strcmp(opt.Par.Constr.RFpeak.Type,'1')  || strcmp(opt.Par.Constr.RFpeak.Type,'c')
        
        if opt.Par.Constr.RFpeak.Power == 2
            Peak = max(c3)+opt.Par.Constr.RFpeak.LimCor;
            Peak = Peak./opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
        else
            Peak = max(abs(rf_long));
            Peak = Peak./opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
        end
    else
        Peak = 0;
    end
    
    con.Peak = Peak;
    
    %% Print
    
    if strcmp(opt.Par.Constr.RFpeak.Type,'1')
        opt.con = con;
        Print_QN_LBFGS([],'con',[],con);
        %             Print_Evaluations([],[],[],'con',[],[],[],[],[],con,[]);
    end
    
    
end
end

function varargout = Rearrange_controls(varargin)
% function varargout = Rearrange_controls(varargin)
%
%   This script rearranges the shape of the control arrays
%

if length(varargin) == 2
    
    
    wx = varargin{1};
    wy = varargin{2};
    
    
    wx = permute(wx,[2,1]);
    wy = permute(wy,[2,1]);
    
    varargout{1} = [wx(:).',wy(:).'];
    
    % Input					< pTx x N >
    % permute				< N x pTx >
    % (:)					< N*pTx x 1 >
    % .'					< 1 x N*pTx >
    
    % | w(n=1,s=1) w(n=2,s=1)       ... w(n=N,s=1)|
    % | w(n=1,s=2) w(n=2,s=2)       ... w(n=N,s=2)|
    % |                       .                   |
    % |                       .                   |
    % |                       .                   |
    % | w(n=1,s=pTx) w(n=2,s=pTx) ... w(n=N,s=pTx)|
    %
    %				to
    %
    % | w(n=1,s=1) w(n=1,s=2) ... w(n=1,s=pTx)|
    % | w(n=2,s=1) w(n=2,s=2) ... w(n=1,s=pTx)|
    % |                   .                   |
    % |                   .                   |
    % |                   .                   |
    % | w(n=N,s=1) w(n=N,s=2) ... w(n=N,s=pTx)|
    %
    %				to
    %
    % | w(n=1,s=1)  |
    % | w(n=2,s=1)  |
    % |      .      |
    % |      .      |
    % |      .      |
    % | w(n=N,s=1)  |
    % | w(n=1,s=2)  |
    % | w(n=2,s=2)  |
    % |      .      |
    % |      .      |
    % |      .      |
    % | w(n=N,s=2)  |
    % |      .      |
    % |      .      |
    % |      .      |
    % |      .      |
    % |      .      |
    % |      .      |
    % | w(n=1,s=pTx)|
    % | w(n=2,s=pTx)|
    % |      .      |
    % |      .      |
    % |      .      |
    % | w(n=N,s=pTx)|
    %
    %				to
    %
    %  | w(n=1,s=1) w(n=2,s=1) ... w(n=N,s=1) w(n=1,s=2) w(n=2,s=2) ... w(n=N,s=2) ... w(n=N,s=pTx)|
    %
    % lastly: output = [wx,wy] in a row array.
    
    if nargout ~= 1
        error('Wrong number of outputs, expected one')
    end
    
elseif length(varargin) == 5
    wx = varargin{1};
    wy = varargin{2};
    gx = varargin{3};
    gy = varargin{4};
    gz = varargin{5};
    
    wx = permute(wx,[2,1]);
    wy = permute(wy,[2,1]);
    if nargout ~= 1
        error('Wrong number of outputs, expected one')
    end
    varargout{1} = [wx(:).',wy(:).',gx(:).',gy(:).',gz(:).'];
    
elseif length(varargin) == 4
    
    x = varargin{1};
    N = varargin{2};
    gchannel = varargin{3};
    redun = varargin{4};
    
    if nargout == gchannel
        
        
        
        
        if length(x) == 3*N
            gx = x(1:N);
            gy = x(N+1:2*N);
            gz = x(2*N+1:3*N);
            
            if nargout ~= 3
                error('Wrong number of outputs, expected three')
            end
            
            varargout{1} = gx(:).';
            varargout{2} = gy(:).';
            varargout{3} = gz(:).';
            
        else
            error('Input has wrong length')
        end
        
    
        
        
        
    end    
    
    
elseif length(varargin) == 3
    
    x = varargin{1};
    N = varargin{2};
    pTx = varargin{3};
    
    
    if nargout > 1
        
        
        
        
        if length(x) == 2*N*pTx
            wx = x(1:N*pTx);
            wy = x(N*pTx+1:end);
            
            if nargout ~= 2
                error('Wrong number of outputs, expected two')
            end
            
            varargout{1} = permute(reshape(wx,N,pTx),[2,1]);
            varargout{2} = permute(reshape(wy,N,pTx),[2,1]);
            
        elseif length(x) == 2*N*pTx+ 3*N;
            if nargout ~= 5
                error('Wrong number of outputs, expected five')
            end
            
            wx = x(1:N*pTx);
            wy = x(N*pTx+1      :2*N*pTx);
            
            
            varargout{1} = permute(reshape(wx,N,pTx),[2,1]);
            varargout{2} = permute(reshape(wy,N,pTx),[2,1]);
            varargout{3} = x(2*N*pTx+1    :2*N*pTx+N);
            varargout{4} = x(2*N*pTx+N+1  :2*N*pTx+2*N);
            varargout{5} = x(2*N*pTx+2*N+1:end);
        elseif length(x) < 2*N*pTx
            
            N = length(x)/2/pTx;
            
            wx = x(1:N*pTx);
            wy = x(N*pTx+1:end);
            
            if nargout ~= 2
                error('Wrong number of outputs, expected two')
            end
            
            varargout{1} = permute(reshape(wx,N,pTx),[2,1]);
            varargout{2} = permute(reshape(wy,N,pTx),[2,1]);
            
        else
            error('Input has wrong length')
        end
        
    else
        
        
        x_ = reshape(x,[pTx,2*N]);
        x = x_(:,1:N).';
        y = x_(:,N+1:end).';
        varargout{1} = [x(:).',y(:).'];
        
        
        
    end
end
end

function Eff = Efficiency_QN_LBFGS(spc,opt,M_T)

if nargin ~= 3
    M_T = opt.M_T;
end

A_x = spc.Md(1:3:end);
A_y = spc.Md(2:3:end);
A_z = spc.Md(3:3:end);

B_x = M_T(1:3:end);
B_y = M_T(2:3:end);
B_z = M_T(3:3:end);
Phi = sum(A_x.*B_x+A_y.*B_y+A_z.*B_z);


Eff = Phi/spc.P;
end


function opt = Arrange_Constraints_4_QN_LBFGS(spc,khr,opt)
% function [opt,Msg] = Arrange_Constraints_and_Penalties(spc,khr,opt)
%
%   This script arranges the various constraints and penalties

Msg = [];


if strcmp(opt.Par.Constr.RFpeak.Type,'1') || strcmp(opt.Par.Constr.RFpeak.Type,'c')
    
    switch opt.Par.Constr.RFpeak.Unit
        
        case 'W'
            
            opt.Par.Constr.RFpeak.Power = 2;
            
            if strcmp(spc.B1_nom_amp_unit,'T/V')
                
                opt.Par.Constr.RFpeak.Peak_unit_cor_factor = 8*50*khr.gamma^2/spc.f_B1_val^2; % (rad/s)^2/W
                
                opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
                opt.Par.Constr.RFpeak.UnitCor = '(rad/s)^2';
                
                
                opt.Par.Constr.RFpeak.Peak_unit_cor_factor_VV_2_W = khr.gamma^2*spc.B1_nom_amp_val^2./opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
            else
                Msg = sprintf('blOCh__opt: Par.Constr.RFpeak.Unit''s value (%s) is incompatible with B1 map unit (%s)',opt.Par.Constr.RFpeak.Unit,spc.B1_nom_amp_unit);
            end
            
        case 'V'
            
            opt.Par.Constr.RFpeak.Power = 1;
            if strcmp(spc.B1_nom_amp_unit,'T/V')
                
                opt.Par.Constr.RFpeak.Peak_unit_cor_factor = khr.gamma/spc.f_B1_val;
                
                opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
                opt.Par.Constr.RFpeak.UnitCor = 'rad/s';
            else
                Msg = sprintf('blOCh__opt: Par.Constr.RFpeak.Unit''s value (%s) is incompatible with B1 map unit (%s)',opt.Par.Constr.RFpeak.Unit,spc.B1_nom_amp_unit);
            end
        case 'rad/s'
            
            opt.Par.Constr.RFpeak.Power = 1;
            
            opt.Par.Constr.RFpeak.Peak_unit_cor_factor = 1;
            
            opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
            opt.Par.Constr.RFpeak.UnitCor = 'rad/s';
        case '(rad/s)^2'
            
            opt.Par.Constr.RFpeak.Power = 2;
            
            opt.Par.Constr.RFpeak.Peak_unit_cor_factor = 1;
            
            opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
            opt.Par.Constr.RFpeak.UnitCor = '(rad/s)^2';
        case 'T'
            
            opt.Par.Constr.RFpeak.Power = 1;
            
            opt.Par.Constr.RFpeak.Peak_unit_cor_factor = khr.gamma;
            opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
            opt.Par.Constr.RFpeak.UnitCor = 'rad/s';
        case 'Hz'
            
            opt.Par.Constr.RFpeak.Power = 1;
            opt.Par.Constr.RFpeak.Peak_unit_cor_factor = 2*pi;
            opt.Par.Constr.RFpeak.LimCor = opt.Par.Constr.RFpeak.Lim*opt.Par.Constr.RFpeak.Peak_unit_cor_factor;
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
%% for fmincon's, fmin's and xlim's

if ~isfield(opt,'w1m')
    opt.w1m = 2*pi*10e3;
    
end


% max_unit_length = max([numel(opt.LocSARLim{2}),])


Nfo = sprintf('    |Par.PeakLim = %1.2e %s = %1.2e %s\n',...
    opt.Par.Constr.RFpeak.Lim,opt.Par.Constr.RFpeak.Unit,opt.Par.Constr.RFpeak.LimCor,opt.Par.Constr.RFpeak.UnitCor);
Display_Message(Nfo,1);


Go = 0;




end

function [stop] = Output_QN_LBFGS(x,optimValues,state,spc,khr,opt)
% function [stop] = Output_QN_LBFGS(x,optimValues,state,spc,khr,opt)
%
%   This script is an output function for the QN_LBFGS functions


global M_T history Duration opt2

% 		disp('output')
stop = false;
% 		optimValues

% to begin with:
% it calls first init
% then it calls iter
% in the next iterations:
% it calls only iter.

switch state
    case 'init'
        %       disp('init')
        history.x = zeros(opt2.MaxIter,length(x)); % This repairs the bug with history.x being stored in matlabs weird background. So now you don't have to "clear all" if you change number of controls between two optimizations
        if isfield(opt2,'maxtest')
            history.x(1,:) = x;
        else
            history.x(1,:) = x;
        end
        history.fval(1)= optimValues.fval;
        
        % 		opt2.k = 1;
        history.safeiter  = optimValues.iteration+1;
        
        [opt2.u,opt2.v] = Rearrange_controls(x,opt2.N,spc.pTx);
        opt2.M_T = M_T;
        history.M_T{1} = M_T;
        
        
        Fun = Efficiency_QN_LBFGS(spc,opt,M_T);
        con = Constr1_QN_LBFGS(spc,khr,opt2);
        Print_QN_LBFGS(opt2,'all',Fun,con);
        %             Print_Evaluations([],[],opt2,'all',Fun,Eff,RMSExy,RMSExyz,Pen,con,pen);
        
        
        history.Fun(1) = Fun;
        
        history.con(1) = con;
        
        %       history.Durations(1) = 0;
        opt2.con = history.con;
        opt2.Fun = history.Fun;
        
        %         opt2 = Plot_Progress(spc,khr,opt2);
        opt2.dFun = [];
        opt2 = Progress_QN_LBFGS(spc,khr,opt2);
        opt2.dFun = 0; % abs(history.fval(iter))-abs(history.fval(iter-1));
        history.Constr.RFpeak.Viol = opt2.Par.Constr.RFpeak.Viol;
        
    case 'iter'
        %       disp('iter')
        if optimValues.iteration > 0
            iter = optimValues.iteration+1;
            
            opt2.k = iter;
            if isfield(opt2,'maxtest')
                history.x(iter,:) = x;
            else
                history.x(iter,:) = x;
            end
            history.fval(iter)= optimValues.fval;
            history.M_T{iter} = M_T;
            %         here = sum(history.M_T{iter})
            history.safeiter  = iter;
            
            
            
            opt2.M_T = M_T;
            
            [opt2.u,opt2.v] = Rearrange_controls(x,opt2.N,spc.pTx);
            
            Fun = Efficiency_QN_LBFGS(spc,opt,M_T);
            con = Constr1_QN_LBFGS(spc,khr,opt2);
            
            Print_QN_LBFGS(opt2,'all',Fun,con);
            %                 Print_Evaluations([],[],opt2,'all',Fun,Eff,RMSExy,RMSExyz,Pen,con,pen);
            
            
            history.Fun(iter) = Fun;
            
            history.con(iter) = con;
            opt2.con = history.con;
            opt2.Fun = history.Fun;
            %             opt2 = Plot_Progress(spc,khr,opt2);
            opt2.dFun = history.Fun(iter)-history.Fun(iter-1);%1; %abs(history.fval(iter))-abs(history.fval(iter-1));
            opt2.xintermed = x;
            opt2 = Progress_QN_LBFGS(spc,khr,opt2);
            history.dFun(iter) = opt2.dFun;
            history.Constr.RFpeak.Viol = opt2.Par.Constr.RFpeak.Viol;
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

function opt = Progress_QN_LBFGS(spc,khr,opt)
% function opt = Progress_QN_LBFGS(spc,khr,opt)
%
%   This script evaluates the overall progress of the optimization based on threshholds,
%   iteration limits etc.





% Constraints
if opt.Par.Constr.Cope > 0
    
    
    
    P_con_Lim_check = opt.Par.Constr.RFpeak.Lim-opt.con(opt.k).Peak;
    
    
    if P_con_Lim_check < 0 && ( strcmp(opt.Par.Constr.RFpeak.Type,'1')||strcmp(opt.Par.Constr.RFpeak.Type,'c'))
        
        opt.Par.Constr.RFpeak.Viol = 1;
        
        Msg = sprintf('Peak constraint violated');
        Display_Message(Msg,1);
    else
        opt.Par.Constr.RFpeak.Viol = 0;
    end
    
    
    
else
    opt.Par.Constr.RFpeak.Viol = 0;
end




switch opt.Par.Constr.Cope
    
    case 0
        
        opt.Go = true;
        opt.ksafe = opt.k;
        opt.k = opt.k+1;
        
    case 1
        
        opt.Go = true;
        opt.ksafe = opt.k;
        opt.k = opt.k+1;
        
    case 2
        
        if opt.Par.Constr.RFpeak.Viol == 0
            
            opt.Go = true;
            opt.ksafe = opt.k;
            opt.k = opt.k+1;
            
            
        else
            opt.Go = false;
            opt.ksafe = opt.k;
        end
        
    case 3
        
        opt.Go = true;
        opt.ksafe = opt.k;
        opt.k = opt.k+1;
        
        
        
end




% save intermediate
if mod(opt.k,opt.ksaveintermediate) == 0
    
    
    tempsave = opt.Save; % take backup
    %
    opt.Save = '000000'; % reset
    opt.Save(2) = tempsave(2); % get initial intention back for data saving
    opt.Save(3) = tempsave(3); % get initial intention back for controls saving
    % thus, if intermediate saving of data and controls were intended it
    % will be saved. But other (bundle, pics, scripts etc.) will not, at
    % least not before the end if intended.
    try opt = Save_Job(spc,khr,opt);
    catch me; Display_Message(['Progress_QN_LBFGS: Save_Job: ',me.message],2); end
    
    opt.Save = tempsave;
    
end
end
