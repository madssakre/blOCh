function varargout = blOCh__spc(Fun,Dim,M0,Md,MPS,varargin)
% function varargout = blOCh__spc(Fun,Dim,M0,Md,MPS,varargin)
%
%   This script develops the spatial dependent struct for blOCh.
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



if isempty(Fun)
    warning off
    spc.Status = 1;
    Msg = [];
    tempprint = 0;
    if spc.Status
        
        try
            valid_Dimensions = {...
                '0+1D',...
                '1DSI',...
                '1DAP',...
                '1DRL',...
                '2DAx',...
                '2DCo',...
                '2DSa',...
                '3D',...
                '1+1DSI',...
                '1+1DAP',...
                '1+1DRL',...
                '2+1DAx',...
                '2+1DCo',...
                '2+1DSa',...
                '3+1D'};
            valid_Mds = {'',[]};
            valid_M0s = {'',[],'MSL','MNI'};
            valid_MPSs = {'',[]};
            
            default_Print = 1;
            valid_Prints = [-1,0,1];
            
            
            default_pTx = 1;
            
            default_v0fac = 1;
            default_B1fac = 1;
            default_T1fac = 1;
            default_T2fac = 1;
            
            default_v0add = 0;
            default_B1add = 0;
            default_T1add = 0;
            default_T2add = 0;
            
            
            valid_Suppressions = {'off','on'};
            
            default_Suppression = valid_Suppressions{1};
            
            valid_Shows = [1,0];
            default_Show = 1;
            
            default_multiband = 1;
            
            default_M0flip = [0,0];
            default_Mdflip = [90,0];
            
            valid_Mpens = {'off','on'};
            default_Mpen = 'off';
            
            default_r0 = [0,0,0];
            default_L0 = [0.5e-2];
            
            switch Dim
                case {'1DSI','1+1DSI'}
                    default_TH = [5e-3];
                    default_c_TH = 0;
                    
                    default_L = [24e-2];
                    default_R = [20];
                    
                    Dimspc = [1,1];
                    
                    
                case {'1DAP','1+1DAP'}
                    default_TH = [5e-3];
                    default_c_TH = 0;
                    default_L = [24e-2];
                    default_R = [20];
                    
                    Dimspc = [1,1];
                    
                    
                    
                case {'1DRL','1+1DRL'}
                    default_TH = [5e-3];
                    default_c_TH = 0;
                    default_L = [24e-2];
                    default_R = [20];
                    
                    Dimspc = [1,1];
                    
                    
                case {'2DAx','2+1DAx'}
                    default_TH = [5e-3,5e-3];
                    default_c_TH = [0,0];
                    
                    
                    default_L = [24e-2,24e-2];
                    default_R = [20,20];
                    
                   Dimspc = [1,2];
                    
                case {'2DSa','2+1DSa'}
                    default_TH = [5e-3,5e-3];
                    default_c_TH = [0,0];
                    default_L = [24e-2,24e-2];
                    default_R = [20,20];
                    
                    Dimspc = [1,2];
                    
                case {'2DCo','2+1DCo'}
                    default_TH = [5e-3,5e-3];
                    default_c_TH = [0,0];
                    
                    default_L = [24e-2,24e-2];
                    default_R = [20,20];
                    
                    Dimspc = [1,2];
                    
                case {'0+1D'}
                    
                    default_TH = [];
                    default_c_TH = [];
                    
                    default_L = [];
                    default_R = [];
                    
                    
                    Dimspc = [0,0];
                    
                case {'3D','3+1D'}
                    
                    default_TH = [5e-3,5e-3,5e-3];
                    default_c_TH = [0,0,0];
                    
                    default_L = [24e-2,24e-2,24e-2];
                    default_R = [20,20,10];
                    
                    
                    Dimspc = [1,3];
            end
            
            if strcmp(Dim(2),'+')
                default_Rv = 20;
                default_BW = 1000;
                default_pBW = 100;
                default_c_pBW = 0;
                default_c_BW = 0;
                Rv_lim(1) = 1;
                Rv_lim(2) = Inf;
            else
                default_Rv = 1;
                default_BW = 0;
                default_pBW = 0;
                default_c_pBW = 0;
                default_c_BW = 0;
                Rv_lim(1) = 0;
                Rv_lim(2) = 2;
            end
            
            switch Dim(1:2)
                case {'1D','1+','0+'}
                    default_Profile = [];
                case {'2D','2+'}
                    default_Profile = {[-1 1 0 0],[0 0 -1 1]};
                case {'3D','3+'}
                    default_Profile = {[-1 1 0 0 0 0],[0 0 -1 1 0 0],[0 0 0 0 -1 1]};
                    
            end
        catch me
            Msg = ['blOCh__spc: ',me.message,' First argument to blOCh__spc must be BUG'];
        end
        
        
        
        p = inputParser;
        spc.Status = 0;
        try p.addRequired('Fun', @(x)isempty(x));
        try p.addRequired('Dim',@(x)any(strcmpi(x,valid_Dimensions)));
            try p.addRequired('M0',@(x)any(strcmpi(x,valid_M0s))||exist(x,'file'));
                try p.addRequired('Md',@(x)any(strcmpi(x,valid_Mds))||exist(x,'file'));
                    try p.addRequired('MPS',@(x)any(strcmpi(x,valid_MPSs))||exist(x,'file'));
                        try p.addParamValue('Show',default_Show,@(x)any(valid_Shows));
                            try p.addParamValue('L',default_L,@(x)validateattributes(x,{'numeric'},{'size',Dimspc,'real','>',0}));
                                try p.addParamValue('R',default_R,@(x)validateattributes(x,{'numeric'},{'size',Dimspc,'integer','nonnegative','finite'}));
                                    try p.addParamValue('r0',default_r0,@(x)validateattributes(x,{'numeric'},{'size',[1,3],'real','finite'}));
                                        try p.addParamValue('L0',default_L0,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>',0}));
                                            try p.addParamValue('BW',default_BW,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>',0}));
                                                try p.addParamValue('Rv',default_Rv,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'integer','>',Rv_lim(1),'<',Rv_lim(2)})); %%
                                                    try p.addParamValue('pBW',default_pBW,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','nonnegative','finite','>',0}));
                                                        try p.addParamValue('c_pBW',default_c_pBW,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','finite'}));
                                                            try p.addParamValue('c_BW',default_c_BW,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','finite'}));
                                                                try p.addParamValue('c_TH',default_c_TH,@(x)validateattributes(x,{'numeric'},{'size',Dimspc,'real','finite'}));
                                                                    try p.addParamValue('TH',default_TH,@(x)validateattributes(x,{'numeric'},{'size',Dimspc,'real','finite','>',0}));
                                                                        try p.addParamValue('Print',default_Print,@(x)any(valid_Prints));
                                                                            try p.addParamValue('pTx',default_pTx,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real','integer','positive'}));
                                                                                try p.addParamValue('Mdflip',default_Mdflip,@(x)validateattributes(x,{'numeric'},{'size',[1,2],'finite','real'}));
                                                                                    try p.addParamValue('Mpen',default_Mpen,@(x)any(strcmpi(x,valid_Mpens)));
                                                                                        
                                                                                        try p.addParamValue('Profiles',default_Profile,@(x)validateattributes(x,{'cell'},{'row'}));
                                                                                            try p.addParamValue('v0fac',default_v0fac,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real'}));
                                                                                                try p.addParamValue('B1fac',default_B1fac,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real'}));
                                                                                                    try p.addParamValue('T1fac',default_T1fac,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real'}));
                                                                                                        try p.addParamValue('T2fac',default_T2fac,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real'}));
                                                                                                            try p.addParamValue('v0add',default_v0add,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real'}));
                                                                                                                try p.addParamValue('B1add',default_B1add,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real'}));
                                                                                                                    try p.addParamValue('T1add',default_T1add,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real'}));
                                                                                                                        try p.addParamValue('T2add',default_T2add,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'finite','real'}));
                                                                                                                            try p.addParamValue('multiband',default_multiband,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>=',1}));
                                                                                                                                try p.addParamValue('M0flip',default_M0flip,@(x)validateattributes(x,{'numeric'},{'size',[1,2],'finite','real','>=',-180,'<=',180}));
                                                                                                                                    try p.addParamValue('Suppression',default_Suppression,@(x)any(strcmpi(x,valid_Suppressions)));
                                                                                                                                        spc.Status = 1;
                                                                                                                                    catch  me;Msg = ['blOCh__spc: ',me.message,' Suppression: on/off']; end
                                                                                                                                catch  me;Msg = ['blOCh__spc: ',me.message,' M0flip: flip-angle and -phase of initial magnetization, [FA(deg.),FP(deg.)]']; end
                                                                                                                            catch  me;Msg = ['blOCh__spc: ',me.message,' multiband: number of bands in multiband excitation']; end
                                                                                                                        catch  me;Msg = ['blOCh__spc: ',me.message,' T2add: T2map offset'];end
                                                                                                                    catch  me;Msg = ['blOCh__spc: ',me.message,' T1add: T1map offset'];end
                                                                                                                catch  me;Msg = ['blOCh__spc: ',me.message,' B1add: B1map offset'];end
                                                                                                            catch  me;Msg = ['blOCh__spc: ',me.message,' v0add: v0map offset'];end
                                                                                                        catch  me;Msg = ['blOCh__spc: ',me.message,' T2fac: T2map factor'];end
                                                                                                    catch  me;Msg = ['blOCh__spc: ',me.message,' T1fac: T1map factor'];end
                                                                                                catch  me;Msg = ['blOCh__spc: ',me.message,' B1fac: B1map factor'];end
                                                                                            catch  me;Msg = ['blOCh__spc: ',me.message,' v0fac: v0map factor'];end
                                                                                        catch  me;Msg = ['blOCh__spc: ',me.message,' Profiles: profile coordinates for displaying'];end
                                                                                    catch  me;Msg = ['blOCh__spc: ',me.message,' Mpen: on/off'];end
                                                                                    
                                                                                catch  me;Msg = ['blOCh__spc: ',me.message,' Mdflip: flip-angle and -phase of desired magnetization, [FA(deg),FP(deg)]'];end
                                                                            catch  me;Msg = ['blOCh__spc: ',me.message,' pTx: number of Tx channels'];end
                                                                        catch  me;Msg = ['blOCh__spc: ',me.message,' Print: flag related to what is being printed'];end
                                                                    catch  me;Msg = ['blOCh__spc: ',me.message, ' TH: slice thickness, eg. in 1D slice selection, or sidelength in 2D (3D) square (box) excitation'];end
                                                                catch  me;Msg = ['blOCh__spc: ',me.message,' c_TH: center of the ROI (m), eg. in 1D slice selection, or 2D (3D) square (box) excitation'];end
                                                            catch  me;Msg = ['blOCh__spc: ',me.message,' c_BW: central frequency of the band width (Hz), relative to resonance 0 Hz'];end
                                                        catch  me;Msg = ['blOCh__spc: ',me.message,' c_pBW: central frequency of the pass band (Hz), relative to resonance 0 Hz'];end
                                                    catch  me;Msg = ['blOCh__spc: ',me.message,' pBW: pass bandwidth (Hz)'];end
                                                catch  me;Msg = ['blOCh__spc: ',me.message,' Rv: spectral resolution (#)'];end
                                            catch  me;Msg = ['blOCh__spc: ',me.message,' BW: spectral bandwidth (Hz)'];end
                                        catch  me;Msg = ['blOCh__spc: ',me.message,' L0: cut-out width (m), eg. 3D data to 1D or 2D spc'];end
                                    catch  me;Msg = ['blOCh__spc: ',me.message,' r0: cut-out position [x0,y0,z0] (m), ie., center of nD FOX data cut-out from mD input data (m>n)'];end
                                    
                                catch  me;
                                    if Dimspc == 1
                                        Msg = ['blOCh__spc: ',me.message,' R: spatial resolution (#)'];
                                    elseif Dimspc == 2
                                        Msg = ['blOCh__spc: ',me.message,' R: spatial resolution, [(#),(#)]'];
                                    elseif Dimspc == 3
                                        Msg = ['blOCh__spc: ',me.message,' R: spatial resolution, [(#),(#),(#)]'];
                                    end
                                end
                            catch  me;
                                if Dimspc == 1
                                    Msg = ['blOCh__spc: ',me.message,' L: FOX length (m)'];
                                elseif Dimspc == 2
                                    Msg = ['blOCh__spc: ',me.message,' L: FOX length, [(m),(m)]'];
                                elseif Dimspc == 3
                                    Msg = ['blOCh__spc: ',me.message,' L: FOX length, [(m),(m),(m)]'];
                                end
                                
                            end
                        catch  me;Msg = ['blOCh__spc: ',me.message,' Show: boolean to show spc or not'];end
                    catch  me;Msg = ['blOCh__spc: ',me.message,' MPS: Maps choice: '''', [], or a filename of a MPS file.'];end
                catch  me;Msg = ['blOCh__spc: ',me.message,' Md: Desired magnetization choice: '''', [], or a filename of a Md file or image file.'];end
            catch  me;Msg = ['blOCh__spc: ',me.message,' M0: Initial Magnetization choice: '''', [], or a filename of a M0 file or image file.'];end
            
        catch  me;Msg = ['blOCh__spc: ',me.message,' Dim: Dimension of OC problem BUG'];end
        catch  me;Display_Message(['blOCh__spc: Fun: ',me.message],2);end
        
        if spc.Status
            
            % -----------------
            try
                p.parse(Fun,Dim,M0,Md,MPS,varargin{:});
                
                spc.Print = p.Results.Print;
                
                if spc.Print == -1
                    tempprint = -1;
                    E = dbstack('-completenames');
                    Name = E(2,1).name;
                    
                    spc.Print = fopen([Name,'.out'],'a');
                else
                    tempprint = 0;
                end
                
            catch me
                Msg = ['blOCh__spc: ',me.message];
                spc.Status = 0;
            end
            
            
            
            if ~isempty(Msg)
                Display_Message(Msg,1); % if spc.Print == 0 and this happens in executable it could crash it.
                spc.Status = 0;
            end
        else
            Display_Message(Msg,1);
        end
    end
    
    if spc.Status
        spc.Print = p.Results.Print;
        spc.pTx = p.Results.pTx;
        spc.Dim = p.Results.Dim;
        
        
        spc.L0 = p.Results.L0;
        spc.r0 = p.Results.r0;
        
        switch Dim
            case {'1DSI','1+1DSI'}
                spc.L = [0,0,p.Results.L];
                spc.R = [1,1,p.Results.R];
                spc.c_TH = [0,0,p.Results.c_TH];
                spc.TH = [0,0,p.Results.TH];
                spc.D = [...
                    spc.r0(1)-spc.L0/2,spc.r0(2)-spc.L0/2,spc.r0(3)-spc.L(3)/2;...
                    spc.r0(1)+spc.L0/2,spc.r0(2)+spc.L0/2,spc.r0(3)+spc.L(3)/2];
            case {'1DAP','1+1DAP'}
                spc.L = [0,p.Results.L,0];
                spc.R = [1,p.Results.R,1];
                spc.c_TH = [0,p.Results.c_TH,0];
                spc.TH = [0,p.Results.TH,0];
                spc.D = [...
                    spc.r0(1)-spc.L0/2,spc.r0(2)-spc.L(2)/2,spc.r0(3)-spc.L0/2;...
                    spc.r0(1)+spc.L0/2,spc.r0(2)+spc.L(2)/2,spc.r0(3)+spc.L0/2];
            case {'1DRL','1+1DRL'}
                spc.L = [p.Results.L,0,0];
                spc.R = [p.Results.R,1,1];
                spc.c_TH = [p.Results.c_TH,0,0];
                spc.TH = [p.Results.TH,0,0];
                spc.D = [...
                    spc.r0(1)-spc.L(1)/2,spc.r0(2)-spc.L0/2,spc.r0(3)-spc.L0/2;...
                    spc.r0(1)+spc.L(1)/2,spc.r0(2)+spc.L0/2,spc.r0(3)+spc.L0/2];
            case {'2DAx','2+1DAx'}
                spc.L = [p.Results.L,0];
                spc.R = [p.Results.R,1];
                spc.c_TH = [p.Results.c_TH,0];
                spc.TH = [p.Results.TH,0];
                spc.D = [...
                    spc.r0(1)-spc.L(1)/2,spc.r0(2)-spc.L(2)/2,spc.r0(3)-spc.L0/2;...
                    spc.r0(1)+spc.L(1)/2,spc.r0(2)+spc.L(2)/2,spc.r0(3)+spc.L0/2];
            case {'2DSa','2+1DSa'}
                spc.L = [0,p.Results.L];
                spc.R = [1,p.Results.R];
                spc.c_TH = [0,p.Results.c_TH];
                spc.TH = [0,p.Results.TH];
                spc.D = [...
                    spc.r0(1)-spc.L0/2,spc.r0(2)-spc.L(2)/2,spc.r0(3)-spc.L(3)/2;...
                    spc.r0(1)+spc.L0/2,spc.r0(2)+spc.L(2)/2,spc.r0(3)+spc.L(3)/2];
            case {'2DCo','2+1DCo'}
                spc.L = [p.Results.L(1),0,p.Results.L(2)];
                spc.R = [p.Results.R(1),1,p.Results.R(2)];
                spc.c_TH = [p.Results.c_TH(1),0,p.Results.c_TH(2)];
                spc.TH = [p.Results.TH(1),0,p.Results.TH(2)];
                spc.D = [...
                    spc.r0(1)-spc.L(1)/2,spc.r0(2)-spc.L0/2,spc.r0(3)-spc.L(3)/2;...
                    spc.r0(1)+spc.L(1)/2,spc.r0(2)+spc.L0/2,spc.r0(3)+spc.L(3)/2];
            case {'3D','3+1D'}
                spc.L = p.Results.L;
                spc.R = p.Results.R;
                spc.TH = p.Results.TH;
                spc.c_TH = p.Results.c_TH;
                spc.D = [spc.r0(1)-spc.L(1)/2,spc.r0(2)-spc.L(2)/2,spc.r0(3)-spc.L(3)/2;...
                    spc.r0(1)+spc.L(1)/2,spc.r0(2)+spc.L(2)/2,spc.r0(3)+spc.L(3)/2];
                
            case {'0+1D'}
                spc.L = 0;
                spc.R = [1,1,1];
                spc.TH = 0;
                spc.c_TH = 0;
                spc.D = [0,0,0;...
                    0,0,0];
                
                
            otherwise
                
        end
        
        if strcmp(Dim(2),'+')
            spc.BW = p.Results.BW;
            spc.c_pBW = p.Results.c_pBW;
            spc.c_BW = p.Results.c_BW;
            spc.pBW = p.Results.pBW;
            spc.Rv = p.Results.Rv;
        else
            if p.Results.BW > 0
                Nfo = ['blOCh__spc: BW is obsolete.'];
                Display_Message(Nfo,spc.Print);
            end
            if p.Results.Rv > 1
                Nfo = ['blOCh__spc: Rv is obsolete.'];
                Display_Message(Nfo,spc.Print);
            end
            if p.Results.c_pBW ~= 0
                Nfo = ['blOCh__spc: c_pBW is obsolete.'];
                Display_Message(Nfo,spc.Print);
            end
            if p.Results.pBW ~= 0
                Nfo = ['blOCh__spc: pBW is obsolete.'];
                Display_Message(Nfo,spc.Print);
            end
            spc.BW = 0;
            spc.c_BW = 0;
            spc.c_pBW = 0;
            spc.pBW = 0;
            spc.Rv = 1;
        end
        
        
        spc.Dv = [spc.c_BW-spc.BW/2;spc.c_BW+spc.BW/2];
        
        spc.v0fac = p.Results.v0fac;
        spc.B1fac = p.Results.B1fac;
        spc.T1fac = p.Results.T1fac;
        spc.T2fac = p.Results.T2fac;
        spc.v0add = p.Results.v0add;
        spc.B1add = p.Results.B1add;
        spc.T1add = p.Results.T1add;
        spc.T2add = p.Results.T2add;
        
        spc.MPS = p.Results.MPS;
        spc.M0 = p.Results.M0;
        spc.Md = p.Results.Md;
        
        spc.Profiles = p.Results.Profiles;
        spc.Mdflip = p.Results.Mdflip;
        spc.Mpen = p.Results.Mpen;
        spc.M0flip = p.Results.M0flip;
        
        spc.Suppression = p.Results.Suppression;
        
        spc.multiband = p.Results.multiband;
        spc.Show = p.Results.Show;
        spc.Profiles = p.Results.Profiles;
        Msg = [];
        
        
        spc.Status = 1;
        
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            spc.Status = 0;
        end
    end
    %% Check M0
    
    if spc.Status
        
        [m0,Msg] = Load_m0(spc);
        
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            spc.Status = 0;
        else
            spc.Status = 1;
        end
    end
    %% Check Md
    if spc.Status
        
        spc.Status = 0;
        [md,Msg] = Load_md(spc);
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            spc.Status = 0;
        else
            spc.Status = 1;
        end
    end
    
    %% Check MPS
    
    if spc.Status
        spc.Status = 0;
        [mps,Msg] = Load_mps(spc);
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            spc.Status = 0;
        else
            spc.Status = 1;
        end
    end
    %% Check orientations
    
    
    if spc.Status
        spc.Status = 0;
        Msg = check_Orientations_Intel(spc,m0,md,mps);
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            spc.Status = 0;
        else
            spc.Status = 1;
        end
    end
    
    %% Check dimensions
    
    if spc.Status
        spc.Status = 0;
        Msg = check_Dimensions(spc.D,md.D,m0.D,mps.D);
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            spc.Status = 0;
        else
            spc.Status = 1;
        end
    end
    
    %% Check other parms
    if spc.Status
        Msg = '';
        spc.Status = 0;
        switch Dimspc(2)
            
            case 1
            case 2
                
                
                
                
                for n = 1:length(spc.Profiles)
                    if ~isempty(spc.Profiles{n})
                        try validateattributes(spc.Profiles{n},{'double'},{'size',[1,4],'<=',1,'>=',-1})
                            
                        catch me
                            Msg = ['Profiles{',num2str(n),'}: ',me.message,' Profiles: numeric 1-by-4 array of profile coordinates e.g. [xi xf yi yf] = [-1 1 0 0.3] (norm. units).'];
                        end
                    end
                end
                
            case 3
                
        end
        
        
        
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            spc.Status = 0;
        else
            
            spc.Status = 1;
            
        end
    end
    
    
    
    
    if spc.Status
        clear default* valid*
        
        spc.Status = 0;
        
        try
            [v0map,Xi,Yi,Zi,Vi] = Morph(mps.Dim,spc.Dim,mps.v0map,mps.D,mps.Dv,spc.L,spc.r0,spc.L0,spc.R,spc.Rv,spc.c_BW,spc.BW);
            [B1map,Xi,Yi,Zi,Vi] = Morph(mps.Dim,spc.Dim,mps.B1map,mps.D,mps.Dv,spc.L,spc.r0,spc.L0,spc.R,spc.Rv,spc.c_BW,spc.BW);
            [T1map,Xi,Yi,Zi,Vi] = Morph(mps.Dim,spc.Dim,mps.T1map,mps.D,mps.Dv,spc.L,spc.r0,spc.L0,spc.R,spc.Rv,spc.c_BW,spc.BW);
            [T2map,Xi,Yi,Zi,Vi] = Morph(mps.Dim,spc.Dim,mps.T2map,mps.D,mps.Dv,spc.L,spc.r0,spc.L0,spc.R,spc.Rv,spc.c_BW,spc.BW);
            [M0,Xi,Yi,Zi,Vi]    = Morph( m0.Dim,spc.Dim,m0.M0,m0.D,m0.Dv,spc.L,spc.r0,spc.L0,spc.R,spc.Rv,spc.c_BW,spc.BW);
            if strcmp(spc.Mpen,'on')
                [Mp1,Xi,Yi,Zi,Vi]    = Morph( md.Dim,spc.Dim,md.MP1,md.D,md.Dv,spc.L,spc.r0,spc.L0,spc.R,spc.Rv,spc.c_BW,spc.BW);
                [Mp2,Xi,Yi,Zi,Vi]    = Morph( md.Dim,spc.Dim,md.MP2,md.D,md.Dv,spc.L,spc.r0,spc.L0,spc.R,spc.Rv,spc.c_BW,spc.BW);
            end
            [Md,Xi,Yi,Zi,Vi]    = Morph( md.Dim,spc.Dim,md.Md,md.D,md.Dv,spc.L,spc.r0,spc.L0,spc.R,spc.Rv,spc.c_BW,spc.BW);
            
            
            
            
            Subject = sum(abs(B1map),5);
            
            v0map = v0map.*spc.v0fac + spc.v0add;
            B1map = B1map.*spc.B1fac + spc.B1add;
            T1map = T1map.*spc.T1fac + spc.T1add;
            T2map = T2map.*spc.T2fac + spc.T2add;
            
            
            idx_Subject = find(Subject);
            
            X_subject = Xi(idx_Subject);
            Y_subject = Yi(idx_Subject);
            Z_subject = Zi(idx_Subject);
            V_subject = Vi(idx_Subject);
            
            M0x = M0(:,:,:,:,1);
            M0y = M0(:,:,:,:,2);
            M0z = M0(:,:,:,:,3);
            
            M0x_subject = M0x(idx_Subject);
            M0y_subject = M0y(idx_Subject);
            M0z_subject = M0z(idx_Subject);
            
            M0 = [M0x_subject(:)';M0y_subject(:)';M0z_subject(:)'];
            
            Mdx = Md(:,:,:,:,1);
            Mdy = Md(:,:,:,:,2);
            Mdz = Md(:,:,:,:,3);
            
            Mdx_subject = Mdx(idx_Subject);
            Mdy_subject = Mdy(idx_Subject);
            Mdz_subject = Mdz(idx_Subject);
            
            Md = [Mdx_subject(:)';Mdy_subject(:)';Mdz_subject(:)'];
            
            
            if strcmp(spc.Mpen,'on')
                Mpx1 = Mp1(:,:,:,:,1);
                Mpy1 = Mp1(:,:,:,:,2);
                Mpz1 = Mp1(:,:,:,:,3);
                
                Mpx_subject = Mpx1(idx_Subject);
                Mpy_subject = Mpy1(idx_Subject);
                Mpz_subject = Mpz1(idx_Subject);
                
                Mp1 = [Mpx_subject(:)';Mpy_subject(:)';Mpz_subject(:)'];
                
                Mpx2 = Mp2(:,:,:,:,1);
                Mpy2 = Mp2(:,:,:,:,2);
                Mpz2 = Mp2(:,:,:,:,3);
                
                Mpx_subject = Mpx2(idx_Subject);
                Mpy_subject = Mpy2(idx_Subject);
                Mpz_subject = Mpz2(idx_Subject);
                
                Mp2 = [Mpx_subject(:)';Mpy_subject(:)';Mpz_subject(:)'];
                
            end
            
            
            
            T1map_subject = T1map(idx_Subject);
            T2map_subject = T2map(idx_Subject);
            v0map_subject = v0map(idx_Subject);
            
            
            
            
            
            B1map_subject = zeros(numel(idx_Subject),size(B1map,5));
            for n = 1:size(B1map,5)
                temp = B1map(:,:,:,:,n);
                temp2 = temp(idx_Subject);
                B1map_subject(:,n) = temp2(:);
            end
            
            %% populate spc struct with listed maps etc.
            
            
            spc.idxtot = idx_Subject;
            
            spc.w0map = v0map_subject(:).*2*pi;
            spc.T1map = T1map_subject(:);
            spc.T2map = T2map_subject(:);
            spc.X = X_subject(:);
            spc.Y = Y_subject(:);
            spc.Z = Z_subject(:);
            spc.V = V_subject(:);
            spc.M0 = M0(:);
            spc.Md = Md(:);
            if strcmp(spc.Mpen,'on')
                spc.Mp1 = Mp1(:);
                spc.Mp2 = Mp2(:);
            end
            spc.P = numel(idx_Subject);
            Mdxyc = complex(spc.Md(1:3:end),spc.Md(2:3:end));
            
            spc.idxinroi = find(abs(Mdxyc));
            spc.idxoutroi = find(~abs(Mdxyc));
            
            Mdxy0 = [spc.Md(1:3:end).';spc.Md(2:3:end).';zeros(1,spc.P)];
            spc.Mdxyc = Mdxyc(:);
            spc.Mdxy0 = Mdxy0(:);
            % about B1map
            
            spc.f_B1_val = 1/mps.B1_nom_amp_val;
            
            spc.B1map = B1map_subject;
            spc.B1_nom_amp_val = mps.B1_nom_amp_val;
            spc.B1_nom_amp_unit = mps.B1_nom_amp_unit;
            spc.B1_type = mps.B1_type;
            switch spc.B1_type
                case 'absolute_valued'
                case 'relative_valued'
                    spc.B1_ref_voltage = mps.B1_ref_voltage;
                    
            end
            %% populate temp spc struct for Show_spc with gridded maps etc.
            if spc.Show > 0
                spctemp = spc;
                spctemp = rmfield(spctemp,'Status');
                spctemp.w0map = v0map.*2*pi;
                
                spctemp.B1map = B1map;
                spctemp.T1map = T1map;
                spctemp.T2map = T2map;
                spctemp.X = Xi;
                spctemp.Y = Yi;
                spctemp.Z = Zi;
                spctemp.V = Vi;
                spctemp.M0 = cat(5,M0x,M0y,M0z);
                spctemp.Md = cat(5,Mdx,Mdy,Mdz);
                
                if strcmp(spc.Mpen,'on')
                    spctemp.Mp1 = cat(5,Mpx1,Mpy1,Mpz1);
                    spctemp.Mp2 = cat(5,Mpx2,Mpy2,Mpz2);
                end
            end
        catch me
            [ST,I] = dbstack('-completenames');
            Msg = ['blOCh__spc: ',me.message,' Something wrong compiling spc'];
        end
        
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            spc.Status = 0;
        else
            
            spc.Status = 1;
        end
        
    else
        spc.Status = 0;
    end
    
    if spc.Status
        if spc.Show > 0
                     try
            spc.fig = Show_spc(spctemp);
            
                    spc = rmfield(spc,'spctemp');
                     catch me
                        Msg = ['blOCh__spc: ',me.message,' Something is wrong with Show'];
                     end
            
        end
    end
    if tempprint == -1
        fclose(spc.Print);
        spc.Print = -1;
    end
    
    varargout{1} = spc;
    
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
            Display_Message(['blOCh__spc: Fun: This switch needs to be expanded'],2)
            
    end
    
    varargout = v;
    
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



function MD = Spect_Expand_MD(MD,MDnoroi,Rv,pBW,c_pBW,BW,c_BW)
% function MD = Spect_Expand_MD(MD,MDnoroi,Rv,pBW,c_pBW,BW)
%
%   This scripts expands the MD array to include a spectral dimension
%
Dv = [c_BW-BW/2;c_BW+BW/2];

fl = (c_pBW-pBW/2);
fr = (c_pBW+pBW/2);
f = linspace(Dv(1),Dv(2),Rv);

[Mnl,idx1] = min(abs(f-fl));
[Mnr,idx2] = min(abs(f-fr));



if idx1 < 2
    idx1 = 2;
end
if idx2 > Rv
    idx2 = Rv;
end
if idx2 < 1
    idx2 = 1;
end
if idx1 > Rv
    idx1 = Rv;
end
temp1 = repmat(MDnoroi,[1,1,1,idx1-1,1]);
temp2 = repmat(MD,[1,1,1,idx2-idx1+1,1]);
temp3 = repmat(MDnoroi,[1,1,1,Rv-idx2,1]);

MD = cat(4,temp1,temp2,temp3);
end




function Out = Assign_FA(MD11,MD22,MD11noroi,MD22noroi,flip)
% function Out = Assign_FA(MD11,MD22,MD11noroi,MD22noroi,flip)
%
%   This script assigns a flip angle to the Md array
%

MD11x = MD11(:,:,:,:,1);
MD11y = MD11(:,:,:,:,2);
MD11z = MD11(:,:,:,:,3);
MD11xnoroi = MD11noroi(:,:,:,:,1);
MD11ynoroi = MD11noroi(:,:,:,:,2);
MD11znoroi = MD11noroi(:,:,:,:,3);

idx11_x = find(MD11x-MD11xnoroi);
idx11_y = find(MD11y-MD11ynoroi);
idx11_z = find(MD11z-MD11znoroi);

MD22x = MD22(:,:,:,:,1);
MD22y = MD22(:,:,:,:,2);
MD22z = MD22(:,:,:,:,3);
MD22xnoroi = MD22noroi(:,:,:,:,1);
MD22ynoroi = MD22noroi(:,:,:,:,2);
MD22znoroi = MD22noroi(:,:,:,:,3);

idx22_x = find(MD22x-MD22xnoroi);
idx22_y = find(MD22y-MD22ynoroi);
idx22_z = find(MD22z-MD22znoroi);


idxmem_x = ismember(idx11_x,idx22_x);
idxmem_y = ismember(idx11_y,idx22_y);
idxmem_z = ismember(idx11_z,idx22_z);

idx_x = idx11_x(idxmem_x);
idx_y = idx11_y(idxmem_y);
idx_z = idx11_z(idxmem_z);

MDx = zeros(size(MD11,1),size(MD11,2),size(MD11,3),1,1);
MDx(idx_x) = sind(flip(1)).*cosd(flip(2));
MDy = zeros(size(MD11,1),size(MD11,2),size(MD11,3),1,1);
MDy(idx_y) = sind(flip(1)).*sind(flip(2));
MDz = ones(size(MD11,1),size(MD11,2),size(MD11,3),1,1);
MDz(idx_z) = cosd(flip(1));

Out = cat(5,MDx,MDy,MDz);

end
function Msg = check_Dimensions(D_spc,D_md,D_m0,D_mps)
% function Msg = check_Dimensions(D_spc,D_md,D_m0,D_mps)
%
%   This script checks the dimensions
%


Msg = {};



test_lb = D_m0(1,:) <= D_spc(1,:);
test_ub = D_m0(2,:) >= D_spc(2,:);
n = 1;
if test_lb(1) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and M0 x-component: Dx(spc) = %1.2e m is less than Dx(M0) = %1.2e m\n',D_spc(1,1),D_m0(1,1));n = n+1;
end
if test_lb(2) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and M0 y-component: Dy(spc) = %1.2e m is less than Dy(M0) = %1.2e m\n',D_spc(1,2),D_m0(1,2));n = n+1;
end
if test_lb(3) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and M0 z-component: Dz(spc) = %1.2e m is less than Dz(M0) = %1.2e m\n',D_spc(1,3),D_m0(1,3));n = n+1;
end
if test_ub(1) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and M0 x-component: Dx(spc) = %1.2e m is more than Dx(M0) = %1.2e m\n',D_spc(2,1),D_m0(2,1));n = n+1;
end
if test_ub(2) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and M0 y-component: Dy(spc) = %1.2e m is more than Dy(M0) = %1.2e m\n',D_spc(2,2),D_m0(2,2));n = n+1;
end
if test_ub(3) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and M0 z-component: Dz(spc) = %1.2e m is more than Dz(M0) = %1.2e m\n',D_spc(2,3),D_m0(2,3));n = n+1;
end
% format long
% test_lb = (D_md(1,:) - D_spc(1,:))
% test_ub = (D_md(2,:) - D_spc(2,:))
test_lb = D_md(1,:) <= D_spc(1,:);
test_ub = D_md(2,:) >= D_spc(2,:);
if test_lb(1) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and Md x-component: Dx(spc) = %1.2e m is less than Dx(Md) = %1.2e m\n',D_spc(1,1),D_md(1,1));n = n+1;
end
if test_lb(2) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and Md y-component: Dy(spc) = %1.2e m is less than Dy(Md) = %1.2e m\n',D_spc(1,2),D_md(1,2));n = n+1;
end
if test_lb(3) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and Md z-component: Dz(spc) = %1.2e m is less than Dz(Md) = %1.2e m\n',D_spc(1,3),D_md(1,3));n = n+1;
end
if test_ub(1) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and Md x-component: Dx(spc) = %1.2e m is more than Dx(Md) = %1.2e m\n',D_spc(2,1),D_md(2,1));n = n+1;
end
if test_ub(2) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and Md y-component: Dy(spc) = %1.2e m is more than Dy(Md) = %1.2e m\n',D_spc(2,2),D_md(2,2));n = n+1;
end
if test_ub(3) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and Md z-component: Dz(spc) = %1.2e m is more than Dz(Md) = %1.2e m\n',D_spc(2,3),D_md(2,3));n = n+1;
end
test_lb = D_mps(1,:) <= D_spc(1,:);
test_ub = D_mps(2,:) >= D_spc(2,:);

if test_lb(1) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and MPS x-component: Dx(spc) = %1.2e m is less than Dx(MPS) = %1.2e m\n',D_spc(1,1),D_mps(1,1));n = n+1;
end
if test_lb(2) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and MPS y-component: Dy(spc) = %1.2e m is less than Dy(MPS) = %1.2e m\n',D_spc(1,2),D_mps(1,2));n = n+1;
end
if test_lb(3) == 0
    Msg{n} = sprintf('blOCh__spc: Lower bound dimension violation in spc and MPS z-component: Dz(spc) = %1.2e m is less than Dz(MPS) = %1.2e m\n',D_spc(1,3),D_mps(1,3));n = n+1;
end
if test_ub(1) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and MPS x-component: Dx(spc) = %1.2e m is more than Dx(MPS) = %1.2e m\n',D_spc(2,1),D_mps(2,1));n = n+1;
end
if test_ub(2) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and MPS y-component: Dy(spc) = %1.2e m is more than Dy(MPS) = %1.2e m\n',D_spc(2,2),D_mps(2,2));n = n+1;
end
if test_ub(3) == 0
    Msg{n} = sprintf('blOCh__spc: Upper bound dimension violation in spc and MPS z-component: Dz(spc) = %1.2e m is more than Dz(MPS) = %1.2e m\n',D_spc(2,3),D_mps(2,3));n = n+1;
end








end

function Msg = check_Orientations_Intel(spc,m0,md,mps)
% function Msg = check_Orientations_Intel(spc,m0,md,mps)
%
%   This script checks the orientations

R_spc = spc.R;
D_spc = spc.D;
Rv_spc = spc.Rv;

R_m0 = m0.R;
D_m0 = m0.D;
Rv_m0 = m0.Rv;

R_md = md.R;
D_md = md.D;
Rv_md = md.Rv;

R_mps = mps.R;
D_mps = mps.D;
Rv_mps = mps.Rv;

Intel = [...
    1	0	0	1	0	0	0	1	1	0	1	1	1	1
    0	1	0	0	1	0	1	1	0	1	1	0	1	1
    0	0	1	0	0	1	1	0	1	1	0	1	1	1
    1	0	0	1	0	0	0	1	1	0	1	1	1	1
    0	1	0	0	1	0	1	1	0	1	1	0	1	1
    0	0	1	0	0	1	1	0	1	1	0	1	1	1
    0	0	0	0	0	0	1	0	0	1	0	0	1	1
    0	0	0	0	0	0	0	1	0	0	1	0	1	1
    0	0	0	0	0	0	0	0	1	0	0	1	1	1
    0	0	0	0	0	0	1	0	0	1	0	0	1	1
    0	0	0	0	0	0	0	1	0	0	1	0	1	1
    0	0	0	0	0	0	0	0	1	0	0	1	1	1
    0	0	0	0	0	0	0	0	0	0	0	0	1	1
    0	0	0	0	0	0	0	0	0	0	0	0	1	1];


[O_spc_R,N_spc] = getOrR(R_spc,Rv_spc);
[O_md_R,N_md] = getOrR(R_md,Rv_md);
[O_m0_R,N_m0] = getOrR(R_m0,Rv_m0);
[O_mps_R,N_mps] = getOrR(R_mps,Rv_mps);

Msg = {};
n= 1;



if ~strcmpi(O_spc_R,O_mps_R)
    
    if Intel(N_spc,N_mps) ~= 1
        Msg{n} = sprintf('blOCh__spc: Orientation mismatch between spc (%s) and mps (%s)\n',O_spc_R,O_mps_R);n = n+1;
    end
end
if ~strcmpi(O_spc_R,O_md_R)
    if Intel(N_spc,N_md) ~= 1
        Msg{n} = sprintf('blOCh__spc: Orientation mismatch between spc (%s) and md (%s)\n',O_spc_R,O_md_R);n = n+1;
    end
end
if ~strcmpi(O_spc_R,O_m0_R)
    if Intel(N_spc,N_m0) ~= 1
        Msg{n} = sprintf('blOCh__spc: Orientation mismatch between spc (%s) and m0 (%s)\n',O_spc_R,O_m0_R);n = n+1;
    end
end

end

function O = Collapse_Dimension(O,Collapse,a_x,a_y,a_z,a_v,b_x,b_y,b_z,b_v)
% function O = Collapse_Dimension(O,Collapse,a_x,a_y,a_z,a_v,b_x,b_y,b_z,b_v)
%
%   This script collapses unnecesary spatial dimensions
%

if Collapse(1)
    O = mean(O(:,a_x:b_x,:,:,:),2);
end
if Collapse(2)
    O = mean(O(a_y:b_y,:,:,:,:),1);
end
if Collapse(3)
    O = mean(O(:,:,a_z:b_z,:,:),3);
end
if Collapse(4)
    O = mean(O(:,:,:,a_v:b_v,:),4);
end
end



function O = Expand_Spectral(I,Expand,R)
% function O = Expand_Spectral(I,Expand,R)
%
%   This script expands I in the spectral dimension
%

if Expand
    
    O = repmat(I,[1,1,1,R,1]);
else
    O = I;
end
end

function [R,idx1,idx2] = find_R_and_idxs(D1,D2)
% function [R,idx1,idx2] = find_R_and_idxs(D1,D2)
%
% Given a dimension span D1, and a another dimension span D2, within D1
% this function calculates a resolution R of a fictitious D2' span which
% holds both D2 and D1 numbers.
% In that new span D2' the D1 numbers will have indices idx1 and idx2.
%
%

temp = [D1(1),D2(1),D2(2),D1(2)];
tempmax = max(abs(temp));

if tempmax > 100
    D = round([D1(1),D2(1),D2(2),D1(2)].*100);
else
    D = round([D1(1),D2(1),D2(2),D1(2)].*100000);
end
G1 = gcd(abs(D(1)),abs(D(2)));
G2 = gcd(abs(D(1)),abs(D(3)));
G3 = gcd(abs(D(1)),abs(D(4)));
G4 = gcd(abs(D(2)),abs(D(3)));
G5 = gcd(abs(D(2)),abs(D(4)));
G6 = gcd(abs(D(3)),abs(D(4)));

G = min([G1,G2,G3,G4,G5,G6]);
if G  == 0
    G = 1;
end
Di = [D(1):G:D(end)];
R = length(Di);

idx1 = find(D(2)==Di);
idx2 = find(D(3)==Di);

end


function [Or,Num] = getOrR(R,Rv)
% function [Or,Num] = getOrR(R,Rv)
%
%   This script gets the orientation from R and Rv
%

if nargin == 1
    Rv = 1;
end


Dim = zeros(1,3);
if R(1) > 1
    Dim(1) = 1;
end
if R(2) > 1
    Dim(2) = 1;
end
if R(3) > 1
    Dim(3) = 1;
end





if Dim(1) == 0 && Dim(2) == 0 && Dim(3) == 1 && Rv == 1
    Or = '1DSI';
    Num = 1;
elseif Dim(1) == 0 && Dim(2) == 1 && Dim(3) == 0 && Rv == 1
    Or = '1DAP';
    Num = 2;
elseif Dim(1) == 1 && Dim(2) == 0 && Dim(3) == 0 && Rv == 1
    Or = '1DRL';
    Num = 3;
elseif Dim(1) == 0 && Dim(2) == 0 && Dim(3) == 1 && Rv > 1
    Or = '1+1DSI';
    Num = 4;
elseif Dim(1) == 0 && Dim(2) == 1 && Dim(3) == 0 && Rv > 1
    Or = '1+1DAP';
    Num = 5;
elseif Dim(1) == 1 && Dim(2) == 0 && Dim(3) == 0 && Rv > 1
    Or = '1+1DRL';
    Num = 6;
elseif Dim(1) == 1 && Dim(2) == 1 && Dim(3) == 0 && Rv == 1
    Or = '2DAx';
    Num = 7;
elseif Dim(2) == 1 && Dim(3) == 1 && Dim(1) == 0 && Rv == 1
    Or = '2DSa';
    Num = 8;
elseif Dim(1) == 1 && Dim(3) == 1 && Dim(2) == 0 && Rv == 1
    Or = '2DCo';
    Num = 9;
elseif Dim(1) == 1 && Dim(2) == 1 && Dim(3) == 0 && Rv > 1
    Or = '2+1DAx';
    Num = 10;
elseif Dim(2) == 1 && Dim(3) == 1 && Dim(1) == 0 && Rv > 1
    Or = '2+1DSa';
    Num = 11;
elseif Dim(1) == 1 && Dim(3) == 1 && Dim(2) == 0 && Rv > 1
    Or = '2+1DCo';
    Num = 12;
elseif Dim(1) == 1 && Dim(2) == 1 && Dim(3) == 1 && Rv == 1
    Or = '3D';
    Num = 13;
elseif Dim(1) == 1 && Dim(2) == 1 && Dim(3) == 1 && Rv > 1
    Or = '3+1D';
    Num = 14;
elseif Dim(1) == 0 && Dim(2) == 0 && Dim(3) == 0 && Rv > 1
    Or = '0+1D';
    Num = 15;
end
end



function [O,Xi,Yi,Zi,Vi] = Interp(Dim,I,D,R,Dv,Rv,Di,Ri,Dvi,Rvi)
% function [O,Xi,Yi,Zi,Vi] = Interp(Dim,I,D,R,Dv,Rv,Di,Ri,Dvi,Rvi)
%
%   This script interpolates spatial arrays
%


[Y,X,Z,V]     = ndgrid(linspace(D(1,2),D(2,2),R(2)),linspace(D(1,1),D(2,1),R(1)),linspace(D(1,3),D(2,3),R(3)),linspace(Dv(1),Dv(2),Rv));
[Yi,Xi,Zi,Vi] = ndgrid(linspace(Di(1,2),Di(2,2),Ri(2)),linspace(Di(1,1),Di(2,1),Ri(1)),linspace(Di(1,3),Di(2,3),Ri(3)),linspace(Dvi(1),Dvi(2),Rvi));

O = zeros(Ri(2),Ri(1),Ri(3),Rvi,size(I,5));
switch Dim
    case '0+1D'
        for n = 1:size(I,5)
            temp = permute(interp1(squeeze(V),squeeze(I(:,:,:,:,n)),squeeze(Vi)),[1,2,3]);
            O(1,1,1,:,n) = temp;
            %          O(:,:,:,:,n) = permute(interp1(squeeze(Z),squeeze(I(:,:,:,:,n)),squeeze(Zi)),[2,3,1]);
        end
    case '1DSI'
        for n = 1:size(I,5)
            temp = permute(interp1(squeeze(Z),squeeze(I(:,:,:,:,n)),squeeze(Zi)),[2,3,1]);
            O(1,1,:,1,n) = temp;
            %          O(:,:,:,:,n) = permute(interp1(squeeze(Z),squeeze(I(:,:,:,:,n)),squeeze(Zi)),[2,3,1]);
        end
    case '1DAP'
        for n = 1:size(I,5)
            O(:,:,:,:,n) = interp1(Y,I(:,:,:,:,n),Yi);
        end
    case '1DRL'
        %       O = interp1(X,squeeze(I),Xi);
        for n = 1:size(I,5)
            O(:,:,:,:,n) = interp1(X,squeeze(I(:,:,:,:,n)),Xi);
        end
    case '1+1DSI'
        for n = 1:size(I,5)
            temp = interp2(permute(squeeze(Z),[2,1]),permute(squeeze(V),[2,1]),permute(squeeze(I(:,:,:,:,n)),[2,1]),permute(squeeze(Zi),[2,1]),permute(squeeze(Vi),[2,1]));
            temp2 = permute(permute(temp,[2,1]),[3,4,1,2]);
            O(1,1,:,:,n) = temp2;
        end
    case '1+1DAP'
        for n = 1:size(I,5)
            temp = interp2(permute(squeeze(Y),[2,1]),permute(squeeze(V),[2,1]),permute(squeeze(I(:,:,:,:,n)),[2,1]),permute(squeeze(Yi),[2,1]),permute(squeeze(Vi),[2,1]));
            temp2 = permute(permute(temp,[2,1]),[1,3,4,2]);
            O(:,1,1,:,n) = temp2;
        end
    case '1+1DRL'
        for n = 1:size(I,5)
            O(:,:,:,:,n) = permute(permute(interp2(permute(squeeze(X),[2,1]),permute(squeeze(V),[2,1]),permute(squeeze(I(:,:,:,:,n)),[2,1]),permute(squeeze(Xi),[2,1]),permute(squeeze(Vi),[2,1])),[2,1]),[3,1,4,2]);
        end
    case '2DAx'
        for n = 1:size(I,5)
            O(:,:,:,:,n) = interp2(X,Y,I(:,:,:,:,n),Xi,Yi);
        end
    case '2DCo'
        for n = 1:size(I,5)
            temp = interp2(permute(squeeze(X),[2,1]),permute(squeeze(Z),[2,1]),permute(squeeze(I(:,:,:,:,n)),[2,1]),permute(squeeze(Xi),[2,1]),permute(squeeze(Zi),[2,1]));
            temp2 = permute(temp,[3,2,1]);
            O(1,:,:,1,n) = temp2;
        end
    case '2DSa'
        for n = 1:size(I,5)
            temp = interp2(permute(squeeze(Y),[2,1]),permute(squeeze(Z),[2,1]),permute(squeeze(I(:,:,:,:,n)),[2,1]),permute(squeeze(Yi),[2,1]),permute(squeeze(Zi),[2,1]));
            temp2 = permute(permute(temp,[2,1]),[1,3,2]);
            O(:,1,:,1,n) = temp2;
        end
    case '3D'
        for n = 1:size(I,5)
            O(:,:,:,:,n) = interp3(X,Y,Z,I(:,:,:,:,n),Xi,Yi,Zi);
        end
    case '2+1DAx'
        for n = 1:size(I,5)
            O(:,:,:,:,n) = permute(interp3(squeeze(X),squeeze(Y),squeeze(V),squeeze(I(:,:,:,:,n)),squeeze(Xi),squeeze(Yi),squeeze(Vi)),[1,2,4,3]);
        end
    case '2+1DCo'
        for n = 1:size(I,5)
            O(:,:,:,:,n) = permute(interp3(permute(squeeze(X),[2,1,3]),permute(squeeze(Z),[2,1,3]),permute(squeeze(V),[2,1,3]),permute(squeeze(I(:,:,:,:,n)),[2,1,3]),permute(squeeze(Xi),[2,1,3]),permute(squeeze(Zi),[2,1,3]),permute(squeeze(Vi),[2,1,3])),[4,2,1,3]);
        end
    case '2+1DSa'
        for n = 1:size(I,5)
            O(:,:,:,:,n) = permute(interp3(permute(squeeze(Y),[2,1,3]),permute(squeeze(Z),[2,1,3]),permute(squeeze(V),[2,1,3]),permute(squeeze(I(:,:,:,:,n)),[2,1,3]),permute(squeeze(Yi),[2,1,3]),permute(squeeze(Zi),[2,1,3]),permute(squeeze(Vi),[2,1,3])),[2,4,1,3]);
        end
    case '3+1D'
        for n = 1:size(I,5)
            O(:,:,:,:,n) = interpn(Y,X,Z,V,I(:,:,:,:,n),Yi,Xi,Zi,Vi);
        end
end
end

function [m0,Msg] = Load_m0(spc)
%
%   This script load the m0 array
%

M0 = spc.M0;
M0flip = spc.M0flip;

% TODOs make MSL general



Status = 0;


if strcmp(M0,'') || isempty(M0)
    
    [M0,D,Dv,BW,L] = Get_default_M0(spc);
    
    
    
elseif strcmp(M0,'MNI')
    if isunix
        
        Sep = '/';
    else
        Sep = '\';
    end
    
    Path = strrep(mfilename('fullpath'),['scripts',Sep,mfilename],['inputdata',Sep,'M0',Sep]);
    
    switch spc.Dim
        case {'2DAx','2+1DAx','3D','3+1D','1DAP','1+1DAP','1DRL','1+1DRL','0+1D'}
            load([Path,'MNIpdAx.m0'],'-mat')
        case {'2DCo','2+1DCo','1DSI','1+1DSI'}
            load([Path,'MNIpdCo.m0'],'-mat')
        case {'2DSa','2+1DSa'}
            load([Path,'MNIpdSa.m0'],'-mat')
    end
    [A,B,C] = size(M0);
    Phant = repmat(M0,[1,1,1,spc.Rv]);
    
    M0 = ones(A,B,C,spc.Rv,3);
    M0(:,:,:,:,3) = Phant.*cosd(M0flip(1));
    M0(:,:,:,:,1) = Phant.*sind(M0flip(1)).*cosd(M0flip(2));
    M0(:,:,:,:,2) = Phant.*sind(M0flip(1)).*sind(M0flip(2));
    
    %    M0 = repmat(M0,[1,1,1,spc.Rv,1]);
    
    %    D = spc.D;
    Dv = spc.Dv;
    BW = spc.BW;
    L = (D(2,:)-D(1,:));
    %    L = spc.L;
    
    
elseif strcmp(M0,'MSL')
    Rmax = max(max(spc.R(2),spc.R(1)),spc.R(3));
    Phant = MSL_3D_phantom([Rmax,Rmax,Rmax]);
    
    M0 = ones(Rmax,Rmax,Rmax,spc.Rv,3);
    M0(:,:,:,:,3) = Phant.*cosd(M0flip(1));
    M0(:,:,:,:,1) = Phant.*sind(M0flip(1)).*cosd(M0flip(2));
    M0(:,:,:,:,2) = Phant.*sind(M0flip(1)).*sind(M0flip(2));
    
    M0 = repmat(M0,[1,1,1,spc.Rv,1]);
    
    D = spc.D;
    Dv = spc.Dv;
    BW = spc.BW;
    L = spc.L;
    
else
    
    try validateattributes(M0,{'char'},{'2d'})
        
        [pathstr, name, extension] = fileparts(M0);
        
        switch extension
            case {'.jpg','.png','.tif','.jpeg','.tiff'}
                
                
                try
                    temp = imread(spc.M0);
                    temp = double(temp(:,:,1));
                    temp = temp./max(temp(:));
                    temp = temp(1:end-10,:);
                    switch spc.Dim(end-1:end)
                        case 'Ax'
                            temp = flipud(temp);
                        case 'Sa'
                            temp = rot90(temp,3);
                        case 'Co'
                            temp = flipud(temp);
                    end
                    
                    switch spc.Dim
                        case '3D'
                            temp = flipud(temp);
                        case '3+1D'
                            temp = flipud(temp);
                    end
                    spc.M0typ  = 1;
                    if spc.M0typ == 0
                        M0z = temp;
                        idxtemp = find(temp);
                        temp2 = ones(size(temp))-temp;
                        M0x = zeros(size(temp));
                        
                        M0y = M0x;
                        M0x(idxtemp) = temp2(idxtemp).*sind(spc.M0flip(1)).*cosd(spc.M0flip(2));
                        M0y(idxtemp) = temp2(idxtemp).*sind(spc.M0flip(1)).*sind(spc.M0flip(2));
                        
                        
                        
                    elseif spc.M0typ == 1
                        M0z = ones(size(temp))-temp;
                        M0x = zeros(size(temp));
                        idxtemp = find(temp);
                        M0y = M0x;
                        M0x(idxtemp) = temp(idxtemp).*sind(spc.M0flip(1)).*cosd(spc.M0flip(2));
                        M0y(idxtemp) = temp(idxtemp).*sind(spc.M0flip(1)).*sind(spc.M0flip(2));
                        
                        M0 = cat(5,M0x,M0y,M0z);
                        
                        
                    end
                    
                    
                    
                    
                    
                    
                    switch spc.Dim
                        case {'1DSI','1DAP','1DRL'}
                            
                            Nfo = ['blOCh__spc: image file for given spc.Dim is not supported. Assigning default M0.'];
                            Display_Message(Nfo,spc.Print);
                            
                            [M0,D,Dv,BW,L] = Get_default_M0(spc);
                            
                        case {'1+1DSI','1+1DAP','1+1DRL','0+1D'}
                            
                            Nfo = ['blOCh__spc: image file for given spc.Dim is not supported. Assigning default M0.'];
                            Display_Message(Nfo,spc.Print);
                            
                            [M0,D,Dv,BW,L] = Get_default_M0(spc);
                            
                            
                        case '2DAx'
                            
                        case '2DCo'
                            
                            M0 = permute(M0,[3,1,2,4,5]);
                        case '2DSa'
                            
                            
                            M0 = permute(M0,[2,3,1,4,5]);
                        case '2+1DAx'
                            M0 = Spect_Expand_M0(M0,spc);
                            
                        case '2+1DCo'
                            
                            M0 = permute(M0,[3,1,2,4,5]);
                            M0 = Spect_Expand_M0(M0,spc);
                            
                        case '2+1DSa'
                            M0 = permute(M0,[2,3,1,4,5]);
                            M0 = Spect_Expand_M0(M0,spc);
                            
                        case {'3D'}
                            M0 = Spat_Expand_M0(M0,spc);
                            
                        case {'3+1D'}
                            M0 = Spat_Expand_M0(M0,spc);
                            
                            M0 = Spect_Expand_M0(M0,spc);
                            
                            
                    end
                    
                    
                    
                    R = [size(M0,2),size(M0,1),size(M0,3)];
                    
                    D = spc.D;
                    Dv = spc.Dv;
                    BW = spc.BW;
                    L = spc.L;
                    
                    
                    
                catch me;Msg = ['blOCh__spc: ',me.message,' Someting is wrong'];end
                
            case '.m0'
                
                try m0 = load(M0,'-mat');
                    try validateattributes(m0.M0,{'double'},{'real','finite'});
                        R = [size(m0.M0,2),size(m0.M0,1),size(m0.M0,3)];
                        try validateattributes(m0.Dv,{'double'},{'size',[2,1],'real','finite'});
                            Dv = m0.Dv;
                            BW = (Dv(2,:)-Dv(1,:));
                            try validateattributes(m0.D,{'double'},{'size',[2,3],'real','finite'});
                                D = m0.D;
                                L = (D(2,:)-D(1,:));
                            catch me;Msg = ['blOCh__spc: ',me.message,' The M0 did not include a D variable or a proper D variable'];end
                        catch me;Msg = ['blOCh__spc: ',me.message,' The M0 did not include a Dv variable or a proper Dv variable'];end
                        
                        if size(m0.M0,5) == 1
                            temp = m0.M0;
                            M0z = temp.*cosd(M0flip(1));
                            M0x = temp.*sind(M0flip(1)).*cosd(M0flip(2));
                            M0y = temp.*sind(M0flip(1)).*sind(M0flip(2));
                            
                            M0 = cat(5,M0x,M0y,M0z);
                            
                        elseif size(m0.M0,5) == 3
                            % maybe apply rotation here
                            M0 = m0.M0;
                        else
                            Msg = ['blOCh__spc: size(M0,5) must be 1 (for M0z) or 3 (for M0x, M0y, M0z)'];
                        end
                        
                    catch me;Msg = ['blOCh__spc: ',me.message,' The M0 must include a M0 variable.'];end
                catch me;Msg = ['blOCh__spc: ',me.message,' Something is wrong.'];end
                
                
                
                
                
            otherwise
                Msg = ['blOCh__spc: ','The M0 file is supposed to be an image file (*.jpg, *.jpeg, *.tif, *.tiff, *.png), a mat-file (*.m0), or be [] or ''''.'];
        end
        
    catch me;Msg = ['blOCh__spc: ',me.message,' This is supposed to be the filename of the 2D M0'];end
    
end
try
    m0.M0 = M0;
    m0.R = [size(M0,2),size(M0,1),size(M0,3)];
    m0.Rv = size(M0,4);
    m0.Dim = getOrR(m0.R,m0.Rv);
    m0.D = D;
    m0.Dv = Dv;
    m0.L = L;
    m0.BW = BW;
    m0.Status = 1;
    Msg = '';
catch me;Msg = ['blOCh__spc: ',me.message,' Compiling M0 went wrong'];end


end

function [M0,D,Dv,BW,L] = Get_default_M0(spc)

M0 = spc.M0;
M0flip = spc.M0flip;

M0 = ones([[spc.R(2),spc.R(1),spc.R(3)],spc.Rv,3]);
M0(:,:,:,:,3) = cosd(M0flip(1));
M0(:,:,:,:,1) = sind(M0flip(1)).*cosd(M0flip(2));
M0(:,:,:,:,2) = sind(M0flip(1)).*sind(M0flip(2));

D = spc.D;
Dv = spc.Dv;
BW = spc.BW;
L = spc.L;
end

function M0 = Spect_Expand_M0(M0,spc)

M0 = repmat(M0,[1,1,1,spc.Rv,1]);
end


function M0 = Spat_Expand_M0(M0,spc)

M0 = repmat(M0,[1,1,spc.R(3),1,1]);
end


function [md,Msg] = Load_md(spc)
% function [md,Msg] = Load_md(spc)
%
%   This script loads the md array
%

Mdflip = spc.Mdflip;


Status = 0;


if strcmp(spc.Md,'') || isempty(spc.Md)
    
    [MD,D,Dv,BW,L] = Get_default_Md(spc);
    
    
    
    
    
else
    
    try validateattributes(spc.Md,{'char'},{'2d'})
        
        [pathstr, name, extension] = fileparts(spc.Md);
        
        switch extension
            
            case {'.jpg','.png','.tif','.jpeg','.tiff'}
                
                try
                    
                    temp = imread(spc.Md);
                    temp = double(temp(:,:,1));
                    temp = temp./max(temp(:));
                    %                temp = temp(1:end-10,:);
                    switch spc.Dim(end-1:end)
                        case 'Ax'
                            temp = flipud(temp);
                        case 'Sa'
                            temp = rot90(temp,3);
                        case 'Co'
                            temp = flipud(temp);
                    end
                    
                    switch spc.Dim
                        case '3D'
                            temp = flipud(temp);
                        case '3+1D'
                            temp = flipud(temp);
                    end
                   
                    Mdz = ones(size(temp));
                    Mdx = zeros(size(temp));
                    idxtemp = find(temp);
                    Mdy = Mdx;
                    Mdx(idxtemp) = temp(idxtemp).*sind(spc.Mdflip(1)).*cosd(spc.Mdflip(2));
                    Mdy(idxtemp) = temp(idxtemp).*sind(spc.Mdflip(1)).*sind(spc.Mdflip(2));
                    Mdz(idxtemp) = cosd(spc.Mdflip(1));
                    MD = cat(5,Mdx,Mdy,Mdz);
                   
                    switch spc.Dim
                        case {'1DSI','1DAP','1DRL'}
                            
                            Nfo = ['blOCh__spc: image file for given spc.Dim is not supported. Assigning default Md.'];
                            Display_Message(Nfo,spc.Print);
                            
                            [MD,D,Dv,BW,L] = Get_default_Md(spc);
                            
                        case {'1+1DSI','1+1DAP','1+1DRL','0+1D'}
                            
                            Nfo = ['blOCh__spc: image file for given spc.Dim is not supported. Assigning default Md.'];
                            Display_Message(Nfo,spc.Print);
                            
                            [MD,D,Dv,BW,L] = Get_default_Md(spc);
                            
                            
                        case '2DAx'
                            
                            
                        case '2DCo'
                            
                            MD = permute(MD,[3,1,2,4,5]);
                            
                        case '2DSa'
                            
                            
                            MD = permute(MD,[2,3,1,4,5]);
                            
                        case '2+1DAx'
                            MDnoroi = zeros(size(MD));
                            MDnoroi(:,:,:,:,3) = 1;
                            MD = Spect_Expand_MD(MD,MDnoroi,spc.Rv,spc.pBW,spc.c_pBW,spc.BW,spc.c_BW);
                            
                        case '2+1DCo'
                            
                            MD = permute(MD,[3,1,2,4,5]);
                            MDnoroi = zeros(size(MD));
                            MDnoroi(:,:,:,:,3) = 1;
                            MD = Spect_Expand_MD(MD,MDnoroi,spc.Rv,spc.pBW,spc.c_pBW,spc.BW,spc.c_BW);
                            
                        case '2+1DSa'
                            MD = permute(MD,[2,3,1,4,5]);
                            MDnoroi = zeros(size(MD));
                            MDnoroi(:,:,:,:,3) = 1;
                            MD = Spect_Expand_MD(MD,MDnoroi,spc.Rv,spc.pBW,spc.c_pBW,spc.BW,spc.c_BW);
                            
                        case {'3D'}
                            MDnoroi = zeros(size(MD));
                            MDnoroi(:,:,:,:,3) = 1;
                            MD = Spat_Expand_MD(MD,MDnoroi,spc.R(3),spc.TH(3),spc.c_TH(3),spc.L(3));
                            
                        case {'3+1D'}
                            MDnoroi = zeros(size(MD));
                            MDnoroi(:,:,:,:,3) = 1;
                            MD = Spat_Expand_MD(MD,MDnoroi,spc.R(3),spc.TH(3),spc.c_TH(3),spc.L(3));
                            MDnoroi = zeros(size(MD));
                            MDnoroi(:,:,:,:,3) = 1;
                            MD = Spect_Expand_MD(MD,MDnoroi,spc.Rv,spc.pBW,spc.c_pBW,spc.BW,spc.c_BW);
                            
                    end
                    
                    
                    
                    R = [size(MD,2),size(MD,1),size(MD,3)];
                    
                    D = spc.D;
                    Dv = spc.Dv;
                    BW = spc.BW;
                    L = spc.L;
                    
                    
                    
                catch me;Msg = ['blOCh__spc: ',me.message,' Someting is wrong'];end
                
            case '.md'
                
                try md = load(spc.Md,'-mat');
                catch me;Msg = ['blOCh__spc: ',me.message,' Something is wrong.'];end
                
                
                try validateattributes(md.Md,{'double'},{'real','finite'});
                    R = [size(md.Md,2),size(md.Md,1),size(md.Md,3)];
                    
                    
                    try validateattributes(md.Dv,{'double'},{'size',[2,1],'real','finite'});
                        Dv = md.Dv;
                        BW = (Dv(2,:)-Dv(1,:));
                        try validateattributes(md.D,{'double'},{'size',[2,3],'real','finite'});
                            
                            D = md.D;
                            L = (D(2,:)-D(1,:));
                        catch me;Msg = ['blOCh__spc: ',me.message,' The Md did not include a D variable or a proper D variable'];end
                        
                    catch me;Msg = ['blOCh__spc: ',me.message,' The Md did not include a Dv variable or a proper Dv variable'];end
                    
                    
                    if size(md.Md,5) == 1
                        temp = md.Md;
                        Mdz = temp.*cosd(Mdflip(1));
                        Mdx = temp.*sind(Mdflip(1)).*cosd(Mdflip(2));
                        Mdy = temp.*sind(Mdflip(1)).*sind(Mdflip(2));
                        
                        MD = cat(5,Mdx,Mdy,Mdz);
                        
                    elseif size(md.Md,5) == 3
                        % maybe apply rotation here
                        MD = md.Md;
                    else
                        Msg = ['blOCh__spc: size(Md,5) must be 1 (for Mdz) or 3 (for Mdx, Mdy, Mdz)'];
                    end
                    
                catch me;Msg = ['blOCh__spc: ',me.message,' The Md must include a Md variable.'];end
                
                
            otherwise
                Msg = ['blOCh__spc: ','The Md file is supposed to be an image file (*.jpg, *.jpeg, *.tif, *.tiff, *.png), a mat-file (*.md), or be [] or ''''.'];
        end
        
    catch me;Msg = ['blOCh__spc: ',me.message,' This is supposed to be the filename of the 2D Md'];end
    
end

if strcmp(spc.Suppression,'on')
    MD = Convert_MD_to_Suppression(MD);
end

if strcmp(spc.Mpen,'on')
    MP1 = Get_MP(MD,spc,90);
    MP2 = Get_MP(MD,spc,-90);
end



try
    md.Md = MD;
    
    if strcmp(spc.Mpen,'on')
        md.MP1 = MP1;
        md.MP2 = MP2;
    end
    md.R = [size(MD,2),size(MD,1),size(MD,3)];
    md.Rv = size(MD,4);
    md.Dim = getOrR(md.R,md.Rv);
    
    md.D = D;
    md.Dv = Dv;
    md.L = L;
    md.BW = BW;
    md.Status = 1;
    Msg = '';
catch me;Msg = ['blOCh__spc: ',me.message,' Compiling Md went wrong'];end


end



function MD = Convert_MD_to_Suppression(MD)


MDx = MD(:,:,:,:,1);
MDy = MD(:,:,:,:,2);
MDz = MD(:,:,:,:,3);

theta  = acosd(MDz);
phi = atan2(MDy,MDx);

theta_max = max(theta(:));
theta_min = min(theta(:));

theta2 = theta_min+theta_max-theta;

MDxn = sind(theta2).*cosd(phi);
MDyn = sind(theta2).*sind(phi);
MDzn = cosd(theta2);


MD = cat(5,MDxn,MDyn,MDzn);


end

function [MP] = Get_MP(MD,spc,theta)

MP = zeros(size(MD));
MP(:,:,:,:,3) = zeros(size(MP(:,:,:,:,3)));
MP(:,:,:,:,1) = ones(size(MP(:,:,:,:,3))).*cosd(spc.Mdflip(2)+theta);
MP(:,:,:,:,2) = ones(size(MP(:,:,:,:,3))).*sind(spc.Mdflip(2)+theta);

end

function [MD,D,Dv,BW,L] = Get_default_Md(spc)

switch spc.Dim
    case {'1DSI','1DAP','1DRL'}
        MD = Make_1D_Md(spc);
        
    case {'1+1DSI','1+1DAP','1+1DRL'}
        
        [MD,MDnoroi] = Make_1D_Md(spc);
        
        MD = Spect_Expand_MD(MD,MDnoroi,spc.Rv,spc.pBW,spc.c_pBW,spc.BW,spc.c_BW);
        
    case {'0+1D'}
        
        MD = zeros([spc.R(2),spc.R(1),spc.R(3),spc.Rv,3]);
        MD(:,:,:,:,3) = 1;
        fl = (spc.c_pBW-spc.pBW/2);
        fr = (spc.c_pBW+spc.pBW/2);
        f = linspace(spc.Dv(1),spc.Dv(2),spc.Rv);
        
        [Mnl,idx1] = min(abs(f-fl));
        [Mnr,idx2] = min(abs(f-fr));
        
        
        if idx1 < 2
            idx1 = 2;
        end
        if idx2 > spc.Rv
            idx2 = spc.Rv;
        end
        if idx2 < 1
            idx2 = 1;
        end
        if idx1 > spc.Rv
            idx1 = spc.Rv;
        end
            MD(1,1,1,idx1:idx2,3) = cosd(spc.Mdflip(1));
            MD(1,1,1,idx1:idx2,1) = sind(spc.Mdflip(1)).*cosd(spc.Mdflip(2));
            MD(1,1,1,idx1:idx2,2) = sind(spc.Mdflip(1)).*sind(spc.Mdflip(2));

        
    case {'2DAx','2DCo','2DSa'}
        
        [MD,MDnoroi] = Make_2D_Md(spc);
        
        
    case {'2+1DAx','2+1DCo','2+1DSa'}
        
        [MD,MDnoroi] = Make_2D_Md(spc);
        
        MD = Spect_Expand_MD(MD,MDnoroi,spc.Rv,spc.pBW,spc.c_pBW,spc.BW,spc.c_BW);
        
    case {'3D'}
        
        [MD,MDnoroi] = Make_3D_Md(spc);
        
    case {'3+1D'}
        
        [MD,MDnoroi] = Make_3D_Md(spc);
        
        MD = Spect_Expand_MD(MD,MDnoroi,spc.Rv,spc.pBW,spc.c_pBW,spc.BW,spc.c_BW);
        
        
        
end

D = spc.D;
Dv = spc.Dv;
BW = spc.BW;
L = spc.L;
end



function [mps,Msg] = Load_mps(spc)
% function [mps,Msg] = Load_mps(spc)
%
%   This script loads the mps data.
%
Msg = [];
mps = [];
MPS = spc.MPS;
% useMPSdim = spc.useMPSdim;

if isunix
    
    Sep = '/';
else
    Sep = '\';
end

if ~strcmpi(MPS,'')
    Msg = [];
    try validateattributes(MPS,{'char'},{'2d'})
        [pathstr, name, extension] = fileparts(MPS);
        
        
        if strcmp(pathstr,'')
            Path = strrep(mfilename('fullpath'),['scripts',Sep,mfilename],['inputdata',Sep,'MPS',Sep]);
        else
            Path = [pathstr,Sep];
        end
        File = [Path, name, extension];
        
        
        
        
        switch extension
            case '.mps'
                
                try mps_temp = load(File,'-mat');
                    
                    try validateattributes(mps_temp.R,{'double'},{'size',[1,3],'integer','finite'});
                        R = mps_temp.R;
                        try validateattributes(mps_temp.Rv,{'double'},{'size',[1,1],'integer','finite'});
                            Rv = mps_temp.Rv;
                            
                            try validateattributes(mps_temp.D,{'double'},{'size',[2,3],'real','finite'});
                                D = mps_temp.D;
                                L = (D(2,:)-D(1,:));
                            catch me;Msg = ['blOCh__spc: ',me.message,' The MPS did not include a D variable or a proper D variable'];
                                
                            end
                            try validateattributes(mps_temp.Dv,{'double'},{'size',[2,1],'real','finite'});
                                Dv = mps_temp.Dv;
                                BW = (Dv(2)-Dv(1));
                            catch me;Msg = ['blOCh__spc: ',me.message,' The MPS did not include a Dv variable or a proper Dv variable'];
                                
                            end
                            try validateattributes(mps_temp.v0map,{'double'},{'size',[R(2),R(1),R(3),Rv],'real','finite'});
                                existv0 = 1;
                            catch me;Nfo = ['blOCh__spc: ',me.message,' The MPS did not include a v0map or a proper v0map. Making default.'];
                                Display_Message(Nfo,spc.Print);existv0 = 0;
                            end
                            
                            try validateattributes(mps_temp.B1map,{'double'},{'size',[R(2),R(1),R(3),Rv,spc.pTx],'finite'});
                                existB1 = 1;
                                pTx = size(mps_temp.B1map,5);
                                Display_Message(sprintf('Assuming the B1 map is normalized to nominal value. The reference voltage specifies what 1 in the B1map corresponds to.'),spc.Print);
                                
                                try validateattributes(mps_temp.B1_type,{'char'},{'nonempty'});
                                    
                                    val_type = {'absolute_valued','relative_valued'};
                                    if ~isempty(validatestring(mps_temp.B1_type,val_type))
                                        B1_type = mps_temp.B1_type;
                                    else
                                        Display_Message(' The MPS did not include a B1_type or a proper B1_type. Assigning ''relative_valued''',spc.Print);
                                        B1_type = 'relative_valued';
                                    end
                                catch me;Nfo = ['blOCh__spc: ',me.message,' The MPS did not include a B1_type or a proper B1_type. Assigning ''absolute_valued'''];
                                    Display_Message(Nfo,spc.Print);
                                    B1_type = 'relative_valued';
                                end
                                
                                
                                B1map = mps_temp.B1map;
                                
                                try validateattributes(mps_temp.B1_nom_amp_val,{'numeric'},{'size',[1,1],'real','positiv','finite'});
                                    
                                    B1_nom_amp_val = mps_temp.B1_nom_amp_val;
                                    
                                catch me;Nfo = ['blOCh__spc: ',me.message,' The MPS did not include a B1_nom_amp_val or a proper B1_nom_amp_val. Assigning unity'];
                                    Display_Message(Nfo,spc.Print);
                                    B1_nom_amp_val = 1;
                                end
                                try validateattributes(mps_temp.B1_nom_amp_unit,{'char'},{'nonempty'});
                                    
                                    val_units = {'T/V','V','relative'};
                                    if ~isempty(validatestring(mps_temp.B1_nom_amp_unit,val_units))
                                        B1_nom_amp_unit = mps_temp.B1_nom_amp_unit;
                                    else
                                        Display_Message(' The MPS did not include a B1_nom_amp_unit or a proper B1_nom_amp_unit. Assigning ''relative''',spc.Print);
                                        B1_nom_amp_unit = 'T/V';
                                    end
                                catch me;Nfo = ['blOCh__spc: ',me.message,' The MPS did not include a B1_nom_amp_unit or a proper B1_nom_amp_unit. Assigning ''relative'''];
                                    Display_Message(Nfo,spc.Print);
                                    B1_nom_amp_unit = 'relative';
                                end
                                
                                switch B1_type
                                    case 'relative_valued'
                                        try validateattributes(mps_temp.B1_ref_voltage,{'numeric'},{'size',[1,1],'real','positiv','finite'});
                                            
                                            B1_ref_voltage = mps_temp.B1_ref_voltage;
                                            
                                        catch me;Nfo = ['blOCh__spc: ',me.message,' The MPS did not include a B1_ref_voltage or a proper B1_ref_voltage. Assigning unity'];
                                            Display_Message(Nfo,spc.Print);
                                            B1_ref_voltage = 1;
                                        end
                                    case 'absolute_valued'
                                        
                                end
                                
                                
                                
                            catch me;Nfo = ['blOCh__spc: ',me.message,' The MPS did not include a B1map or a proper B1map. Making default.'];
                                Display_Message(Nfo,spc.Print);existB1 = 0;
                            end
                            if ~existB1
                                pTx = spc.pTx;
                                % make this better
                                if pTx > 1
                                    Dphasor = 2*pi/pTx;
                                    phasor = [0:Dphasor:2*pi-Dphasor];
                                    B1map = zeros(R(2),R(1),R(3),Rv,pTx);
                                    for p = 1:pTx
                                        B1map(:,:,:,:,p) = ones(R(2),R(1),R(3),Rv)./pTx.*exp(1i.*phasor(p).*ones(R(2),R(1),R(3),Rv));
                                    end
                                    %                            B1map = ones(R(2),R(1),R(3),Rv);
                                    %                            B1map = repmat(B1map,[1,1,1,1,pTx])./pTx;
                                    %                           B1map = repmat(B1map,[1,1,1,1,pTx]);
                                    
                                    %                            B1map = complex(B1map,zeros(size(B1map)));
                                else
                                    B1map = complex(ones(R(2),R(1),R(3),Rv),zeros(R(2),R(1),R(3),Rv));
                                    
                                end
                                B1_ref_voltage = 1;
                                B1_nom_amp_unit = 'relative';
                                B1_nom_amp_val = 1;
                                B1_type = 'relative_valued';
                            end
                            
                            
                            
                            try validateattributes(mps_temp.T1map,{'double'},{'size',[R(2),R(1),R(3),Rv],'real','finite','>=',0});
                                existT1 = 1;
                            catch me;Nfo = ['blOCh__spc: ',me.message,' The MPS did not include a T1map or a proper T1map. Making default.'];
                                Display_Message(Nfo,spc.Print);existT1 = 0;
                            end
                            
                            try validateattributes(mps_temp.T2map,{'double'},{'size',[R(2),R(1),R(3),Rv],'real','finite','>=',0});
                                existT2 = 1;
                            catch me;Nfo = ['blOCh__spc: ',me.message,' The MPS did not include a T2map or a proper T2map. Making default.'];
                                Display_Message(Nfo,spc.Print);existT2 = 0;
                            end
                            
                            if ~existv0
                                v0map = zeros(R(2),R(1),R(3),Rv);
                            else
                                v0map = mps_temp.v0map;
                            end
                            
                            
                            
                            if ~existT1
                                T1map = ones(R(2),R(1),R(3),Rv).*1e16;
                            else
                                T1map = mps_temp.T1map;
                            end
                            if ~existT2
                                T2map = ones(R(2),R(1),R(3),Rv).*1e16;
                            else
                                T2map = mps_temp.T2map;
                            end
                            
                            
                            
                        catch me;Msg = ['blOCh__spc: ',me.message,' The MPS must include Rv variable: [Rv] (#)'];end
                        
                        
                    catch me;Msg = ['blOCh__spc: ',me.message,' The MPS must include R variable: [Rx,Ry,Rz] (#)'];end
                catch me;Msg = ['blOCh__spc: ',me.message,' Something is wrong'];end
                
            otherwise
                Msg = ['blOCh__spc: ','The MPS file is supposed to be a "*.mps"-file (mat-file)'];
        end
    catch me; Msg = ['blOCh__spc: ',me.message,' This is supposed to be the filename of the 1D MPS, i.e., "*.m"-file (mat-file)'];end
    
    if ~isempty(Msg)
        Display_Message(Msg,spc.Print);
        Status = 0;
    else
        
        Status = 1;
        
    end
else
    R = spc.R;
    Rv = spc.Rv;
    D = spc.D;
    Dv = spc.Dv;
    
    L = spc.L;
    BW = spc.BW;
    v0map = zeros(R(2),R(1),R(3),Rv);
    
    pTx = spc.pTx;
    % make this better
    if pTx > 1
        Dphasor = 2*pi/pTx;
        phasor = [0:Dphasor:2*pi-Dphasor];
        B1map = zeros(R(2),R(1),R(3),Rv,pTx);
        for p = 1:pTx
            B1map(:,:,:,:,p) = ones(R(2),R(1),R(3),Rv)./pTx.*exp(1i.*phasor(p).*ones(R(2),R(1),R(3),Rv));
        end
        %                            B1map = ones(R(2),R(1),R(3),Rv);
        %                            B1map = repmat(B1map,[1,1,1,1,pTx])./pTx;
        %                           B1map = repmat(B1map,[1,1,1,1,pTx]);
        
        %                            B1map = complex(B1map,zeros(size(B1map)));
    else
        B1map = complex(ones(R(2),R(1),R(3),Rv),zeros(R(2),R(1),R(3),Rv));
        
    end
    B1_ref_voltage = 1;
    B1_nom_amp_unit = 'sensitivity';
    B1_nom_amp_val = 1;
    B1_type = 'relative_valued';
    T1map = ones(R(2),R(1),R(3),Rv).*1e16;
    
    T2map = ones(R(2),R(1),R(3),Rv).*1e16;
    Status = 1;
end

if Status
    try
        mps.v0map = v0map;
        mps.B1map = B1map;
        if exist('B1_ref_voltage')
            mps.B1_ref_voltage = B1_ref_voltage;
        end
        mps.B1_nom_amp_unit = B1_nom_amp_unit;
        mps.B1_nom_amp_val = B1_nom_amp_val;
        mps.B1_type = B1_type;
        mps.T1map = T1map;
        mps.T2map = T2map;
        mps.R = R;
        mps.Rv = Rv;
        mps.Dim = getOrR(R,Rv);
        mps.L = L;
        mps.BW = BW;
        mps.D = D;
        mps.Dv = Dv;
        mps.pTx = pTx;
        mps.Status = 1;
        Msg = '';
    catch me;Msg = ['blOCh__spc: ',me.message,' Compiling MPS went wrong'];end
end


end


function [MD,MDnoroi] = Make_1D_Md(spc)
% function [MD,MDnoroi] = Make_1D_Md(spc)
%
%   This script makes a 1D Md array
%

MD = zeros([spc.R(2),spc.R(1),spc.R(3),1,3]);
MD(:,:,:,:,3) = 1;
MDnoroi = MD;
switch spc.Dim
    case {'1DAP','1+1DAP'}
        %       idx1 = round((spc.c_TH(2)-spc.TH(2)/2)./spc.L(2).*spc.R(2)+spc.R(2)/2);
        %       idx2 = round((spc.c_TH(2)+spc.TH(2)/2)./spc.L(2).*spc.R(2)+spc.R(2)/2);
        
        fl = (spc.c_TH(2)-spc.TH(2)/2);
        fr = (spc.c_TH(2)+spc.TH(2)/2);
        f = linspace(spc.D(1,2),spc.D(2,2),spc.R(2));
        
        [Mnl,idx1] = min(abs(f-fl));
        [Mnr,idx2] = min(abs(f-fr));
        
        
        if idx1 < 1
            idx1 = 1;
        end
        if idx2 > spc.R(2)
            idx2 = spc.R(2);
        end
        if idx2 < 1
            idx2 = 1;
        end
        if idx1 > spc.R(2)
            idx1 = spc.R(2);
        end
        MD(idx1:idx2,1,1,1,3) = cosd(spc.Mdflip(1));
        MD(idx1:idx2,1,1,1,1) = sind(spc.Mdflip(1)).*cosd(spc.Mdflip(2));
        MD(idx1:idx2,1,1,1,2) = sind(spc.Mdflip(1)).*sind(spc.Mdflip(2));
    case {'1DRL','1+1DRL'}
        fl = (spc.c_TH(1)-spc.TH(1)/2);
        fr = (spc.c_TH(1)+spc.TH(1)/2);
        f = linspace(spc.D(1,1),spc.D(2,1),spc.R(1));
        
        [Mnl,idx1] = min(abs(f-fl));
        [Mnr,idx2] = min(abs(f-fr));
        if idx1 < 1
            idx1 = 1;
        end
        if idx2 > spc.R(1)
            idx2 = spc.R(1);
        end
        if idx2 < 1
            idx2 = 1;
        end
        if idx1 > spc.R(1)
            idx1 = spc.R(1);
        end
        MD(1,idx1:idx2,1,1,3) = cosd(spc.Mdflip(1));
        MD(1,idx1:idx2,1,1,1) = sind(spc.Mdflip(1)).*cosd(spc.Mdflip(2));
        MD(1,idx1:idx2,1,1,2) = sind(spc.Mdflip(1)).*sind(spc.Mdflip(2));
    case {'1DSI','1+1DSI'}
        
        fl = (spc.c_TH(3)-spc.TH(3)/2);
        fr = (spc.c_TH(3)+spc.TH(3)/2);
        f = linspace(spc.D(1,3),spc.D(2,3),spc.R(3));
        
        [Mnl,idx1] = min(abs(f-fl));
        [Mnr,idx2] = min(abs(f-fr));
        if idx1 < 1
            idx1 = 1;
        end
        if idx2 > spc.R(3)
            idx2 = spc.R(3);
        end
        if idx2 < 1
            idx2 = 1;
        end
        if idx1 > spc.R(3)
            idx1 = spc.R(3);
        end
        MD(1,1,idx1:idx2,1,3) = cosd(spc.Mdflip(1));
        MD(1,1,idx1:idx2,1,1) = sind(spc.Mdflip(1)).*cosd(spc.Mdflip(2));
        MD(1,1,idx1:idx2,1,2) = sind(spc.Mdflip(1)).*sind(spc.Mdflip(2));
        
end

end

function [Md,Mdnoroi] = Make_2D_Md(spc)
% function [Md,Mdnoroi] = Make_2D_Md(spc)
%
%   This script makes a 2D Md array
%


switch spc.Dim
    
    case {'2DAx','2+1DAx'}
        temp = spc;
        temp.Dim = '1DAP';
        temp.R = [1,spc.R(2),1];
        [MD1,MD1noroi] = Make_1D_Md(temp);
        temp = spc;
        temp.Dim = '1DRL';
        temp.R = [spc.R(1),1,1];
        [MD2,MD2noroi] = Make_1D_Md(temp);
        
        MD11 = repmat(MD1,[1,spc.R(1),1,1,1]);
        MD22 = repmat(MD2,[spc.R(2),1,1,1,1]);
        
        MD11noroi = repmat(MD1noroi,[1,spc.R(1),1,1,1]);
        MD22noroi = repmat(MD2noroi,[spc.R(2),1,1,1,1]);
        
        Md = Assign_FA(MD11,MD22,MD11noroi,MD22noroi,spc.Mdflip);
        Mdnoroi = MD11noroi;
        
        
        
        
        % d=1
    case {'2DCo','2+1DCo'}
        
        temp = spc;
        
        temp.Dim = '1DSI';
        temp.R = [1,1,spc.R(3)];
        [MD1,MD1noroi] = Make_1D_Md(temp);
        
        temp = spc;
        
        temp.Dim = '1DRL';
        temp.R = [spc.R(1),1,1];
        [MD2,MD2noroi] = Make_1D_Md(temp);
        
        MD11 = repmat(MD1,[1,spc.R(1),1,1,1]);
        MD22 = repmat(MD2,[1,1,spc.R(3),1,1]);
        MD11noroi = repmat(MD1noroi,[1,spc.R(1),1,1,1]);
        MD22noroi = repmat(MD2noroi,[1,1,spc.R(3),1,1]);
        
        Md = Assign_FA(MD11,MD22,MD11noroi,MD22noroi,spc.Mdflip);
        Mdnoroi = MD11noroi;
        
        
        
    case {'2DSa','2+1DSa'}
        
        temp = spc;
        temp.Dim = '1DSI';
        temp.R = [1,1,spc.R(3)];
        [MD1,MD1noroi] = Make_1D_Md(temp);
        temp = spc;
        temp.Dim = '1DAP';
        temp.R = [1,spc.R(2),1];
        [MD2,MD2noroi] = Make_1D_Md(temp);
        
        MD11 = repmat(MD1,[spc.R(2),1,1,1,1]);
        MD22 = repmat(MD2,[1,1,spc.R(3),1,1]);
        
        MD11noroi = repmat(MD1noroi,[spc.R(2),1,1,1,1]);
        MD22noroi = repmat(MD2noroi,[1,1,spc.R(3),1,1]);
        
        Md = Assign_FA(MD11,MD22,MD11noroi,MD22noroi,spc.Mdflip);
        Mdnoroi = MD11noroi;
        
end

end



function [Md,Mdnoroi] = Make_3D_Md(spc)
% function [Md,Mdnoroi] = Make_3D_Md(spc)
%
%   This script makes a 3D Md array
%

temp = spc;
temp.Dim = '2DAx';
temp.R = [temp.R(1),temp.R(2),1];
[MD1,MD1noroi] = Make_2D_Md(temp);

temp = spc;
temp.Dim = '2DSa';
temp.R = [1,temp.R(2),temp.R(3)];
[MD2,MD2noroi] = Make_2D_Md(temp);


MD11 = repmat(MD1,[1,1,spc.R(3),1,1]);
MD22 = repmat(MD2,[1,spc.R(1),1,1]);

MD11noroi = repmat(MD1noroi,[1,1,spc.R(3),1,1]);
MD22noroi = repmat(MD2noroi,[1,spc.R(1),1,1]);

Md = Assign_FA(MD11,MD22,MD11noroi,MD22noroi,spc.Mdflip);
Mdnoroi = MD11noroi;
end



function [O,Xi,Yi,Zi,Vi] = Morph(Dim_i,Dim_o,I,D_i,Dv_i,L,r0,L0,R_o,Rv_o,c_BW,BW)
% function [O,Xi,Yi,Zi,Vi] = Morph(Dim_i,Dim_o,I,D_i,Dv_i,L,r0,L0,R_o,Rv_o,c_BW,BW)
%
%   This script morphs the spatial data struct.
%

R_i = [size(I,2),size(I,1),size(I,3)];
Rv_i = size(I,4);
switch Dim_o
    case '1DSI'
        
        D_o = [...
            r0(1)-L0/2,r0(2)-L0/2,r0(3)-L(3)/2;...
            r0(1)+L0/2,r0(2)+L0/2,r0(3)+L(3)/2];
        
    case '1DAP'
        
        D_o = [...
            r0(1)-L0/2,r0(2)-L(2)/2,r0(3)-L0/2;...
            r0(1)+L0/2,r0(2)+L(2)/2,r0(3)+L0/2];
        
    case '1DRL'
        
        D_o = [...
            r0(1)-L(1)/2,r0(2)-L0/2,r0(3)-L0/2;...
            r0(1)+L(1)/2,r0(2)+L0/2,r0(3)+L0/2];
        
    case '1+1DSI'
        
        
        D_o = [...
            r0(1)-L0/2,r0(2)-L0/2,r0(3)-L(3)/2;...
            r0(1)+L0/2,r0(2)+L0/2,r0(3)+L(3)/2];
    case '1+1DAP'
        
        
        D_o = [...
            r0(1)-L0/2,r0(2)-L(2)/2,r0(3)-L0/2;...
            r0(1)+L0/2,r0(2)+L(2)/2,r0(3)+L0/2];
        
    case '1+1DRL'
        
        
        D_o = [...
            r0(1)-L(1)/2,r0(2)-L0/2,r0(3)-L0/2;...
            r0(1)+L(1)/2,r0(2)+L0/2,r0(3)+L0/2];
        
    case '2DAx'
        
        D_o = [...
            r0(1)-L(1)/2,r0(2)-L(2)/2,r0(3)-L0/2;...
            r0(1)+L(1)/2,r0(2)+L(2)/2,r0(3)+L0/2];
        
    case '2DSa'
        
        
        D_o = [...
            r0(1)-L0/2,r0(2)-L(2)/2,r0(3)-L(3)/2;...
            r0(1)+L0/2,r0(2)+L(2)/2,r0(3)+L(3)/2];
        
    case '2DCo'
        
        D_o = [...
            r0(1)-L(1)/2,r0(2)-L0/2,r0(3)-L(3)/2;...
            r0(1)+L(1)/2,r0(2)+L0/2,r0(3)+L(3)/2];
        
    case '2+1DAx'
        
        D_o = [...
            r0(1)-L(1)/2,r0(2)-L(2)/2,r0(3)-L0/2;...
            r0(1)+L(1)/2,r0(2)+L(2)/2,r0(3)+L0/2];
        
    case '2+1DSa'
        
        D_o = [...
            r0(1)-L0/2,r0(2)-L(2)/2,r0(3)-L(3)/2;...
            r0(1)+L0/2,r0(2)+L(2)/2,r0(3)+L(3)/2];
        
    case '2+1DCo'
        
        D_o = [...
            r0(1)-L(1)/2,r0(2)-L0/2,r0(3)-L(3)/2;...
            r0(1)+L(1)/2,r0(2)+L0/2,r0(3)+L(3)/2];
        
    case '3D'
        
        D_o = [...
            r0(1)-L(1)/2,r0(2)-L(2)/2,r0(3)-L(3)/2;...
            r0(1)+L(1)/2,r0(2)+L(2)/2,r0(3)+L(3)/2];
        
    case '3+1D'
        
        D_o = [...
            r0(1)-L(1)/2,r0(2)-L(2)/2,r0(3)-L(3)/2;...
            r0(1)+L(1)/2,r0(2)+L(2)/2,r0(3)+L(3)/2];
        
    case '0+1D'
        
        D_o = [...
            0,0,0;0,0,0];
        
        
    otherwise
        
end
Dv_o = [c_BW-BW/2;c_BW+BW/2];


%%

[R_x,a_x,b_x] = find_R_and_idxs(D_i(:,1),D_o(:,1));
[R_y,a_y,b_y] = find_R_and_idxs(D_i(:,2),D_o(:,2));
[R_z,a_z,b_z] = find_R_and_idxs(D_i(:,3),D_o(:,3));
[R_v,a_v,b_v] = find_R_and_idxs(Dv_i,Dv_o);

Morph = [Dim_i,'_2_',Dim_o];
Status = 1;
switch Morph
    
        
    case {'2+1DAx_2_2DAx','2+1DSa_2_2DSa','2+1DCo_2_2DCo'}
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,1];
        Expand = 0;
    case {'2DAx_2_2+1DAx','2DSa_2_2+1DSa','2DCo_2_2+1DCo'}
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i; % before R_v
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,0];
        Expand = 1;
        
    case {'3D_2_3D','2DAx_2_2DAx','2DSa_2_2DSa','2DCo_2_2DCo'}
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,0];
        Expand = 0;
        % 3+1D --> 3D
    case '3+1D_2_3D'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,1];
        Expand = 0;
        
        
        % 3D --> 3+1D
    case '3D_2_3+1D'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,0];
        Expand = 1;
        % 3+1D --> 3+1D
    case {'3+1D_2_3+1D','2+1DAx_2_2+1DAx','2+1DSa_2_2+1DSa','2+1DCo_2_2+1DCo','0+1D_2_0+1D'}
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,0];
        Expand = 0;
        
        
        % 3D --> 2D
    case '3D_2_2DAx'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1:2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1:2),1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1:2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 0;
        
    case '3D_2_2DSa'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2:3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2:3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2:3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 0;
        
        
    case '3D_2_2DCo'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,R_i(3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 0;
        
        
        % 3D --> 2+1D
        
    case '3D_2_2+1DAx'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1:2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1:2),1];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [D_o(:,1:2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 1;
        
    case '3D_2_2+1DCo'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,R_i(3)];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 1;
        
    case '3D_2_2+1DSa'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2:3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2:3)];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2:3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 1;
        
        
        % 3+1D --> 2+1D
        
    case '3+1D_2_2+1DAx'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1:2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1:2),1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1:2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 0;
        
        
    case '3+1D_2_2+1DSa'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2:3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2:3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2:3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 0;
        
        
    case '3+1D_2_2+1DCo'
        
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,R_i(3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 0;
        
        % 3+1D --> 2D
    case '3+1D_2_2DAx'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1:2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1:2),1];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1:2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,1];
        Expand = 0;
        
        
    case '3+1D_2_2DCo'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1:2),1];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,1];
        Expand = 0;
        
        
        
    case '3+1D_2_2DSa'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2:3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2:3)];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2:3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,1];
        Expand = 0;
        
        
        
        
        % 3D --> 1D
    case '3D_2_1DRL'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,1,0];
        Expand = 0;
        
    case '3D_2_1DAP'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,1,0];
        Expand = 0;
        
    case '3D_2_1DSI'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,1,0,0];
        Expand = 0;
        
        
        
        
        % 3+1D --> 1D
    case '3+1D_2_1DRL'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,1,1];
        Expand = 0;
        
    case '3+1D_2_1DAP'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,1,1];
        Expand = 0;
        
    case '3+1D_2_1DSI'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,1,0,1];
        Expand = 0;
        
        
        
        % 3D --> 1+1D
    case '3D_2_1+1DRL'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,1,0];
        Expand = 1;
        
    case '3D_2_1+1DAP'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,1,0];
        Expand = 1;
        
    case '3D_2_1+1DSI'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,1,0,0];
        Expand = 1;
        
        % 3+1D --> 1+1D
    case '3+1D_2_1+1DRL'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,1,0];
        Expand = 0;
        
    case '3+1D_2_1+1DAP'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,1,0];
        Expand = 0;
        
    case '3+1D_2_1+1DSI'
        
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,1,0];
        Expand = 0;
        
        % 1D --> 1D
        
    case {'1DAP_2_1DAP','1DSI_2_1DSI','1DRL_2_1DRL'}
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,0];
        Expand = 0;
        
        % 1+1D --> 1D
        
    case {'1+1DAP_2_1DAP','1+1DSI_2_1DSI','1+1DRL_2_1DRL'}
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,1];
        Expand = 0;
        % 1D --> 1+1D
        
    case {'1DAP_2_1+1DAP','1DSI_2_1+1DSI','1DRL_2_1+1DRL'}
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,0];
        Expand = 1;
        
        % 1+1D --> 1+1D
        
    case {'1+1DAP_2_1+1DAP','1+1DSI_2_1+1DSI','1+1DRL_2_1+1DRL'}
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = R_i;
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = R_i;
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = D_o;
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,0,0];
        Expand = 0;
        
        %% 2D --> 1D
        
    case '2DAx_2_1DAP'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2),1];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 0;
        
        
    case '2DAx_2_1DRL'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,1];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 0;
        
    case '2DCo_2_1DSI'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,1,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 0;
        
    case '2DCo_2_1DRL'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),1,R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 0;
        
    case '2DSa_2_1DAP'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [1,R_i(2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 0;
        
    case '2DSa_2_1DSI'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [1,R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 0;
        %% 2+1D --> 1D
    case '2+1DAx_2_1DAP'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2),1];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,1];
        Expand = 0;
    case '2+1DAx_2_1DRL'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,1];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,1];
        Expand = 0;
    case '2+1DCo_2_1DSI'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,1,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,1];
        Expand = 0;
    case '2+1DCo_2_1DRL'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),1,R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,1];
        Expand = 0;
        
    case '2+1DSa_2_1DAP'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [1,R_i(2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,1];
        Expand = 0;
    case '2+1DSa_2_1DSI'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [1,R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = R_v;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = [c_pBW;c_pBW];
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = [c_pBW;c_pBW];
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,1];
        Expand = 0;
        %% 2D --> 1+1D
    case '2DAx_2_1+1DAP'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2),1];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 1;
        
        
    case '2DAx_2_1+1DRL'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,1];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 1;
        
    case '2DCo_2_1+1DSI'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,1,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 1;
        
    case '2DCo_2_1+1DRL'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),1,R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 1;
        
    case '2DSa_2_1+1DAP'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [1,R_i(2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 1;
        
    case '2DSa_2_1+1DSI'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [1,R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_o;
        Rv_i_intp2 = Rv_o;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 1;
        
        
        
        %% 2+1D --> 1+1D
        
    case '2+1DAx_2_1+1DAP'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,R_i(2),1];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 0;
        
        
    case '2+1DAx_2_1+1DRL'
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),R_y,1];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 0;
        
    case '2+1DCo_2_1+1DSI'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_x,1,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [1,0,0,0];
        Expand = 0;
        
    case '2+1DCo_2_1+1DRL'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [R_i(1),1,R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [R_i(1),1,1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [D_o(:,1),[r0(2);r0(2)],[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 0;
        
    case '2+1DSa_2_1+1DAP'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [1,R_i(2),R_z];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,R_i(2),1];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],D_o(:,2),[r0(3);r0(3)]];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,0,1,0];
        Expand = 0;
        
    case '2+1DSa_2_1+1DSI'
        
        D_i_intp1 = D_i;
        R_i_intp1 = R_i;
        Dv_i_intp1 = Dv_i;
        Rv_i_intp1 = Rv_i;
        
        D_o_intp1 = D_i;
        R_o_intp1 = [1,R_y,R_i(3)];
        Dv_o_intp1 = Dv_i;
        Rv_o_intp1 = Rv_i;
        
        D_i_intp2 = D_i;
        R_i_intp2 = [1,1,R_i(3)];
        Dv_i_intp2 = Dv_i;
        Rv_i_intp2 = Rv_i;
        
        D_o_intp2 = [[r0(1);r0(1)],[r0(2);r0(2)],D_o(:,3)];
        R_o_intp2 = R_o;
        Dv_o_intp2 = Dv_o;
        Rv_o_intp2 = Rv_o;
        
        
        Collapse = [0,1,0,0];
        Expand = 0;
        
    otherwise
        Status = 0;
end

if Status
%% Do first interpolation
if strcmp(Dim_i,'2DAx') && strcmp(Dim_o,'2+1DAx')
    d = 0;
end
O = Interp(Dim_i,I,D_i_intp1,R_i_intp1,Dv_i_intp1,Rv_i_intp1,D_o_intp1,R_o_intp1,Dv_o_intp1,Rv_o_intp1);

%% Collapse dimensions

O = Collapse_Dimension(O,Collapse,a_x,a_y,a_z,a_v,b_x,b_y,b_z,b_v);

%% Expand spectral

O = Expand_Spectral(O,Expand,Rv_o);

%% Do last Interpolation
% tic
[O,Xi,Yi,Zi,Vi] = Interp(Dim_o,O,D_i_intp2,R_i_intp2,Dv_i_intp2,Rv_i_intp2,D_o_intp2,R_o_intp2,Dv_o_intp2,Rv_o_intp2);
% toc
else
    error('Morph: currently doesn''t support %s as a choice',Morph)
end

end



function MD = Spat_Expand_MD(MD,MDnoroi,R,TH,c_TH,L)
% function MD = Spat_Expand_MD(MD,MDnoroi,R,TH,c_TH,L)
%
%   This script expands the MD array in a spatial dimension
%

idx1 = round((c_TH-TH/2)./L.*R+R/2);

idx2 = round((c_TH+TH/2)./L.*R+R/2);
if idx1 < 1
    idx1 = 1;
end
if idx2 > R
    idx2 = R;
end
if idx2 < 1
    idx2 = 1;
end
if idx1 > R
    idx1 = R;
end
temp1 = repmat(MDnoroi,[1,1,idx1-1,1,1]);
temp2 = repmat(MD,[1,1,idx2-idx1+1,1,1]);
temp3 = repmat(MDnoroi,[1,1,R-idx2,1,1]);

MD = cat(3,temp1,temp2,temp3);
end

function grid = list2grid(list,R,idx,Dim)
% function grid = list2grid(list,R,idx,Dim)
%
%
%

%%

pTx = size(list,2);

if Dim ~= 1 && Dim ~= 3 && Dim ~= 2
    error('list2grid: received a list of points with unidentifiable anatomy.')
end
if Dim == 3 && pTx > 1
    error('list2grid: received a list of points with unidentifiable anatomy.')
end

if numel(R) ~= 4
    error('list2grid: received an improper resolution. It must be [Rx,Ry,Rz,R(4)]')
end







if isempty(idx)
    
    if Dim == 1
        
        grid = zeros(R(2),R(1),R(3),R(4),pTx);
        for s = 1:pTx
            g1 = reshape(list(:,s),R);
            grid(:,:,:,:,s) = g1;
        end
        
        
        
    elseif Dim == 3
        g1 = reshape(list(1:3:end),[R(2),R(1),R(3),R(4)]);
        g2 = reshape(list(2:3:end),[R(2),R(1),R(3),R(4)]);
        g3 = reshape(list(3:3:end),[R(2),R(1),R(3),R(4)]);
        grid = cat(5,g1,g2,g3);
    end
    
    
else
    
    if Dim == 1
        
        grid = zeros(R(2),R(1),R(3),R(4),pTx);
        
        
        if pTx>1
            g1 = zeros(R(2),R(1),R(3),R(4));
            for s = 1:pTx
                g1(idx) = list(:,s);
                grid(:,:,:,:,s) = g1;
            end
        else
            grid(idx) = list;
        end
    elseif Dim == 2
        
        % This is for gridding a list of complex Mxy
        
        
        l1 = real(list);
        l2 = imag(list);
        
        g1_  = zeros(R(2),R(1),R(3),R(4));
        g2_  = zeros(R(2),R(1),R(3),R(4));
        g3_  = zeros(R(2),R(1),R(3),R(4));
        
        
        [sub1,sub2] = ind2sub([R(2),R(1),R(3),R(4)],idx);
        for n = 1:length(idx)
            
            g1_(sub1(n),sub2(n)) = l1(n);
            g2_(sub1(n),sub2(n)) = l2(n);
%             g3_(sub1(n),sub2(n)) = l3(n);
        end
        
        grid = cat(5,g1_,g2_,g3_);
       
        
    elseif Dim == 3
        
        
        
        
        l1 = list(1:3:end);
        l2 = list(2:3:end);
        l3 = list(3:3:end);
        
        g1_  = zeros(R(2),R(1),R(3),R(4));
        g2_  = zeros(R(2),R(1),R(3),R(4));
        g3_  = zeros(R(2),R(1),R(3),R(4));
        
        
        [sub1,sub2] = ind2sub([R(2),R(1),R(3),R(4)],idx);
        for n = 1:length(idx)
            
            g1_(sub1(n),sub2(n)) = l1(n);
            g2_(sub1(n),sub2(n)) = l2(n);
            g3_(sub1(n),sub2(n)) = l3(n);
        end
        
        grid = cat(5,g1_,g2_,g3_);
    end
    
    
end
end


%% Show Functions

function fig = Show_spc(spc)
global hspc

%% Position stuff


    fig = figure;




hspc.spc = spc;
mp = get(0, 'MonitorPositions');
if size(mp,1) > 1
    mp = mp(1,:);
end

set(gcf,'Units','pixels')
set(gcf,'Position',[1,1,mp(3).*0.9,mp(4).*0.9])


% set(gcf,'Units','points')
WHratio = get(gcf,'Position');
WHratio = WHratio(3)/WHratio(4);
set(gcf,'Units','normalized')

Sliderheight = 0.04;
Textheight = 0.03;
Textwidth = 0.04;
Spacer  = 0.005;


hspc.listbox1= uicontrol('Style','listbox','Callback',@listbox1_Callback);
hspc.popupmenu1= uicontrol('Style','popupmenu','Callback',@popupmenu1_Callback);

hspc.axes1 = axes('Parent',gcf);
hspc.axes2 = axes('Parent',gcf);
% hspc.axes3 = axes('Parent',gcf);

hspc.slider1= uicontrol('Style','slider','Callback',@slider1_Callback,'Visible','off');
hspc.slider2= uicontrol('Style','slider','Callback',@slider2_Callback);
hspc.slider3= uicontrol('Style','slider','Callback',@slider3_Callback);

hspc.text1= uicontrol('Style','text','Callback',@text1_Callback);
hspc.text2= uicontrol('Style','text','Callback',@text2_Callback);
hspc.text3= uicontrol('Style','text','Callback',@text3_Callback);
hspc.text4= uicontrol('Style','text','Callback',@text4_Callback);
hspc.text5= uicontrol('Style','text','Callback',@text5_Callback);
hspc.text6= uicontrol('Style','text','Callback',@text6_Callback);
hspc.text7= uicontrol('Style','text','Callback',@text7_Callback);
hspc.text8= uicontrol('Style','text','Callback',@text8_Callback);
hspc.text9= uicontrol('Style','text','Callback',@text9_Callback);


set(hspc.slider1,'Units','normalized')
set(hspc.slider2,'Units','normalized')
set(hspc.slider3,'Units','normalized')

set(hspc.text1,'Units','normalized')
set(hspc.text2,'Units','normalized')
set(hspc.text3,'Units','normalized')
set(hspc.text4,'Units','normalized')
set(hspc.text5,'Units','normalized')
set(hspc.text6,'Units','normalized')
set(hspc.text7,'Units','normalized')
set(hspc.text9,'Units','normalized')
set(hspc.axes1,'Units','normalized')
set(hspc.axes2,'Units','normalized')
% set(hspc.axes3,'Units','normalized')

set(hspc.listbox1,'Units','normalized')
set(hspc.popupmenu1,'Units','normalized')

r.x_hspc.slider1 = 0.65;
r.x_hspc.slider2 = 0.65;
r.x_hspc.slider3 = 0.65;

r.y_hspc.slider1 = 0.3-Spacer-Textheight-Spacer-Sliderheight;
r.y_hspc.slider2 = 0.3-Spacer-Textheight-Spacer-Sliderheight-Spacer-Textheight-Spacer-Sliderheight;
r.y_hspc.slider3 = 0.3-Spacer-Textheight-Spacer-Sliderheight-Spacer-Textheight-Spacer-Sliderheight-Spacer-Textheight-Spacer-Sliderheight;

r.w_hspc.slider1 = 0.3;
r.w_hspc.slider2 = 0.3;
r.w_hspc.slider3 = 0.3;

r.h_hspc.slider1 = Sliderheight;
r.h_hspc.slider2 = Sliderheight;
r.h_hspc.slider3 = Sliderheight;
% r.x_axes1 = 0.05;
r.x_axes1 = 0.1;
r.xo_axes1 = 0.05;
r.x_axes2 = 0.05;
% r.x_axes3 = 0.05;

% r.y_axes1 = 0.3;
r.yo_axes1 = 0.1;
r.y_axes1 = 0.2;
% r.y_axes2 = 0.25-Spacer-0.05;
r.y_axes2 = 0.05;
% r.y_axes3 = 0.25-Spacer-0.05-Spacer-Sliderheight-0.05;

% r.w_axes1 = 0.5;
r.w_axes2 = 0.5;
% r.w_axes3 = 0.5;

% r.h_axes1 = 0.5;
r.ho_axes1 = 0.8;
r.h_axes1 = 0.7;
r.wo_axes1 = 0.8/WHratio;

r.w_axes1 = 0.7/WHratio;
% r.h_axes2 = 0.05;
r.h_axes2 = 0.025;
% r.h_axes3 = 0.05;


r.x_text1 = 0.65;
r.x_text2 = 0.65;
r.x_text3 = 0.65;
r.x_text4 = 0;
r.x_text5 = 0.05+0.5+Spacer;

r.x_text6 = 0;
r.x_text7 = 0.05+0.5+Spacer;

r.x_text9 = 0.1;
r.y_text1 = 0.3-Textheight-Spacer;
r.y_text2 = 0.3-Textheight-Spacer-Sliderheight-Spacer-Textheight-Spacer;
r.y_text3 = 0.3-Textheight-Spacer-Sliderheight-Spacer-Textheight-Spacer-Sliderheight-Spacer-Textheight-Spacer;
% r.y_text4 = 0.25-Textheight-Spacer;
% r.y_text5 = 0.25-Textheight-Spacer;
r.y_text4 = 0.05;
r.y_text5 = 0.05;
r.y_text9 = 0.0;
r.y_text6 = 0.25-Spacer-Textheight-Spacer-Sliderheight-Spacer-Textheight;
r.y_text7 = 0.25-Spacer-Textheight-Spacer-Sliderheight-Spacer-Textheight;
r.w_text1 = 0.3;
r.w_text2 = 0.3;
r.w_text3 = 0.3;
r.w_text4 = Textwidth;
r.w_text5 = Textwidth;
r.w_text6 = Textwidth;
r.w_text7 = Textwidth;
r.w_text9 = 1;
r.h_text1 = Textheight;
r.h_text2 = Textheight;
r.h_text3 = Textheight;
r.h_text4 = Textheight;
r.h_text5 = Textheight;
r.h_text6 = Textheight;
r.h_text7 = Textheight;
r.h_text9 = Textheight;

r.x_popupmenu1 = 0.05;
r.y_popupmenu1 = 0.92;
r.w_popupmenu1 = 0.5;
r.h_popupmenu1 = 0.05;

r.x_listbox1 = 0.65;
r.y_listbox1 = 0.3;
r.w_listbox1 = 0.3;
r.h_listbox1 = 0.65;

set(hspc.slider1,'Position',[r.x_hspc.slider1,r.y_hspc.slider1,r.w_hspc.slider1,r.h_hspc.slider1])
set(hspc.slider2,'Position',[r.x_hspc.slider2,r.y_hspc.slider2,r.w_hspc.slider2,r.h_hspc.slider2])
set(hspc.slider3,'Position',[r.x_hspc.slider3,r.y_hspc.slider3,r.w_hspc.slider3,r.h_hspc.slider3])

set(hspc.text1,'Position',[r.x_text1,r.y_text1,r.w_text1,r.h_text1])
set(hspc.text2,'Position',[r.x_text2,r.y_text2,r.w_text2,r.h_text2])
set(hspc.text3,'Position',[r.x_text3,r.y_text3,r.w_text3,r.h_text3])
set(hspc.text4,'Position',[r.x_text4,r.y_text4,r.w_text4,r.h_text4])
set(hspc.text5,'Position',[r.x_text5,r.y_text5,r.w_text5,r.h_text5])
set(hspc.text6,'Position',[r.x_text6,r.y_text6,r.w_text6,r.h_text6])
set(hspc.text7,'Position',[r.x_text7,r.y_text7,r.w_text7,r.h_text7])
set(hspc.text9,'Position',[r.x_text9,r.y_text9,r.w_text9,r.h_text9])
set(hspc.text9,'String','blOCh uses the LPS coordinate system (R --x--> L, A --y--> P, I --z--> S, (i,j,k) = (y,x,z))')

set(hspc.axes1,'OuterPosition',[r.xo_axes1,r.yo_axes1,r.wo_axes1,r.ho_axes1])
set(hspc.axes1,'Position',[r.x_axes1,r.y_axes1,r.w_axes1,r.h_axes1])
set(hspc.axes2,'Position',[r.x_axes2,r.y_axes2,r.w_axes2,r.h_axes2])
% set(hspc.axes3,'Position',[r.x_axes3,r.y_axes3,r.w_axes3,r.h_axes3])
% get(hspc.text1,'Fontsize')
set(hspc.listbox1,'Position',[r.x_listbox1,r.y_listbox1,r.w_listbox1,r.h_listbox1])
set(hspc.popupmenu1,'Position',[r.x_popupmenu1,r.y_popupmenu1,r.w_popupmenu1,r.h_popupmenu1])

% set(hspc.axes3,'Visible','off')
set(hspc.text6,'Visible','off')
set(hspc.text7,'Visible','off')


if hspc.spc.Dim(1:2) == '1D'
    
    set(hspc.text4,'Visible','off')
    set(hspc.text5,'Visible','off')
    set(hspc.axes2,'Visible','off')
elseif hspc.spc.Dim(1:2) == '0+'
    set(hspc.axes2,'Visible','off')
    set(hspc.text4,'Visible','off')
    set(hspc.text5,'Visible','off')
end


if hspc.spc.pTx > 1
    set(hspc.slider1,'Min',1)
    set(hspc.slider1,'Max',hspc.spc.pTx)
    ChannelSliderStep = [1, 1] / (hspc.spc.pTx - 1);
    set(hspc.slider1,'SliderStep',ChannelSliderStep)
    set(hspc.slider1,'Visible','off','Value',1)
    set(hspc.text1,'Visible','off','String',sprintf('Tx %i of %i',get(hspc.slider1,'Value'),hspc.spc.pTx))
    hspc.pTxNo = 1;

else
    set(hspc.slider1,'Visible','off','Value',1)
    set(hspc.text1,'Visible','off');
    hspc.pTxNo = 1;
end

if hspc.spc.Dim(1) == '3'
    set(hspc.slider3,'Min',1)
    set(hspc.slider3,'Max',hspc.spc.R(3))
    SliceSliderStep = [1, 1] / (hspc.spc.R(3) - 1);
    set(hspc.slider3,'SliderStep',SliceSliderStep)
    set(hspc.slider3,'Visible','on','Value',1)
    
    hspc.SlNo = 1;
    
    hspc.Z = linspace(hspc.spc.D(1,3),hspc.spc.D(2,3),hspc.spc.R(3));
    set(hspc.text3,'String',sprintf('Slice %i of %i (%f m)',hspc.SlNo,hspc.spc.R(3),hspc.Z(hspc.SlNo)))
    set(hspc.text3,'Visible','on')
    

else
    set(hspc.slider3,'Visible','off','Value',1)
    set(hspc.text3,'Visible','off');
    hspc.SlNo = 1;
end

if hspc.spc.Dim(2) == '+'
    if hspc.spc.Dim(1) == '1'
        set(hspc.slider2,'Visible','off','Value',1)
        set(hspc.text2,'Visible','off');
        hspc.freqNo = 1;
        hspc.freq = linspace(hspc.spc.Dv(1),hspc.spc.Dv(2),hspc.spc.Rv);
    elseif hspc.spc.Dim(1) == '0'
        set(hspc.slider2,'Visible','off','Value',1)
        set(hspc.text2,'Visible','off');
        hspc.freqNo = 1;
        hspc.freq = linspace(hspc.spc.Dv(1),hspc.spc.Dv(2),hspc.spc.Rv);
    else
        hspc.freq = linspace(hspc.spc.Dv(1),hspc.spc.Dv(2),hspc.spc.Rv);
        
        set(hspc.slider2,'Min',1)
        set(hspc.slider2,'Max',hspc.spc.Rv)
        SliceSliderStep = [1, 1] / (hspc.spc.Rv - 1);
        set(hspc.slider2,'SliderStep',SliceSliderStep)
        hspc.freqNo = round(hspc.spc.Rv/2);
        set(hspc.slider2,'Visible','on','Value',hspc.freqNo)
        set(hspc.text2,'Visible','on','String',sprintf('Frequency: %.2f [Hz]',hspc.freq(hspc.freqNo)))
        

    end
else
    set(hspc.slider2,'Visible','off','Value',1)
    set(hspc.text2,'Visible','off');
    hspc.freqNo = 1;
end

Populate_Listbox;

String = cell(1,21);
if isfield(hspc.spc,'Md')
    String{1} = 'Desired Magnetization, |Mxy| [M0]';
    String{2} = 'Desired Magnetization, Mx [M0]';
    String{3} = 'Desired Magnetization, My [M0]';
    String{4} = 'Desired Magnetization, Mz [M0]';
    String{5} = 'Desired Flipangle, [deg.]';
    String{6} = 'Desired Flipphase, [deg.]';
    
end
if isfield(hspc.spc,'M0')
    String{7} = 'Initial Flipangle, [deg.]';
    String{8} = 'Initial Flipphase, [deg.]';
    String{9} = 'Initial Magnetization, |Mxy| [M0]';
    String{10} = 'Initial Magnetization, Mx [M0]';
    String{11} = 'Initial Magnetization, My [M0]';
    String{12} = 'Initial Magnetization, Mz [M0]';
end
if isfield(hspc.spc,'w0map')
    String{13} = 'B0 map, v0 [Hz]';
end
if isfield(hspc.spc,'B1map')
    
    String{14} = 'B1 sensitivity map, |B1|  [norm.]';
    String{15} = 'B1 sensitivity map, arg(B1)  [rad.]';
    String{16} = 'B1 sensitivity map, Re(B1)  [norm.]';
    String{17} = 'B1 sensitivity map, Im(B1)  [norm.]';
end
if isfield(hspc.spc,'T1map')
    
    String{18} = 'T1 map, [s]';
    String{19} = 'R1 map, [1/s]';
end
if isfield(hspc.spc,'T2map')
    
    String{20} = 'T2 map, [s]';
    String{21} = 'R2 map, [1/s]';
end
if isfield(hspc.spc,'Mpen')
String{22} = 'Penalty Magnetization 1, Mx [M0]';
String{23} = 'Penalty Magnetization 1, My [M0]';
String{24} = 'Penalty Magnetization 1, Mz [M0]';

String{25} = 'Penalty Magnetization 2, Mx [M0]';
String{26} = 'Penalty Magnetization 2, My [M0]';
String{27} = 'Penalty Magnetization 2, Mz [M0]';
end


set(hspc.popupmenu1, 'String', String);
set(hspc.popupmenu1, 'Value',1);




[hspc.colmap.Jet,hspc.colmap.Gray] = GetColormaps;

hspc.Colormapno = 1;



hObj.hspc = hspc;


if iscolumn(hspc.spc.X) % Is this right?? What happens in 1D AP/ LR case
    hspc.spc.M0 = blOCh__3_7__list2grid(hspc.spc.M0,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,3);
    hspc.spc.Md = blOCh__3_7__list2grid(hspc.spc.Md,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,3);
    hspc.spc.X = blOCh__3_7__list2grid(hspc.spc.X,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,1);
    hspc.spc.Y = blOCh__3_7__list2grid(hspc.spc.Y,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,1);
    hspc.spc.Z = blOCh__3_7__list2grid(hspc.spc.Z,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,1);
    hspc.spc.V = blOCh__3_7__list2grid(hspc.spc.V,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,1);
    hspc.spc.B1map = blOCh__3_7__list2grid(hspc.spc.B1map,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,1);
    hspc.spc.w0map = blOCh__3_7__list2grid(hspc.spc.w0map,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,1);
    hspc.spc.T1map = blOCh__3_7__list2grid(hspc.spc.T1map,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,1);
    hspc.spc.T2map = blOCh__3_7__list2grid(hspc.spc.T2map,[hspc.spc.R,hspc.spc.Rv],hspc.spc.idxtot,1);
    
else
    Plotting;
end


guidata(fig,hspc);




end

function popupmenu1_Callback(hObj, evd)
global hspc
hspc.popupmenu1select = get(hObj,'Value');
Plotting;
hspc.output = hObj;
guidata(hObj);

end

function listbox1_Callback(hObj, evd)
global hspc
end

function slider1_Callback(hObj, evd)
global hspc
hspc.pTxNo = round(get(hObj, 'Value'));
set(hspc.text1,'String',sprintf('Tx %i of %i',hspc.pTxNo,hspc.spc.pTx))
Plotting;
hspc.output = hObj;
guidata(hObj);
end

function slider2_Callback(hObj, evd)
global hspc
hspc.freqNo = round(get(hObj, 'Value'));
set(hObj, 'Value',hspc.freqNo);
set(hspc.text2,'Visible','on','String',sprintf('Frequency: %.2f [Hz]',hspc.freq(hspc.freqNo)))
Plotting;
hspc.output = hObj;
guidata(hObj);
end
function slider3_Callback(hObj, evd)
global hspc
hspc.SlNo = round(get(hObj, 'Value'));

set(hspc.text3,'String',sprintf('Slice %i of %i (%f m)',hspc.SlNo,hspc.spc.R(3),hspc.Z(hspc.SlNo)))
Plotting;
hspc.output = hObj;
guidata(hObj);
end


function Populate_Listbox
global hspc

names = fieldnames(hspc.spc);
vals = cell(length(names),1);
for n = 1:length(names)
    
    if ischar(getfield(hspc.spc,names{n}))
        vals{n} = getfield(hspc.spc,names{n});
    elseif isnumeric(getfield(hspc.spc,names{n}))
        [A,B,C,D,E] = size(getfield(hspc.spc,names{n}));
        % 		[A,B,C,D,E]
        
        if E == 1
            
            if D == 1 && C == 1 && B == 1 && A == 1
                vals{n} = num2str(getfield(hspc.spc,names{n}));
            elseif D == 1 && C == 1 && B ~= 1 && A == 1
                if B < 6
                    temp = getfield(hspc.spc,names{n});
                    
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
                    temp = getfield(hspc.spc,names{n});
                    
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
                    temp = getfield(hspc.spc,names{n});
                    
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
        
        
    elseif iscell(getfield(hspc.spc,names{n}))
        [A,B,C,D,E] = size(getfield(hspc.spc,names{n}));
        if E == 1
            
            if D == 1 && C == 1 && B == 1 && A == 1
                vals{n} = num2str(getfield(hspc.spc,names{n}));
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

    elseif isstruct(getfield(hspc.spc,names{n}))
        [A,B] = size(getfield(hspc.spc,names{n}));
        vals{n} = sprintf('<%ix%i struct>',A,B);
    end
end

List = cell(length(names),1);

for n = 1:length(names)
    List{n} = sprintf('%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s',names{n},vals{n});
    
end

set(hspc.listbox1,'String',List)





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

function Plotting
global hspc
what2plot;


switch hspc.spc.Dim(1:2)
    

    case {'1D','0+'}
        
        
        axes(hspc.axes1)
        plot(hspc.ordinate)
        set(hspc.axes1,'XTick',hspc.Dim1axis)
        set(hspc.axes1,'XTickLabel',hspc.Dim1ticklabel)
        xlabel(hspc.Dim2label)
        %       ylabel(Ylabel)
        
    case {'1+','2D','2+','3D','3+'}
        axes(hspc.axes1)
        imagesc(hspc.ordinate,[hspc.ordinate_mn,hspc.ordinate_mx])
        %       caxis([hspc.ordinate_mn_win,hspc.ordinate_mx_win])
        
        set(hspc.axes1,'XTick',hspc.Dim2axis,'YTick',hspc.Dim1axis)
        set(hspc.axes1,'XTickLabel',hspc.Dim2ticklabel,'YTickLabel',hspc.Dim1ticklabel)
        %       set(hspc.axes1,'Ydir','normal')
        xlabel(hspc.Dim2label)
        ylabel(hspc.Dim1label)
        
        
        if hspc.Colormapno == 1
            Colmap = hspc.colmap.Gray;
            colormap gray(512);
        else
            Colmap = hspc.colmap.Jet;
            colormap jet(512);
        end
        set(gca,'YDir','normal')
        axes(hspc.axes2)
        imagesc(permute(Colmap,[2,1,3])), axis off
end



end

function what2plot
global hspc
% get(popupmenu1,'Value')
set(hspc.slider1,'Visible','off')
set(hspc.text1,'Visible','off')
switch get(hspc.popupmenu1,'Value')
    case 1
        
        hspc.Colormapno = 1;
        hspc.ordinate = abs(complex(hspc.spc.Md(:,:,:,:,1),hspc.spc.Md(:,:,:,:,2)));
        hspc.ordinate_mn = 0;
        hspc.ordinate_mx = 1;
        %       'Desired Magnetization, |Mxy| [M0]';
    case 2
        hspc.Colormapno = 1;
        %       'Desired Magnetization, Mx [M0]';
        hspc.ordinate = hspc.spc.Md(:,:,:,:,1);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 3
        hspc.Colormapno = 1;
        %       'Desired Magnetization, My [M0]';
        hspc.ordinate = hspc.spc.Md(:,:,:,:,2);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 4
        %       'Desired Magnetization, Mz [M0]';
        hspc.ordinate = hspc.spc.Md(:,:,:,:,3);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 5
        hspc.Colormapno = 1;
        %       'Desired Flipangle, [deg.]';
        hspc.ordinate = ones(hspc.spc.R(2),hspc.spc.R(1),hspc.spc.R(3),hspc.spc.Rv).*hspc.spc.Mdflip(1);
        hspc.ordinate_mn = -180;
        hspc.ordinate_mx = 180;
    case 6
        hspc.Colormapno = 1;
        %       'Desired Flipphase, [deg.]';
        
        hspc.ordinate = ones(hspc.spc.R(2),hspc.spc.R(1),hspc.spc.R(3),hspc.spc.Rv).*hspc.spc.Mdflip(2);
        hspc.ordinate_mn = -180;
        hspc.ordinate_mx = 180;
    case 7
        hspc.Colormapno = 1;
        %       'Initial Flipangle, [deg.]';
        hspc.ordinate = ones(hspc.spc.R(2),hspc.spc.R(1),hspc.spc.R(3),hspc.spc.Rv).*hspc.spc.M0flip(1);
        hspc.ordinate_mn = -180;
        hspc.ordinate_mx = 180;
    case 8
        hspc.Colormapno = 1;
        %       'Initial Flipphase, [deg.]';
        hspc.ordinate = ones(hspc.spc.R(2),hspc.spc.R(1),hspc.spc.R(3),hspc.spc.Rv).*hspc.spc.M0flip(2);
        hspc.ordinate_mn = -180;
        hspc.ordinate_mx = 180;
    case 9
        hspc.Colormapno = 1;
        %       'Initial Magnetization, |Mxy| [M0]';
        hspc.ordinate = abs(complex(hspc.spc.M0(:,:,:,:,1),hspc.spc.M0(:,:,:,:,2)));
        hspc.ordinate_mn = 0;
        hspc.ordinate_mx = 1;
    case 10
        hspc.Colormapno = 1;
        %       'Initial Magnetization, Mx [M0]';
        hspc.ordinate = hspc.spc.M0(:,:,:,:,1);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 11
        hspc.Colormapno = 1;
        %       'Initial Magnetization, My [M0]';
        hspc.ordinate = hspc.spc.M0(:,:,:,:,2);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 12
        hspc.Colormapno = 1;
        %       'Initial Magnetization, Mz [M0]';
        hspc.ordinate = hspc.spc.M0(:,:,:,:,3);
        %       hspc.ordinate_mn = -1;
        %       hspc.ordinate_mx = 1;
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = min(hspc.ordinate(:));
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mn = hspc.ordinate_mn-eps;
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 13
        hspc.Colormapno = 2;
        %       'B0 map, v0 [Hz]';
        hspc.ordinate = hspc.spc.w0map./2/pi;
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = min(hspc.ordinate(:));
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mn = hspc.ordinate_mn-eps;
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 14
        hspc.Colormapno = 2;
        %       'B1 sensitivity map, |B1|  [norm.]';
        if hspc.spc.pTx > 1
            set(hspc.slider1,'Visible','on')
            set(hspc.text1,'Visible','on')
        else
            set(hspc.slider1,'Visible','off')
            set(hspc.text1,'Visible','off')
        end
        hspc.ordinate = abs(hspc.spc.B1map(:,:,:,:,hspc.pTxNo));
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = 0;
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 15
        hspc.Colormapno = 2;
        %       'B1 sensitivity map, arg(B1)  [rad.]';
        if hspc.spc.pTx > 1
            set(hspc.slider1,'Visible','on')
            set(hspc.text1,'Visible','on')
        else
            set(hspc.slider1,'Visible','off')
            set(hspc.text1,'Visible','off')
        end
        hspc.ordinate = angle(hspc.spc.B1map(:,:,:,:,hspc.pTxNo));
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = min(hspc.ordinate(:));
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mn = hspc.ordinate_mn-eps;
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
        %       hspc.ordinate_mn = -pi;
        %       hspc.ordinate_mx = pi;
    case 16
        hspc.Colormapno = 2;
        %       B1 sensitivity map, Re(B1)  [norm.]';
        if hspc.spc.pTx > 1
            set(hspc.slider1,'Visible','on')
            set(hspc.text1,'Visible','on')
        else
            set(hspc.slider1,'Visible','off')
            set(hspc.text1,'Visible','off')
        end
        hspc.ordinate = real(hspc.spc.B1map(:,:,:,:,hspc.pTxNo));
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = min(hspc.ordinate(:));
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mn = hspc.ordinate_mn-eps;
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 17
        hspc.Colormapno = 2;
        %       B1 sensitivity map, Im(B1)  [norm.]';
        if hspc.spc.pTx > 1
            set(hspc.slider1,'Visible','on')
            set(hspc.text1,'Visible','on')
        else
            set(hspc.slider1,'Visible','off')
            set(hspc.text1,'Visible','off')
        end
        hspc.ordinate = imag(hspc.spc.B1map(:,:,:,:,hspc.pTxNo));
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = min(hspc.ordinate(:));
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mn = hspc.ordinate_mn-eps;
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 18
        hspc.Colormapno = 2;
        %       'T1 map, [s]';
        hspc.ordinate = hspc.spc.T1map;
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = 0;
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 19
        hspc.Colormapno = 2;
        %       'R1 map, [1/s]';
        hspc.ordinate = 1./hspc.spc.T1map;
        hspc.ordinate(isnan(hspc.ordinate)) = 0;
        hspc.ordinate(isinf(hspc.ordinate)) = 0;
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = 0;
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 20
        hspc.Colormapno = 2;
        %       'T2 map, [s]';
        hspc.ordinate = hspc.spc.T2map;
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = 0;
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 21
        hspc.Colormapno = 2;
        %       'R2 map, [1/s]';
        hspc.ordinate = 1./hspc.spc.T2map;
        hspc.ordinate(isnan(hspc.ordinate)) = 0;
        hspc.ordinate(isinf(hspc.ordinate)) = 0;
        hspc.ordinate_mx = max(hspc.ordinate(:));
        hspc.ordinate_mn = 0;
        if hspc.ordinate_mn == hspc.ordinate_mx
            hspc.ordinate_mx = hspc.ordinate_mx+eps;
        end
    case 22
        hspc.Colormapno = 1;
        %       'Penalty Magnetization, Mx [M0]';
        hspc.ordinate = hspc.spc.Mp1(:,:,:,:,1);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 23
        hspc.Colormapno = 1;
        %       'Penalty Magnetization, My [M0]';
        hspc.ordinate = hspc.spc.Mp1(:,:,:,:,2);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 24
        %       'Penalty Magnetization, Mz [M0]';
        hspc.ordinate = hspc.spc.Mp1(:,:,:,:,3);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 25
        hspc.Colormapno = 1;
        %       'Penalty Magnetization, Mx [M0]';
        hspc.ordinate = hspc.spc.Mp2(:,:,:,:,1);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 26
        hspc.Colormapno = 1;
        %       'Penalty Magnetization, My [M0]';
        hspc.ordinate = hspc.spc.Mp2(:,:,:,:,2);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;
    case 27
        %       'Penalty Magnetization, Mz [M0]';
        hspc.ordinate = hspc.spc.Mp2(:,:,:,:,3);
        hspc.ordinate_mn = -1;
        hspc.ordinate_mx = 1;    
        
        
end

hspc.ordinate = hspc.ordinate;
set(hspc.text4,'String',sprintf('%1.1e',hspc.ordinate_mn));
set(hspc.text5,'String',sprintf('%1.1e',hspc.ordinate_mx));

switch hspc.spc.Dim
    
    case '0+1D'
        hspc.Dim1ticklabel = linspace(hspc.spc.Dv(1),hspc.spc.Dv(2),5);
        hspc.Dim1axis = linspace(1,hspc.spc.Rv,5);
        hspc.Dim2label = ' v [Hz]';
        hspc.ordinate = squeeze(hspc.ordinate);
        
    case '1DSI'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,3),hspc.spc.D(2,3),5);
        %       hspc.Dim1ticklabel = linspace(hspc.ordinate_mn,hspc.ordinate_mx,5);
        
        hspc.Dim1axis = linspace(1,hspc.spc.R(3),5);
        hspc.Dim2label = ' z [m]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '1DAP'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,2),hspc.spc.D(2,2),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(2),5);
        hspc.Dim2label = ' y [m]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '1DRL'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,1),hspc.spc.D(2,1),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(1),5);
        hspc.Dim2label = ' x [m]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '1+1DSI'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,3),hspc.spc.D(2,3),5);
        hspc.Dim2ticklabel = linspace(hspc.spc.Dv(1),hspc.spc.Dv(2),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(3),5);
        hspc.Dim2axis = linspace(1,hspc.spc.Rv,5);
        hspc.Dim1label = ' z [m]';
        hspc.Dim2label = ' v [Hz]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '1+1DAP'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,3),hspc.spc.D(2,3),5);
        hspc.Dim2ticklabel = linspace(hspc.spc.Dv(1),hspc.spc.Dv(2),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(2),5);
        hspc.Dim2axis = linspace(1,hspc.spc.Rv,5);
        hspc.Dim1label = ' y [m]';
        hspc.Dim2label = ' v [Hz]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '1+1DRL'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,1),hspc.spc.D(2,1),5);
        hspc.Dim2ticklabel = linspace(hspc.spc.Dv(1),hspc.spc.Dv(2),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(1),5);
        hspc.Dim2axis = linspace(1,hspc.spc.Rv,5);
        hspc.Dim1label = ' x [m]';
        hspc.Dim2label = ' v [Hz]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '2DAx'
        hspc.Dim2ticklabel = linspace(hspc.spc.D(1,1),hspc.spc.D(2,1),5);
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,2),hspc.spc.D(2,2),5);
        hspc.Dim2axis = linspace(1,hspc.spc.R(1),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(2),5);
        hspc.Dim1label = ' y [m]';
        hspc.Dim2label = ' x [m]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '2DCo'
        hspc.Dim2ticklabel = linspace(hspc.spc.D(1,1),hspc.spc.D(2,1),5);
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,3),hspc.spc.D(2,3),5);
        hspc.Dim2axis = linspace(1,hspc.spc.R(1),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(3),5);
        
        hspc.Dim2label = ' x [m]';
        hspc.Dim1label = ' z [m]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '2DSa'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,2),hspc.spc.D(2,2),5);
        hspc.Dim2ticklabel = linspace(hspc.spc.D(1,3),hspc.spc.D(2,3),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(2),5);
        hspc.Dim2axis = linspace(1,hspc.spc.R(3),5);
        hspc.Dim2label = ' y [m]';
        hspc.Dim1label = ' z [m]';
        hspc.ordinate = squeeze(hspc.ordinate);
    case '2+1DAx'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,1),hspc.spc.D(2,1),5);
        hspc.Dim2ticklabel = linspace(hspc.spc.D(1,2),hspc.spc.D(2,2),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(1),5);
        hspc.Dim2axis = linspace(1,hspc.spc.R(2),5);
        hspc.Dim2label = ' y [m]';
        hspc.Dim1label = ' x [m]';
        hspc.ordinate = squeeze(hspc.ordinate(:,:,:,hspc.freqNo));
    case '2+1DCo'
        hspc.Dim2ticklabel = linspace(hspc.spc.D(1,1),hspc.spc.D(2,1),5);
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,3),hspc.spc.D(2,3),5);
        hspc.Dim2axis = linspace(1,hspc.spc.R(1),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(3),5);
        
        hspc.Dim2label = ' x [m]';
        hspc.Dim1label = ' z [m]';
        hspc.ordinate = squeeze(hspc.ordinate(:,:,:,hspc.freqNo));
    case '2+1DSa'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,2),hspc.spc.D(2,2),5);
        hspc.Dim2ticklabel = linspace(hspc.spc.D(1,3),hspc.spc.D(2,3),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(2),5);
        hspc.Dim2axis = linspace(1,hspc.spc.R(3),5);
        hspc.Dim2label = ' y [m]';
        hspc.Dim1label = ' z [m]';
        hspc.ordinate = squeeze(hspc.ordinate(:,:,:,hspc.freqNo));
    case '3D'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,1),hspc.spc.D(2,1),5);
        hspc.Dim2ticklabel = linspace(hspc.spc.D(1,2),hspc.spc.D(2,2),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(1),5);
        hspc.Dim2axis = linspace(1,hspc.spc.R(2),5);
        hspc.Dim2label = ' x [m]';
        hspc.Dim1label = ' y [m]';
        hspc.ordinate = squeeze(hspc.ordinate(:,:,hspc.SlNo,1));
    case '3+1D'
        hspc.Dim1ticklabel = linspace(hspc.spc.D(1,1),hspc.spc.D(2,1),5);
        hspc.Dim2ticklabel = linspace(hspc.spc.D(1,2),hspc.spc.D(2,2),5);
        hspc.Dim1axis = linspace(1,hspc.spc.R(1),5);
        hspc.Dim2axis = linspace(1,hspc.spc.R(2),5);
        hspc.Dim2label = ' x [m]';
        hspc.Dim1label = ' y [m]';
        hspc.ordinate = squeeze(hspc.ordinate(:,:,hspc.SlNo,hspc.freqNo));
end
end