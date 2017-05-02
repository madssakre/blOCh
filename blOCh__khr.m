function varargout = blOCh__khr(Fun,spc,Traj,Par,varargin)
% function varargout = blOCh__khr(Fun,spc,Traj,varargin)
%
%   This script makes the time dependent input data structures for blOCh
%
%
%     If setting up a new scanner system:
%
%       Search this file for the comment "SYS" and do as told there
%
%     If implementing a new trajectory:
%
%       Search this file for the comment "IMP" and do as told there.
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
    if spc.Status
        Msg = '';
        khr.Status = 1;
        khr.Print = 1;%spc.Print;
        if isunix
            
            Sep = '/';
        else
            Sep = '\';
        end
    else
        khr = [];
        khr.Status = 0;
    end
    
    if khr.Status
        khr.Status = 0;

    valid_Trajs   = {'none4spect','none4spect2','1Dplat','1Dsl','1Dslref','2Dspii','2Dspio','1+1DTrapOsc','Random'};
        

        
        
        % SYS: Put the name of the system into this cell array
        
        valid_Systems{1} = 'Unity';
        valid_Systems{2} = '3T';
        valid_Systems{3} = '7T';
        
        % SYS: If that system is supposed to be default. Change the number
        % here accordingly.
        
        default_System    = valid_Systems{2};
        
        % SYS: Input the max gradient amplitude of the hardware [T/m]
        
        actual_GmHW(1) = 1;
        actual_GmHW(2) = 40e-3;
        actual_GmHW(3) = 40e-3;
        
        
        % SYS: Input the max gradient slewrate of the hardware [T/m/s]
        
        actual_SmHW(1) = 1;
        actual_SmHW(2) = 150;
        actual_SmHW(3) = 170;
        
        % SYS: Input the max RF amplitude of the hardware [Hz]
        actual_RFmHW(1) = 1;
        actual_RFmHW(2) = 700;
        actual_RFmHW(3) = 700;
        
        
        % SYS: If you need another nuclei, input here
        
        valid_Nuclei      = {'1H','13C','19F'};
        
        default_Nuc = valid_Nuclei{1};
        
        switch spc.Dim
            case {'1DSI','1DAP','1DRL'}
                Dimkhr = [1,1];
            case {'2DAx','2DCo','2DSa'}
                Dimkhr = [1,2];
            case {'3D'}
                Dimkhr = [1,3];
            case {'1+1DSI','1+1DAP','1+1DRL'}
                Dimkhr = [1,1];
            case {'2+1DAx','2+1DCo','2+1DSa'}
                Dimkhr = [1,2];
            case '3+1D'
                Dimkhr = [1,3];
            case '0+1D'
                Dimkhr = [0,0];
        end
        
        
        
        
        
        %% GENERAL
        
        
        valid_MaskOnlys = [0,1];
        
        default_R = ceil(spc.R/2);
        default_Rv = ceil(spc.Rv/2);
        default_L = spc.L;
        
        default_Gmpct     = 100;
        default_Smpct     = 100;
        default_GmHW      = actual_GmHW(1);
        default_SmHW      = actual_SmHW(1);
        default_RFmHW     = actual_RFmHW(1);
        default_Show      = 1;
        default_dt        = 10e-6;
        default_MaskOnly        = valid_MaskOnlys(1);
        
        
        
        if ~isempty(Msg)
            Display_Message(Msg,khr.Print);
        else
            khr.Status = 1;
        end
    end
    
    if khr.Status
        khr.Status = 0;
        p = inputParser;
        try p.addRequired('Fun', @(x)isempty(x));
            try p.addRequired('spc', @(x)isstruct(x));
                try p.addRequired('Traj',@(x)Validate_Traj(x,valid_Trajs));
                    try p.addRequired('Par', @(x)isstruct(x)||isempty(x));
                        try p.addParamValue('L',default_L,@(x)validateattributes(x,{'numeric'},{'size',Dimkhr,'real','>',0}));
                            try p.addParamValue('R',default_R,@(x)validateattributes(x,{'numeric'},{'size',Dimkhr,'integer','nonnegative','finite'}));
                                try p.addParamValue('Rv',default_Rv,@(x)validateattributes(x,{'numeric'},{'size',Dimkhr,'integer','nonnegative','finite'}));
                                    
                                    try p.addParamValue('Nuc',default_Nuc,@(x)any(strcmpi(x,valid_Nuclei)));
                                        try p.addParamValue('Sys',default_System,@(x)any(strcmpi(x,valid_Systems)));
                                            try p.addParamValue('Gmpct',default_Gmpct,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>',0,'<=',1000}));
                                                try p.addParamValue('Smpct',default_Smpct,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>',0,'<=',1000}));
                                                    try p.addParamValue('GmHW',default_GmHW,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>=',0}));
                                                        try p.addParamValue('SmHW',default_SmHW,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>=',0}));
                                                            try p.addParamValue('RFmHW',default_RFmHW,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>=',0}));
                                                                try p.addParamValue('dt',default_dt,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','>',0,'finite'}));
                                                                    try p.addParamValue('Show',default_Show,@(x)validateattributes(x,{'numeric'},{'size',[1,1],'real','integer','nonnegative'}));
                                                                        try p.addParamValue('MaskOnly',default_MaskOnly,@(x)any(valid_MaskOnlys));
                                                                            catch me;Msg = ['blOCh__khr: ',me.message];end
                                                                    catch me;Msg = ['blOCh__khr: ',me.message];end
                                                                catch me;Msg = ['blOCh__khr: ',me.message];end
                                                            catch  me;Msg = ['blOCh__khr: ',me.message];end
                                                        catch  me; Msg = ['blOCh__khr: ',me.message];end
                                                    catch  me;Msg = ['blOCh__khr: ',me.message];end
                                                catch  me;Msg = ['blOCh__khr: ',me.message];end
                                            catch  me;Msg = ['blOCh__khr: ',me.message];end
                                        catch  me;Msg = ['blOCh__khr: ',me.message];end
                                    catch  me;Msg = ['blOCh__khr: ',me.message];end
                                catch  me;Msg = ['blOCh__khr: ',me.message];end
                            catch  me;
                                if Dimkhr == 1
                                    Msg = ['blOCh__khr: ',me.message,' R: spatial resolution (#)'];
                                elseif Dimkhr == 2
                                    Msg = ['blOCh__khr: ',me.message,' R: spatial resolution, [(#),(#)]'];
                                elseif Dimkhr == 3
                                    Msg = ['blOCh__khr: ',me.message,' R: spatial resolution, [(#),(#),(#)]'];
                                end
                            end
                        catch  me;
                            if Dimkhr == 1
                                Msg = ['blOCh__khr: ',me.message,' L: FOX length (m)'];
                            elseif Dimkhr == 2
                                Msg = ['blOCh__khr: ',me.message,' L: FOX length, [(m),(m)]'];
                            elseif Dimkhr == 3
                                Msg = ['blOCh__khr: ',me.message,' L: FOX length, [(m),(m),(m)]'];
                            end
                            
                        end
                    catch  me;Msg = ['blOCh__khr: ',me.message];end
                catch  me;Msg = ['blOCh__khr: ',me.message];end
            catch  me;Msg = ['blOCh__khr: ',me.message];end
        catch  MException;Display_Message(['blOCh__khr: Fun: ',MException.message],2);end
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
            khr.Status = 0;
        else
            khr.Status = 1;
        end
    end
    
    if khr.Status
        khr.Status = 0;
        try p.parse(Fun,spc,Traj,Par,varargin{:});
        catch me;Msg = ['blOCh__khr: ',me.message];end
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
        else
            khr.Status = 1;
        end
    end
    
    if khr.Status
        khr.Status = 0;
        khr.Traj = p.Results.Traj;
        tempPar = p.Results.Par;
        
        khr.dum.valid_Systems = valid_Systems;
        khr.dum.actual_GmHW = actual_GmHW;
        khr.dum.actual_SmHW = actual_SmHW;
        khr.dum.actual_RFmHW = actual_RFmHW;
        khr.MaskOnly = p.Results.MaskOnly;
        khr.Show = p.Results.Show;
        [pathstr, name, extension] = fileparts(khr.Traj);
        
        if any(strcmpi(khr.Traj,valid_Trajs))
            
            GoTo = 1;
            
        else
            switch extension
                
                case {'gra','.gra'}
                    GoTo = 3;
                case {'txt','.txt'}
                    GoTo = 4;
                case ''
                    fileID = fopen(khr.Traj);
                    if fileID == -1
                        fileID = fopen([khr.Traj,'.m']);
                    end
                    if fileID ~= -1
                        GoTo = 2;
                        
                        
                    end
                    
                otherwise
                    GoTo = 0;
            end
            
        end
        khr.Par = p.Results.Par;
        switch GoTo 
            
            case {1,2}
           idx = find(strcmpi(p.Results.Sys,valid_Systems));
            
            khr.L = p.Results.L;
            khr.R = p.Results.R;
            khr.Rv = p.Results.Rv;
            
            khr.GmHW = actual_GmHW(idx);
            khr.SmHW = actual_SmHW(idx);
            khr.RFmHW = actual_RFmHW(idx);
            khr.Gmpct = p.Results.Gmpct;
            khr.Smpct = p.Results.Smpct;
            khr.Nuc = p.Results.Nuc;
            khr.Sys = p.Results.Sys;
            
            khr.dt = p.Results.dt;
            
            % SYS: if you need another nuclei input the gyromagnetic ratio here
            % accordingly.
            
            switch khr.Nuc
                case '1H'
                    khr.gamma = 26.7522128e7;
                case '13C'
                    khr.gamma = 6.7282840e7;
                case '19F'
                    khr.gamma = 25.1814800e7;
            end
            khr.gammabar = khr.gamma/2/pi;
                    % IMP: Go to this script and change accordingly
                    
                    [temp,Msg] = Make_khr(spc,khr);

            
            case 3
                [temp,Msg] = Load_mat_khr(khr);
            case 4
            khr.Nuc = p.Results.Nuc;
                    khr.Sys = p.Results.Sys;
                    
                    khr.dt = p.Results.dt;
                    
                    [temp,Msg] = Load_ascii_khr(khr);
        end

           if 0 
            
                % IMP: Make a case for your trajectory. In that case make a
                % default par-struct. You assign special parameters to the
                % input struct Par, and all variables not given in this struct will be
                % read from the default par-struct.
                
                switch khr.Traj
                    
                    case 'Random'
                        Par.Gipct     = 10;
                    case 'spokes'
                        
                    case '1+1DTrapOsc'
                        Par.TBW = 8;
                        Par.T = 3.2e-3;
                    case {'1Dsl','1Dslref'}
                        Par.TBW = 8;
                        Par.T = 3.2e-3;
                    case '2DrandT'
                        Par.TBW = 8;
                        Par.T = 3.2e-3;
                        
                        
                        
                        
                        
                    case '3Dsphstckvdspii'
                        Par.Lvd = [max(spc.L(1:2)),max(spc.L(1:2))];
                        Par.rvd = [0.33 0.67];
                    case '3Dstckspi_rep'
                        Par.Rv = ceil(spc.Rv./2);

             
                    otherwise
                        Par = [];
                end
                
                
                % IMP: Here you can add a special validation function in a
                % case. Otherwise parameters to go on will simply those of
                % khr.Par. And if any missing the default will be loaded. From
                % this point on, any errors associated with ill parameters,
                % say, wrong array types etc. is in your hands.

            
            
        end
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
        else
            khr.Status = 1;
        end
    end
    
    if khr.Status
        khr.Status = 0;
        try
            
            khr.dt = temp.dt;
            
            
            
            
            if isfield(temp,'Tsub')
                khr.Tsub = temp.Tsub;
                khr.BWnonpadded = temp.BWnonpadded;
                khr.BWpadded = temp.BWpadded;
                khr.Npad = temp.Npad;
                khr.Npadded = temp.Npadded;
            end
            khr.Gmpct = temp.Gmpct;
            if isfield(temp,'Gipct')
                
                khr.Gipct = temp.Gipct;
            else
                khr.Gipct = 10;
            end
            khr.gamma = temp.gamma;
            khr.Smpct = temp.Smpct;
            khr.GmHW = temp.GmHW;
            khr.SmHW = temp.SmHW;
            khr.RFmHW = temp.RFmHW;
            khr.Sys = temp.Sys;
            
            
            khr.Dim = spc.Dim;
            if khr.MaskOnly
                
                khr.k = temp.k(:,temp.m_idx_from:temp.m_idx_to);
                khr.g = temp.g(:,temp.m_idx_from:temp.m_idx_to);
                khr.s = temp.s(:,temp.m_idx_from:temp.m_idx_to);
                khr.t = temp.t(temp.m_idx_from:temp.m_idx_to);
                khr.m = temp.m(temp.m_idx_from:temp.m_idx_to);
                khr.N = length(khr.t);
                khr.m_idx_from =1;
                khr.m_idx_to = khr.N;
                khr.T = khr.dt.*khr.N;
            else
                khr.k = temp.k;
                khr.g = temp.g;
                khr.s = temp.s;
                khr.t = temp.t;
                khr.m = temp.m;
                khr.m_idx_from =temp.m_idx_from;
                khr.m_idx_to = temp.m_idx_to;
                khr.T = temp.T;
                khr.N = temp.N;
            end
            
        catch  me;Msg = ['blOCh__khr: ',me.message];end
        
        
        
        
        if ~isempty(Msg)
            Display_Message(Msg,spc.Print);
        else
            khr.Status = 1;
        end
        
        if khr.Status
            
            if khr.Show > 0
                
                khr.fig = Show_khr(khr);
            end
            
        end
        
    end
    varargout{1} = khr;
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

function test = Validate_Traj(x,ValidTrajs)
if ischar(x)
    if ismember(x,ValidTrajs)
        test = true;
        Display_Message(sprintf('Validate_Traj: Traj "%s" is OK',x),1);
    else
        A = exist(x);
        %         B = str2func(x)
        if A == 1 || A == 7 || A == 3 || A == 6 || A == 4 || A == 8 || A == 5
            test = false;
            Display_Message(sprintf('Validate_Traj: Something is wrong with Traj "%s"',x),2);
        elseif A == 2
            
            [pathstr, name, extension] = fileparts(x);
%             khr.Show = p.Results.Show;
            switch extension
                
                case {'gra','.gra','txt','.txt'}

                    Gofurther = 1;
                    fileID = fopen(x);
                case ''
                    fileID = fopen(x);
                    if fileID == -1
                        fileID = fopen([x,'.m']);
                    end
                    if fileID ~= -1
                        Gofurther = 1;
                    end
                    
                otherwise
                    Gofurther = 0;
            end
            
            
            if Gofurther
                line = fgetl(fileID);
                fclose(fileID);
                
                if strcmp(line(1:8),'function')
                    
                    test = true;
                    Display_Message(sprintf('Validate_Traj: File for Traj "%s" exists and is presumably a function for generating a trajectory',x),1);
                elseif strcmp(line(1:6),'MATLAB')
                    Display_Message(sprintf('Validate_Traj: File for Traj "%s" exists and is presumably a mat-file for generating a trajectory',x),1);
                    test = true;
                else
                    idx = strfind(line,' ');
                    if ~isempty(idx)
                         temp = line(1:idx(1)-1);
                         
                         temp2 = str2double(temp);
                         
                         if isnumeric(temp2)
                             Display_Message(sprintf('Validate_Traj: File for Traj "%s" exists and is presumably an ascii-file with a trajectory',x),1);
                             test = true;
                         else
                             Display_Message(sprintf('Validate_Traj: Something is wrong with Traj "%s"',x),2);
            
                         end
                         
                    else
                        Display_Message(sprintf('Validate_Traj: Something is wrong with Traj "%s"',x),2);
            
                    end
                    
                    
                end
            else
              Display_Message(sprintf('Validate_Traj: Something is wrong with Traj "%s"',x),2);
              
            end
        else
            Display_Message(sprintf('Validate_Traj: Something is wrong with Traj "%s"',x),2);
                
        end

    end
else
    Display_Message(sprintf('Validate_Traj: The Traj must be a literal string'),2);
    test = false;
end

end

function test = Validate_Par(x)
if isstruct(x)
    test = true;
else
    Display_Message(sprintf('Validate_Par: The Par input must be a struct'),2);
    test = false;
end

end

function [khr,Msg] = Make_khr(spc,khr)

Msg = '';
switch khr.Traj
    
    
    
    case {'1+1DTrapOsc'}
        
        switch spc.Dim
            
            case '1+1DRL'
                
                [k,g,s,t,m,khr.TH] = Traj_Trapesoidal_Oscillation(khr.Par.TBW,khr.gammabar,khr.Gmpct,khr.GmHW,khr.Smpct,khr.SmHW,spc.BW,khr.dt,khr.Rv,spc.TH(1));
                
                
                khr.N = length(k);
                khr.k = zeros(3,khr.N);
                khr.g = zeros(3,khr.N);
                khr.s = zeros(3,khr.N);
                khr.t = t;
                khr.m = m;
                
                khr.k(1,:) = k;
                khr.g(1,:) = g;
                khr.s(1,:) = s;
            case '1+1DAP'
                
                [k,g,s,t,m,khr.TH] = Traj_Trapesoidal_Oscillation(khr.Par.TBW,khr.gammabar,khr.Gmpct,khr.GmHW,khr.Smpct,khr.SmHW,spc.BW,khr.dt,khr.Rv,spc.TH(2));
                
                
                khr.N = length(k);
                khr.k = zeros(3,khr.N);
                khr.g = zeros(3,khr.N);
                khr.s = zeros(3,khr.N);
                khr.t = t;
                khr.m = m;
                
                khr.k(2,:) = k;
                khr.g(2,:) = g;
                khr.s(2,:) = s;
            case '1+1DSI'
                
                [k,g,s,t,m,khr.TH] = Traj_Trapesoidal_Oscillation(khr.Par.TBW,khr.gammabar,khr.Gmpct,khr.GmHW,khr.Smpct,khr.SmHW,spc.BW,khr.dt,khr.Rv,spc.TH(3));
                
                
                khr.N = length(k);
                khr.k = zeros(3,khr.N);
                khr.g = zeros(3,khr.N);
                khr.s = zeros(3,khr.N);
                khr.t = t;
                khr.m = m;
                
                khr.k(3,:) = k;
                khr.g(3,:) = g;
                khr.s(3,:) = s;
        end
        
        
        
        
    case '2Dspii'
        khr.R = max(khr.R(1:2));
        khr.L = max(khr.L(1:2));
        
        Par.T = 3.2e-3;
        Par.Acc = [1,1,1,1];
                        
        khr.Par = Get_NewPar(khr.Par,Par);
        
        
        [khr.k,khr.g,khr.s,khr.t,khr.m,khr.N,khr.Ntot] = Traj_2Dspii(khr.GmHW,khr.SmHW,khr.Gmpct,khr.Smpct,khr.gamma,khr.dt,khr.L./max(khr.Par.Acc(1:2)),khr.R./max(khr.Par.Acc(1:2)),'In');
        
        
        khr.N = length(khr.g);
        
        khr.Dim = khr.Traj(1:2);
        
        temp = khr;
        switch spc.Dim
            case '2DSa'
                khr.k(2,:) = temp.k(1,:);
                khr.k(3,:) = temp.k(2,:);
                khr.k(1,:) = temp.k(3,:);
                
                khr.g(2,:) = temp.g(1,:);
                khr.g(3,:) = temp.g(2,:);
                khr.g(1,:) = temp.g(3,:);
                
                khr.s(2,:) = temp.s(1,:);
                khr.s(3,:) = temp.s(2,:);
                khr.s(1,:) = temp.s(3,:);
                
                
            otherwise
                
        end
        
    case '2Dspio'
        khr.R = max(khr.R(1:2));
        khr.L = max(khr.L(1:2));
        
        Par.T = 3.2e-3;
        Par.Acc = [1,1,1,1];
                        
        khr.Par = Get_NewPar(khr.Par,Par);
        
        [khr.k,khr.g,khr.s,khr.t,khr.m,khr.N,khr.Ntot] = Traj_2Dspii(khr.GmHW,khr.SmHW,khr.Gmpct,khr.Smpct,khr.gamma,khr.dt,khr.L./max(khr.Acc(1:2)),khr.R./max(khr.Acc(1:2)),'In');
        
        
        khr.N = length(khr.g);
       
        
        khr.Dim = khr.Traj(1:2);
        
        temp = khr;
        switch spc.Dim
            case '2DSa'
                khr.k(2,:) = temp.k(1,:);
                khr.k(3,:) = temp.k(2,:);
                khr.k(1,:) = temp.k(3,:);
                
                khr.g(2,:) = temp.g(1,:);
                khr.g(3,:) = temp.g(2,:);
                khr.g(1,:) = temp.g(3,:);
                
                khr.s(2,:) = temp.s(1,:);
                khr.s(3,:) = temp.s(2,:);
                khr.s(1,:) = temp.s(3,:);
                
                
            otherwise
                
        end
        khr.k = khr.k(:,end:-1:1);
        khr.g = khr.g(:,end:-1:1);
        khr.s = khr.s(:,end:-1:1);
        khr.m = khr.m(:,end:-1:1);
        
        
        
    case 'Random'
        
        Par.T = 3.2e-3;
        Par.Gipct = 50;     
        khr.Par = Get_NewPar(khr.Par,Par);
        
        khr.N = ceil(khr.Par.T/khr.dt);
        
        switch spc.Dim
            
            case '1DRL'
                khr.g = Traj_Random(khr.N,[1],khr.GmHW*Par.Gipct/100);
                
            case '1DAP'
                khr.g = Traj_Random(khr.N,[2],khr.GmHW*Par.Gipct/100);
            case '1DSI'
                khr.g = Traj_Random(khr.N,[3],khr.GmHW*Par.Gipct/100);
            case '1+1DRL'
                khr.g = Traj_Random(khr.N,[1],khr.GmHW*Par.Gipct/100);
            case '1+1DAP'
                khr.g = Traj_Random(khr.N,[2],khr.GmHW*Par.Gipct/100);
            case '1+1DSI'
                khr.g = Traj_Random(khr.N,[3],khr.GmHW*Par.Gipct/100);
            case '2DAx'
                khr.g = Traj_Random(khr.N,[1,2],khr.GmHW*Par.Gipct/100);
            case '2DCo'
                khr.g = Traj_Random(khr.N,[1,3],khr.GmHW*Par.Gipct/100);
            case '2DSa'
                khr.g = Traj_Random(khr.N,[2,3],khr.GmHW*Par.Gipct/100);
            case '2+1DAx'
                khr.g = Traj_Random(khr.N,[1,2],khr.GmHW*Par.Gipct/100);
            case '2+1DCo'
                khr.g = Traj_Random(khr.N,[1,3],khr.GmHW*Par.Gipct/100);
            case '2+1DSa'
                khr.g = Traj_Random(khr.N,[2,3],khr.GmHW*Par.Gipct/100);
            case '3D'
                khr.g = Traj_Random(khr.N,[1,2,3],khr.GmHW*Par.Gipct/100);
            case '3+1D'
                khr.g = Traj_Random(khr.N,[1,2,3],khr.GmHW*Par.Gipct/100);
        end
        
        khr.k = cumsum(khr.g)*khr.gamma*khr.dt;
        khr.s = diff([[0;0;0],khr.g])./khr.dt;
        khr.m = true(1,khr.N);
        khr.t = [0:khr.N-1]*khr.dt;
        
        
    case '1Dsl'
        Par.TBW = 4;
        Par.T = 3.2e-3;
        khr.Par = Get_NewPar(khr.Par,Par);
        switch spc.Dim
            
            case '1DRL'
                
                [k,g,s,khr.t,khr.m,khr.N] = Traj_1Dsl(khr.Par.TBW,khr.Par.T,khr.gamma,spc.TH(1),khr.dt);
                khr.k = zeros(3,khr.N);
                khr.k(1,:) = k;
                khr.g = zeros(3,khr.N);
                khr.g(1,:) = g;
                khr.s = zeros(3,khr.N);
                khr.s(1,:) = s;
            case '1DAP'
                [k,g,s,khr.t,khr.m,khr.N] = Traj_1Dsl(khr.Par.TBW,khr.Par.T,khr.gamma,spc.TH(2),khr.dt);
                khr.k = zeros(3,khr.N);
                khr.k(2,:) = k;
                khr.g = zeros(3,khr.N);
                khr.g(2,:) = g;
                khr.s = zeros(3,khr.N);
                khr.s(2,:) = s;
            case '1DSI'
                [k,g,s,khr.t,khr.m,khr.N] = Traj_1Dsl(khr.Par.TBW,khr.Par.T,khr.gamma,spc.TH(3),khr.dt);
                khr.k = zeros(3,khr.N);
                khr.k(3,:) = k;
                khr.g = zeros(3,khr.N);
                khr.g(3,:) = g;
                khr.s = zeros(3,khr.N);
                khr.s(3,:) = s;
        end
        
        
         khr.Dim = khr.Traj(1:2);
        
    case '1Dslref'
        Par.TBW = 4;
        Par.T = 3.2e-3;
        khr.Par = Get_NewPar(khr.Par,Par);
        switch spc.Dim
            
            case '1DRL'
                
                [k,g,s,khr.t,khr.m,khr.N] = Traj_1Dslref(khr.Par.TBW,khr.Par.T,khr.gamma,spc.TH(1),khr.dt);
                khr.k = zeros(3,khr.N);
                khr.k(1,:) = k;
                khr.g = zeros(3,khr.N);
                khr.g(1,:) = g;
                khr.s = zeros(3,khr.N);
                khr.s(1,:) = s;
            case '1DAP'
                [k,g,s,khr.t,khr.m,khr.N] = Traj_1Dslref(khr.Par.TBW,khr.Par.T,khr.gamma,spc.TH(2),khr.dt);
                khr.k = zeros(3,khr.N);
                khr.k(2,:) = k;
                khr.g = zeros(3,khr.N);
                khr.g(2,:) = g;
                khr.s = zeros(3,khr.N);
                khr.s(2,:) = s;
            case '1DSI'
                [k,g,s,khr.t,khr.m,khr.N] = Traj_1Dslref(khr.Par.TBW,khr.Par.T,khr.gamma,spc.TH(3),khr.dt);
                khr.k = zeros(3,khr.N);
                khr.k(3,:) = k;
                khr.g = zeros(3,khr.N);
                khr.g(3,:) = g;
                khr.s = zeros(3,khr.N);
                khr.s(3,:) = s;
        end
        
        
        khr.Dim = khr.Traj(1:2);
    case '1Dplat'
        khr.Dim = khr.Traj(1:2);
        [khr.k,khr.g,khr.s,khr.t,khr.ktot,khr.gtot,khr.stot,khr.ttot,khr.N,khr.Ntot] = Traj_1Dplat(khr.GmHW,khr.SmHW,khr.Gmpct,khr.Smpct,khr.gamma,khr.dt);
        
        
    case 'none4spect'
        Par.BW = spc.BW;
        
        Par.dt = 10e-6;
        
        Par.T = 1/spc.pBW*2;
        Par.Rv = Par.T/Par.dt;
        
        khr.Par = Get_NewPar(khr.Par,Par);
        
        [khr.k,khr.g,khr.s,khr.t,khr.m,khr.N] = Traj_none4spect(khr.Par.BW,khr.Par.T,khr.Par.dt);
    case 'none4spect2'
        % default
        Par.NOC = 1000;
        Par.Tactive = 1e-3;
        Par.Tdead = 200e-6;
        
        % compare default with user input
        khr.Par = Get_NewPar(khr.Par,Par);
        
        % correct values to agree with dt and number of controls.
        khr.Par.dt = khr.Par.Tactive/khr.Par.NOC;
        khr.Par.Tactive= khr.Par.NOC*khr.Par.dt;
        khr.Par.Tdead = ceil(khr.Par.Tdead/khr.Par.dt)*khr.Par.dt;
        
        [khr.k,khr.g,khr.s,khr.t,khr.m,khr.N] = Traj_none4spect2(khr.Par.Tdead,khr.Par.Tactive,khr.Par.dt);        
    otherwise
        
        fun = str2func(khr.Traj);
        khr = fun(spc,khr);
        
        
        
        
end



khr.T = khr.dt.*khr.N;

khr.m_idx_from = find(diff([0,khr.m])>0);
khr.m_idx_to = find(diff([khr.m,0])<0);

khr.Dim = spc.Dim;
if ~isempty(Msg)
    Display_Message(Msg,spc.Print);
    Status = 0;
end






end

function [ting,Msg] = Load_mat_khr(khr)

Msg = [];
try Load = load(khr.Traj,'-mat');
    try validateattributes(Load.g,{'double'},{'size',[3,NaN],'real','finite'});
        ting.g = Load.g;
        ting.N = size(ting.g,2);
        try validateattributes(Load.dt,{'double'},{'size',[1,1],'real','finite','>',0});
            ting.dt = Load.dt;
            try validateattributes(Load.Nuc,{'char'},{'row'});
                ting.Nuc = Load.Nuc;
                switch ting.Nuc
                    case '1H'
                        ting.gamma = 2.6751289e8;
                    case '13C'
                        ting.gamma = 67.262e6;
                    case '19F'
                        ting.gamma = 251.662e6;
                    otherwise
                        Msg = ['Load_mat_khr: Nuc: must be 1H, 13C or 19F'];
                end
                % @(x)any(strcmpi(x,valid_Dimensions))
                try validateattributes(Load.Sys,{'char'},{'row'});
                    ting.Sys = Load.Sys;
                    % I need to implement this line:
                    %                @(x)any(strcmpi(x,khr.dum.valid_Systems))
                    
                    
                    
                    idx = find(strcmpi(ting.Sys,khr.dum.valid_Systems));
                    
                    ting.GmHW = khr.dum.actual_GmHW(idx);
                    ting.SmHW = khr.dum.actual_SmHW(idx);
                    ting.RFmHW = khr.dum.actual_RFmHW(idx);
                    try validateattributes(Load.m,{'logical'},{'size',[1,ting.N]});
                        ting.m = Load.m;
                        ting.m_idx_from = find(diff([0,ting.m])>0);
                        ting.m_idx_to = find(diff([ting.m,0])<0);
                        existm = 1;
                    catch me;Nfo = ['Load_mat_khr: ',me.message,' m: making default'];
                        Display_Message(Nfo,khr.Print);
                        existm = 0;
                        
                    end
                    
                    try validateattributes(Load.k,{'double'},{'size',[3,ting.N],'real','finite'});
                        ting.k = Load.k;
                        existk = 1;
                    catch me;Nfo = ['Load_mat_khr: ',me.message,' k: making default'];
                        Display_Message(Nfo,khr.Print);
                        existk = 0;
                        
                    end
                    try validateattributes(Load.s,{'double'},{'size',[3,ting.N+1],'real','finite'});
                        ting.s = Load.s;
                        exists = 1;
                    catch me;Nfo = ['Load_mat_khr: ',me.message,' s: making default'];
                        Display_Message(Nfo,khr.Print);
                        exists = 0;
                    end
                    try validateattributes(Load.t,{'double'},{'size',[3,ting.N],'real','finite'});
                        ting.t = Load.t;
                        existt = 1;
                    catch me;Nfo = ['Load_mat_khr: ',me.message,' t: making default'];
                        Display_Message(Nfo,khr.Print);
                        existt = 0;
                    end
                    if ~existm
                        ting.m = true(1,ting.N);
                        ting.m_idx_from = 1;
                        ting.m_idx_to = ting.N;
                        
                    end
                    if ~existk
                        ting.k = ting.gamma.*trapz(ting.g,2).*ting.dt;
                    end
                    if ~exists
                        ting.s = diff([zeros(3,1),ting.g,zeros(3,1)],2,2)./ting.dt;
                    end
                    if ~existt
                        ting.t = [0:ting.N-1].*ting.dt;
                        %                   ting.t = -ting.t(end:-1:1);
                    end
                    ting.T = ting.N.*ting.dt;
                    Gabs = sqrt(ting.g(1,:).^2+ting.g(2,:).^2+ting.g(3,:).^2);
                    Gmax = Gabs(:);
                    ting.Gmpct = Gmax./ting.GmHW*100;
                    if ting.Gmpct > 100
                        Msg = ['Load_mat_khr: Gmpct: is violated'];
                    end
                    Sabs = sqrt(ting.s(1,:).^2+ting.s(2,:).^2+ting.s(3,:).^2);
                    Smax = Sabs(:);
                    ting.Smpct = Smax./ting.SmHW*100;
                    if ting.Smpct > 100
                        Msg = ['Load_mat_khr: Smpct: is violated'];
                    end
                    
                catch me;Msg = ['Load_mat_khr: ',me.message,' Sys: was not included or proper'];end
            catch me;Msg = ['Load_mat_khr: ',me.message,' Nuc: must be included'];end
        catch me;Msg = ['Load_mat_khr: ',me.message,' The File must include a dt variable.'];end
    catch me;Msg = ['Load_mat_khr: ',me.message,' g: was not included or proper'];end
catch me;Msg = ['Load_mat_khr: ',me.message,' The khr mat is faulty'];end

end

function [ting,Msg] = Load_ascii_khr(khr)
% function [ting,Msg] = Load_ascii_khr(khr)
%
%   This script loads khr data from an ascii file

Msg = [];
try g = load(khr.Traj,'-ascii');
    try validateattributes(g,{'double'},{'size',[3,NaN],'real','finite'});
        ting.Sys = khr.Sys;
        ting.Nuc = khr.Nuc;
        switch ting.Nuc
            case '1H'
                ting.gamma = 2.6751289e8;
            case '13C'
                ting.gamma = 67.262e6;
            case '19F'
                ting.gamma = 251.662e6;
            otherwise
                Msg = ['Load_ascii_khr:: Nuc: must be 1H, 13C or 19F'];
        end
        
        idx = find(strcmpi(ting.Sys,khr.dum.valid_Systems));
        
        ting.GmHW = khr.dum.actual_GmHW(idx);
        ting.SmHW = khr.dum.actual_SmHW(idx);
        ting.RFmHW = khr.dum.actual_RFmHW(idx);
        
        ting.g = g;
        ting.N = size(ting.g,2);
        ting.dt = khr.dt;
        
        
        
        
        
        ting.m = true(1,ting.N);
        ting.m_idx_from = 1;
        ting.m_idx_to = ting.N;
        
        ting.t = [0:ting.N-1].*ting.dt;
        
        ting.k = ting.gamma.*cumtrapz(ting.t,ting.g,2).*ting.dt-sum(g,2).*ting.dt*ting.gamma;

        ting.s = diff([zeros(3,1),ting.g],1,2)./ting.dt;

        ting.T = ting.N*ting.dt;
        
        Gabs = sqrt(ting.g(1,:).^2+ting.g(2,:).^2+ting.g(3,:).^2);
        Gmax = Gabs(:);
        ting.Gmpct = Gmax./ting.GmHW*100;
        if ting.Gmpct > 100
            Msg = ['Load_ascii_khr:: Gmpct: is violated'];
        end
        Sabs = sqrt(ting.s(1,:).^2+ting.s(2,:).^2+ting.s(3,:).^2);
        Smax = Sabs(:);
        ting.Smpct = Smax./ting.SmHW*100;
        if ting.Smpct > 100
            Msg = ['Load_ascii_khr:: Smpct: is violated'];
        end
        
    catch me;Msg = ['Load_ascii_khr:: ',me.message,' g: was not proper'];end
catch me;Msg = ['Load_ascii_khr:: ',me.message,' The ascii khr is faulty'];end

end

function [k,g,s,t,m,TH] = Traj_Trapesoidal_Oscillation(TBW,gammabar,Gmpct,GmHW,Smpct,SmHW,BW,dt,Rv,TH)



Gm = Gmpct*GmHW/100;
Sm = Smpct*SmHW/100;
T_osc = 1/(2*100);
[N_osc,N_pla,N_plahalf,N_sub,N_rmp] = GetN('Tosc',T_osc,dt,Gm,Sm);

[T_osc,T_pla,T_plahalf,T_sub,T_rmp] = GetT(dt,N_osc,N_pla,N_plahalf,N_sub,N_rmp);



Gamp = TBW/(gammabar*TH*T_pla);

Gapct = Gamp/GmHW*100;

if abs(Gapct) > 100
    Gapct = 100;
    Gamp = GmHW;
    
    TH = (TBW)/(gammabar*Gamp*T_pla);
    Nfo = sprintf('Traj_Trapesoidal_Oscillation: TH was adjusted to %1.2e m\n',TH);
    Display_Message(Nfo,spc.Print);
    
else
end
[g,N_burp,scale] = Getg(Gamp,N_rmp,N_pla,N_plahalf,Rv);

N_tot = length(g);


T_tot = N_tot*dt;

k = cumsum(g)*dt*gammabar-sum(g,2).*dt*gamma;
t = [0:dt:T_tot-dt];
m = true(1,N_tot);



s = diff([0,g])./dt;




end

function [k,g,s,t,ktot,gtot,stot,ttot,N,Ntot] = Traj_1Dplat(GmHW,SmHW,Gmpct,Smpct,gamma,dt)
a = SmHW*Smpct; % bgu
N = 500;
Ntot = N;

s = zeros(3,N);


t = [-N*dt+dt:dt:0];

g = zeros(3,N);
g(3,:) = ones(1,N)*GmHW*Gmpct/100;



k = cumsum(g,2).*dt*gamma-sum(g,2).*dt*gamma;

ktot =k;
gtot = g;
stot = s;
ttot = t;

end

function [k,g,s,t,m,N] = Traj_1Dsl(TBW,T,gamma,TH,dt)

BW = TBW/T;

Gsl = 2*pi/gamma*BW/TH;


N = round(T/dt);

g = ones(1,N)*Gsl;

k = cumsum(g,2).*dt*gamma-sum(g,2).*dt*gamma;
s = zeros(1,N);
m = true(1,N);

t = [-N*dt+dt:dt:0];




end

function [k,g,s,t,m,N] = Traj_1Dslref(TBW,T,gamma,TH,dt)
% with refocusing loop. No slewrate considered.
BW = TBW/T;

Gsl = 2*pi/gamma*BW/TH;


N = round(T/dt);



g = [0,ones(1,N)*Gsl,-ones(1,round(N/2))*Gsl,0];

k = cumsum(g,2).*dt*gamma-sum(g,2).*dt*gamma;
s = [0,diff(g)];
m = false(size(g));
m(g > 0) = true;

t= [0:dt:length(g)*dt-dt];
t = t-round(N/2)*dt;

N = length(g);

end

function [k,g,s,t,m,N] = Traj_none4spect(BW,T,dt)

dt2 = 1/BW;

t = 0:dt:T-dt;
N=length(t);
k = zeros(3,N);
g = zeros(3,N);
s = zeros(3,N);
m = true(1,N);
end
function [k,g,s,t,m,N] = Traj_none4spect2(Tdead,Tactive,dt)
% includes a deadtime after the pulse


t = 0:dt:Tactive+Tdead-dt;
Ndead = ceil(Tdead/dt);
N=length(t);
k = zeros(3,N);
g = zeros(3,N);
s = zeros(3,N);
m = true(1,N);
m(end-Ndead+1:end) = false;


end

function g = Traj_Random(N,Dim,Gi)

g = zeros(3,N);
for n = 1:numel(Dim)
    rand('seed',sum(100*clock()));
    g_ = (rand(1,N)-0.5)*2.*Gi;
    
    if Dim(n) == 1
        g(1,:) = g_;
    elseif Dim(n) == 2
        g(2,:) = g_;
    elseif Dim(n) == 3
        g(3,:) = g_;
    end
    pause(1)
end




end

function [ktot,gtot,stot,ttot,m,N,Ntot,theta] = Traj_2Dspii(GmHW,SmHW,Gmpct,Smpct,gamma,dt,L,R,Direction)
% function [ktot,gtot,stot,ttot,m,N,Ntot,theta] = Make_2Dspii(GmHW,SmHW,Gmpct,Smpct,gamma,dt,L,R,Direction)
%
%   This script makes a 2d inward spiral based on the Glover'1999 approach


Smax = SmHW*Smpct/100;
Gmax = GmHW*Gmpct/100;

gammabar = gamma/2/pi;
dx = L/R;
lambda = 1/(2*pi*L);
beta = gammabar*Smax/lambda;
kmax = R/(2*L);
dr = 1./(2*kmax);
a2 = (9*beta/4)^(1/3);
Lambda = 10;
thetamax = kmax/lambda;
Ts = (3*gamma*Gmax/(4*pi*lambda*a2^2))^3;
thetas = (1/2*beta.*Ts.^2)./(Lambda+beta/(2*a2)*Ts.^(4/3));
Tg = pi*lambda*(thetamax^2-thetas^2)/(gamma*Gmax);
Ns = round(Ts/dt);
Ng = round(Tg/dt);

if thetamax > thetas
    Msg = sprintf(1,'Make_2Dspii: It contains a slewrate limited portion and a amplitude limited portion');
    Display_Message(Msg,spc.Print);
    
    Nt = Ns+Ng;
    Tacq = Ts + Tg;
    
    ts = linspace(0,Ts,Ns);
    tg = linspace(Ts+dt,Tacq,Ng);
    ttot = [ts,tg];
    
    theta1 = (beta/2.*ts.^2)./(Lambda+beta/(2*a2)*ts.^(4/3));
    theta1dot = -2/3*a2*beta.*ts.*(beta.*ts.^(4/3)-6*a2*Lambda)./(beta.*ts.^(4/3)+6*a2*Lambda).^2;
    
    
    theta2 = sqrt(thetas.^2+gamma/(pi*lambda)*Gmax*(tg-Ts));
    
    theta2dot = Gmax*gamma./(2*lambda*sqrt(pi/lambda*(Gmax*gamma.*tg-Gmax*gamma*Ts+lambda*pi*thetas^2)));
    
    k1 = lambda.*theta1.*complex(cos(theta1),sin(theta1));
    k2 = lambda.*theta2.*complex(cos(theta2),sin(theta2));
    k_ = [k1,k2];
    theta = [theta1,theta2];
else
    
    
    
    Tacq = 2*pi*L/(3*1)*sqrt(pi/(gamma*Smax*dx^3));
    Nt = round(Tacq/dt);
    ts = linspace(0,Tacq,Nt);
    ttot = [ts];
    theta1 = (beta/2.*ts.^2)./(Lambda+beta/(2*a2)*ts.^(4/3));
    %     theta1dot = -2/3*a2*beta.*ts.*(beta.*ts.^(4/3)-6*a2*Lambda)./(beta.*ts.^(4/3)+6*a2*Lambda).^2;
    theta = theta1;
    
    k_ = lambda.*theta1.*complex(cos(theta1),sin(theta1));
    
end


pha = -pi/4;

alpha = Smax/1.5;

k_ = k_.*exp(-1i*(angle(k_(end))+pha));
g_ = [0, diff(k_,1,2)./(dt*gammabar)];

ge = real(g_(end));

trmp = abs(ge)/alpha;
tplateau = trmp/10;

nrmp = round(trmp/dt);
nplateau = round(tplateau/dt);

grmp = [ones(1,nplateau)*abs(ge),linspace(abs(ge),0,nrmp)];
grmp = complex(-grmp,grmp);

grmp = grmp(end:-1:1);

g_ = g_(end:-1:1);

nramp = length(grmp);

k_ = gammabar*cumsum([grmp,g_])*dt;

ke = k_(end);

kex = real(ke);
key = imag(ke);
sx = sign(kex);
sy = sign(key);

[kmx,idxmx] = max([abs(kex),abs(key)]);


t_tip = sqrt(kmx/gammabar/alpha)*2;


ntipup = ceil(t_tip/2/dt);
ntipdown = floor(t_tip/2/dt);

if ntipup == ntipdown
    ntipdown = ntipdown-1;
end

ntip = ntipdown+ntipup;

tip = [linspace(0,1,ntipup),linspace(1,0,ntipdown)];
tip = tip./sum(tip);

g_tip = complex(-sx.*tip.*kex,-sy.*tip.*key)/gammabar/dt;

g_f = [g_tip,grmp,g_];


gprep = [g_tip,grmp];
gtraj = g_;

gtot = [gprep,gtraj];

ntot = length(gtot);
nprep = length(gprep); % MSV efter 3 okt 2012

ntraj = length(gtraj);



ktot_ = gammabar*cumsum(gtot)*dt;

k = ktot_(nprep+1:end);








switch Direction
    
    case 'In'
        
        stot_ = [0,diff(gtot,1,2)/dt];
        k_ = ktot_(nprep+1:end);
        g_ = gtot(nprep+1:end);
        gtot_ = gtot;
        s_ = stot_(nprep+1:end);
        ttot_ = [-ntot*dt+dt:dt:0];
        t_ = ttot_(nprep+1:end);
        
        Ntot = length(ttot_);
        N = length(k_);
        
        m = false(1,ntot);
        m(nprep+1:end) = true;
        
    case 'Out'
        
        ktot_ = ktot_(end:-1:1);
        gtot_ = gtot(end:-1:1);
        stot_ = [0,diff(gtot_,1,2)/dt];
        k_ = ktot_(1:ntraj);
        g_ = gtot_(1:ntraj);
        s_ = stot_(1:ntraj);
        
        
        ttot_ = [0:dt:ntot*dt-dt];
        t_ = ttot(1:ntraj);
        
        Ntot = length(ttot_);
        N = length(k_);
        
        m = true(1,ntot);
        m(end-nprep+1:end) = false;
        
    otherwise
        
        error('Make_2Dspii: Specify Direction: In/Out')
        
end


k = [real(k_);imag(k_);zeros(1,N)];
g = [real(g_);imag(g_);zeros(1,N)];
s = [real(s_);imag(s_);zeros(1,N)];
t = t_;

ktot = [real(ktot_);imag(ktot_);zeros(1,Ntot)];
gtot = [real(gtot_);imag(gtot_);zeros(1,Ntot)];
stot = [real(stot_);imag(stot_);zeros(1,Ntot)];
ttot = ttot_;






end

function [k,g,s,t,m] = Traj_1p1Dsl_rep(TBW,T,gamma,TH,dt)

BW = TBW/T;

Gsl = 2*pi/gamma*BW/(TH);


N = round(T/dt);

g = ones(1,N)*Gsl;

k = cumsum(g).*dt*gamma;
s = zeros(1,N);


t = [-N*dt+dt:dt:0];

m = true(size(k));


end

function z_sc = Get_z_sc(Nz,vord)


z_sc1 = linspace(-1,1,Nz);

switch vord
    case 'sym'
        [temp,idx] = sort(abs(z_sc1));
        z_sc2 = z_sc1(idx);
        z_sc = z_sc2(end:-1:1);
    case 'seq'
        
        z_sc = z_sc1;
end




end


function [N_osc,N_pla,N_plahalf,N_sub,N_rmp] = GetN(Type,varargin)

switch Type
    
    case 'Tosc'
        % varargin{1} = T_tot;
        % varargin{2} = dt;
        % varargin{3} = Gmax;
        % varargin{4} = Smax;
        
        
        N_sub = round(varargin{1}./2./varargin{2});
        if mod(N_sub,2)
            N_sub = N_sub+1;
        end
        
        N_rmp = ceil(varargin{3}./varargin{4}./varargin{2});
        
        N_pla = N_sub - 2*N_rmp;
        
        if N_pla < 2
            N_pla = 2;
        end
        if mod(N_pla,2)
            N_pla = N_pla+1;
        end
        
        N_sub = 2*N_rmp + N_pla;
        N_osc = 2*N_sub;
        N_plahalf = N_pla/2;
        
    case 'Npla'
        % varargin{1} = T_pla;
        % varargin{2} = dt;
        % varargin{3}  = N_rmp;
        T_pla = varargin{1};
        dt = varargin{2};
        
        N_pla = ceil(T_pla/dt);
        
        
        N_rmp = varargin{3};
        if mod(N_pla,2)
            N_pla = N_pla+1;
            
        end
        N_sub = 2*N_rmp + N_pla;
        N_osc = 2*N_sub;
        N_plahalf = N_pla/2;
        
end
end

function [g,N_burp,scale] = Getg(Gamp,N_rmp,N_pla,N_plahalf,Ncyc)

scale = (2*N_plahalf+N_rmp)/(2*(N_plahalf+N_rmp));

gburp = [linspace(0,1-1/N_rmp,N_rmp),ones(1,N_plahalf),linspace(1-1/N_rmp,0,N_rmp)].*scale.*Gamp;
%  _
% / \
g1 = [linspace(0,-1+1/N_rmp,N_rmp),-ones(1,N_pla),linspace(-1+1/N_rmp,1-1/N_rmp,2*N_rmp),ones(1,N_plahalf)].*Gamp;
%      _
%     /
% \__/

g2 = [ones(1,N_plahalf),linspace(1-1/N_rmp,-1+1/N_rmp,2*N_rmp),-ones(1,N_pla),linspace(-1+1/N_rmp,1-1/N_rmp,2*N_rmp),ones(1,N_plahalf)].*Gamp;
% _      _
%  \    /
%   \__/
N_burp = length(gburp);



gref = [gburp,g1];
g = gref;

for i=2:(Ncyc-1)
    g = [g g2];
    
end
g = [g gref(end:-1:1)];

end

function [varargout] = GetT(dt,varargin)

for n = 1:length(varargin)
    varargout{n} = varargin{n}*dt;
end
end

function nP = Get_NewPar(uP,dP)

if isempty(uP)
    nP = dP;
else
dPnames = recursivefieldnames(dP);
uPnames = recursivefieldnames(uP);
nP = struct;
N = length(dPnames);
for n = 1:N
    curName = dPnames{n};
    test = find(strcmp(curName,uPnames));
    if isempty(test)
       eval(['nP.',curName,'=dP.',curName,';'])
    else
        eval(['nP.',curName,'=uP.',curName,';'])
    end
    
end
end

end

function O = recursivefieldnames(I)
f=0; q=1; U={};
O=cellfun(@(x)strcat('I.',x),fieldnames(I),'UniformOutput',0);

while q~=0
    H={}; q=length(O);
    for i=1:length(O) 
        if isstruct(eval(O{i}))==1
            if f~=1
                A=fieldnames(eval(O{i}));
                A=cellfun(@(x)strcat(sprintf('%s.',O{i}),x),A,'UniformOutput',0);
                H=cat(1,H,A);
            elseif f==1
                H=cat(1,H,O{i});
                q=q-1;
            end
            U=cat(1,U,O{i});
        else
            H=cat(1,H,O{i}); q=q-1;
        end
    end
    O = H;
end
O=cellfun(@(x)sprintf('%s',x(3:end)),O,'UniformOutput',0);
end












function fig = Show_khr(khr)
global h
h.khr = khr;

    fig = figure;
    
mp = get(0, 'MonitorPositions');
if size(mp,1) > 1
    mp = mp(1,:);
end

set(gcf,'Units','pixels')
set(gcf,'Position',[100,1,mp(3).*0.9,mp(4).*0.9])

h.thisfig = gcf;

h.listbox1= uicontrol('Style','listbox','Callback',@listbox1_Callback);
set(h.listbox1,'Units','normalized')
h.x_listbox1 = 0;
h.y_listbox1 = 0;
h.w_listbox1 = 0.1;
h.h_listbox1 = 0.5;
set(h.listbox1,'Position',[h.x_listbox1,h.y_listbox1,h.w_listbox1,h.h_listbox1])
bg1 = uibuttongroup('Visible','off',...
    'Position',[0 0.9 .1 0.1],...
    'SelectionChangedFcn',@P1_SelectionChangeFcn);

h.Show_k = uicontrol(bg1,'Style',...
    'radiobutton',...
    'String','k',...
    'Units','normalized',...
    'Position',[0.1 0.7 0.9,0.33],...
    'HandleVisibility','off');

h.Show_g = uicontrol(bg1,'Style','radiobutton',...
    'String','g',...
    'Units','normalized',...
    'Position',[0.1 0.4 0.9,0.33],...
    'HandleVisibility','off');
%
h.Show_s = uicontrol(bg1,'Style','radiobutton',...
    'String','s',...
    'Units','normalized',...
    'Position',[0.1 0.1 0.9,0.33],...
    'HandleVisibility','off');

bg1.Visible = 'on';

bg2 = uibuttongroup('Visible','off',...
    'Position',[0 0.5 .1 0.4],...
    'SelectionChangedFcn',@P2_SelectionChangeFcn);

h.Show_3D = uicontrol(bg2,'Style',...
    'radiobutton',...
    'String','(x,y,z)',...
    'Units','normalized',...
    'Position',[0.1 0.8 0.9,0.125],...
    'HandleVisibility','off');

h.Show_2Dxy = uicontrol(bg2,'Style','radiobutton',...
    'String','(x,y)',...
    'Units','normalized',...
    'Position',[0.1 0.7 0.9,0.125],...
    'HandleVisibility','off');

h.Show_2Dxz = uicontrol(bg2,'Style','radiobutton',...
    'String','(x,z)',...
    'Units','normalized',...
    'Position',[0.1 0.6 0.9,0.125],...
    'HandleVisibility','off');

h.Show_2Dyz = uicontrol(bg2,'Style','radiobutton',...
    'String','(y,z)',...
    'Units','normalized',...
    'Position',[0.1 0.5 0.9,0.125],...
    'HandleVisibility','off');

h.Show_Xt = uicontrol(bg2,'Style','radiobutton',...
    'String','x(t)',...
    'Units','normalized',...
    'Position',[0.1 0.4 0.9,0.125],...
    'HandleVisibility','off');
h.Show_Yt = uicontrol(bg2,'Style','radiobutton',...
    'String','y(t)',...
    'Units','normalized',...
    'Position',[0.1 0.3 0.9,0.125],...
    'HandleVisibility','off');
h.Show_Zt = uicontrol(bg2,'Style','radiobutton',...
    'String','z(t)',...
    'Units','normalized',...
    'Position',[0.1 0.2 0.9,0.125],...
    'HandleVisibility','off');

h.Show_XtYtZt = uicontrol(bg2,'Style','radiobutton',...
    'String','x(t),y(t),z(t)',...
    'Units','normalized',...
    'Position',[0.1 0.1 0.9,0.125],...
    'HandleVisibility','on');
h.Show_Magnitude = uicontrol(bg2,'Style','checkbox',...
    'String','| * |',...
    'Units','normalized',...
    'Position',[0.1 0.0 0.9,0.125],...
    'HandleVisibility','on','Callback',@Show_Magnitude_Callback);



bg2.Visible = 'on';




if ~isfield(h.khr,'t')
    h.khr.t = [0:h.khr.dt:h.khr.N*h.khr.dt-h.khr.dt];
end
if ~isfield(h.khr,'s')
    set(h.Show_s,'Visible','off')
end
if ~isfield(h.khr,'g') && isfield(h.khr,'k')
    set(h.Show_g,'Visible','off')
    h.Type1 = 'k';
    h.Type2 = 'XtYt';
end
if ~isfield(h.khr,'k') && isfield(h.khr,'g')
    set(h.Show_k,'Visible','off')
    h.Type1 = 'g';
    h.Type2 = 'XtYt';
elseif ~isfield(h.khr,'k') && ~isfield(h.khr,'g')
    error('There must be at least a field named g or k')
else
    h.Type1 = 'g';
end

switch h.Type1
    case 'k'
        set(h.Show_k,'Value',1)
        set(h.Show_g,'Value',0)
        set(h.Show_s,'Value',0)
        
    case 'g'
        set(h.Show_k,'Value',0)
        set(h.Show_g,'Value',1)
        set(h.Show_s,'Value',0)
    case 's'
        set(h.Show_k,'Value',0)
        set(h.Show_g,'Value',0)
        set(h.Show_s,'Value',1)
end



h = Populate_Listbox(h);

h.ShowMagn = 0;


h.DispAll = axes('Parent',gcf);
h.DispX = axes('Parent',gcf);
h.DispY = axes('Parent',gcf);
h.DispZ = axes('Parent',gcf);

set(h.DispAll,'Units','normalized')
set(h.DispX,'Units','normalized')
set(h.DispY,'Units','normalized')
set(h.DispZ,'Units','normalized')

set(h.DispX,'Position'  ,[0.15 0.7 0.8 0.25])
set(h.DispY,'Position'  ,[0.15 0.375 0.8 0.25])
set(h.DispZ,'Position'  ,[0.15 0.05 0.8 0.25])

set(h.DispAll,'Position',[0.15 0.05 0.8 0.9])



switch h.khr.Dim
    
    case {'1DAP','1DRL','1DSI','1+1DAP','1+1DRL','1+1DSI'}
        
        switch h.khr.Dim
            case {'1DRL','1+1DRL'}
                h.Type2 = 'Xt';
            case {'1DAP','1+1DAP'}
                h.Type2 = 'Yt';
            case {'1DSI','1+1DSI'}
                h.Type2 = 'Zt';
        end
        
        set(h.DispAll,'Visible','on')
        set(h.DispX,'Visible','off')
        set(h.DispY,'Visible','off')
        set(h.DispZ,'Visible','off')
        
        set(h.Show_3D,'Visible','off','Value',0)
        set(h.Show_2Dxy,'Visible','off','Value',0)
        set(h.Show_2Dxz,'Visible','off','Value',0)
        set(h.Show_2Dyz,'Visible','off','Value',0)
        set(h.Show_Xt,'Visible','off','Value',0)
        set(h.Show_Yt,'Visible','off','Value',0)
        set(h.Show_XtYtZt,'Visible','off','Value',0)
        set(h.Show_Zt,'Visible','on','Value',1)
        
        set(h.Show_Magnitude,'Visible','on','Value',0)
        
        
        
        
    case {'2DAx','2DCo','2DSa','2+1DAx','2+1DCo','2+1DSa'}
        
        h.Type2 = 'XtYtZt';
        set(h.DispAll,'Visible','off')
        set(h.DispX,'Visible','on')
        set(h.DispY,'Visible','on')
        set(h.DispZ,'Visible','on')
        
        
        set(h.Show_3D,'Visible','off','Value',0)
        set(h.Show_2Dxy,'Visible','on','Value',0)
        set(h.Show_2Dxz,'Visible','on','Value',0)
        set(h.Show_2Dyz,'Visible','on','Value',0)
        set(h.Show_Xt,'Visible','on','Value',0)
        set(h.Show_Yt,'Visible','on','Value',0)
        set(h.Show_XtYtZt,'Visible','on','Value',1)
        set(h.Show_Zt,'Visible','on','Value',0)
        
        set(h.Show_Magnitude,'Visible','on','Value',0)
        
        
        
        
        
    case {'3D','3+1D','0+1D'}
        
        h.Type2 = 'XtYtZt';
        set(h.DispAll,'Visible','off')
        set(h.DispX,'Visible','on')
        set(h.DispY,'Visible','on')
        set(h.DispZ,'Visible','on')
        
        set(h.Show_3D,'Visible','on','Value',0)
        set(h.Show_2Dxy,'Visible','on','Value',0)
        set(h.Show_2Dxz,'Visible','on','Value',0)
        set(h.Show_2Dyz,'Visible','on','Value',0)
        set(h.Show_Xt,'Visible','on','Value',0)
        set(h.Show_Yt,'Visible','on','Value',0)
        set(h.Show_XtYtZt,'Visible','on','Value',1)
        set(h.Show_Zt,'Visible','on','Value',0)
        
        set(h.Show_Magnitude,'Visible','on','Value',0)
        
end

h = Plotting(h);

end

function Show_Magnitude_Callback(hObj,evd)
global h

h.ShowMagn = get(h.Show_Magnitude,'Value');
h = Plotting(h);

end

function P2_SelectionChangeFcn(hObj,evd)
global h
set(h.Show_Magnitude,'Value',0)
h.ShowMagn = get(h.Show_Magnitude,'Value');
newButton=get(evd.NewValue,'String');
% disp('Panel2')
switch newButton
    case 'x(t)'
        h.Type2 = 'Xt';
        
        
    case 'y(t)'
        h.Type2 = 'Yt';
        
        
    case 'z(t)'
        h.Type2 = 'Zt';
        
        
        
        
        
    case 'x(t),y(t),z(t)'
        
        h.Type2 = 'XtYtZt';
    case '(x,y,z)'
        
        
        
        h.Type2 = '(x,y,z)';
    case '(x,y)'
        
        h.Type2 = '2Dxy';
    case '(x,z)'
        
        h.Type2 = '2Dxz';
    case '(y,z)'
        
        
        h.Type2 = '2Dyz';
end

h = Plotting(h);


end


function P1_SelectionChangeFcn(hObj,evd)
global h
set(h.Show_Magnitude,'Value',0)
h.ShowMagn = get(h.Show_Magnitude,'Value');

newButton=get(evd.NewValue,'String');

switch newButton
    case 'k'
        h.Type1 = 'k';
    case 'g'
        h.Type1 = 'g';
    case 's'
        h.Type1 = 's';
        
end

h = Plotting(h);

end

function h = Plotting(h)



switch h.Type1
    case 'k'
        
        
        Array = h.khr.k;
    case 'g'
        Array = h.khr.g;
        
    case 's'
        Array = h.khr.s;
end

switch h.Type2
    case 'Xt'
        cla(h.DispX)
        cla(h.DispY)
        cla(h.DispZ)
        cla(h.DispAll)
        set(h.DispX,'Visible','off')
        set(h.DispY,'Visible','off')
        set(h.DispZ,'Visible','off')
        set(h.DispAll,'Visible','on')
        
        y = Array(1,:);
        x = h.khr.t;
        
        if h.ShowMagn
            y = abs(y);
        end

        axes(h.DispAll)
        plot(x,y,'r','linewidth',2)
                axis([min(x)-eps,max(x)+eps,min(y)-eps,max(y)+eps])

    case 'Yt'
        cla(h.DispX)
        cla(h.DispY)
        cla(h.DispZ)
        cla(h.DispAll)
        set(h.DispX,'Visible','off')
        set(h.DispY,'Visible','off')
        set(h.DispZ,'Visible','off')
        set(h.DispAll,'Visible','on')
        y = Array(2,:);
        x = h.khr.t;
        
        if h.ShowMagn
            y = abs(y);
        end

        axes(h.DispAll)
        plot(x,y,'b','linewidth',2)
        %         axis([Xmin,Xmax,Ymin1,Ymax1])
        axis([min(x)-eps,max(x)+eps,min(y)-eps,max(y)+eps])
        hand = h.DispAll;

    case 'Zt'
        
        cla(h.DispX)
        cla(h.DispY)
        cla(h.DispZ)
        cla(h.DispAll)
        set(h.DispX,'Visible','off')
        set(h.DispY,'Visible','off')
        set(h.DispZ,'Visible','off')
        set(h.DispAll,'Visible','on')
        
        y = Array(3,:);
        x = h.khr.t;
        
        if h.ShowMagn
            y = abs(y);
        end

        axes(h.DispAll)
        plot(x,y,'g','linewidth',2)
        axis([min(x)-eps,max(x)+eps,min(y)-eps,max(y)+eps])

    case 'XtYtZt'
        
        
        if h.ShowMagn
            cla(h.DispX)
            cla(h.DispY)
            cla(h.DispZ)
            cla(h.DispAll)
            set(h.DispX,'Visible','off')
            set(h.DispY,'Visible','off')
            set(h.DispZ,'Visible','off')
            set(h.DispAll,'Visible','on')
            
            y1 = Array(1,:);
            y2 = Array(2,:);
            y3 = Array(3,:);
            x = h.khr.t;
            
            y = sqrt(y1.^2+y2.^2+y3.^2);

            axes(h.DispAll)
            plot(x,y,'k','linewidth',2)
            axis([min(x)-eps,max(x)+eps,min(y)-eps,max(y)+eps])

        else
            cla(h.DispX)
            cla(h.DispY)
            cla(h.DispZ)
            cla(h.DispAll)
            set(h.DispX,'Visible','on')
            set(h.DispY,'Visible','on')
            set(h.DispZ,'Visible','on')
            set(h.DispAll,'Visible','off')
            
            y1 = Array(1,:);
            y2 = Array(2,:);
            y3 = Array(3,:);
            x = h.khr.t;

            axes(h.DispX)
            plot(x,y1,'r','linewidth',2)
axis([min(x)-eps,max(x)+eps,min(y1)-eps,max(y1)+eps])
            axes(h.DispY)
            plot(x,y2,'b','linewidth',2)
            axis([min(x)-eps,max(x)+eps,min(y2)-eps,max(y2)+eps])
            axes(h.DispZ)
            plot(x,y3,'g','linewidth',2)
            axis([min(x)-eps,max(x)+eps,min(y3)-eps,max(y3)+eps])
        end
        
    case '3D'
        
        cla(h.DispX)
        cla(h.DispY)
        cla(h.DispZ)
        cla(h.DispAll)
        cla(h.DispPos)
        set(h.DispX,'Visible','off')
        set(h.DispY,'Visible','off')
        set(h.DispZ,'Visible','off')
        set(h.DispAll,'Visible','on')
        set(h.DispPos,'Visible','off')
        x = Array(1,:);
        y = Array(2,:);
        z = Array(3,:);
        
        
        
        axes(h.DispAll)
        plot3(x,y,z,'k','linewidth',2)
        axis([min(x)-eps,max(x)+eps,min(y)-eps,max(y)+eps,min(z)-eps,max(z)+eps])
        
    case '2Dxy'
        cla(h.DispX)
        cla(h.DispY)
        cla(h.DispZ)
        cla(h.DispAll)
        
        set(h.DispX,'Visible','off')
        set(h.DispY,'Visible','off')
        set(h.DispZ,'Visible','off')
        set(h.DispAll,'Visible','on')
        x = Array(1,:);
        y = Array(2,:);
        
        axes(h.DispAll)
        plot(x,y,'m','linewidth',2)
        axis([min(x)-eps,max(x)+eps,min(y)-eps,max(y)+eps])
        
    case '2Dxz'
        cla(h.DispX)
        cla(h.DispY)
        cla(h.DispZ)
        cla(h.DispAll)
        
        set(h.DispX,'Visible','off')
        set(h.DispY,'Visible','off')
        set(h.DispZ,'Visible','off')
        set(h.DispAll,'Visible','on')
        x = Array(1,:);
        z = Array(3,:);
        axes(h.DispAll)
        plot(x,z,'linewidth',2,'Color',[210	105	 30]./255)
        axis([min(x)-eps,max(x)+eps,min(z)-eps,max(z)+eps])
    case '2Dyz'
        cla(h.DispX)
        cla(h.DispY)
        cla(h.DispZ)
        cla(h.DispAll)
        
        set(h.DispX,'Visible','off')
        set(h.DispY,'Visible','off')
        set(h.DispZ,'Visible','off')
        set(h.DispAll,'Visible','on')
        y = Array(2,:);
        z = Array(3,:);
        axes(h.DispAll)
        plot(y,z,'c','linewidth',2)
        axis([min(y)-eps,max(y)+eps,min(z)-eps,max(z)+eps])
end

drawnow


end

function listbox1_Callback(hObj,evd)

end

function h = Populate_Listbox(h)


names = fieldnames(h.khr);
vals = cell(length(names),1);
for n = 1:length(names)
    
    if ischar(getfield(h.khr,names{n}))
        vals{n} = getfield(h.khr,names{n});
    elseif isnumeric(getfield(h.khr,names{n}))
        [A,B] = size(getfield(h.khr,names{n}));
        
        if A == 1 && B == 1
            vals{n} = num2str(getfield(h.khr,names{n}));
        else
            if B < 6
                temp = getfield(h.khr,names{n});
                
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
                vals{n} = sprintf('<%ix%i double>',A,B);
            end
        end
    elseif iscell(getfield(h.khr,names{n}))
        [A,B] = size(getfield(h.khr,names{n}));
        vals{n} = sprintf('<%ix%i cell>',A,B);
    elseif isstruct(getfield(h.khr,names{n}))
        [A,B] = size(getfield(h.khr,names{n}));
        vals{n} = sprintf('<%ix%i struct>',A,B);
    end
end

List = cell(length(names),1);

for n = 1:length(names)
    List{n} = sprintf('%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s',names{n},vals{n});
    
end

set(h.listbox1,'String',List)


end


