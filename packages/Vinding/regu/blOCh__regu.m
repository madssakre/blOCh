function opt = blOCh__regu(spc,khr,opt)
% function opt = blOCh__regu(spc,khr,opt)
%
%  This is a wrapper for parts of the regularization tool box by Per Christian
%  Hansen who has this license attached:
%
%   
% Copyright (c) 2015, Per Christian Hansen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of DTU Compute nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
%
%
%   This wrapper script is licensed as:
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
%
%   Down to business:
%
%   This version:
%   Tikhonov regularization and best lambda finding by L-curve of PCHansen
%   regu method.
%
%   Output is uo and vo and more. Since the regu does not give u and v per
%   iteration uo and vo are those corresponding to the best lambda value.
%
%   Crude clipping is availabe as well as scaling from STA to LTA.
%   See further limitations so far on clipping (units etc.)
%
%   When using Scale, it is important to set Mdflip even if you load a Md
%   from file where it might ignore Mdflip because Md is defined already.
%   This is important because the Scale is calculated from Mdflip

tic_Opt_Wall_Clock_Time = tic;

Par.lambda = [1e-3,1e-4,1e-5,1e-6,1e-7,1e-8];
Par.lambda = logspace(-16,16,25); % shifts

Par.Constr.peak_uv.Lim = 1000;
Par.Constr.peak_uv.Unit = 'Hz';

Par.Constr.peak_a.Lim = 1000;
Par.Constr.peak_a.Unit = 'Hz';

Par.Constr.peak_uv.Cope = 'Ignore';
Par.Constr.peak_a.Cope = 'Constr';

Par.Show = 1;

Par.ScaleX2Flip = -1 ;
Par.Temporary_STA_Flip = -1;


Par.Method = 'Tikhonov';

opt.Par = blOCh__khr('Get_NewPar',[],[],[],opt.Par,Par);

opt.Par.Num_lambda = numel(opt.Par.lambda);


if opt.Par.Temporary_STA_Flip > 0
    
    B1map = blOCh__spc('list2grid',[],[],[],[],spc.B1map,[spc.R,spc.Rv],spc.idxtot,1);
    Subject = sum(abs(B1map),5);
    idx_Subject = find(Subject);
    
    Md = blOCh__spc('list2grid',[],[],[],[],spc.Md,[spc.R,spc.Rv],spc.idxtot,3);
    
    
    Mdz = Md(:,:,:,:,3);
    
    IdxIn = find(Mdz~=1);
    
    Mdx = Md(:,:,:,:,1);
    Mdy = Md(:,:,:,:,2);
    
    Mdx(IdxIn) = sind(opt.Par.Temporary_STA_Flip);
    Mdz(IdxIn) = cosd(opt.Par.Temporary_STA_Flip);
    
    
    Mdx_subject = Mdx(idx_Subject);
    Mdy_subject = Mdy(idx_Subject);
    Mdz_subject = Mdz(idx_Subject);
    
    Md = [Mdx_subject(:)';Mdy_subject(:)';Mdz_subject(:)'];
    
    Md = Md(:);
    
else
    Md = spc.Md;
end


%%


if opt.Par.Num_lambda > 1
    
    
    % this will pick the x that best matches spc.Md and has the minimum
    % ||x||^2. There will be no x clipping in this stage
    % This is the normal procedure
    
    b = complex(Md(1:3:end),Md(2:3:end));
    A = Get_System_Matrix(spc,khr);
    tic_Opt_Only_Time = tic;
    switch opt.Par.Method
        case 'Tikhonov'
            [U,s,V] = csvd (A);
            [x,norm_Axb_vs_lam,norm_x_vs_lam]= tikhonov (U,s,V,b,opt.Par.lambda);
            
            [reg_corner_,rho_,eta_,reg_param_] = l_curve (U,s,b);
        otherwise
            error('Only Tikhonov has been imlemented at the time of writing')
    end
    opt.Opt_Only_Time = toc(tic_Opt_Only_Time);
    opt.x = x; % store for external use
    Idx = find(abs(opt.Par.lambda-reg_corner_)==min(abs(opt.Par.lambda-reg_corner_)));
    opt.Par.Best_lambda_Number = Idx;
    opt.Par.Best_lambda = opt.Par.lambda(Idx);
    
    
    
    if opt.Show
        figure
        loglog(rho_,eta_)
        hold on
        loglog(norm_Axb_vs_lam(Idx),norm_x_vs_lam(Idx),'*')
    end
    
    
else
    
    
    % this will induce no lambda optimization since only one lambda
    % value is submitted. There will be no x clipping in this stage
    % NOT tested thoroughly
    
    b = complex(spc.Md(1:3:end),spc.Md(2:3:end));
    A = Get_System_Matrix(spc,khr);
    switch opt.Par.Method
        case 'Tikhonov'
            [U,s,V] = csvd (A);
            [x,norm_Axb_vs_lam,norm_x_vs_lam]= tikhonov (U,s,V,b,opt.Par.lambda);
            
            %                 [reg_corner_,rho_,eta_,reg_param_] = l_curve (U,s,b);
        otherwise
            error('Only Tikhonov has been imlemented at the time of writing')
    end
    opt.x = x; % store for external use
    Idx = 1;
    opt.Par.Best_lambda_Number = 1;
    opt.Par.Best_lambda = opt.Par.lambda;
    
    
end
% this will scale x according to
% this will be no x clipping in this stage

if opt.Par.Temporary_STA_Flip > 0
    
    CurrentFlip = opt.Par.Temporary_STA_Flip;
    Scale = opt.Par.ScaleX2Flip./CurrentFlip;
    
else
    if opt.Par.ScaleX2Flip > 0
        
        CurrentFlip = spc.Mdflip(1);
        Scale = opt.Par.ScaleX2Flip./CurrentFlip;
        
        
        
    else
        
        
        Scale = 1;
        
    end
end



%%
if spc.pTx > 1
    opt.uo = zeros(spc.pTx,opt.N);
    opt.vo = zeros(spc.pTx,opt.N);
    g = 1;
    for s = 1:spc.pTx
        tempuv = x(:,Idx).*Scale;
        tempuv = tempuv(:).';
        
        opt.uo(s,:) = real(tempuv(g:g+khr.N-1))*khr.gamma;
        opt.vo(s,:) = imag(tempuv(g:g+khr.N-1))*khr.gamma;
        g = g + opt.N;
    end
    
else
    tempuv = x(:,Idx).*Scale;
    tempuv = tempuv(:).';
    opt.uo = real(tempuv)*khr.gamma;
    opt.vo = imag(tempuv)*khr.gamma;
end




%%

if strcmp(opt.Par.Constr.peak_uv.Cope,'Constr')
    
    switch opt.Par.Constr.peak_uv.Unit
        case 'Hz'
            
            Idx1 = find(opt.uo/2/pi>opt.Par.Constr.peak_uv.Lim);
            Idx2 = find(opt.uo/2/pi<-opt.Par.Constr.peak_uv.Lim);
            Idx3 = find(opt.vo/2/pi>opt.Par.Constr.peak_uv.Lim);
            Idx4 = find(opt.vo/2/pi<-opt.Par.Constr.peak_uv.Lim);
            
            if ~isempty(Idx1)
            opt.uo(Idx1) = opt.Par.Constr.peak_uv.Lim*2*pi;
            end
            if ~isempty(Idx2)
            opt.uo(Idx2) = -opt.Par.Constr.peak_uv.Lim*2*pi;
            end
            if ~isempty(Idx3)
            opt.vo(Idx3) = opt.Par.Constr.peak_uv.Lim*2*pi;
            end
            if ~isempty(Idx4)
            opt.vo(Idx4) = -opt.Par.Constr.peak_uv.Lim*2*pi;
            end
            if ~isempty(Idx1) || ~isempty(Idx2) || ~isempty(Idx3) || ~isempty(Idx4)
            disp('Clipping uo and/ or vo')
            opt.Par.Clipped_uv = 1;
            else
                opt.Par.Clipped_uv = 0;
            end
            
        otherwise
            error('Only clipping with units in Hz has been imlemented at the time of writing')
    end
end

opt.ao = abs(complex(opt.uo,opt.vo));
opt.po = angle(complex(opt.uo,opt.vo));

if strcmp(opt.Par.Constr.peak_a.Cope,'Constr')
    
    switch opt.Par.Constr.peak_uv.Unit
        case 'Hz'
            
            Idx = find(opt.ao/2/pi>opt.Par.Constr.peak_a.Lim);
            
            if ~isempty(Idx)
            opt.ao(Idx) = opt.Par.Constr.peak_a.Lim*2*pi;
            
            disp('Clipping ao')
            opt.Par.Clipped_a = 1;
            else
                opt.Par.Clipped_a = 0;
            end
            
            
        otherwise
            error('Only clipping with units in Hz has been imlemented at the time of writing')
    end
    
end





opt.uo = real(opt.ao.*exp(1i.*opt.po));
opt.vo = imag(opt.ao.*exp(1i.*opt.po));

opt.fo  = diff(cat(2,opt.po,zeros(spc.pTx,1)),1,2)./opt.dt;

opt.gxo = opt.g(1,:);
opt.gyo = opt.g(2,:);
opt.gzo = opt.g(3,:);

opt.k = -1;
opt.Fun = norm_Axb_vs_lam;
opt.Xnorm = norm_x_vs_lam;

opt.Opt_Wall_Clock_Time = toc(tic_Opt_Wall_Clock_Time);
opt.Go = false;
end

function Afull = Get_System_Matrix(spc,khr)
D = cell(spc.pTx,1);
for s = 1:spc.pTx
    D{s} = diag(spc.B1map(:,s));
end

A1 = 1i.*khr.gamma*khr.dt;
A2 = exp(1i.*(spc.w0map+2*pi*spc.V).*(khr.t-khr.T));
A3 = exp(1i.*2*pi*(spc.X*khr.k(1,:)+spc.Y*khr.k(2,:)+spc.Z*khr.k(3,:)));
A = A1.*A2.*A3;

% Assemble
% tic
Afull = [];
for s = 1:spc.pTx
    Afull = [Afull,D{s}*A];
end
% toc
% tic
% Afull = complex(zeros(size(A,1),size(A,2)*spc.pTx),zeros(size(A,1),size(A,2)*spc.pTx));
% k = 1;
% for s = 1:spc.pTx
%     Afull(:,k:k+size(A,2)-1) = D{s}*A;
%     k = k+size(A,2);
% end
% toc
% max(max(abs(Afull-Afull1)))
end



function x_prog = Rescalex(x,opt)
if strcmp(opt.Par.optRF,'on') && strcmp(opt.Par.optGRA,'on')
    x_prog = [x(opt.idx_u(1):opt.idx_v(2)).*opt.maxRF,x(opt.idx_gx(1):opt.idx_gz(2)).*opt.maxG];
elseif strcmp(opt.Par.optRF,'off') && strcmp(opt.Par.optGRA,'on')
    x_prog = x.*opt.maxG;
elseif strcmp(opt.Par.optRF,'on') && strcmp(opt.Par.optGRA,'off')
    x_prog = x.*opt.maxRF;
end

end





%%

function [reg_corner,rho,eta,reg_param] = l_curve(U,sm,b,method,L,V)
%L_CURVE Plot the L-curve and find its "corner".
%
% Copyright (c) 2015, Per Christian Hansen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of DTU Compute nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% [reg_corner,rho,eta,reg_param] =
%                  l_curve(U,s,b,method)
%                  l_curve(U,sm,b,method)  ,  sm = [sigma,mu]
%                  l_curve(U,s,b,method,L,V)
%
% Plots the L-shaped curve of eta, the solution norm || x || or
% semi-norm || L x ||, as a function of rho, the residual norm
% || A x - b ||, for the following methods:
%    method = 'Tikh'  : Tikhonov regularization   (solid line )
%    method = 'tsvd'  : truncated SVD or GSVD     (o markers  )
%    method = 'dsvd'  : damped SVD or GSVD        (dotted line)
%    method = 'mtsvd' : modified TSVD             (x markers  )
% The corresponding reg. parameters are returned in reg_param.  If no
% method is specified then 'Tikh' is default.  For other methods use plot_lc.
%
% Note that 'Tikh', 'tsvd' and 'dsvd' require either U and s (standard-
% form regularization) computed by the function csvd, or U and sm (general-
% form regularization) computed by the function cgsvd, while 'mtvsd'
% requires U and s as well as L and V computed by the function csvd.
%
% If any output arguments are specified, then the corner of the L-curve
% is identified and the corresponding reg. parameter reg_corner is
% returned.  Use routine l_corner if an upper bound on eta is required.

% Reference: P. C. Hansen & D. P. O'Leary, "The use of the L-curve in
% the regularization of discrete ill-posed problems",  SIAM J. Sci.
% Comput. 14 (1993), pp. 1487-1503.

% Per Christian Hansen, DTU Compute, October 27, 2010.

% Set defaults.
if (nargin==3), method='Tikh'; end  % Tikhonov reg. is default.
npoints = 200;  % Number of points on the L-curve for Tikh and dsvd.
smin_ratio = 16*eps;  % Smallest regularization parameter.

% Initialization.
[m,n] = size(U); [p,ps] = size(sm);
if (nargout > 0), locate = 1; else locate = 0; end
beta = U'*b; beta2 = norm(b)^2 - norm(beta)^2;
if (ps==1)
    s = sm; beta = beta(1:p);
else
    s = sm(p:-1:1,1)./sm(p:-1:1,2); beta = beta(p:-1:1);
end
xi = beta(1:p)./s;
xi( isinf(xi) ) = 0;

if (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4))
    
    eta = zeros(npoints,1); rho = eta; reg_param = eta; s2 = s.^2;
    reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
    ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
    for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
    for i=1:npoints
        f = s2./(s2 + reg_param(i)^2);
        eta(i) = norm(f.*xi);
        rho(i) = norm((1-f).*beta(1:p));
    end
    if (m > n & beta2 > 0), rho = sqrt(rho.^2 + beta2); end
    marker = '-'; txt = 'Tikh.';
    
elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4))
    
    eta = zeros(p,1); rho = eta;
    eta(1) = abs(xi(1))^2;
    for k=2:p, eta(k) = eta(k-1) + abs(xi(k))^2; end
    eta = sqrt(eta);
    if (m > n)
        if (beta2 > 0), rho(p) = beta2; else rho(p) = eps^2; end
    else
        rho(p) = eps^2;
    end
    for k=p-1:-1:1, rho(k) = rho(k+1) + abs(beta(k+1))^2; end
    rho = sqrt(rho);
    reg_param = (1:p)'; marker = 'o';
    if (ps==1)
        U = U(:,1:p); txt = 'TSVD';
    else
        U = U(:,1:p); txt = 'TGSVD';
    end
    
elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))
    
    eta = zeros(npoints,1); rho = eta; reg_param = eta;
    reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
    ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
    for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
    for i=1:npoints
        f = s./(s + reg_param(i));
        eta(i) = norm(f.*xi);
        rho(i) = norm((1-f).*beta(1:p));
    end
    if (m > n & beta2 > 0), rho = sqrt(rho.^2 + beta2); end
    marker = ':';
    if (ps==1), txt = 'DSVD'; else txt = 'DGSVD'; end
    
elseif (strncmp(method,'mtsv',4))
    
    if (nargin~=6)
        error('The matrices L and V must also be specified')
    end
    [p,n] = size(L); rho = zeros(p,1); eta = rho;
    [Q,R] = qr(L*V(:,n:-1:n-p),0);
    for i=1:p
        k = n-p+i;
        Lxk = L*V(:,1:k)*xi(1:k);
        zk = R(1:n-k,1:n-k)\(Q(:,1:n-k)'*Lxk); zk = zk(n-k:-1:1);
        eta(i) = norm(Q(:,n-k+1:p)'*Lxk);
        if (i < p)
            rho(i) = norm(beta(k+1:n) + s(k+1:n).*zk);
        else
            rho(i) = eps;
        end
    end
    if (m > n & beta2 > 0), rho = sqrt(rho.^2 + beta2); end
    reg_param = (n-p+1:n)'; txt = 'MTSVD';
    U = U(:,reg_param); sm = sm(reg_param);
    marker = 'x'; ps = 2;  % General form regularization.
    
else
    error('Illegal method')
end

% Locate the "corner" of the L-curve, if required.
if (locate)
    [reg_corner,rho_c,eta_c] = l_corner(rho,eta,reg_param,U,sm,b,method);
end
if 0 % msv maybe in the future
    % Make plot.
    plot_lc(rho,eta,marker,ps,reg_param);
    if locate
        ax = axis;
        HoldState = ishold; hold on;
        loglog([min(rho)/100,rho_c],[eta_c,eta_c],':r',...
            [rho_c,rho_c],[min(eta)/100,eta_c],':r')
        title(['L-curve, ',txt,' corner at ',num2str(reg_corner)]);
        axis(ax)
        if (~HoldState), hold off; end
    end
end
end


function [x_lambda,rho,eta] = tikhonov(U,s,V,b,lambda,x_0)
%TIKHONOV Tikhonov regularization.
%
% Copyright (c) 2015, Per Christian Hansen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of DTU Compute nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% [x_lambda,rho,eta] = tikhonov(U,s,V,b,lambda,x_0)
% [x_lambda,rho,eta] = tikhonov(U,sm,X,b,lambda,x_0) ,  sm = [sigma,mu]
%
% Computes the Tikhonov regularized solution x_lambda, given the SVD or
% GSVD as computed via csvd or cgsvd, respectively.  If the SVD is used,
% i.e. if U, s, and V are specified, then standard-form regularization
% is applied:
%    min { || A x - b ||^2 + lambda^2 || x - x_0 ||^2 } .
% If, on the other hand, the GSVD is used, i.e. if U, sm, and X are
% specified, then general-form regularization is applied:
%    min { || A x - b ||^2 + lambda^2 || L (x - x_0) ||^2 } .
%
% If an initial estimate x_0 is not specified, then x_0 = 0 is used.
%
% Note that x_0 cannot be used if A is underdetermined and L ~= I.
%
% If lambda is a vector, then x_lambda is a matrix such that
%    x_lambda = [ x_lambda(1), x_lambda(2), ... ] .
%
% The solution norm (standard-form case) or seminorm (general-form
% case) and the residual norm are returned in eta and rho.

% Per Christian Hansen, DTU Compute, April 14, 2003.

% Reference: A. N. Tikhonov & V. Y. Arsenin, "Solutions of Ill-Posed
% Problems", Wiley, 1977.

% Initialization.
if (min(lambda)<0)
    error('Illegal regularization parameter lambda')
end
m = size(U,1);
n = size(V,1);
[p,ps] = size(s);
beta = U(:,1:p)'*b;
zeta = s(:,1).*beta;
ll = length(lambda); x_lambda = zeros(n,ll);
rho = zeros(ll,1); eta = zeros(ll,1);

% Treat each lambda separately.
if (ps==1)
    
    % The standard-form case.
    if (nargin==6), omega = V'*x_0; end
    for i=1:ll
        if (nargin==5)
            x_lambda(:,i) = V(:,1:p)*(zeta./(s.^2 + lambda(i)^2));
            rho(i) = lambda(i)^2*norm(beta./(s.^2 + lambda(i)^2));
        else
            x_lambda(:,i) = V(:,1:p)*...
                ((zeta + lambda(i)^2*omega)./(s.^2 + lambda(i)^2));
            rho(i) = lambda(i)^2*norm((beta - s.*omega)./(s.^2 + lambda(i)^2));
        end
        eta(i) = norm(x_lambda(:,i));
    end
    if (nargout > 1 & size(U,1) > p)
        rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
    end
    
elseif (m>=n)
    
    % The overdetermined or square general-form case.
    gamma2 = (s(:,1)./s(:,2)).^2;
    if (nargin==6), omega = V\x_0; omega = omega(1:p); end
    if (p==n)
        x0 = zeros(n,1);
    else
        x0 = V(:,p+1:n)*U(:,p+1:n)'*b;
    end
    for i=1:ll
        if (nargin==5)
            xi = zeta./(s(:,1).^2 + lambda(i)^2*s(:,2).^2);
            x_lambda(:,i) = V(:,1:p)*xi + x0;
            rho(i) = lambda(i)^2*norm(beta./(gamma2 + lambda(i)^2));
        else
            xi = (zeta + lambda(i)^2*(s(:,2).^2).*omega)./...
                (s(:,1).^2 + lambda(i)^2*s(:,2).^2);
            x_lambda(:,i) = V(:,1:p)*xi + x0;
            rho(i) = lambda(i)^2*norm((beta - s(:,1).*omega)./...
                (gamma2 + lambda(i)^2));
        end
        eta(i) = norm(s(:,2).*xi);
    end
    if (nargout > 1 & size(U,1) > p)
        rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
    end
    
else
    
    % The underdetermined general-form case.
    gamma2 = (s(:,1)./s(:,2)).^2;
    if (nargin==6), error('x_0 not allowed'), end
    if (p==m)
        x0 = zeros(n,1);
    else
        x0 = V(:,p+1:m)*U(:,p+1:m)'*b;
    end
    for i=1:ll
        xi = zeta./(s(:,1).^2 + lambda(i)^2*s(:,2).^2);
        x_lambda(:,i) = V(:,1:p)*xi + x0;
        rho(i) = lambda(i)^2*norm(beta./(gamma2 + lambda(i)^2));
        eta(i) = norm(s(:,2).*xi);
    end
    
end

end


function [U,s,V] = csvd(A,tst)
%CSVD Compact singular value decomposition.
%
% Copyright (c) 2015, Per Christian Hansen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of DTU Compute nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% s = csvd(A)
% [U,s,V] = csvd(A)
% [U,s,V] = csvd(A,'full')
%
% Computes the compact form of the SVD of A:
%    A = U*diag(s)*V',
% where
%    U  is  m-by-min(m,n)
%    s  is  min(m,n)-by-1
%    V  is  n-by-min(m,n).
%
% If a second argument is present, the full U and V are returned.
%
% Per Christian Hansen, IMM, 06/22/93.
%


if (nargin==1)
    if (nargout > 1)
        [m,n] = size(A);
        if (m >= n)
            [U,s,V] = svd(full(A),0); s = diag(s);
        else
            [V,s,U] = svd(full(A)',0); s = diag(s);
        end
    else
        U = svd(full(A));
    end
else
    if (nargout > 1)
        [U,s,V] = svd(full(A)); s = diag(s);
    else
        U = svd(full(A));
    end
end

end


function [reg_c,rho_c,eta_c] = l_corner(rho,eta,reg_param,U,s,b,method,M)
%L_CORNER Locate the "corner" of the L-curve.
%
%
% Copyright (c) 2015, Per Christian Hansen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of DTU Compute nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% [reg_c,rho_c,eta_c] =
%        l_corner(rho,eta,reg_param)
%        l_corner(rho,eta,reg_param,U,s,b,method,M)
%        l_corner(rho,eta,reg_param,U,sm,b,method,M) ,  sm = [sigma,mu]
%
% Locates the "corner" of the L-curve in log-log scale.
%
% It is assumed that corresponding values of || A x - b ||, || L x ||,
% and the regularization parameter are stored in the arrays rho, eta,
% and reg_param, respectively (such as the output from routine l_curve).
%
% If nargin = 3, then no particular method is assumed, and if
% nargin = 2 then it is issumed that reg_param = 1:length(rho).
%
% If nargin >= 6, then the following methods are allowed:
%    method = 'Tikh'  : Tikhonov regularization
%    method = 'tsvd'  : truncated SVD or GSVD
%    method = 'dsvd'  : damped SVD or GSVD
%    method = 'mtsvd' : modified TSVD,
% and if no method is specified, 'Tikh' is default.  If the Spline Toolbox
% is not available, then only 'Tikh' and 'dsvd' can be used.
%
% An eighth argument M specifies an upper bound for eta, below which
% the corner should be found.
%
% Per Christian Hansen, DTU Compute, January 31, 2015.
%
% Ensure that rho and eta are column vectors.
rho = rho(:); eta = eta(:);

% Set default regularization method.
if (nargin <= 3)
    method = 'none';
    if (nargin==2), reg_param = (1:length(rho))'; end
else
    if (nargin==6), method = 'Tikh'; end
end

% Set this logical variable to 1 (true) if the corner algorithm
% should always be used, even if the Spline Toolbox is available.
alwayscorner = 0;

% Set threshold for skipping very small singular values in the
% analysis of a discrete L-curve.
s_thr = eps;  % Neglect singular values less than s_thr.

% Set default parameters for treatment of discrete L-curve.
deg   = 2;  % Degree of local smooting polynomial.
q     = 2;  % Half-width of local smoothing interval.
order = 4;  % Order of fitting 2-D spline curve.

% Initialization.
if (length(rho) < order)
    error('Too few data points for L-curve analysis')
end
if (nargin > 3)
    [p,ps] = size(s); [m,n] = size(U);
    beta = U'*b; b0 = b - U*beta;
    if (ps==2)
        s = s(p:-1:1,1)./s(p:-1:1,2);
        beta = beta(p:-1:1);
    end
    xi = beta./s;
    if (m>n)  % Take of the least-squares residual.
        beta = [beta;norm(b0)];
    end
end

% Restrict the analysis of the L-curve according to M (if specified).
if (nargin==8)
    index = find(eta < M);
    rho = rho(index); eta = eta(index); reg_param = reg_param(index);
end

if (strncmp(method,'Tikh',4) || strncmp(method,'tikh',4))
    
    % The L-curve is differentiable; computation of curvature in
    % log-log scale is easy.
    
    % Compute g = - curvature of L-curve.
    g = blOCh__lcfun(reg_param,s,beta,xi);
    
    % Locate the corner.  If the curvature is negative everywhere,
    % then define the leftmost point of the L-curve as the corner.
    [~,gi] = min(g);
    reg_c = fminbnd('blOCh__lcfun',...
        reg_param(min(gi+1,length(g))),reg_param(max(gi-1,1)),...
        optimset('Display','off'),s,beta,xi); % Minimizer.
    kappa_max = - blOCh__lcfun(reg_c,s,beta,xi); % Maximum curvature.
    
    if (kappa_max < 0)
        lr = length(rho);
        reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
    else
        f = (s.^2)./(s.^2 + reg_c^2);
        eta_c = norm(f.*xi);
        rho_c = norm((1-f).*beta(1:length(f)));
        if (m>n), rho_c = sqrt(rho_c^2 + norm(b0)^2); end
    end
    
elseif (strncmp(method,'tsvd',4) || strncmp(method,'tgsv',4) || ...
        strncmp(method,'mtsv',4) || strncmp(method,'none',4))
    
    % Use the adaptive pruning algorithm to find the corner, if the
    % Spline Toolbox is not available.
    if ~exist('splines','dir') || alwayscorner
        %error('The Spline Toolbox in not available so l_corner cannot be used')
        reg_c = corner(rho,eta);
        rho_c = rho(reg_c);
        eta_c = eta(reg_c);
        return
    end
    
    % Otherwise use local smoothing followed by fitting a 2-D spline curve
    % to the smoothed discrete L-curve. Restrict the analysis of the L-curve
    % according to s_thr.
    if (nargin > 3)
        if (nargin==8)       % In case the bound M is in action.
            s = s(index,:);
        end
        index = find(s > s_thr);
        rho = rho(index); eta = eta(index); reg_param = reg_param(index);
    end
    
    % Convert to logarithms.
    lr = length(rho);
    lrho = log(rho); leta = log(eta); slrho = lrho; sleta = leta;
    
    % For all interior points k = q+1:length(rho)-q-1 on the discrete
    % L-curve, perform local smoothing with a polynomial of degree deg
    % to the points k-q:k+q.
    v = (-q:q)'; A = zeros(2*q+1,deg+1); A(:,1) = ones(length(v),1);
    for j = 2:deg+1, A(:,j) = A(:,j-1).*v; end
    for k = q+1:lr-q-1
        cr = A\lrho(k+v); slrho(k) = cr(1);
        ce = A\leta(k+v); sleta(k) = ce(1);
    end
    
    % Fit a 2-D spline curve to the smoothed discrete L-curve.
    sp = spmak((1:lr+order),[slrho';sleta']);
    pp = ppbrk(sp2pp(sp),[4,lr+1]);
    
    % Extract abscissa and ordinate splines and differentiate them.
    % Compute as many function values as default in spleval.
    P     = spleval(pp);  dpp   = fnder(pp);
    D     = spleval(dpp); ddpp  = fnder(pp,2);
    DD    = spleval(ddpp);
    ppx   = P(1,:);       ppy   = P(2,:);
    dppx  = D(1,:);       dppy  = D(2,:);
    ddppx = DD(1,:);      ddppy = DD(2,:);
    
    % Compute the corner of the discretized .spline curve via max. curvature.
    % No need to refine this corner, since the final regularization
    % parameter is discrete anyway.
    % Define curvature = 0 where both dppx and dppy are zero.
    k1    = dppx.*ddppy - ddppx.*dppy;
    k2    = (dppx.^2 + dppy.^2).^(1.5);
    I_nz  = find(k2 ~= 0);
    kappa = zeros(1,length(dppx));
    kappa(I_nz) = -k1(I_nz)./k2(I_nz);
    [kmax,ikmax] = max(kappa);
    x_corner = ppx(ikmax); y_corner = ppy(ikmax);
    
    % Locate the point on the discrete L-curve which is closest to the
    % corner of the spline curve.  Prefer a point below and to the
    % left of the corner.  If the curvature is negative everywhere,
    % then define the leftmost point of the L-curve as the corner.
    if (kmax < 0)
        reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
    else
        index = find(lrho < x_corner & leta < y_corner);
        if ~isempty(index)
            [~,rpi] = min((lrho(index)-x_corner).^2 + (leta(index)-y_corner).^2);
            rpi = index(rpi);
        else
            [~,rpi] = min((lrho-x_corner).^2 + (leta-y_corner).^2);
        end
        reg_c = reg_param(rpi); rho_c = rho(rpi); eta_c = eta(rpi);
    end
    
elseif (strncmp(method,'dsvd',4) || strncmp(method,'dgsv',4))
    
    % The L-curve is differentiable; computation of curvature in
    % log-log scale is easy.
    
    % Compute g = - curvature of L-curve.
    g = blOCh__lcfun(reg_param,s,beta,xi,1);
    
    % Locate the corner.  If the curvature is negative everywhere,
    % then define the leftmost point of the L-curve as the corner.
    [~,gi] = min(g);
    reg_c = fminbnd('blOCh__lcfun',...
        reg_param(min(gi+1,length(g))),reg_param(max(gi-1,1)),...
        optimset('Display','off'),s,beta,xi,1); % Minimizer.
    kappa_max = - blOCh__lcfun(reg_c,s,beta,xi,1); % Maximum curvature.
    
    if (kappa_max < 0)
        lr = length(rho);
        reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
    else
        f = s./(s + reg_c);
        eta_c = norm(f.*xi);
        rho_c = norm((1-f).*beta(1:length(f)));
        if (m>n), rho_c = sqrt(rho_c^2 + norm(b0)^2); end
    end
    
else
    error('Illegal method')
end

end


