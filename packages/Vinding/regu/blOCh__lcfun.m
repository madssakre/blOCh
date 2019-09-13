function g = blOCh__lcfun(lambda,s,beta,xi,fifth)
% function g = blOCh__lcfun(lambda,s,beta,xi,fifth)
%
%  This is a part of the regularization tool box by Per Christian
%  Hansen who has this license:
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
%   Down to business:
%
% Auxiliary routine for l_corner; computes the NEGATIVE of the curvature.
% Note: lambda may be a vector.  PCH, DTU Compute, Jan. 31, 2015.
%
% Due to some matlab issue this function cannot be embedded in the
% blOCh__5__regu file because fminbnd calls it and it needs it in a
% separate file. Mads 180507

% Initialization.
phi = zeros(size(lambda)); dphi = phi; psi = phi; dpsi = phi;
eta = phi; rho = phi;
if length(beta) > length(s)  % A possible least squares residual.
    LS = true;
    rhoLS2 = beta(end)^2;
    beta = beta(1:end-1);
else
    LS = false;
end

% Compute some intermediate quantities.
for i = 1:length(lambda)
    if (nargin==4)
        f  = (s.^2)./(s.^2 + lambda(i)^2);
    else
        f  = s./(s + lambda(i));
    end
    cf = 1 - f;
    eta(i) = norm(f.*xi);
    rho(i) = norm(cf.*beta);
    f1 = -2*f.*cf/lambda(i);
    f2 = -f1.*(3-4*f)/lambda(i);
    phi(i)  = sum(f.*f1.*abs(xi).^2);
    psi(i)  = sum(cf.*f1.*abs(beta).^2);
    dphi(i) = sum((f1.^2 + f.*f2).*abs(xi).^2);
    dpsi(i) = sum((-f1.^2 + cf.*f2).*abs(beta).^2);
end
if LS  % Take care of a possible least squares residual.
    rho = sqrt(rho.^2 + rhoLS2);
end

% Now compute the first and second derivatives of eta and rho
% with respect to lambda;
deta  =  phi./eta;
drho  = -psi./rho;
ddeta =  dphi./eta - deta.*(deta./eta);
ddrho = -dpsi./rho - drho.*(drho./rho);

% Convert to derivatives of log(eta) and log(rho).
dlogeta  = deta./eta;
dlogrho  = drho./rho;
ddlogeta = ddeta./eta - (dlogeta).^2;
ddlogrho = ddrho./rho - (dlogrho).^2;

% Let g = curvature.
g = - (dlogrho.*ddlogeta - ddlogrho.*dlogeta)./...
    (dlogrho.^2 + dlogeta.^2).^(1.5);

end