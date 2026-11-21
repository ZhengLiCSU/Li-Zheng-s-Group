function [F,DH,DK,diagnostic] = soil_cvfe_flux(h,K,G,conductivityLocation)
%SOIL_CVFE_FLUX Physical Darcy flow through median-dual triangular subfaces.
% [F,DH,DK,diagnostic]=soil_cvfe_flux(h,K,G) takes pressure head h and
% conductivity K at the ns LOCAL soil nodes and averages K per triangle.
% soil_cvfe_flux(h,Ke,G,'element') instead takes element conductivities,
% e.g. obtained by positive quadrature of a nonlinear constitutive law.
%
% DH is the sparse partial derivative dF/dh with conductivity held fixed.
% DK is dF/dKnodal in default 'nodal' mode, or dF/dKe in 'element' mode.
% Constitutive chain terms belong to the caller: DH+DK*dK/dh, DK*dK/dT.
% No flux clipping, conductivity floor or algebraic cotangent upwinding is
% used. F>0 denotes physical flow from G.from to G.to. The three subfaces
% per triangle must remain separate when constructing advective heat flux.
% H=h-G.z follows positive-down depth. K is scalar isotropic conductivity.

if nargin<4 || isempty(conductivityLocation)
    conductivityLocation='nodal';
end
location=validatestring(conductivityLocation,{'nodal','element'},mfilename, ...
    'conductivityLocation',4);
validateattributes(h,{'numeric'},{'vector','real','finite','numel',G.ns}, ...
    mfilename,'h',1);
validateattributes(K,{'numeric'},{'vector','real','finite','nonnegative'}, ...
    mfilename,'K',2);
h=h(:); K=K(:);
switch location
    case 'nodal'
        if numel(K)~=G.ns
            error('soil_cvfe_flux:ConductivitySize','Nodal K must have ns entries.');
        end
        Ke=G.meanK*K;
    case 'element'
        if numel(K)~=G.ne
            error('soil_cvfe_flux:ConductivitySize','Element Ke must have ne entries.');
        end
        Ke=K;
end
H=h-G.z; unitFlow=G.D0*H; faceK=Ke(G.elementIndex);
F=faceK.*unitFlow;
if nargout>=2
    DH=spdiags(faceK,0,G.nf,G.nf)*G.D0;
end
if nargout>=3
    DK=spdiags(unitFlow,0,G.nf,G.nf)*G.elementSelect;
    if strcmp(location,'nodal')
        DK=DK*G.meanK;
    end
end
if nargout>=4
    % MATLAB vector-by-vector indexing keeps the source orientation when
    % ne=1; reshape explicitly so a single triangle is still 1-by-3.
    elementH=reshape(H(G.localTri(:)),G.ne,3);
    gradH=[sum(G.gradX.*elementH,2),sum(G.gradY.*elementH,2)];
    diagnostic=struct('head',H,'elementK',Ke,'gradientHead',gradH, ...
        'darcyElement',-Ke.*gradH,'unitFlow',unitFlow, ...
        'conductivityLocation',location);
end
end
