% CCA Canonical correlation analysis between data and targets
%
%    [W1,W2,RHO] = CCA(A,N)
%    [W1,W2,RHO] = CCA(A,FRAC)
%
%  INPUT
%    A           Dataset
%    N  or FRAC  Number of dimensions (>= 1) or minimal correlation (< 1)
%                to retain. Default: N = inf.
%
%  OUTPUT
%    W1          Affine CCA mapping for data.
%    W2          Affine CCA mapping for targets.
%    RHO         Squared canonical correlations.
%
%  DESCRIPTION
%  This routine performs a canonical correlation analysis (CCA).
%
%  SEE ALSO
%  MAPPINGS, DATASETS, PCLDC, KLLDC, PCA, KLM, FISHERM

function [w1,w2,rho] = cca (a,frac)

        prtrace(mfilename);

        truefrac = [];

        % Default: preserve all dimensions (identity mapping).
        if (nargin < 2) | (isempty(frac))
                frac = inf;
                prwarning (3,'no dimensionality given, only ordering dimensions');
        end

       	mapname = 'Canonical Correlation Analysis';

        % Empty mapping: return straightaway.
        if (nargin < 1) | (isempty(a))
                w = mapping(type,frac);
                w = setname(w,mapname);
                return
        end

        a = prdataset(a);   % make sure we have a dataset
        isvaldset(a,1);   % at least 1 object per class

        [m,k,c] = getsize(a); 
        x = getdata(a); y = gettargets(a); 

        % Shift mean of x and y to origin.
        xb = +(x*scalem(prdataset(x))); yb = +(y*scalem(prdataset(y)));

        % Always perform a REDUCM, as (especially) Y may live in a subspace.
        xu = reducm(prdataset(xb)); xb = +(xb*xu); 
				yu = reducm(prdataset(yb)); yb = +(yb*yu);
        korg = k; [m,k] = size(xb);
        nc = min(size(xu,2),size(yu,2));

        % Calculate covariance matrix of the data and the targets,
				% and the cross-covariance matrix.

				Gxx = (xb'*xb)/(m-1); Gyy = (yb'*yb)/(m-1); Gxy = (xb'*yb)/(m-1);

        [Fx,Vx] = eig(inv(Gxx)*Gxy*inv(Gyy)*Gxy');
        [Fy,Vy] = eig(inv(Gyy)*Gxy'*inv(Gxx)*Gxy);

        % v = V(I) contains the sorted eigenvalues:
        [vx,Ix] = sort(-abs(diag(Vx))); [vy,Iy] = sort(-abs(diag(Vy)));
				vx = -vx(1:nc); vy = -vy(1:nc); Ix = Ix(1:nc); Iy = Iy(1:nc);

        if (frac == inf)                        % Return all dimensions, decorrelated and ordered.
                n = nc; rho = vx;
        elseif (frac == 0)                      % Just return cumulative retained variance.
                w1 = cumsum(vx)/sum(vx);
                w2 = cumsum(vy)/sum(vy); 
								rho = [];
                return
        elseif (frac >= 1)                      % Return FRAC dimensions.
                n = abs(frac); if (n > nc), error('illegal dimensionality requested'); end
                Ix = Ix(1:n); Iy = Iy(1:n); sv = sum(vx); rho = vx(1:n);
        elseif (frac > 0)                       % Return the N dimensions that retain at 
																								% least FRAC squared correlation.
                J = find(cumsum(vx)/sum(vx) > frac);
                n = J(1); Ix = Ix(1:n); Iy = Iy(1:n); rho = vx(1:n);
        end

        % If needed, apply pre-calculated projection to (M-1) dimensional subspace.
        if (~isempty(xu))
                xrot = xu.data.rot*Fx(:,Ix); xoff = xu.data.offset*Fx(:,Ix);
                yrot = yu.data.rot*Fy(:,Iy); yoff = yu.data.offset*Fy(:,Iy);
        else
                xrot = Fx(:,Ix);             xoff = -mean(x*Fx(:,Ix));
                yrot = Fy(:,Iy);             yoff = -mean(y*Fy(:,Iy));
        end

        % Construct affine mapping.
        w1 = affine(xrot,xoff,x); w1 = setname(w1,mapname);
        w2 = affine(yrot,yoff,y); w2 = setname(w2,mapname);

return

