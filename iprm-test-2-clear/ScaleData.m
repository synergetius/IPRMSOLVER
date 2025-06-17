function [cscale,rscale] = ScaleData(A,iprint,scltol)



  if iprint > 0
    fprintf('\ngmscale: Geometric-Mean scaling of matrix')
    fprintf('\n-------\n                 Max col ratio')
  end

  [m,n]   = size(A);
  A       = abs(A);    % Work with |Aij|
  maxpass = 10;
  aratio  = 1e+50;
  damp    = 1e-4;
  small   = 1e-8;
  rscale  = ones(m,1);
  cscale  = ones(n,1);

%---------------------------------------------------------------
% Main loop.
%---------------------------------------------------------------
  for npass = 0:maxpass

    % Find the largest column ratio.
    % Also set new column scales (except on pass 0).

    rscale(rscale==0) = 1;
    Rinv    = diag(sparse(1./rscale));
    SA      = Rinv*A;
    [I,J,V] = find(SA);
    invSA   = sparse(I,J,1./V,m,n);
    cmax    = full(max(   SA))';   % column vector
    cmin    = full(max(invSA))';   % column vector
    cmin    = 1./(cmin + eps);
    sratio  = max( cmax./cmin );   % Max col ratio
    if npass > 0
      cscale = sqrt( max(cmin, damp*cmax) .* cmax );
    end

    if iprint > 0
      fprintf('\n  After %2g %19.2f', npass, sratio)
    end

    if npass >= 2 && sratio >= aratio*scltol, break; end
    if npass == maxpass, break; end
    aratio  = sratio;

    % Set new row scales for the next pass.

    cscale(cscale==0) = 1;
    Cinv    = diag(sparse(1./cscale));
    %disp("A");
    %disp(size(A));
    SA      = A*Cinv;                  % Scaled A
    [I,J,V] = find(SA);
    invSA   = sparse(I,J,1./V,m,n);
    rmax    = full(max(   SA,[],2));   % column vector
    rmin    = full(max(invSA,[],2));   % column vector
    rmin    = 1./(rmin + eps);
    rscale  = sqrt( max(rmin, damp*rmax) .* rmax );
  end
%---------------------------------------------------------------
% End of main loop.
%---------------------------------------------------------------

% Reset column scales so the biggest element
% in each scaled column will be 1.
% Again, allow for empty rows and columns.

  rscale(rscale==0) = 1;
  Rinv    = diag(sparse(1./rscale));
  SA      = Rinv*A;
  [I,J,V] = find(SA);
  cscale  = full(max(SA))';   % column vector
  cscale(cscale==0) = 1;

% Find the min and max scales.

  if iprint>0
    [rmin,imin] = min(rscale);
    [rmax,imax] = max(rscale);
    [cmin,jmin] = min(cscale);
    [cmax,jmax] = max(cscale);

    fprintf('\n\n  Min scale               Max scale')
    fprintf('\n  Row %6g %9.1e    Row %6g %9.1e', imin, rmin, imax, rmax)
    fprintf('\n  Col %6g %9.1e    Col %6g %9.1e', jmin, cmin, jmax, cmax)
  end

% end of gmscale
