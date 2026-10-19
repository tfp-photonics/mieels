function sol = solve( obj, q )
%  SOLVE - Solve Mie equations.
%
%  Usage for obj = miebeam.miesolver :
%    sol = solve( obj, q )
%  Input
%    q    :  inhomogeneity of incoming fields
%  Output
%    sol  :  solution vector

%  transition matrix
t = transition( obj, q.k0 );
%  solve Mie equations
n = size( t, 3 );
Q = reshape( cat( 1, q.e, q.h ), 4, 1, n, [] );
vec = squeeze( pagemtimes( t, Q ) );

%  set solution vector
sol.b1 = reshape( vec( 1, :, : ), n, [] );
sol.c1 = reshape( vec( 2, :, : ), n, [] );
sol.b2 = reshape( vec( 3, :, : ), n, [] );
sol.c2 = reshape( vec( 4, :, : ), n, [] );    
%  wavenumber of light in vacuum
sol.k0 = q.k0;
