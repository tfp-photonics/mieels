function P = scattering( obj, sol )
%  SCATTERING - Scattered power for Mie solution.
%
%  Usage for obj = miebeam.miesolver :
%    P = scattering( obj, sol )
%  Input
%    sol  :  solution vector
%  Output
%    P    :  scattered power

%  material properties
mat = obj.mat2;
[ k, Z ] = deal( mat.k( sol.k0 ), mat.Z( sol.k0 ) );
%  scattered power
P = 0.5 / ( Z * k ^ 2 ) * sum( abs( sol.b2 ) .^ 2 + abs( sol.c2 ) .^ 2 ); 
