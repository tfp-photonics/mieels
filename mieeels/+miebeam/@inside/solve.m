function sol = solve( obj, b, k0 )
%  SOLVE - Solve Mie equation for inside beam.
%
%  Usage for obj = miebeam.inside :
%    sol = solve( obj, b, k0 )
%  Input
%    b        :  impact parameter
%    k0       :  wavenumber of light in vacuum
%  Output
%    sol      :  Mie solution

if isscalar( b )
  %  solve Mie equations
  sol = solve( obj.mie, qinc( obj, b, k0 ) ); 
  %  add impact parameter to solution
  sol.b = b;
else
  sol = arrayfun( @( b ) solve( obj, b, k0 ), b, 'uniform', 1 );
  %  concatenate solutions
  b1 = horzcat( sol.b1 );  c1 = horzcat( sol.c1 );
  b2 = horzcat( sol.b2 );  c2 = horzcat( sol.c2 );
  %  set output
  sol = struct( 'b1', b1, 'c1', c1, 'b2', b2, 'c2', c2, 'b', b, 'k0', k0 );
end
