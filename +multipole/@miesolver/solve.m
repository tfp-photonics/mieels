function sol = solve( obj, q, varargin )
%  SOLVE - Solve Mie equation for incoming multipole coefficients.
%
%  Usage for obj = multipole.miesolver :
%    sol = solve( obj, q, PropertyPairs )
%  Input
%    q            :  incoming multipole coefficients
%  PropertyName
%    excitation   :  'out' or 'in' for excitation from outside or inside
%  Output
%    sol          :  multipole solution at sphere outside

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'excitation', 'out' );
%  parse input
parse( p, varargin{ : } );

%  Mie coefficients
k0 = q.k0;
[ a, b ] = miecoefficients( obj, k0, obj.tab.l, varargin{ : } );
%  scattered Mie coefficients for outside fields
a = - a .* q.a;  
b = - b .* q.b;
%  set output
switch p.Results.excitation
  case 'out'
    sol = multipole.solution( obj.lmax, obj.embedding, k0, a, b, q.a, q.b );
  otherwise
    sol = multipole.solution( obj.lmax, obj.embedding, k0, a, b );
end
