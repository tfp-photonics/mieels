function obj = init( obj, mat1, mat2, diameter, varargin )
%  INIT - Initialize Mie solver.
%
%  Usage for obj = multipole.miesolver :
%    obj = init( obj, mat1, mat2, diameter, lmax )
%  Input
%    mat1     :  
%    lmax     :  maximal degree for multipole expansion

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'lmax', 20 );
%  parse input
parse( p, varargin{ : } );

%  save input
[ obj.mat1, obj.mat2, obj.diameter ] = deal( mat1, mat2, diameter );
%  angular degrees and orders
obj.lmax = p.Results.lmax;
[ obj.tab.l, obj.tab.m ] = multipole.sphtable( p.Results.lmax );
