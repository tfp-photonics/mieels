function data = work( obj, sol )
%  WORK - Work performed by electron against electric field.
%
%  Usage for obj = miebeam.inside :
%    data = work( obj, sol )
%  Input
%    sol        :  Mie solution 
%  Output
%    data.work  :  work integrand
%    data.z     :  propagation distances
%    data.w     :  integration weights

%  radius of sphere and impact parameter
a = 0.5 * obj.mie.diameter;
b = sol.b;
%  crossing point
z0 = sqrt( a ^ 2 - b ^ 2 );

%  integration positions and weights for outside region
[ z1, w1 ] = misc.lgwt( obj.nz, z0, z0 + a );
[ z2, w2 ] = misc.lgwt( obj.nz, z0 + a, obj.zlim * a );
%  integration positions and weights for inside region
[ z3, w3 ] = misc.lgwt( obj.nz, -z0, z0 );
%  concatenate positions and weights
z = vertcat( -z1, -z2, z1, z2, z3 );
w = vertcat(  w1,  w2, w1, w2, w3 );
%  sort positions and weights
[ data.z, ind ] = sort( z );  data.w = w( ind );

%  electromagnetic fields
pos = [ b, 0, 0 ] + data.z * [ 0, 0, 1 ];
e = fields( obj.mie, sol, pos );
%  work integrand
fac = exp( - 1i * sol.k0 / obj.vel * data.z );
data.work = real( fac .* e( :, 3 ) );
