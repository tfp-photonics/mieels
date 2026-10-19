function sol = select( ~, sol, ind )
%  SELECT - Select specific elements.
%
%  Usage for obj = miebeam.miesolver :
%    sol = select( obj, sol, ind )
%  Input
%    sol    :  Mie solution
%    ind    :  index to selcted elements
%  Output
%    sol    :  solution to selected elements

for name = [ "b1", "c1", "b2", "c2" ]
  sol.( name ) = sol.( name )( :, ind );
end
