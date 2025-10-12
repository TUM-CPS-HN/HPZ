% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the Minkowski sum of two zonotopic sets, Z = X + Y
%   Syntax:
%       Z = plus(X,Y) = X + Y
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object) or n x 1 vector
%       Y - zonotopic set in R^n (hybZono, conZono, or zono object) or n x 1 vector
%   Outputs:
%       Z - zonotopic set in R^n (hybZono, conZono, or zono object)
%   Notes:
%       Overloaded '+' operator
%       If X is a hybZono and Y is a conZono, then Z is a hybZono
%       If X is a conZono and Y is a zono, then Z is a conZono
%       If X and Y are zono, then Z is a zono
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = hpz_plus(obj1,obj2)


Gc = [obj1.Gc obj2.Gc];
Gb = [obj1.Gb obj2.Gb];
c = obj1.c + obj2.c;
Ac = blkdiag(obj1.Ac,obj2.Ac);
Ab = blkdiag(obj1.Ab,obj2.Ab);
b = [obj1.b; obj2.b];
E = blkdiag(obj1.E,obj2.E);
EC = blkdiag(obj1.EC,obj2.EC);
% out = hybZono(Gc,Gb,c,Ac,Ab,b);

out = hybPolyZono(c,Gc,Gb,E,Ac,Ab,b,EC);

end