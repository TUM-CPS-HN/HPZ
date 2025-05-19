function res = my_and_(cPZ,S)

% compute lower bound
H=S.A;
f=S.b;
c=cPZ.c;
d_max = f-H*c - sum(abs(H*cPZ.G),2);
l = supportFunc_(cPZ,S.A','lower','interval',[]);
l=d_max;
% add additional constraints
A = [S.A * cPZ.G, -0.5*(S.b-l)];
b = 0.5*(S.b+l) - S.A*cPZ.c;
EC = blkdiag(cPZ.E,1);
G = [cPZ.G,zeros(size(cPZ.G,1),size(cPZ.G,2)),zeros(size(cPZ.G,1),1)];
E = [cPZ.E; zeros(1,size(cPZ.E,2))];
E=[E,zeros(size(E,1),size(E,2)),zeros(size(E,1),1)];
% cPZ.id = [cPZ.id;max(cPZ.id) + 1];
A = blkdiag(cPZ.A, A);
b = [cPZ.b; b];
EC = [[cPZ.EC;zeros(1,size(cPZ.EC,2))], EC];
res = conPolyZono(c,G,E,A,b,EC);


end

