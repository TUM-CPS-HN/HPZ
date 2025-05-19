
c = [0;0];
G = [1 0 -0; 0 1 -1];
E = [1 0 0; 0 1 1; 0 0 1];

A = [1 1 -1];
b = 1;
EC = [0 1 2; 1 0 0; 0 1 0];

cPZ = conPolyZono(c, G,E, A, b,EC);


P1 = polytope([0 1],1);
res1 = my_and_(cPZ,P1);

figure; hold on;
xlim([-5,5]); ylim([-5,5]);
plot(P1,[1,2],'r');
num_split=16;
plot(res1,[1,2],'b','LineWidth',2,'Splits',num_split);
plot(cPZ,[1,2],'g','LineWidth',2,'Splits',num_split);


