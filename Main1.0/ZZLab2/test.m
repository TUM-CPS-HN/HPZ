close all;
clear all;
Gc = [3 -3 1; 2 1 -2];          % Continuous generator matrix
Gb = [6 -6 2; 4 2 -4];          % Binary generator matrix
c = [0;0];                      % Center
Ac = [1 1 1];                   % Continuous constraint matrix
Ab = [1 1 1];                   % Binary constraint matrix
b = 1;                          % Constraint offset vector
E = [2 0 1; 1 2 0; 1 1 1];
EC = [1 0 0; 0 2 0; 0 0 1];

EC2 = [1 0 0; 1 2 1; 1 0 1];
HPZ = hybPolyZono(c,Gc,Gb,E,Ac,Ab,b,EC,[],[1;2;3]);
HPZ2 = hybPolyZono(c,2*Gc,Gb,E,2*Ac,Ab,0.1*b,EC2,[],[1;5;6]);
% HPZ = hybPolyZono(c,Gc,Gb,E,Ac,Ab,b,EC,[]);
% HPZ2 = hybPolyZono(c,2*Gc,Gb,E,2*Ac,Ab,0.1*b,EC2,[]);
% figure;
% plotHPZ(HPZ);
% res=plotAsConPolyZono(HPZ); 
% hold on;

% HPZ3 = hpz_plus(HPZ,HPZ);
HPZ3 = Exact_hpz_plus(HPZ,HPZ2);


   % c = [0;0];
   % G = [2 0 1;0 2 1];
   % E = [1 0 3;0 1 1];
   % A = [1 -1];
   % b = 0;
   % EC = [2 0; 0 1];
   % cPZ1 = conPolyZono(c,G,E,A,b,EC);
   % 
   % M = 0.1 * [3 1;2 4];
   % cPZ2 = M * cPZ1;
   % 
   % res = exactPlus(cPZ1,cPZ2);



% Z = hybZono(Gc,Gb,c,Ac,Ab,b);
% figure;
% plot(Z, 'b'); hold on;


% CPZ = conPolyZono(c,Gc,E,Ac,b,EC);
% figure;
% plot(CPZ); hold on;
