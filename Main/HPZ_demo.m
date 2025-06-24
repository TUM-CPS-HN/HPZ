close all;
clear all;
figure;
hold on;


% c = [0;0];
% G = [1 0 1.5; 0 1 2];
% E = [1 0 0; 0 1 0; 0 0 1];
% A = [1 2 0.5];
% b = 0;
% EC = [1 0 0; 0 1 0; 0 0 3];
% 
% cPZ = conPolyZono(c,G,E,A,b,EC);
% plot(cPZ,[1,2],'FaceColor','r','Splits',8);

c1 = [0;0];
G1 = [1 0 1.5 0.5; 0 1 2 -2];
E1 = [1 0 0 0; 0 1 0 2; 0 0 1 1];
A1 = [1 2 0.5];
b1 = 0;
EC1 = [1 0 0; 0 1 0; 0 0 3];


cPZ1 = conPolyZono(c1,G1,E1,A1,b1,EC1);
handleX0=plot(cPZ1,[1,2],'Linestyle','--','Splits',8);



% figure
% c=c1;
% Gc=G1;
% Gb=[3 1;1 3];
% E=E1;
% Ac=A1;
% Ab=[0 0 ];
% b=b1;
% EC=EC1;
% HPZ = hybPolyZono(c,Gc,Gb,E,Ac,Ab,b,EC);
% 
% hold on;
% plotHPZ(HPZ);

c=c1;
Gc=G1;
Gb=[3 3 2;1 3 4];
Gb=[3 1 4;1 3 -4];
E=E1;
Ac=A1;
Ab=[0 0 0]; %8
% Ab=[0 0 3.5]; %8 different but same
% Ab=3.5/3*[1 1 1];%8 different 
Ab=1.5*[1 1 1];% 6
b=b1;
EC=EC1;
HPZ = hybPolyZono(c,Gc,Gb,E,Ac,Ab,b,EC);

hold on;

Ab=[0 0 0]; %8
plotHPZ(HPZ);

HPZ = hybPolyZono(c,Gc,Gb,E,Ac,Ab,b,EC);

hold on;
plotHPZ(HPZ);


HZ=hybZono(Gc,Gb,c,[],[],[]);
plot(HZ,plotOptions('FaceColor', 'none', ...      % 填充色    
     'EdgeColor', 'g', ...     
     'EdgeAlpha', 1, ...
     'LineWidth', 1));


xlabel('$z_1$','Interpreter', 'latex');
ylabel('$z_2$','Interpreter', 'latex');



% handleX0=plot([100;100]+cPZ1,[1,2],'Linestyle','--','Splits',8,'Color','r');
% handleX1=plot([100;100]+cPZ1,[1,2],'Splits',8,'FaceColor','r', 'EdgeColor','k', 'FaceAlpha',0.3);
% handleX2=plot([100;100]+cPZ1,[1,2],'Splits',8,'FaceColor','b', 'EdgeColor','k', 'FaceAlpha',0.3);
% legend([handleX0,handleX1,handleX2],{'${\mathcal{CPZ}}$','${\mathcal{HPZ}_1}$','${\mathcal{HPZ}_2}$'}, ...
% 'Interpreter', 'latex', 'Location', 'best');
% legend(handleX1,'${\mathcal{HPZ}_1}$', ...
% 'Interpreter', 'latex', 'Location', 'best');
% legend(handleX2,'${\mathcal{HPZ}_2}$', ...
% 'Interpreter', 'latex', 'Location', 'best');
