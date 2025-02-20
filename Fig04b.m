% Stress
% Using phases elastic left:  C_1*eps1,  right: \bar{C}*eps0
%   eps1 = a_c1*eps0,  matrix modulus C is used to compute a_c1

clear
E0 = 1.0;       %  Young's modulus of matrix  [0.8  2.0]
v0 = 0.25;    %  Poisson's ratio of matrix
v1 = 0.20;    %  Poisson's ratio of particle
c = 0.8;
c2 = 0.6;
% E1 = 5.0;      %  Young's modulus of particle   [5.3  10.0]
E1 =  [1.0 3.0 5.0 7.0 9.0 11.0];


al  = 0.9;
for i = 1:6
    Et(i) = al*E0 + (1.0-al)*E1(i);
    a_c1(i) = Et(i)/(Et(i)+(E1(i)-Et(i))*(1.0-c));
end

% legend({'Reuss/Sachs model','Eshelby method','Mori-Tanaka method', ...
%     'Self-consistent method', 'Voigt/Taylor model'},'FontSize',10);

h = 0.02;
xl = [-1.0:h:0.0-h];
x0 = 0;
xr = [h:h:1.0];
Nl = size(xl,2);
for k = 1:6
    for i = 1:Nl
        a_c1s = a_c1(k).*ones(1,Nl);
        yl(k,i) = 0.5*E1(k)*a_c1s(i)*a_c1s(i); 
    end
end

for k = 1:6
    Ebb(k) = E0 + c*(E1(k) - E0)*a_c1(k);
end

Nr = size(xr,2);
for  k = 1:6
    a_crs = ones(1,Nr);
    for i = 1:Nr
       yr(k,i) = 0.5*Ebb(k)*a_crs(i)*a_crs(i);  
    end
    y0(k) = 0.5*(yl(k,Nl)+yr(k,1));
end

figure;    axis([-1.0 1.0 0.0 3.0]);  hold on;
% grid on; 
% yyaxis left

plot([xl(1:7:Nl) x0 xr(1:7:Nr)],[yl(1,1:7:Nl) y0(1) yr(1,1:7:Nr)], ...
    'o', 'markersize', 4,'markerfacecolor',[0.6 0.6 0.9],'markeredgecolor',[0.6 0.6 0.9], ...
    'color',[0.6 0.6 0.9], 'linewidth',2.0,'linestyle','-'); hold on;
plot([xl(1:7:Nl) x0 xr(1:7:Nr)],[yl(2,1:7:Nl) y0(2) yr(2,1:7:Nr)], ...
    's', 'markersize', 5,'markerfacecolor',[0.5 0.5 0.9],'markeredgecolor',[0.5 0.5 0.9], ...
    'color',[0.5 0.5 0.9], 'linewidth',2.0,'linestyle','-'); hold on;
plot([xl(1:7:Nl) x0 xr(1:7:Nr)],[yl(3,1:7:Nl) y0(3) yr(3,1:7:Nr)], ...
    'd', 'markersize', 5,'markerfacecolor',[0.4 0.4 0.9],'markeredgecolor',[0.4 0.4 0.9], ...
    'color',[0.4 0.4 0.9], 'linewidth',2.0,'linestyle','-'); hold on;
plot([xl(1:7:Nl) x0 xr(1:7:Nr)],[yl(4,1:7:Nl) y0(4) yr(4,1:7:Nr)], ...
    '>', 'markersize', 4,'markerfacecolor',[0.3 0.3 0.9],'markeredgecolor',[0.3 0.3 0.9], ...
    'color',[0.3 0.3 0.9], 'linewidth',2.0,'linestyle','-'); hold on;
plot([xl(1:7:Nl) x0 xr(1:7:Nr)],[yl(5,1:7:Nl) y0(5) yr(5,1:7:Nr)], ...
    '^', 'markersize', 4,'markerfacecolor',[0.2 0.2 0.9],'markeredgecolor',[0.2 0.2 0.9], ...
    'color',[0.2 0.2 0.9], 'linewidth',2.0,'linestyle','-'); hold on;
plot([xl(1:7:Nl) x0 xr(1:7:Nr)],[yl(6,1:7:Nl) y0(6) yr(6,1:7:Nr)], ...
    '+', 'markersize', 5,'markerfacecolor',[0.1 0.1 0.9],'markeredgecolor',[0.1 0.1 0.9], ...
    'color',[0.1 0.1 0.9], 'linewidth',2.0,'linestyle','-'); hold on;

arrow([-0.10 2.3],[-0.10 2.6],'color','r', 'linewidth',2); hold on;
text(-0.25,2.5, ['$\mathcal{C}_{1}$'],'Interpreter','latex', 'FontSize',16, ...
    'Fontweight','bold','color','black'); hold on;


plot(0,0.05, 'v', 'markersize', 6,'markerfacecolor','r','markeredgecolor','r',...
    'linestyle','none'); hold on;

% legend({'Reuss/Sachs model','Eshelby method','Mori-Tanaka method', ...
%     'Self-consistent method', 'Voigt/Taylor model','This study'},'FontSize',10);
text(-0.65,0.2, 'Before', 'FontSize',12,'Fontweight','bold','color','black'); hold on;
text(0.3,0.2, 'After', 'FontSize',12,'Fontweight','bold','color','black'); hold on;

text(-0.9, 2.85, '(B)','FontSize',16,'color','black'); hold on;


strx={['$c$']};
str1={['$\bar{E}/E_{0}$']};
str2={['$\bar{\nu}/\nu_{0}$']};
xlabel('Homogenization','Interpreter','latex','FontSize',14,'Fontweight',...
     'bold','color','black'); hold on;
ylabel('Surface energy','Interpreter','latex','FontSize',14,'Fontweight',...
     'bold','color','black'); hold on;
title(['Tough inclusion ($\mathcal{C}_{1} > \mathcal{C}_{0}$)'],'Interpreter','latex','FontSize',12,'color','black');
set(gca,'XTickLabel',[],'YTickLabel',[],'FontSize',14);







 
 


% axis equal off