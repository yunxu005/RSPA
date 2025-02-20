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

e_esh = inline('1.0 + 5.0*c*(p0-1)/(2*p0+3)','c','p0');
e_mt  = inline('1.0 + 5.0*c*(p0-1)/(2.0*(1.0-c)*p0+3.0+2.0*c)','c','p0');
e_sc  = inline('-(5.0*c*(1.0-p0)+2.0*p0-3.0)/6.0+sqrt((5.0*c*(1.0-p0)+2.0*p0-3)^2+24.0*p0)/6.0','c','p0');
e_reu = inline('p0/((1.0-c)*p0+c)','c','p0');
e_vt = inline('(1.0-c)+c*p0','c','p0');

% % Dependence on Young's moduli
% p0 = E1/E0;
% y_esh = e_esh(c,p0)*E0;
% y_mt = e_mt(c,p0)*E0;
% y_sc = e_sc(c,p0)*E0;
% y_reu = e_reu(c,p0)*E0;
% y_vt = e_vt(c,p0)*E0;

% Eb =[y_reu  y_esh  y_mt  y_sc  y_vt];     %  [ y_reu  y_esh   y_mt  y_sc  y_vt]
% for k = 1:5
%     a_c1(k) = E0/(E0+(E1-E0)*k*1.0/5.0);
% %     a_c1(k) = Eb(k)/(Eb(k)+(E1-Eb(k))*(1-c));
% end

% a_c1(1) = E0/(E0+(E1-E0)*(1.0-c));
% a_c1(2) = E0/(E0+(E1-E0)*(1.0-c2));
% a_c1(3) = E0/(E0+(E1-E0)*(1-c2))/((1-c)+c*a_c1(2));
% a_c1(4) =y_sc/(y_sc+(E1-y_sc)*(1.0-c2));
% a_c1(5) = 1.0;
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
        yl(k,i) = E1(k)*a_c1s(i); 
    end
end

for k = 1:6
    Ebb(k) = E0 + c*(E1(k) - E0)*a_c1(k);
end

Nr = size(xr,2);
for  k = 1:6
    a_crs = ones(1,Nr);
    for i = 1:Nr
       yr(k,i) = Ebb(k)*a_crs(i);  
    end
    y0(k) = 0.5*(yl(k,Nl)+yr(k,1));
end

figure;    axis([-1.0 1.0 0.0 7.0]);  hold on;
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

arrow([0.1 5.7],[0.1 6.5],'color','r', 'linewidth',2); hold on;
text(0.18,6.2, ['$\mathcal{C}_{1}$'],'Interpreter','latex', 'FontSize',16, ...
    'Fontweight','bold','color','black'); hold on;

% plot([xl x0 xr],[yl(1,:) y0(1) yr(1,:)],'color',[0.7 0.7 0.9], 'linewidth',2.0,'linestyle','-'); hold on;
% plot([xl x0 xr],[yl(2,:) y0(2) yr(2,:)],'color',[0.5 0.5 0.9], 'linewidth',2.0,'linestyle','-'); hold on;
% plot([xl x0 xr],[yl(3,:) y0(3) yr(3,:)],'color',[0.3 0.3 0.9], 'linewidth',2.0,'linestyle','-'); hold on;
% plot([xl x0 xr],[yl(4,:) y0(4) yr(4,:)],'color',[0.1 0.1 0.9], 'linewidth',2.0,'linestyle','-'); hold on;

plot(0,0.15, 'v', 'markersize', 6,'markerfacecolor','r','markeredgecolor','r',...
    'linestyle','none'); hold on;

% legend({'Reuss/Sachs model','Eshelby method','Mori-Tanaka method', ...
%     'Self-consistent method', 'Voigt/Taylor model','This study'},'FontSize',10);
text(-0.65,0.5, 'Before', 'FontSize',12,'Fontweight','bold','color','black'); hold on;
text(0.3,0.5, 'After', 'FontSize',12,'Fontweight','bold','color','black'); hold on;

text(-0.9, 6.7, '(A)','FontSize',16,'color','black'); hold on;


strx={['$c$']};
str1={['$\bar{E}/E_{0}$']};
str2={['$\bar{\nu}/\nu_{0}$']};
xlabel('Homogenization','Interpreter','latex','FontSize',14,'Fontweight',...
     'bold','color','black'); hold on;
ylabel('Stress','Interpreter','latex','FontSize',14,'Fontweight',...
     'bold','color','black'); hold on;
title(['Tough inclusion ($\mathcal{C}_{1} > \mathcal{C}_{0}$)'],'Interpreter','latex','FontSize',12,'color','black');
set(gca,'XTickLabel',[],'YTickLabel',[],'FontSize',14);







 
 


% axis equal off