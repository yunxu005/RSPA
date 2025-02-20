clear         %  Standard code
%    0.057  3.72  200.8;  0.10  4.7  320.6;
%    0.15  5.5  439.1;  0.20   6.4  594.5;
E0 = 55;       %  Young's modulus of matrix
E1 = 68.9;      %  Young's modulus of particle  5.0
v0 = 0.20;    %  Poisson's ratio of matrix
v1 = 0.30;    %  Poisson's ratio of particle  0.20
G0 = 26.2*1.5;     %  Toughness of matrix  1.2
G1 = 3265.6;    %  Toughness of particle    15
c = 0.20;

p0 = E1/E0;
q0 = v1/v0;
p1 = E0/E1;
q1 = v0/v1;
xi = G1/G0;      %  Contrast of Toughness

e0_esh = inline('1.0 + 5.0*c*(p0-1)/(2*p0+3)','c','p0');
e0_mt  = inline('1.0 + 5.0*c*(p0-1)/(2.0*(1.0-c)*p0+3.0+2.0*c)','c','p0');
e0_sc  = inline('-(5.0*c*(1.0-p0)+2.0*p0-3.0)/6.0+sqrt((5.0*c*(1.0-p0)+2.0*p0-3)^2+24.0*p0)/6.0','c','p0');
e0_reu = inline('p0/((1.0-c)*p0+c)','c','p0');

v0_esh = inline('1.0 + 5.0*c*(q0-1)/(2*p0+3)','c','p0','q0');
v0_mt  = inline('1.0 + 5.0*c*(q0-1)/(2.0*(1-c)*p0+3.0+2.0*c)','c','p0','q0');
v0_sc  = inline('1.0 + 5.0*c*(q0-1)*e_sc/(2.0*p0+3.0*e_sc)','c','p0','q0','e_sc');
v0_reu = inline('1.0 + c*(q0-1.0)/((1.0-c)*p0+c)','c','p0','q0');

vb0_esh = inline('(1.0+v0_esh*v0)/(1.0+v0)*(1.0-(1.0+v0_esh)*v0)/(1.0-2*v0)',...
    'v0_esh','v0');
vb0_mt  = inline('(1.0+v0_mt*v0)/(1.0+v0)*(1.0-(1.0+v0_mt)*v0)/(1.0-2*v0)',...
    'v0_mt','v0');
vb0_sc  = inline('(1.0+v0_sc*v0)/(1.0+v0)*(1.0-(1.0+v0_sc)*v0)/(1.0-2*v0)','v0_sc','v0');
vb0_reu = inline('(1.0+v0_reu*v0)/(1.0+v0)*(1.0-(1.0+v0_reu)*v0)/(1.0-2*v0)',...
    'v0_reu','v0');

phi0_esh = inline('e0_esh*vb0_esh','e0_esh','vb0_esh');
phi0_mt  = inline('e0_mt*vb0_mt','e0_mt','vb0_mt');
phi0_sc  = inline('e0_sc*vb0_sc','e0_sc','vb0_sc');
phi0_reu = inline('e0_reu*vb0_reu','e0_reu','vb0_reu');


e1_esh = inline('p1 + 5.0*c*(1-p1)/(2*p0+3)','c','p0','p1');
e1_mt  = inline('p1 + 5.0*c*(1-p1)/(2.0*(1.0-c)*p0+3.0+2.0*c)','c','p0','p1');
e1_sc  = inline('(5.0*c*(1.0-p1)+3.0*p1-2.0)/6.0+sqrt((5.0*c*(1.0-p1)+3.0*p1-2.0)^2+24.0*p1)/6.0','c','p1');
e1_reu = inline('p1 + c*(1.0-p1)/((1.0-c)*p0+c)','c','p0','p1');

v1_esh = inline('q1 + 5.0*c*(1-q1)/(2*p0+3)','c','p0','q1');
v1_mt  = inline('q1 + 5.0*c*(1-q1)/(2.0*(1-c)*p0+3.0+2.0*c)','c','p0','q1');
v1_sc  = inline('q1 + 5.0*c*(1-q1)*e1_sc/(2.0+3.0*e1_sc)','c','q1','e1_sc');
v1_reu = inline('q1 + c*(1.0-q1)/((1.0-c)*p0+c)','c','p0','q1');

vb1_esh = inline('(1.0+v1_esh*v1)/(1.0+v1)*(1.0-(1.0+v1_esh)*v1)/(1.0-2*v1)',...
    'v1_esh','v1');
vb1_mt  = inline('(1.0+v1_mt*v1)/(1.0+v1)*(1.0-(1.0+v1_mt)*v1)/(1.0-2*v1)',...
    'v1_mt','v1');
vb1_sc  = inline('(1.0+v1_sc*v1)/(1.0+v1)*(1.0-(1.0+v1_sc)*v1)/(1.0-2*v1)','v1_sc','v1');
vb1_reu  = inline('(1.0+v1_reu*v1)/(1.0+v1)*(1.0-(1.0+v1_reu)*v1)/(1.0-2*v1)','v1_reu','v1');

phi1_esh = inline('e1_esh*vb1_esh','e1_esh','vb1_esh');
phi1_mt  = inline('e1_mt*vb1_mt','e1_mt','vb1_mt');
phi1_sc  = inline('e1_sc*vb1_sc','e1_sc','vb1_sc');
phi1_reu  = inline('e1_reu*vb1_reu','e1_reu','vb1_reu');

Ps0_reu = phi0_reu(e0_reu(c,p0),vb0_reu(v0_reu(c,p0,q0),v0));
Ps0_mt  = phi0_mt(e0_mt(c,p0),vb0_mt(v0_mt(c,p0,q0),v0));
Ps0_sc  = phi0_sc(e0_sc(c,p0),vb0_sc(v0_sc(c,p0,q0,e0_sc(c,p0)),v0));
    
Ps1_reu = phi1_reu(e1_reu(c,p0,p1),vb1_reu(v1_reu(c,p0,q1),v1));
Ps1_mt  = phi1_mt(e1_mt(c,p0,p1),vb1_mt(v1_mt(c,p0,q1),v1));
Ps1_sc  = phi1_sc(e1_sc(c,p1),vb1_sc(v1_sc(c,q1,e1_sc(c,p1)),v1));
    
Gb_reu = (1.0-c)*Ps0_reu*G0 + c*Ps1_reu*G1;
Gb_mt  = (1.0-c)*Ps0_mt*G0 + c*Ps1_mt*G1;
Gb_sc  = (1.0-c)*Ps0_sc*G0 + c*Ps1_sc*G1;
Gb_avg = (Gb_reu+Gb_mt+Gb_sc)/3.0;
e0_avg = (e0_reu(c,p0)+e0_mt(c,p0)+e0_sc(c,p0))/3.0;

[Gb_reu  Gb_mt  Gb_sc  Gb_avg]
[e0_reu(c,p0) e0_mt(c,p0) e0_sc(c,p0) e0_avg]*E0







