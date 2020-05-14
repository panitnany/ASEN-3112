%% Equivalent beam model
%Part1: Moment of Inertia

%Define constants (outer and inner radii, length of bar, elastic modulus)
r=3/16*0.0254; %in to m
r_i=1/8*0.0254; %in to m
Lb=0.25; %m
E=70*10^9; %Pa

%Using geometry the lever arm distance, cross sectional area of a single
%member, and area moment of inertia are calculated. 

d=sqrt(2)*Lb;
A=pi*(r^2-r_i^2);
Ic=pi/32*(r^4-r_i^4);
I=4*(Ic+A*d^2);

%Finding deflection, start by defining midspan load, total truss length,
%and support reactions.
Load=222.4; %N
L=16*Lb; %m
Ra=Load/2;

%Find the moment as a function of x (piecewise), then integrate twice to
%find deflection (Second Order Method). Use BC's: M(0)=M(L)=0, deflection(0)=deflection(L)=0, so
%no constants are needed. 
syms x
M1=Ra*x;
M2=Ra*x-Load*(x-L/2);
Moment= piecewise( 0 < x < L/2, M1, L>x>=L/2, M2);

% Ode45 from x=[0:L] and with initial deflection v(0)=0 and initial slope
% v'(0)=-L^3/16
[x,v]=ode45(@fun,[0 L], [0; (-L^3)/16]);
deflection = 1000/(E*I)*v(:,1);
plot(x,deflection)
title('Deflection');
xlabel('X (m)');
ylabel('Deflection (mm)');
midpointdeflection=subs(deflection,x,L/2);

function dvdx=fun(x,v)
    if x<=2
        dvdx=[v(2); 111.2*x] %Moment
    else 
        dvdx=[v(2); 111.2*x-222.4*(x-2)]
    end
end