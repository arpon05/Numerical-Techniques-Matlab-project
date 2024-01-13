clc 
close all
clear

%Given that
R=0.5;     %ohm
C=0.001;   %F
L=0.001;   %H

n=1000;
L1=0.01-0; %length
delx=L1/(n+1);
%Boundary condition
v0=0;
v_n=-0.001816;
%For coefficient values.
a=(2*R*C*L)-(delx*L);          
b=(2*R*(delx^2))-(4*R*C*L);  
c=(2*R*C*L)+(delx*L);           
d=0;
%For tridiagonal matrix
A = diag(b*ones(1,n)) + diag(c*ones(1,n-1),1) + diag(a*ones(1,n-1),-1);
%For n x 1 matrix
m=zeros(n-2,1);
m(1:n-2)=d;

B=[d-a*v0
   m(1:n-2)
   d-c*v_n];

v=inv(A)*B   %V
%Boundary condition
t=linspace(0,0.01,n+2) 
%overall matrix of 'v'
v_all=[v0
       v
       v_n]

figure(1)
plot(t,v_all,'linewidth',2);hold on
set(gca,'fontsize',12);
xlabel('t')
ylabel('v and vexact')
grid on

%For ic vs t
v2=v_all;
t2=t;

for k=2:n+1
    
    ic(k)= C*((v2(k+1)-v2(k-1))/ (2*(t(k+1)-t(k))))  %Using central difference
    
end

figure(2)
ic(1)=NaN;
ic(n+2)= NaN;
plot(t,ic,'linewidth',2)
set(gca,'fontsize',12);
xlabel('t')
ylabel('ic')
grid on

%For PR vs t

for j=1:n+2
    
    PR(j)=((v2(j)^2)/R)    %Instantaneous power PR=(V^2)/R
    
end
figure(3)
plot(t,PR,'linewidth',2);hold on
set(gca,'fontsize',12);
xlabel('t')
ylabel('PR and pe')
grid on

%For Pc vs t

for u=1:n+2
    
    Pc(u)=(ic(u))*v2(u)   %Instantaneous power absorbed by capacitor 
    
end

figure(4)
plot(t,Pc,'linewidth',2)
set(gca,'fontsize',12);
xlabel('t')
ylabel('pc')
grid on

%Calculating exact value
alp=1/(2*R*C);
s=-alp;
Ae=-4/C;
%Voltage exact value equation
ve=Ae.*t.*exp(s.*t);
q=zeros(n+2,1);
q(1:n+2)=ve;           %converting ve array (1*n) to (n*1)
%For exact power
pe=(ve.^2)./R;
%Calculating Error
ei=(q-v_all);
%The amount of error En
En=sqrt(sum(ei.^2)./(n+2));
%V exact plot
figure(1)
plot(t,q,'--','linewidth',2)
%P exact plot
figure(3)
plot(t,pe,'--','linewidth',2)

%Find and plot En vs n
n1 = [75:25:1000 2000 3500 5000];     

for i=1:n1
    
    En1=sqrt(sum(ei.^2)./(n1+2));
    
end
figure(5)
plot(n1,En1,'o-','linewidth',2)
xlabel('n')
ylabel('En')
set(gca,'fontsize',12);
grid on
