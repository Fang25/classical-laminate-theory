
%%
% Line 283 - invert local stress plot

%%
clear all;

%%  input 
% force and moment
F=[500 0 0  0  0 0    ]  % [Nx Ny Nxy Mx My Mxy] (N/mm & Nmm/mm)

E1=    121                %Young's moduli-1[GPa]
E2=  8.6                %Ypung's moduli-2[Gpa]
G12=  4.7               %Shear modules   [Gpa]
miu12=    0.27             %Poisson¡¯s ratio
t=    5/4;                 %thickness       [mm]
Dir=[ 0 45 -45 90    ];%directions      [degrees]
fail=[1.04 0.57 0.035 0.114 0.072 0.021 0.011 0.002 0.0064 0.038]; 
%  [thegam[GP]1t 1c 2t 2c tao12 eps1t 1c 2t 2c gamma12]

Nx=F(1,1);
Ny=F(1,2);
Nxy=F(1,3);
Mx=F(1,4);
My=F(1,5);
Mxy=F(1,6);

Xt=fail(1,1);
Xc=fail(1,2);
Yt=fail(1,3);
Yc=fail(1,4);
St=fail(1,5);
Xet=fail(1,6);
Xec=fail(1,7);
Yet=fail(1,8);
Yec=fail(1,9);
Gt=fail(1,10);
F1=1/Xt-1/Xc;
F11=1/(Xt*Xc);
F12=-0.5;
F2=1/Yt-1/Yc;
F22=1/(Yt*Yc);
F6=0;
F66=1/(St^2);
%%
n=size(Dir,2);           %the number of the layers
T(1,1:n)=t;              
h=sum(T);                %total thickness of the laminate
mech=[E1,E2,G12,miu12];   %mechanical properities
%%
A=zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
z0=-h/2;
Q0=layer(mech,0);
for i=1:n  % Qi [GPa]
    Qi{1,i}=layer(mech,Dir(1,i));
   z(1,i)=z0+sum(T(1,1:i)); 
   Ti{1,i}=Tr(Dir(1,i));
end  
for i=1:n;
   if i==1;
       A=Qi{1,i}*(z(1,i)-z0)*1000;
       B=Qi{1,i}*1/2*(z(1,i)^2-z0^2)*1000; 
       D=Qi{1,i}*1/3*(z(1,i)^3-z0^3)*1000;
   else
       A=A+Qi{1,i}*(z(1,i)-z(1,i-1))*1000;
       B=B+Qi{1,i}*1/2*(z(1,i)^2-z(1,i-1)^2)*1000;
       D=D+Qi{1,i}*1/3*(z(1,i)^3-z(1,i-1)^3)*1000;
   end
end  
ABD=[A B;B D];% A[N/mm] /B[N]  D[Nmm]
Exlam=1/h*(A(1,1)-(A(1,2)^2/A(2,2)));
Eylam=1/h*(A(2,2)-(A(1,2)^2/A(1,1)));
Gxylam=A(3,3)/h; 
miuxylam=A(1,2)/A(1,1);
miuyxlam=miuxylam*Eylam/Exlam;
%% deformation 
zk1=[z0,z];
for i=1:n;
   zk(1,i)=(zk1(1,i)+zk1(1,i+1))/2; 
end

def=(inv(ABD)*F')'; 
ex=def(1,1);
ey=def(1,2);
ez=def(1,3);
kx=def(1,4);
ky=def(1,5);
kz=def(1,6);
def1=[ex;ey;ez];

for i=1:n
   defi{1,i}=def1+zk(1,i)*[kx;ky;kz] ;
end
for i=1:n+1;
   defi2{1,i}=def1+zk1(1,i)*[kx;ky;kz];
end
strain=cell2mat(defi2);
T0=Tr(0);
thegam0=Q0*T0'*def1;
% Tim=cell2mat(Ti);
% Qim=cell2mat(Qi);
% thegam1= cell2mat(Ti(1,1))*cell2mat(Qi(1,1))*(strain(1,1)+(strain(1,1)-strain(1,2))/(zk1(1,1)-zk1(1,2))*(-0.20-zk1(1,1)));
% tx=linspace(zk1(1,1),zk1(1,2),10);
%%
for i=1:n;  %thegam={}(GPa)
%thegam{1,i}=Ti{1,i}*inv(Ti{1,i})*Q0*(inv(Ti{1,i}))'*def1;
thegam{1,i}=Ti{1,i}*Qi{1,i}*defi{1,i};% def1
e0i{1,i}=(inv(Ti{1,i}))'*defi{1,i};
e12{1,i}=(inv(Ti{1,i}))'*defi2{1,i};
e12{2,i}=(inv(Ti{1,i}))'*defi2{1,i+1};
end
 for i=1:n;    % for assignment3 q2 xy
    thegam2{1,i}=Qi{1,i}*defi2{1,i};
    thegam2{2,i}=Qi{1,i}*defi2{1,i+1};
 end
  for i=1:n;   % for assignment3 q3 12
    thegam3{1,i}=Ti{1,i}*Qi{1,i}*defi2{1,i};
    thegam3{2,i}=Ti{1,i}*Qi{1,i}*defi2{1,i+1};  
 end
 stress=cell2mat(thegam2);
 stress2=cell2mat(thegam3);
%% failure criterion
%Maximum Stress
for i=1:n;
    for e=1:3;
   if thegam3{1,i}(e,1)>=0 | (e==3);
       fmss{1,i}(e,1)=abs(thegam3{1,i}(e,1)/fail(1,2*e-1));
      
   else
       fmss{1,i}(e,1)=abs(thegam3{1,i}(e,1)/fail(1,2*e));
       
   end
          if thegam3{2,i}(e,1)>=0 | (e==3);
       fmss{2,i}(e,1)=abs(thegam3{2,i}(e,1)/fail(1,2*e-1));
      
   else
       fmss{2,i}(e,1)=abs(thegam3{2,i}(e,1)/fail(1,2*e));
       
   end   
    end
end
%Maximum Strain
for i=1:n;
   for e=1:3;
       if e12{1,i}(e,1)>=0 | e==3;
           fmsn{1,i}(e,1)=abs(e12{1,i}(e,1)/fail(1,5+2*e-1));             
                else
           fmsn{1,i}(e,1)=abs(e12{1,i}(e,1)/fail(1,5+2*e));
       end
              if e12{2,i}(e,1)>=0 | e==3;
           fmsn{2,i}(e,1)=abs(e12{2,i}(e,1)/fail(1,5+2*e-1));             
                else
           fmsn{2,i}(e,1)=abs(e12{2,i}(e,1)/fail(1,5+2*e));
       end
   end
end
%Tsai-Hill
for i=1:n;
for e=0:3;
   if e==1;
       if thegam{1,i}(e,1)>=0;
           X=Xt;
       else 
           X=Xc;
       end
   end
   if e==2;
       if thegam{1,i}(e,1)>=0;
           Y=Yt;
       else 
           Y=Yc;
       end
   end
   if e==3;
       S=St;
   end
end
a=thegam{1,i};
fth{1,i}=a(1,1)^2/X^2-a(1,1)*a(2,1)/X^2+a(2,1)^2/Y^2+a(3,1)^2/S^2;
end
%Tsai-Wu
for i=1:n;
    s1=thegam3{1,i}(1);
    s2=thegam3{1,i}(2);
    s3=thegam3{1,i}(3);
    s1u=thegam3{2,i}(1);
    s2u=thegam3{2,i}(2);
    s3u=thegam3{2,i}(3);
    TW(1,i)=abs(F1*s1+F2*s2+F6*s3+F11*s1^2+F22*s2^2+F66*s3^2+2*F12*s1*s2);
    TW(2,i)=abs(F1*s1u+F2*s2u+F6*s3u+F11*s1u^2+F22*s2u^2+F66*s3u^2+2*F12*s1u*s2u);
end
%%
fpl_1=cell2mat(fmss);
failstress=max(fpl_1(:));
[M1,I1]=max(fpl_1(:));
[I1_row, I1_col] = ind2sub(size(fpl_1),I1)
maxstress=Nx*1/failstress;
fpl_2=cell2mat(fmsn);
failstrain=max(fpl_2(:));
[M2,I2]=max(fpl_2(:));
[I2_row, I2_col] = ind2sub(size(fpl_2),I2)
maxstrain=Nx*1/failstrain

%% figures
figure;
plot1(defi2,h,n,zk1);
figure;
plot2(stress,zk1,n);
figure;
plot3(stress2,zk1,n);
%% plot function -- strain
function plot1(defi2,h,n,zk1);  % plot stain
strain = cell2mat(defi2);
yh=linspace(-h/2,h/2);
xs1=linspace(strain(1,1),strain(1,n+1));
xs2=linspace(strain(2,1),strain(2,n+1));
xs3=linspace(strain(3,1),strain(3,n+1));
plot(xs1,yh,xs2,yh,xs3,yh);

grid on;
hold on;
for i =1:n-1;
   ydot(1,i)=zk1(1,i+1);
end

for i=1:n-1;
p=plot([min(strain(:)) max(strain(:))],[ydot(1,i) ydot(1,i)],"--");
hold on;
p.Color='k';
end

 legend('\epsilon_x','\epsilon_y','\epsilon_{xy}');


title('Strain distribution through thickness');
xlabel('strain');
ylabel('thickness (mm)');
end
%% plot function -- stress  xy
    function plot2(stress,zk1,n)
for i=1:n
 yhi(i,:)=linspace(zk1(1,i),zk1(1,i+1));
 xs1(i,:)=linspace(stress(1,i),stress(4,i));
 xs2(i,:)=linspace(stress(2,i),stress(5,i));
 xs3(i,:)=linspace(stress(3,i),stress(6,i));
end
h1=plot(xs1',yhi','b');
hold on;
h2=plot(xs2',yhi','r');
h3=plot(xs3',yhi','g');
grid on;

for i=1:n-1;
p=plot([min(stress(:)) max(stress(:))],[zk1(1,i+1) zk1(1,i+1)],"--");
hold on;
p.Color='k';
end
legend([h1(1) h2(1) h3(1)],'\sigma_x','\sigma_y','\sigma_{xy}');
title('Stress distribution through thickness');
xlabel('stress (GPa)');
ylabel('thickness (mm)');
    end
%% plot function -- stress 12
    function plot3(stress,zk1,n)
for i=1:n
 yhi(i,:)=linspace(zk1(1,i),zk1(1,i+1));
 xs1(i,:)=linspace(stress(1,i),stress(4,i));
 xs2(i,:)=linspace(stress(2,i),stress(5,i));
 xs3(i,:)=linspace(stress(3,i),stress(6,i));
end
h1=plot(xs1',yhi','b');
hold on;
h2=plot(xs2',yhi','r');
h3=plot(xs3',yhi','g');
grid on;

for i=1:n-1;
p=plot([min(stress(:)) max(stress(:))],[zk1(1,i+1) zk1(1,i+1)],"--");
hold on;
p.Color='k';
end
legend([h1(1) h2(1) h3(1)],'\sigma_1','\sigma_2','\sigma_{12}');
title('Stress distribution through thickness');
xlabel('stress (GPa)');
ylabel('thickness (mm)');
%view(90,-90);
end
%% function (assignment 1)
function [Qbar,Ex,Ey,Gxy,miuxy]=layer(data,a);
%%
E1=data(1,1);
E2=data(1,2);
G12=data(1,3);
miu12=data(1,4);
miu21=miu12*E2/E1;
S11=1/E1;
S12=-(miu21/E2);
S21=-(miu12/E1);
S22=1/E2;
S66=1/G12;
Sl=[S11 S12 0;S21 S22 0;0 0 S66];

m=cosd(a);
n=sind(a);
D=[m^4,      n^4 ,      2*m^2*n^2,      m^2*n^2;
   n^4,      m^4,       2*m^2*n^2,      m^2*n^2;
   4*m^2*n^2,4*m^2*n^2, -8*m^2*n^2,    (m^2-n^2)^2;
   m^2*n^2,   m^2*n^2,   m^4+n^4,     -m^2*n^2;
   2*m^3*n,   -2*m*n^3, 2*(m*n^3-m^3*n),m*n^3-m^3*n;
   2*m*n^3, -2*m^3*n,   2*(m^3*n-m*n^3), m^3*n-m*n^3;];
Sb=(D*[S11 S22 S12 S66]')';
Sbar=[Sb(1,1) Sb(1,4) Sb(1,5);
     Sb(1,4) Sb(1,2) Sb(1,6);
     Sb(1,5) Sb(1,6) Sb(1,3);] ;
 Ex=1/Sbar(1,1);
 Ey=1/Sbar(2,2);
 Gxy=1/Sbar(3,3);
 miuxy=-Ex*Sbar(1,2);
 miuyx=miuxy*Ey/Ex;
 Ql=inv(Sl);
 Qbar=inv(Sbar);
end




function [T]=Tr(a)   % transformation and angle
m=cosd(a);
n=sind(a);
T=[m^2 n^2 2*m*n;
    n^2 m^2 -2*m*n;
    -m*n m*n m^2-n^2];
end
 %%