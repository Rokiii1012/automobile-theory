

%----------------------------------------------------------------------------
clear
% 图1发动机功率外特性和转矩外特性曲线
% 参数


%空载质量      hg0           L        a0     
m0 =3980;     hg0 = 0.8;    L=3.95;   a0=2.2;
%满载质量      hg
m1 =9000;     hg1 = 1.17;             a1 = 2.95; 

F_u1 = linspace(0,35000,100000);
F_u11(1,:) = F_u1;
F_u22 = F_u1.* ((1-0.38)/0.38); 
G = m0*9.8;
for i=1:100000
%    F_u2(1,i)=(G*(sqrt((L-a0)*(L-a0)+4*hg0*L*F_u1(1,i)/G))/hg0-(G*(L-a0)/hg0+2*F_u1(1,i)))/2;
   
   F_u2(1,i) = (G/hg0 * sqrt((L-a0)*(L-a0)+ 4 * hg0 * L * F_u1(1,i) /G) -(G*(L-a0)/hg0+2*F_u1(1,i)) )/2;
   fi_x(1,i) =(F_u1(1,i)+F_u2(1,i))/G;
end



fi = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];


for i=1:10
   
    for j= 1:100000
        
    if(F_u2(1,j)>= F_u1(1,j)*(L-hg0*fi(1,i))/((fi(1,i)*hg0)) - (G*(L-a0))/hg0)
    Fxb2(i,j)=F_u1(1,j)*(L-hg0*fi(1,i))/((fi(1,i)*hg0)) - (G*(L-a0))/hg0;
    else
    Fxb2(i,j) = F_u2(1,j);
    end
    
    if(F_u2(1,j)<=(-fi(1,i)*hg0)/(L+fi(1,i)*hg0)*F_u1(1,j)+ fi(1,i)*G*a0/(L+fi(1,i)*hg0))
    Fxb2(10+i,j) = (-fi(1,i)*hg0)/(L+fi(1,i)*hg0)*F_u1(1,j)+ fi(1,i)*G*a0/(L+fi(1,i)*hg0);
    else
    Fxb2(i+10,j) = F_u2(1,j);
    end
    end
end


F_u3(1,:)=F_u2;
figure
plot(F_u1,F_u2,F_u1,F_u22,F_u1,Fxb2(1,:),"r",F_u1,Fxb2(2,:),"r",F_u1,Fxb2(3,:),"r",F_u1,Fxb2(4,:),"r",F_u1,Fxb2(5,:),"r",F_u1,Fxb2(6,:),"r",F_u1,Fxb2(7,:),"r",F_u1,Fxb2(8,:),"r",F_u1,Fxb2(9,:),"r",F_u1,Fxb2(10,:),"r")
ylim([0,35000])
hold on
plot(F_u1,Fxb2(11,:),"k",F_u1,Fxb2(12,:),"k",F_u1,Fxb2(13,:),"k",F_u1,Fxb2(14,:),"k",F_u1,Fxb2(15,:),"k",F_u1,Fxb2(16,:),"k",F_u1,Fxb2(17,:),"k",F_u1,Fxb2(18,:),"k",F_u1,Fxb2(19,:),"k",F_u1,Fxb2(20,:),"k")
% title('空载')
 gtext('I曲线');
 gtext('β曲线（空载）');
 gtext('f线组');
 gtext('r线组');

F_u1 = linspace(0,85000,100000);
F_u11(2,:) = F_u1;
F_u22 = F_u1.* ((1-0.38)/0.38); 
G = m1*9.8;
for i=1:100000

   F_u2(1,i) = (G/hg1 * sqrt((L-a1)*(L-a1)+ 4 * hg1 * L * F_u1(1,i) /G) -(G*(L-a1)/hg1+2*F_u1(1,i)) )/2;
   
   fi_x(2,i) =(F_u1(1,i)+F_u2(1,i))/G;
end




fi = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];

for i=1:10
   
    for j= 1:100000
    if(F_u2(1,j)>=F_u1(1,j)*(L-hg1*fi(1,i))/((fi(1,i)*hg1)) - (G*(L-a1))/hg1)
    Fxb2(i,j)=F_u1(1,j)*(L-hg1*fi(1,i))/((fi(1,i)*hg1)) - (G*(L-a1))/hg1;
    else
    Fxb2(i,j) = F_u2(1,j);
    end
    
    if(F_u2(1,j)>=F_u1(1,j)*(L-hg1*fi(1,i))/((fi(1,i)*hg1)) - (G*(L-a1))/hg1)
    Fxb2(10+i,j) = (-fi(1,i)*hg1)/(L+fi(1,i)*hg1)*F_u1(1,j)+ fi(1,i)*G*a1/(L+fi(1,i)*hg1);
    else
    Fxb2(10+i,j) = F_u2(1,j);
    end
    
    end
end

fi0=(14334.5+23399.2)/(m1*9.8);
F_u222 = (-1).*F_u1+fi0*(m1*9.8);


F_u3(2,:)=F_u2;
figure
plot(F_u1,F_u2,F_u1,F_u22,F_u1,Fxb2(1,:),F_u1,Fxb2(2,:),F_u1,Fxb2(3,:),F_u1,Fxb2(4,:),F_u1,Fxb2(5,:),F_u1,Fxb2(6,:),F_u1,Fxb2(7,:),F_u1,Fxb2(8,:),F_u1,Fxb2(9,:),F_u1,Fxb2(10,:))
ylim([0,85000])
hold on
plot(F_u1,Fxb2(11,:),F_u1,Fxb2(12,:),F_u1,Fxb2(13,:),F_u1,Fxb2(14,:),F_u1,Fxb2(15,:),F_u1,Fxb2(16,:),F_u1,Fxb2(17,:),F_u1,Fxb2(18,:),F_u1,Fxb2(19,:),F_u1,Fxb2(20,:))
hold on
plot(F_u1,F_u222,"--");

% title('满载')
 gtext('I曲线');
 gtext('β曲线（满载）');
 gtext('f线组');
 gtext('r线组');
  gtext('φ0=0.4278');
 
 
fi0=(14334.5+23399.2)/(m1*9.8);
F_u222 = (-1).*F_u1+fi0*(m1*9.8);
F_u333 = (-1).*F_u1+1.0*(m0*9.8);

figure
% plot(F_u11(1,:),F_u3(1,:),F_u11(2,:),F_u3(2,:));
% plot(F_u1,F_u333,'--');
% hold on
plot(F_u1,F_u22,F_u11(1,1:72213),F_u3(1,1:72213),F_u11(2,:),F_u3(2,:),F_u1,F_u222,"--");
% plot(F_u1,F_u22)
xlabel('Fu1/n');
ylabel('Fu2/n');
 xlim([0 35000])
 ylim([0 34000])
 gtext('I曲线');
  gtext('β曲线（满载）');
    gtext('β曲线（空载）');
      gtext('φ0=0.4278');




%--------------------------------------------
beta =0.38;
%利用附着系数
fi_f = [];
fi_r = [];
z = linspace(0.0001,1,10000);
for i= 1:10000
    
 fi_f(1,i)= beta * z(1,i)*L/((L-a0)+z(1,i)*hg0);
 fi_f(2,i)= beta * z(1,i)*L/((L-a1)+z(1,i)*hg1);
 fi_r(1,i)= (1-beta)*z(1,i)*L/(a0-z(1,i)*hg0);
 fi_r(2,i)= (1-beta)*z(1,i)*L/(a1-z(1,i)*hg1);
 
 E_f(1,i) = z(1,i)/fi_f(1,i).*100;
 E_f(2,i) = z(1,i)/fi_f(2,i).*100;
 E_r(1,i) = z(1,i)/fi_r(1,i).*100;
 E_r(2,i) = z(1,i)/fi_r(2,i).*100;

    
    
end
figure
plot(z,z,"--",z,fi_f(1,:),"k",z,fi_f(2,:),"k",z,fi_r(1,:),"b-.",z,fi_r(2,:),"b-.")
title('利用附着系数关系图');
xlabel('制动强度z/g');
ylabel('利用附着系数');
xlim([0 1]);
ylim([0 1]);
hold on
u1=0.2:0.01:0.8;
z3=0.1+0.85*(u1-0.2);
plot(z3,u1,'k','linewidth',1.5);
hold on
 
%对于最大总质量大于3.5t的货车
z1=0.15:0.01:0.3;
u1=z1+0.08;
u2=z1-0.08;
plot(z1,u1,'k',z1,u2,'k','linewidth',1.5)
hold on
 
%当制动强度z≥0.3时
z2=0.3:0.01:1;
u2=((z2-0.3)/0.74)+0.38;
plot(z2,u2,'k','linewidth',1.5)
hold on
 
x=0:0.01:1;
y=x;
plot(x,y,'k--','linewidth',1.5)
grid on
% title('ECE法规规定的最大总质量超过3.5t货车的制动力分配曲线')
axis([0,0.8,0,0.8])
xlabel('制动强度z/g');
ylabel('利用附着系数φ');
% legend('φ=(z+0.07)/0.85','φ=z+0.08','φ=z-0.08','φ=(z-0.02)/0.74','φ=z')
% plot(z(1,1500:3048),z(1,1500:3048)-0.08,"k",z(1,1500:3048),z(1,1500:3048)+0.08,"k",'linewidth',1);
% plot(z(1,3048:3113),z(1,3048:3113)./0.74-0.02/0.74,"k",z(1,1200:3113),z(1,1200:3113)+0.07/0.85,"k",'linewidth',1);
gtext('φr（空载）');
gtext('φr（满载）');
gtext('φf（空载）');
gtext('φf（满载）');
% 
% 
% gtext('φ=(z+0.07)/0.85');
%  gtext('φ=z+0.08');
%  gtext('φ=z-0.08');
% gtext('φ=(z-0.02)/0.74');
% gtext('φ=z');
gtext('ECE法规');

 
figure
plot(z,E_f(1,:),'k--',z,E_f(2,:),'k--',z,E_r(1,:),'b',z,E_r(2,:),'b');
title('前后轴制动效率曲线');
xlabel('附着系数');
ylabel('制动效率%');
ylim([0 100]);
gtext('Ef');
 gtext('Er');
  gtext('Er');
 gtext('空载');
gtext('满载');



