clear

a = [-3.85445e-12,40.874e-9,-165.44e-6,0.29527,-19.313];
% 发动机转矩
T_tq = [];
% 功率
P_e = [];


%驱动力行驶阻力平衡图
%车轮半径 m    %轴距     %质心至前轴距离(满载)   质心高（满载）
r = 0.367;    L = 3.2;    a1 = 1.947;            h_g = 0.9;
%装载质量 kg    %整车整备质量    %总质量
m1 = 2000;       m2 = 1800;     m3 = 3880;


% 传动系机械效率
yita_T = 0.85;
% 滚动阻力系数
f = 0.013;
%空气阻力系数*迎风面积 m2
C_DA = 2.77;
%主减速器传动比  %变速器传动比
i_0 = 6.33;   
%飞轮转动惯量 kg*m2  %二前轮转动惯量  %四后轮转动惯量
I_f = 0.218;         I_w1 = 1.798;   I_w2 = 3.598;
% 挡位 1档 2档 3档 4档 5档
i_g4 = [6.09,3.09,1.71,1.00];
i_g5 = [5.56,2.769,1.644,1.00,0.793];
%重力

%速度(与转速之间关系为0.377*r*n/i_0/i_g)    % 坡度    %汽车旋转质量换算系数
u_a =[];                                    slope = 0;     delta = 0;

%行驶方程式
F_t = [];



for n_tmp = 600:4000
    %1档
    u_a1(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,1);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    P_e1 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
    F_t1(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,1) * yita_T / r;
    
end


for n_tmp = 600:4000
    %2档
    u_a2(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,2);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    P_e2 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
    F_t2(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,2) * yita_T / r;
    
end

for n_tmp = 600:4000
    %3档
    u_a3(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,3);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    P_e3 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
    F_t3(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,3) * yita_T / r;
    
end


for n_tmp = 600:4000
    %4档
    u_a4(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
    F_t4(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,4) * yita_T / r;
    
end


for n_tmp = 600:4000
    %5档
    u_a5(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
    F_t5(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,5) * yita_T / r;
    
end

for v =1:110
u_a(1,v) = v;
P_load(1,v) = v*f*m3*9.8/3600 + C_DA*v*v*v/76140;
bc (1,v)= P_load(1,v)/yita_T;

    
end

% 绘图
figure


plot(u_a1(1,600:end),P_e1(1,600:end),u_a2(1,600:end),P_e2(1,600:end),u_a3(1,600:end),P_e3(1,600:end),u_a4(1,600:end),P_e4(1,600:end),u_a5(1,600:end),P_e5(1,600:end),u_a,bc);

xlabel('u_a /（km/ h）')
ylabel('Pe / (kW)')
gtext('1档');
gtext('2档');
gtext('3档');
gtext('4档');
gtext('5档');


%-----------------------------------------------------------------------------
%系数矩阵
for P=1:62
P_e(1,P) = P;
c = [0.17768,-5.8629,72.379,-416.46,1326.8];
b1(1,P)= polyval(c,P);

c = [0.043072,-2.0553,36.657,-303.98,1354.7];
b2(1,P)= polyval(c,P);

c = [0.006164,-0.51184,14.524,-189.75,1284.4];
b3(1,P)= polyval(c,P);

c = [0.0018555,-0.18517,7.0035,-121.59,1122.9];
b4(1,P)= polyval(c,P);

c = [0.00068906,-0.091077,4.4763,-98.893,1141.0];
b5(1,P)= polyval(c,P);

c = [0.00035052,-0.05138,2.8593,-73.714,1051.2];
b6(1,P)= polyval(c,P);

c = [0.00028230,-0.047449,2.9788,-84.478,1233.9];
b7(1,P)= polyval(c,P);

 c= [-0.000038568,-0.00075215,0.71113,-45.291,1129.7];
b8(1,P)= polyval(c,P);
end


plot(P_e(1,1:10),b1(1,1:10),P_e(1,1:20),b2(1,1:20),P_e(1,1:end),b3(1,1:end),P_e(1,1:40),b4(1,1:40),P_e(1,1:40),b5(1,1:40),P_e(1,1:50),b5(1,1:50),P_e(1,1:55),b6(1,1:55),P_e(1,1:60),b7(1,1:60),P_e(1,1:62),b8(1,1:62))






%-----------------------------------------------------------------

vector_n = [815,1207,1614,2012,2603,3006,3403,3804];
vector_b =[];

 n_tmp = vector_n(1,1);
 vector_b(1,1) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(2,1) = b1(1,round(P_e5 (1,n_tmp)));
 
 n_tmp = vector_n(1,2);
 vector_b(1,2) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(2,2) = b2(1,round(P_e5 (1,n_tmp)));
 
 n_tmp = vector_n(1,3);
 vector_b(1,3) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(2,3) = b3(1,round(P_e5 (1,n_tmp)));
 
 n_tmp = vector_n(1,4);
 vector_b(1,4) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(2,4) = b4(1,round(P_e5 (1,n_tmp)));
    
 n_tmp = vector_n(1,5);
 vector_b(1,5) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(2,5) = b5(1,round(P_e5 (1,n_tmp)));

  n_tmp = vector_n(1,6);
 vector_b(1,6) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(2,6) = b6(1,round(P_e5 (1,n_tmp)));
 
  n_tmp = vector_n(1,7);
 vector_b(1,7) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(2,7) = b7(1,round(P_e5 (1,n_tmp)));
 
 n_tmp = vector_n(1,8);
 vector_b(1,8) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e5 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(2,8) = b8(1,round(P_e5 (1,n_tmp)));
%   
%   vector_b(:,1) = [];
%   vector_b(:,2) = [];
%   vector_b(:,1) = [];

 
 
 
 
 
 n_tmp = vector_n(1,1);
 vector_b(3,1) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(4,1) = b1(1,round(P_e4 (1,n_tmp)));
 
 n_tmp = vector_n(1,2);
 vector_b(3,2) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(4,2) = b2(1,round(P_e4 (1,n_tmp)));
 
 n_tmp = vector_n(1,3);
 vector_b(3,3) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(4,3) = b3(1,round(P_e4 (1,n_tmp)));
 
 n_tmp = vector_n(1,4);
 vector_b(3,4) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(4,4) = b4(1,round(P_e4 (1,n_tmp)));
    
 n_tmp = vector_n(1,5);
 vector_b(3,5) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(4,5) = b5(1,round(P_e4 (1,n_tmp)));

  n_tmp = vector_n(1,6);
 vector_b(3,6) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(4,6) = b6(1,round(P_e4 (1,n_tmp)));
 
  n_tmp = vector_n(1,7);
 vector_b(3,7) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(4,7) = b7(1,round(P_e4 (1,n_tmp)));
 
 n_tmp = vector_n(1,8);
 vector_b(3,8) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
 T_tq(1,n_tmp)= polyval(a,n_tmp);
 P_e4 (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
 vector_b(4,8) = b7(1,round(P_e4 (1,n_tmp)));
%  
 vector_b(:,1)= [];
 vector_b(:,2)= [];
 
 
 v_tmp =linspace(0,120);
 p1 = polyfit(vector_b(1,:),vector_b(2,:),3);
 f1 = polyval(p1,v_tmp);
 
 p2 = polyfit(vector_b(3,:),vector_b(4,:),3);
 f2 = polyval(p2,v_tmp);
%  figure
%   plot(v_tmp,f1,vector_b(1,1:end),vector_b(2,1:end),'o')
  
 figure
 plot(v_tmp,f1,vector_b(1,1:end),vector_b(2,1:end),'o',v_tmp,f2,vector_b(3,1:end),vector_b(4,1:end),'*')
%  title('最高档和次高档等速燃油消耗率')
 xlabel('u_a - km/h')
 ylabel('燃油消耗率b[g/(kw*h)]')
 xlim([0,90]);
 ylim([0,1000]);
 gtext('5档');
  gtext('4档');


pg = 7.00;
for i = 1:6

    vector_p(1,i) = (vector_b(1,i)*((3880*9.8)*f+C_DA*vector_b(1,i)*vector_b(1,i)/21.15))/(3600*yita_T);
    vector_p(2,i) = vector_p(1,i)*vector_b(2,i)/(1.02*vector_b(1,i)*pg);
    
    
    vector_p(3,i) = (vector_b(3,i)*((3880*9.8)*f+C_DA*vector_b(3,i)*vector_b(3,i)/21.15))/(3600*yita_T);
    vector_p(4,i) = vector_p(3,i)*vector_b(4,i)/(1.02*vector_b(3,i)*pg);
    
end


u_a =[];
u_a = linspace(10,120,6);
 p3 = polyfit(u_a(1,:),vector_p(2,:),3);
 f3 = polyval(p3,v_tmp);
%  figure
%  plot(v_tmp,f3,u_a,vector_p(2,:),'*')
%  title('最高档汽车百公里油耗曲线')
%  xlabel('u_a - km/h')
%  ylabel('Qs/（L/100km）')
% 
%  
 

 u_a =[];
u_a = linspace(10,120,6);
 p4 = polyfit(u_a(1,:),vector_p(4,:),3);
 f4 = polyval(p4,v_tmp);
 figure
 plot(v_tmp,f4,u_a,vector_p(4,:),'o',v_tmp,f3,u_a,vector_p(2,:),'*')
%  axis(0,130,0,100)
%  title('最高档汽车百公里油耗曲线与次高档汽车百公里油耗曲线')
 xlabel('u_a - km/h')
 ylabel('Qs/（L/100km）')
 xlim([0 130])
 ylim([0 30])
  gtext('5档');
  gtext('4档');
 
 
 
 
 
 