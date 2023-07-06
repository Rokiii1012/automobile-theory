clear
% 汽车动力性计算


%----------------------------------------------------------------------------
% 图1发动机功率外特性和转矩外特性曲线
% 参数


% 转速，单位r/min
n = [];
%车轮半径，单位m

% 发动机转矩曲线四次多项式的系数矩阵 a4,a3,a2,a1,a0
a = [-3.85445e-12,40.874e-9,-165.44e-6,0.29527,-19.313];
% 发动机转矩
T_tq = [];
% 功率
P_e = [];


for n_tmp = 600:4000
    n(1,n_tmp) = n_tmp;
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    P_e (1,n_tmp)= T_tq(1,n_tmp)*n_tmp/9550;
end

% 绘图
figure
yyaxis left
plot(n(1,600:end),T_tq(1,600:end));
title('发动机功率外特性和转矩外特性曲线')
xlabel('n /（r / min）')
ylabel('Ttq / (N * m)')

yyaxis right
plot(n(1,600:end),P_e(1,600:end));
ylim([0 60])
ylabel('Pe / (kW)')


%----------------------------------------------------------------------------------------------------
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

% F_w = 
% F_j ;

for n_tmp = 600:4000
    %1档
    u_a1(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,1);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t1(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,1) * yita_T / r;
    
end


for n_tmp = 600:4000
    %2档
    u_a2(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,2);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t2(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,2) * yita_T / r;
    
end

for n_tmp = 600:4000
    %3档
    u_a3(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,3);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t3(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,3) * yita_T / r;
    
end


for n_tmp = 600:4000
    %4档
    u_a4(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t4(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,4) * yita_T / r;
    
end


for n_tmp = 600:4000
    %5档
    u_a5(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t5(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,5) * yita_T / r;
    
end

for n_tmp = 1:110

u_a(1,n_tmp) = n_tmp-1;
G = m3 * 9.8;
F_f = G * f;
F_w (1,n_tmp) =  C_DA * u_a(1,n_tmp) * u_a(1,n_tmp) /21.15;
sum (1,n_tmp) = F_f + F_w(1,n_tmp);    
    
end


% 绘图
figure
 plot(u_a1(1,1000:4000),F_t1(1,1000:4000),u_a2(1,1000:4000),F_t2(1,1000:4000),u_a3(1,1000:4000),F_t3(1,1000:4000),u_a4(1,1000:4000),F_t4(1,1000:4000),u_a5(1,1000:4000),F_t5(1,1000:4000),u_a,sum);
% plot(u_a5(1,1000:4000),F_t5(1,1000:4000),u_a,sum);
% title('汽车驱动力-行驶阻力平衡图')
xlabel('u_a /（km / h）')
ylabel('F_t / N ,(F_f + F_w) / N')



%-----------------------------------------------------------------------------
%动力特性图

for n_tmp = 600:4000
    %1档
    u_a1(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,1);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t1(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,1) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a1(1,n_tmp) * u_a1(1,n_tmp) /21.15;
    D1(1,n_tmp) = (F_t1(1,n_tmp)-F_w(1,n_tmp))/G;
    
end


for n_tmp = 600:4000
    %2档
    u_a2(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,2);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t2(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,2) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a2(1,n_tmp) * u_a2(1,n_tmp) /21.15;
    D2(1,n_tmp) = (F_t2(1,n_tmp)-F_w(1,n_tmp))/G;
    
end

for n_tmp = 600:4000
    %3档
    u_a3(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,3);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t3(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,3) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a3(1,n_tmp) * u_a3(1,n_tmp) /21.15;
    D3(1,n_tmp) = (F_t3(1,n_tmp)-F_w(1,n_tmp))/G;
    
end


for n_tmp = 600:4000
    %4档
    u_a4(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t4(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,4) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a4(1,n_tmp) * u_a4(1,n_tmp) /21.15;
    D4(1,n_tmp) = (F_t4(1,n_tmp)-F_w(1,n_tmp))/G;
    
end


for n_tmp = 600:4000
    %5档
    u_a5(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t5(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,5) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a5(1,n_tmp) * u_a5(1,n_tmp) /21.15;
    D5(1,n_tmp) = (F_t5(1,n_tmp)-F_w(1,n_tmp))/G;
    
end

% 绘图
figure
plot(u_a1(1,1000:4000),D1(1,1000:4000),u_a2(1,1000:4000),D2(1,1000:4000),u_a3(1,1000:4000),D3(1,1000:4000),u_a4(1,1000:4000),D4(1,1000:4000),u_a5(1,1000:4000),D5(1,1000:4000));

xlabel('u_a /（km / h）')
ylabel('D')



%符合率------------------------------------------------------------------------------
for n_tmp = 600:4000
    %1档
    u_a1(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,1);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t1(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,1) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a1(1,n_tmp) * u_a1(1,n_tmp) /21.15;
    D1(1,n_tmp) = (F_t1(1,n_tmp)-F_w(1,n_tmp))/G;
    
%     load_rate1(1,n_tmp) =  F_t1(1,n_tmp)/max(F_t1);
    load_rate1(1,n_tmp) = (F_w(1,n_tmp)+F_f)/F_t1(1,n_tmp);
  
    
end


for n_tmp = 600:4000
    %2档
    u_a2(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,2);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t2(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,2) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a2(1,n_tmp) * u_a2(1,n_tmp) /21.15;
    D2(1,n_tmp) = (F_t2(1,n_tmp)-F_w(1,n_tmp))/G;
%     load_rate2(1,n_tmp) =  F_t2(1,n_tmp)/max(F_t2);
   load_rate2(1,n_tmp) = (F_w(1,n_tmp)+F_f)/F_t2(1,n_tmp);
end

for n_tmp = 600:4000
    %3档
    u_a3(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,3);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t3(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,3) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a3(1,n_tmp) * u_a3(1,n_tmp) /21.15;
    D3(1,n_tmp) = (F_t3(1,n_tmp)-F_w(1,n_tmp))/G;
%     load_rate3(1,n_tmp) =  F_t3(1,n_tmp)/max(F_t3);
      load_rate3(1,n_tmp) = (F_w(1,n_tmp)+F_f)/F_t3(1,n_tmp);
    
end


for n_tmp = 600:4000
    %4档
    u_a4(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t4(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,4) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a4(1,n_tmp) * u_a4(1,n_tmp) /21.15;
    D4(1,n_tmp) = (F_t4(1,n_tmp)-F_w(1,n_tmp))/G;
%     load_rate4(1,n_tmp) =  F_t4(1,n_tmp)/max(F_t4);
    load_rate4(1,n_tmp) = (F_w(1,n_tmp)+F_f)/F_t4(1,n_tmp);
    
end


for n_tmp = 600:4000
    %5档
    u_a5(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t5(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,5) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a5(1,n_tmp) * u_a5(1,n_tmp) /21.15;
    D5(1,n_tmp) = (F_t5(1,n_tmp)-F_w(1,n_tmp))/G;
%     load_rate5(1,n_tmp) =  F_t5(1,n_tmp)/max(F_t5);
load_rate5(1,n_tmp) = (F_w(1,n_tmp)+F_f)/F_t5(1,n_tmp);
    
end

% 绘图
figure
plot(u_a1(1,1000:4000),load_rate1(1,1000:4000).*100,u_a2(1,1000:4000),load_rate2(1,1000:4000).*100,u_a3(1,1000:4000),load_rate3(1,1000:4000).*100,u_a4(1,1000:4000),load_rate4(1,1000:4000).*100,u_a5(1,1000:3500),load_rate5(1,1000:3500).*100);

xlabel('u_a /（km / h）')
ylabel('负荷率（%）')
ylim([0 100]);







%------------------------------------------------------------------------------------------
 %汽车旋转质量换算系数
 
for n_tmp = 600:4000
    %1档
    delta = 1 + (I_w1+I_w2)/((r*r)*m3) + (I_f*i_0*i_0*i_g5(1,1)*i_g5(1,1)*yita_T)/((r*r)*m3);
    u_a1(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,1);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t1(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,1) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a1(1,n_tmp) * u_a1(1,n_tmp) /21.15;
    D1(1,n_tmp) = (F_t1(1,n_tmp)-F_w(1,n_tmp))/G;
    acc1(1,n_tmp) = ((F_t1(1,n_tmp) - F_f - F_w(1,n_tmp))/(m3 * delta));
    acc1_r(1,n_tmp) =  1/acc1(1,n_tmp);
end

u_a2 = [];
for n_tmp = 600:4000
    %2档
    
    delta = 1 + (I_w1+I_w2)/((r*r)*m3) + (I_f*i_0*i_0*i_g5(1,2)*i_g5(1,2)*yita_T)/((r*r)*m3);
    u_a2(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,2);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t2(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,2) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a2(1,n_tmp) * u_a2(1,n_tmp) /21.15;
    D2(1,n_tmp) = (F_t2(1,n_tmp)-F_w(1,n_tmp))/G;
    acc2(1,n_tmp) = ((F_t2(1,n_tmp) - F_f - F_w(1,n_tmp))/(m3 * delta));
    acc2_r(1,n_tmp) =  1/acc2(1,n_tmp);
    
end
u_a3 = [];
for n_tmp = 600:4000
    %3档
    delta = 1 + (I_w1+I_w2)/((r*r)*m3) + (I_f*i_0*i_0*i_g5(1,3)*i_g5(1,3)*yita_T)/((r*r)*m3);
    u_a3(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,3);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t3(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,3) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a3(1,n_tmp) * u_a3(1,n_tmp) /21.15;
    D3(1,n_tmp) = (F_t3(1,n_tmp)-F_w(1,n_tmp))/G;
    acc3(1,n_tmp) = ((F_t3(1,n_tmp) - F_f - F_w(1,n_tmp))/(m3 * delta));
    acc3_r(1,n_tmp) =  1/acc3(1,n_tmp);
    
end

u_a4 = [];
for n_tmp = 600:4000
    %4档
    delta = 1 + (I_w1+I_w2)/((r*r)*m3) + (I_f*i_0*i_0*i_g5(1,4)*i_g5(1,4)*yita_T)/((r*r)*m3);
    u_a4(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,4);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t4(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,4) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a4(1,n_tmp) * u_a4(1,n_tmp) /21.15;
    D4(1,n_tmp) = (F_t4(1,n_tmp)-F_w(1,n_tmp))/G;
    
    acc4(1,n_tmp) = ((F_t4(1,n_tmp) - F_f - F_w(1,n_tmp))/(m3 * delta));
    acc4_r(1,n_tmp) =  1/acc4(1,n_tmp);
    
end

u_a5 = [];
for n_tmp = 600:4000
    %5档
    delta = 1 + (I_w1+I_w2)/((r*r)*m3) + (I_f*i_0*i_0*i_g5(1,5)*i_g5(1,5)*yita_T)/((r*r)*m3);
    u_a5(1,n_tmp) = 0.377 * r * n_tmp / i_0 / i_g5(1,5);
    T_tq(1,n_tmp)= polyval(a,n_tmp);
    F_t5(1,n_tmp) = T_tq(1,n_tmp) * i_0 * i_g5(1,5) * yita_T / r;
    F_w (1,n_tmp) =  C_DA * u_a5(1,n_tmp) * u_a5(1,n_tmp) /21.15;
    D5(1,n_tmp) = (F_t5(1,n_tmp)-F_w(1,n_tmp))/G;
    acc5(1,n_tmp) = ((F_t5(1,n_tmp) - F_f - F_w(1,n_tmp))/(m3 * delta));
    acc5_r(1,n_tmp) =  1/acc5(1,n_tmp);
    
end


figure
plot(u_a1(1,1000:4000),acc1(1,1000:4000),u_a2(1,1000:4000),acc2(1,1000:4000),u_a3(1,1000:4000),acc3(1,1000:4000),u_a4(1,1000:3500),acc4(1,1000:3500),u_a5(1,1000:3000),acc5(1,1000:3000));
xlabel('u_a /（km / h）')
ylabel('a/(m/s2)')


figure
plot(u_a1(1,600:4000),acc1_r(1,600:4000),u_a2(1,600:4000),acc2_r(1,600:4000),u_a3(1,600:4000),acc3_r(1,600:4000),u_a4(1,600:4000),acc4_r(1,600:4000),u_a5(1,1000:3000),acc5_r(1,1000:3000));
xlabel('u_a ')
ylabel('1/a')



% 
% t1 = cumtrapz(acc1_r(1,600:4000)).*((max(u_a1)-min(u_a1(1,600:end)))/3600);
% t2 = cumtrapz(acc2_r(1,600:4000)).*((max(u_a2)-min(u_a2(1,600:end)))/3600);
% t3 = cumtrapz(acc3_r(1,600:4000)).*((max(u_a3)-min(u_a3(1,600:end)))/3600);
% t4 = cumtrapz(acc4_r(1,600:4000)).*((max(u_a4)-min(u_a4(1,600:end)))/3600);
% t5 = cumtrapz(acc5_r(1,600:3500)).*((max(u_a5(1,600:3500))-min(u_a5(1,600:3500)))/2900);


t1 = cumtrapz((max(u_a1)-min(u_a1(1,600:end)))/3600,acc1_r(1,600:4000)./3.6);
t2 = cumtrapz((max(u_a2)-min(u_a2(1,600:end)))/3600,acc2_r(1,600:4000)./3.6);
t3 = cumtrapz((max(u_a3)-min(u_a3(1,600:end)))/3600,acc3_r(1,600:4000)./3.6);
t4 = cumtrapz((max(u_a4)-min(u_a4(1,600:end)))/3600,acc4_r(1,600:4000)./3.6);
t5 = cumtrapz((max(u_a5(1,600:3000))-min(u_a5(1,600:3000)))/2900,acc5_r(1,600:3000)./3.6);

figure
% plot(t1(1,600:610),u_a1(1,600:610),t2(1,600:610),u_a2(1,600:610),t3(1,600:610),u_a3(1,600:610),t4(1,600:610),u_a4(1,600:610),t5(1,600:610),u_a5(1,600:610))
% plot(t1(1,1:1193),u_a1(1,600:1792),t2(1,296:803),u_a2(1,895:1402),t3(1,233:end),u_a3(1,832:4000))
% xlim([0 100])
% plot(t1(1,1:end),u_a1(1,600:4000),t2(1,1:end),u_a2(1,600:4000),t3(1,1:end),u_a3(1,600:4000),t4(1,1:end),u_a4(1,600:4000),t5(1,1:end),u_a5(1,600:3500))
vector_acc = [];
vector_acc(1,1)=0;
vector_acc(2,1)=0;

vector_acc(1,2)=max(t1);
vector_acc(2,2)=max(u_a1);

vector_acc(1,3)=max(t2);
vector_acc(2,3)=max(u_a2);

vector_acc(1,4)=max(t3);
vector_acc(2,4)=max(u_a3);

vector_acc(1,5)=max(t4);
vector_acc(2,5)=max(u_a4);
plot(vector_acc(1,1),vector_acc(2,1),vector_acc(1,2),vector_acc(2,2),vector_acc(1,3),vector_acc(2,3),vector_acc(1,4),vector_acc(2,4),vector_acc(1,5),vector_acc(2,5));
line([vector_acc(1,1),vector_acc(1,2),vector_acc(1,3),vector_acc(1,4),vector_acc(1,5)],[vector_acc(2,1),vector_acc(2,2),vector_acc(2,3),vector_acc(2,4),vector_acc(2,5)])
% title('汽车加速时间曲线')
xlabel('t-s ')
ylabel('ua-km/h')
xlim([0 60])
hold on
% figure
% plot(t4,u_a4(1,600:4000),t5,u_a5(1,600:3000));
plot(t5,u_a5(1,600:3000)+43.5);
% title('最高档和次高档加速时间')
% xlabel('t-s ')
% ylabel('ua-km/h')
ylim([0 95])
% gtext('汽车最高档超车加速时间曲线')
% gtext('汽车原地起步连续换挡加速时间曲线')

