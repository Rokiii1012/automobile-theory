clear
m=1818.2;
lz=3885;
L=3.048;
a=1.463;
b=1.585;
k1=-62618;
k2=-110185;
i=20;
g=9.8;
%% 计算稳定车速稳定性因数 k，特征车速uch m/s
K = (m*(a/k2-b/k1))/(L*L);
uch = sqrt(1/K);
%% 稳态横摆角速度增益曲线
ua = 0:0.1:150;
wr = (ua./L)./(1+K.*(ua.^2));%%横摆角速度增益
wr1 = ua./L;
figure
plot(ua,wr,LineWidth=1,Color=[0,0,0]);
xlabel('ua/(km/h)');
ylabel('稳态横摆增益');
% title('稳态横摆角速度增益曲线');
xlim([0,20]);
us = 22.35;%m/s
%% 车速u=22.35m/s时的转向灵敏度
wrs = (us/L)/(1+K*(us.^2));
%% 3
ay = 0.4*g;%侧向力

f_r_wheel_side_angle_diff = L*K*ay;%前后轮侧偏角之差
%% 4  R/R0
rb = 1 + K.*ua.^2;
rb0 = 1 + 0.^ua.^2;
figure
plot(ua.^2,rb,LineWidth=1,Color=[0.85,0.33,0.10]);
hold on;
plot(ua.^2,rb0,LineWidth=1,Color='k');
set(gca,'ytick',0:1:2);
set(gca,'xtick',0:0:0);
ylim([0,2]);
xlim([0,100]);
xlabel('u^2');
ylabel('R/R0');
gtext('中性转向 K=0')
gtext('转向不足 K>0')
% title('转向半径之比曲线');
%% 5 静态储备系数sm
SM = k2/(k1+k2)-a/L;