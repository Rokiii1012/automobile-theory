%
zat=0.25;%阻尼比
gama=9;%刚度比y
mju=10;%质量比u
fs=3;zats=0.25;g=9.8;a0=10^(-6);f0=1.5;
ua=20;Gqn0=2.56*10^(-4);n0=0.1;detaf=0.2;N=180;
f=detaf*[0:N];
lamta=f/f0;
lamtas=f/fs;
Wf=0*f;
deta=((1-lamta.^2).*(1+gama-1/mju*lamta.^2)-1).^2+4*zat^2*lamta.^2.*(gama-(1/mju+1)*lamta.^2).^2;
z1_q=gama*sqrt(((1-lamta.^2).^2+4*zat^2*lamta.^2)./deta);
z2_z1=sqrt((1+4*zat^2*lamta.^2)./((1-lamta.^2).^2+4*zat^2*lamta.^2));
p_z2=sqrt((1+(2*zats*lamtas).^2)./((1-lamtas.^2).^2+(2*zats*lamtas).^2));
z2_q=gama*sqrt((1+4*zat^2*lamta.^2)./deta);
p_q=p_z2.*z2_q;
jfg_Gqddf=4*pi^2*sqrt(Gqn0*n0^2*ua)*f;
jfg_Gqdd1f=z1_q.*jfg_Gqddf;
jfg_Gqdd2f=z2_q.*jfg_Gqddf;
jfg_Gaf=p_q.*jfg_Gqddf;
sigmaqdd=sqrt(trapz(f,jfg_Gqddf.^2));%路面不平度加速度均方根值
sigmaqdd1=sqrt(trapz(f,jfg_Gqdd1f.^2));%车轮加速度均方根值
sigmaqdd2=sqrt(trapz(f,jfg_Gqdd2f.^2));%车身加速度均方根值
sigmaa=sqrt(trapz(f,jfg_Gaf.^2));%人体加速度均方根值
for i=1:(N+1)
if f(i)<=2
Wf(i)=0.5;
elseif f(i)<=4
Wf(i)=f(i)/4;
elseif f(i)<=12.5
Wf(i)=1;
else
Wf(i)=12.5/f(i);
end
end
kk=Wf.^2.*jfg_Gaf.^2;
aw=sqrt(trapz(f,kk));%加权加速度均方根值
Law=20*log10(aw/a0);%加权振级
disp('路面不平度加速的均方根值为');disp(sigmaqdd);
disp('车轮加速度均方根值为');disp(sigmaqdd1);
disp('车身加速度均方根值为');disp(sigmaqdd2);
disp('人体加速度均方根值为');disp(sigmaa);
disp('加权加速度均方根值为');disp(aw);
disp('加权振级');disp(Law);
figure(1)
plot(f,z1_q,LineWidth=1),xlabel('激振频率'),ylabel('|z1/q|'),title('幅频特性')
figure(2)
plot(f,z2_z1,LineWidth=1),xlabel('激振频率'),ylabel('|z2/z1|'),title('幅频特性')
figure(3)
plot(f,p_z2,LineWidth=1),xlabel('激振频率'),ylabel('|p/z2|'),title('幅频特性')
figure(4)
plot(f,jfg_Gqdd1f,LineWidth=1);xlabel('激振频率'),ylabel('|(Gz1(f))^1/2|'),title('车轮加速度均方根值')
figure(5)
plot(f,jfg_Gqdd2f,LineWidth=1);xlabel('激振频率'),ylabel('|(Gz2(f))^1/2|'),title('车身加速度均方根值'),
figure(6)
plot(f,jfg_Gaf,LineWidth=1);xlabel('激振频率'),ylabel('|(Ga(f))^1/2|'),title('人体加速度均方根值'),

%6.5 第三问
%-------------------------------------------------------------------------------------------
fs=3;zat_s=0.25;g=9.8;
ua=20;Gqn0=2.56*10^(-4);n0=0.1;detaf=0.2;N=180;
f0=1.5;zat=0.25;gama=9;mju=10;
ff0=[0.25:0.05:3];
sigmaz2=0*ff0;
sigmafd=0*ff0;
sigmaFd_G=0*ff0;
M=length(ff0);

for i=1:M
f0=ff0(i);
f=detaf*[0:N];lamta=f/f0;lamtas=f/fs;
deta=((1-lamta.^2).*(1+gama-1/mju*lamta.^2)-1).^2+4*zat^2*lamta.^2.*(gama-(1/mju+1)*lamta.^2).^2;
z2_qdot=2*pi*f*gama.*sqrt((1+4*zat^2*lamta.^2)./deta);
fd_qdot=gama*lamta.^2./(2*pi*f+eps)./sqrt(deta);
Fd_Gqdot=2*pi*f*gama/g.*sqrt(((lamta.^2/(mju+1)-1).^2+4*zat^2*lamta.^2)./deta);
Gq_dotf=4*pi.^2*Gqn0*n0^2*ua;
Gz2f=(z2_qdot).^2*Gq_dotf;
Gfd_qf=(fd_qdot).^2*Gq_dotf;
GFd_Gf=(Fd_Gqdot).^2*Gq_dotf;
sigmaz2(i)=sqrt(trapz(f,Gz2f));%积分然后开平方
sigmafd(i)=sqrt(trapz(f,Gfd_qf));
sigmaFd_G(i)=sqrt(trapz(f,GFd_Gf));
if f0==1.5
sgmz2=sigmaz2(i);
sgmfd=sigmafd(i);
sgmFd_G=sigmaFd_G(i);
end
end
sz2=20*log10(sigmaz2/sgmz2);
sfd=20*log10(sigmafd/sgmfd);
sFd_G=20*log10(sigmaFd_G/sgmFd_G);
figure
%z2 红色实线    fd 蓝色点划      FD/G 黑色虚线
plot(ff0,sz2,'r',ff0,sfd,'b-.',ff0,sFd_G,'k--',LineWidth=1);
axis([0.25 3 -25 15]);
gtext('σq')
gtext('σfd')
gtext('σFd/G')
xlabel('车身部分固有频率f0/Hz'),ylabel('oz2/dB,ofd/dB,oFd/G/dB'),title('三个响应量均方根值随固有频率f0变化的曲线'),;

%-----------------------------------------------------------------------------------------

fs=3;zat_s=0.25;g=9.8;
ua=20;Gqn0=2.56*10^(-4);n0=0.1;detaf=0.2;N=180;
f0=1.5;zat=0.25;gama=9;mju=10;
c=(0.5-0.125)/180;
zatO=[0.125:c:0.5];
sigmaz2=0.*zatO;
sigmafd=0.*zatO;
sigmaFd_G=0.*zatO;
M=length(zatO);

for i=1:M
zat=zatO(i);
f=detaf*[0:N];lamta=f/f0;lamtas=f/fs;
deta=((1-lamta.^2).*(1+gama-1/mju*lamta.^2)-1).^2+4*zat^2*lamta.^2.*(gama-(1/mju+1)*lamta.^2).^2;
z2_qdot=2*pi*f*gama.*sqrt((1+4*zat^2*lamta.^2)./deta);
fd_qdot=gama*lamta.^2./(2*pi*f+eps)./sqrt(deta);
Fd_Gqdot=2*pi*f*gama/g.*sqrt(((lamta.^2/(mju+1)-1).^2+4*zat^2*lamta.^2)./deta);
Gq_dotf=4*pi.^2*Gqn0*n0^2*ua;
Gz2f=(z2_qdot).^2*Gq_dotf;
Gfd_qf=(fd_qdot).^2*Gq_dotf;
GFd_Gf=(Fd_Gqdot).^2*Gq_dotf;
sigmaz2(i)=sqrt(trapz(f,Gz2f));
sigmafd(i)=sqrt(trapz(f,Gfd_qf));
sigmaFd_G(i)=sqrt(trapz(f,GFd_Gf));
if zat==0.25
sgmz2=sigmaz2(i);
sgmfd=sigmafd(i);
sgmFd_G=sigmaFd_G(i);
end
end
sz2=20*log10(sigmaz2/sgmz2);
sfd=20*log10(sigmafd/sgmfd);
sFd_G=20*log10(sigmaFd_G/sgmFd_G);
figure
plot(zatO,sz2,'r',zatO,sfd,'b-.',zatO,sFd_G,'k--',LineWidth=1);%
axis([0.125 0.5 -4 4]);
gtext('σq')
gtext('σfd')
gtext('σFd/G')
xlabel('车身部分阻尼比ζ（截塔）'),ylabel('oz2/dB,ofd/dB,oFd/G/dB'),title('三个响应量均方根值随阻尼比ζ（截塔）变化的曲线');

%------------------------------------------------------------------------



fs=3;zat_s=0.25;g=9.8;
ua=20;Gqn0=2.56*10^(-4);n0=0.1;detaf=0.2;N=180;
f0=1.5;zat=0.25;gama=9;mju=10;
c=(0.5-0.125)/180;
gamaO=[4:0.1:19];
sigmaz2=0.*gamaO;
sigmafd=0.*gamaO;
sigmaFd_G=0.*gamaO;
M=length(gamaO);

for i=1:M
gama=gamaO(i);
f=detaf*[0:N];lamta=f/f0;lamtas=f/fs;
deta=((1-lamta.^2).*(1+gama-1/mju*lamta.^2)-1).^2+4*zat^2*lamta.^2.*(gama-(1/mju+1)*lamta.^2).^2;
z2_qdot=2*pi*f*gama.*sqrt((1+4*zat^2*lamta.^2)./deta);
fd_qdot=gama*lamta.^2./(2*pi*f+eps)./sqrt(deta);
Fd_Gqdot=2*pi*f*gama/g.*sqrt(((lamta.^2/(mju+1)-1).^2+4*zat^2*lamta.^2)./deta);
Gq_dotf=4*pi.^2*Gqn0*n0^2*ua;
Gz2f=(z2_qdot).^2*Gq_dotf;
Gfd_qf=(fd_qdot).^2*Gq_dotf;
GFd_Gf=(Fd_Gqdot).^2*Gq_dotf;
sigmaz2(i)=sqrt(trapz(f,Gz2f));
sigmafd(i)=sqrt(trapz(f,Gfd_qf));
sigmaFd_G(i)=sqrt(trapz(f,GFd_Gf));
if gama==9
sgmz2=sigmaz2(i);
sgmfd=sigmafd(i);
sgmFd_G=sigmaFd_G(i);
end
end
sz2=20*log10(sigmaz2/sgmz2);
sfd=20*log10(sigmafd/sgmfd);
sFd_G=20*log10(sigmaFd_G/sgmFd_G);
figure
plot(gamaO,sz2,'r',gamaO,sfd,'b-.',gamaO,sFd_G,'k--',LineWidth=1);
axis([4 18 -5 6]);
gtext('σq')
gtext('σfd')
gtext('σFd/G')
xlabel('悬架与轮胎的刚度比γ（伽马）'),ylabel('oz2/dB,ofd/dB,oFd/G/dB'),title('三个响应量均方根值随刚度比γ（伽马）变化的曲线');
%---------------------------------------------------------------------------------------------------


fs=3;zat_s=0.25;g=9.8;
ua=20;Gqn0=2.56*10^(-4);n0=0.1;detaf=0.2;N=180;
f0=1.5;zat=0.25;gama=9;mju=10;
c=(0.5-0.125)/180;
mjuO=[5:0.1:20];
sigmaz2=0*mjuO;
sigmafd=0*mjuO;
sigmaFd_G=0*mjuO;
M=length(mjuO);

for i=1:M
mju=mjuO(i);
f=detaf*[0:N];lamta=f/f0;lamtas=f/fs;
deta=((1-lamta.^2).*(1+gama-1/mju*lamta.^2)-1).^2+4*zat^2*lamta.^2.*(gama-(1/mju+1)*lamta.^2).^2;
z2_qdot=2*pi*f*gama.*sqrt((1+4*zat^2*lamta.^2)./deta);
fd_qdot=gama*lamta.^2./(2*pi*f+eps)./sqrt(deta);
Fd_Gqdot=2*pi*f*gama/g.*sqrt(((lamta.^2/(mju+1)-1).^2+4*zat^2*lamta.^2)./deta);
Gq_dotf=4*pi.^2*Gqn0*n0^2*ua;
Gz2f=(z2_qdot).^2*Gq_dotf;
Gfd_qf=(fd_qdot).^2*Gq_dotf;
GFd_Gf=(Fd_Gqdot).^2*Gq_dotf;
sigmaz2(i)=sqrt(trapz(f,Gz2f));
sigmafd(i)=sqrt(trapz(f,Gfd_qf));
sigmaFd_G(i)=sqrt(trapz(f,GFd_Gf));
if mju==10
sgmz2=sigmaz2(i);
sgmfd=sigmafd(i);
sgmFd_G=sigmaFd_G(i);
end
end
sz2=20*log10(sigmaz2/sgmz2);
sfd=20*log10(sigmafd/sgmfd);
sFd_G=20*log10(sigmaFd_G/sgmFd_G);
figure
plot(mjuO,sz2,'r',mjuO,sfd,'b-.',mjuO,sFd_G,'k--',LineWidth=1);
axis([5 20 -2 2]);
gtext('σq')
gtext('σfd')
gtext('σFd/G')
xlabel('车身与车轮部分质量比μ （缪）'),ylabel('oz2/dB,ofd/dB,oFd/G/dB'),title('三个响应量均方根值随质量比μ（缪）变化的曲线');