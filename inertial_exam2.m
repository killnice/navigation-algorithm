% 方案一 利用角增量构建更新四元数
x = -0.10*pi/180;
y = -0.16*pi/180;
z = -10.01*pi/180; %弧度单位制
%注意弧度单位制 注意四元数向量
q0 = [cos(x/2)*cos(y/2)*cos(z/2) + sin(x/2)*sin(y/2)*sin(z/2)
      sin(x/2)*cos(y/2)*cos(z/2) - cos(x/2)*sin(y/2)*sin(z/2)
      cos(x/2)*sin(y/2)*cos(z/2) + sin(x/2)*cos(y/2)*sin(z/2)
      cos(x/2)*cos(y/2)*sin(z/2) - sin(x/2)*sin(y/2)*cos(z/2)];
q0 = reshape(q0,[1,4]);
A=importdata('w.csv'); %csv数据存储文件
b=A.textdata;
time_temp = char(b(:,1));% 提取时间戳字符串 -->cell转char
time_list = [];%存储时间戳 t
dt_list = []; %存储时间戳的差分 dt
w_list = A.data(:,7:9)*pi/180 ;%存储陀螺仪角速度计测量值 w
for i = 2:length(time_temp) 
    time_list = [time_list ; str2double(time_temp(i,8:end))];
end

for i = 2:length(time_list)
    dt_list = [dt_list ;time_list(i)- time_list(i-1)];
end
Q= q0;
angle_list_one =[];
for i = 1:length(w_list)-1 %四元数更新
    
    dtheta = w_list(i,:)*dt_list(i);
    dtheta_mod = norm(dtheta);
    rk_real = cos(dtheta_mod/2);
    if dtheta_mod==0
        rk_img =[0,0,0];
    else 
        rk_img =dtheta/dtheta_mod*sin(dtheta_mod/2);
    end
    rk = [rk_real rk_img];
    Q = quatmultiply(Q,rk);
    x = atan2(2*(Q(1)*Q(2)+Q(3)*Q(4)),1 - 2*(Q(2)^2+Q(3)^2));
    y = asin(2*(Q(1)*Q(3) - Q(2)*Q(4)));
    z = atan2(2*(Q(1)*Q(4)+Q(2)*Q(3)),1 - 2*(Q(3)^2+Q(4)^2));
    angle = [x,y,z]*180/pi;
    angle_list_one = [angle_list_one ; angle];
end
% figure(1)



% 方案二 根据四元数微分方程直接进行数值求解


I = eye(4);
Q_one = q0.';
Q_two = q0.';
Q_three = q0.';
angle_list_two = [];
for i = 1:length(w_list)-1 %四元数更新
    
    wx = w_list(i,1);
    wy = w_list(i,2);
    wz = w_list(i,3);
    
    dtheta = [0, -wx ,-wy , -wz;
              wx,  0 , wz , -wy;
              wy, -wz,  0 ,  wx;
              wz,  wy, -wx,   0;]*dt_list(i);
          
    Q_one = (I + 0.5 * dtheta) * Q_one;
    Q_two = (I + 0.5 * dtheta + 0.5*(0.5 * dtheta)^2) * Q_two;
    Q_three = (I + 0.5 * dtheta + 0.5*(0.5 * dtheta)^2 + 1/6 * (0.5 * dtheta)^3) * Q_three;
    Q = Q_three;
    x = atan2(2*(Q(1)*Q(2)+Q(3)*Q(4)),1 - 2*(Q(2)^2+Q(3)^2));
    y = asin(2*(Q(1)*Q(3) - Q(2)*Q(4)));
    z = atan2(2*(Q(1)*Q(4)+Q(2)*Q(3)),1 - 2*(Q(3)^2+Q(4)^2));
    angle = [x,y,z]*180/pi;
    angle_list_two = [angle_list_two ;angle];
    
     
end
figure(1)
plot(angle_list_one)
title('利用角增量构建更新四元数的曲线')
figure(2)
plot(angle_list_two)
title('四元数微分方程直接进行数值求解')


% 
% plot(angle_list_one(:,3),'r')
% 
% plot(angle_list_two(:,3) - angle_list_one(:,3),'b')
% 
% legend('第一种方案','第二种方案')
% title('两种方案对比')

