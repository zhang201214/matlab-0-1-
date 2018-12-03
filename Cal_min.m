function [min_v]=Cal_min(r1,r2,theta,tf)
% r1=1;
% r2=1;
% theta=pi/3;
mu=4*pi*pi;
% Nmax=fix(tfinal/Tg);
r1_vec = [r1 ;0 ;0];    %%服务卫星P1出发点向量
vc1_vec = [0;sqrt(mu./r1);0];  %%起始点速度
T1 = (2*pi*r1^(1.5))/sqrt(mu);      %%服务卫星轨道周期
T2 = (2*pi*r2^(1.5))/sqrt(mu);      %%目标卫星轨道周期（用于计算最终交汇点）
theta_true = theta + tf*2*pi/T2;            %%实际真近点角差值，备用 
theta= rem(theta_true,2*pi);  
tol = 1*10e-8;
if theta==0 || theta== 2*pi || theta== pi 
     theta = rem(theta+tol,2*pi);
end  
r2_vec = [r2*cos(theta) ;r2*sin(theta); 0];  %%服务卫星P1与目标卫星P2期望交汇点向量         
vc2_vec = [sqrt(mu/r2)*cos(theta+pi/2); sqrt(mu/r2)*sin(theta+pi/2); 0];  %%终点速度      
     
lw = (theta<pi)*0+(theta>=pi)*1;  % 计算lw的值,theta<180, lw = 0; theta>=180,lw = 1 
[Nmax,gt] = ComputeNmax(r1, r2, theta, tf,mu);  %计算允许的最大转移圈数  gt =1 表示 tf>tmNmax ；gt =0 表示 tf<=tmNmax
        
deta_v = 50*ones(2*Nmax+1,1); %用于保存2*Nmax+1个解

for k = 0 : Nmax      
            if(k==0)  %0圈只有一个解
                if(Nmax==0) 
                    [v1_jie,v2_jie] = solve_lambert(r1_vec,r2_vec,tf,mu,lw,0,gt); % Nmax=0 ，当tf>tmNmax,0圈Lambert的解取上肢；当tf<=tmNmax,0圈Lambert的解取下肢；
                    delta_v1(2*k+1) = norm(v1_jie-vc1_vec);
                    delta_v2(2*k+1) = norm(vc2_vec-v2_jie);
                    deta_v(2*k+1) = delta_v1(2*k+1) + delta_v2(2*k+1);          % 计算速度增量
                else 
                    [v1_jie,v2_jie] = solve_lambert(r1_vec,r2_vec,tf,mu,lw,0,1); % Nmax>0时，0圈Lambert的解取上肢
                    delta_v1(2*k+1) = norm(v1_jie-vc1_vec);
                    delta_v2(2*k+1) = norm(vc2_vec-v2_jie);
                    deta_v(2*k+1) = delta_v1(2*k+1) + delta_v2(2*k+1);          % 计算速度增量
                end 
            else   %多圈每一个k对应两个解
                [v1_jie,v2_jie] = solve_lambert(r1_vec,r2_vec,tf,mu,lw,k,1);  %解在上肢 
                delta_v1(2*k) = norm(v1_jie-vc1_vec);
                delta_v2(2*k) = norm(vc2_vec-v2_jie);
                deta_v(2*k) = delta_v1(2*k) + delta_v2(2*k);          % 计算速度增量
                
                [v1_jie,v2_jie] = solve_lambert(r1_vec,r2_vec,tf,mu,lw,k,0);  %解在下肢 
                delta_v1(2*k+1) = norm(v1_jie-vc1_vec);
                delta_v2(2*k+1) = norm(vc2_vec-v2_jie);
                deta_v(2*k+1) = delta_v1(2*k+1) + delta_v2(2*k+1);          % 计算速度增量
            end
end 
min_v=min(deta_v);
end