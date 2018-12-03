function [min_v]=Cal_min(r1,r2,theta,tf)
% r1=1;
% r2=1;
% theta=pi/3;
mu=4*pi*pi;
% Nmax=fix(tfinal/Tg);
r1_vec = [r1 ;0 ;0];    %%��������P1����������
vc1_vec = [0;sqrt(mu./r1);0];  %%��ʼ���ٶ�
T1 = (2*pi*r1^(1.5))/sqrt(mu);      %%�������ǹ������
T2 = (2*pi*r2^(1.5))/sqrt(mu);      %%Ŀ�����ǹ�����ڣ����ڼ������ս���㣩
theta_true = theta + tf*2*pi/T2;            %%ʵ�������ǲ�ֵ������ 
theta= rem(theta_true,2*pi);  
tol = 1*10e-8;
if theta==0 || theta== 2*pi || theta== pi 
     theta = rem(theta+tol,2*pi);
end  
r2_vec = [r2*cos(theta) ;r2*sin(theta); 0];  %%��������P1��Ŀ������P2�������������         
vc2_vec = [sqrt(mu/r2)*cos(theta+pi/2); sqrt(mu/r2)*sin(theta+pi/2); 0];  %%�յ��ٶ�      
     
lw = (theta<pi)*0+(theta>=pi)*1;  % ����lw��ֵ,theta<180, lw = 0; theta>=180,lw = 1 
[Nmax,gt] = ComputeNmax(r1, r2, theta, tf,mu);  %������������ת��Ȧ��  gt =1 ��ʾ tf>tmNmax ��gt =0 ��ʾ tf<=tmNmax
        
deta_v = 50*ones(2*Nmax+1,1); %���ڱ���2*Nmax+1����

for k = 0 : Nmax      
            if(k==0)  %0Ȧֻ��һ����
                if(Nmax==0) 
                    [v1_jie,v2_jie] = solve_lambert(r1_vec,r2_vec,tf,mu,lw,0,gt); % Nmax=0 ����tf>tmNmax,0ȦLambert�Ľ�ȡ��֫����tf<=tmNmax,0ȦLambert�Ľ�ȡ��֫��
                    delta_v1(2*k+1) = norm(v1_jie-vc1_vec);
                    delta_v2(2*k+1) = norm(vc2_vec-v2_jie);
                    deta_v(2*k+1) = delta_v1(2*k+1) + delta_v2(2*k+1);          % �����ٶ�����
                else 
                    [v1_jie,v2_jie] = solve_lambert(r1_vec,r2_vec,tf,mu,lw,0,1); % Nmax>0ʱ��0ȦLambert�Ľ�ȡ��֫
                    delta_v1(2*k+1) = norm(v1_jie-vc1_vec);
                    delta_v2(2*k+1) = norm(vc2_vec-v2_jie);
                    deta_v(2*k+1) = delta_v1(2*k+1) + delta_v2(2*k+1);          % �����ٶ�����
                end 
            else   %��Ȧÿһ��k��Ӧ������
                [v1_jie,v2_jie] = solve_lambert(r1_vec,r2_vec,tf,mu,lw,k,1);  %������֫ 
                delta_v1(2*k) = norm(v1_jie-vc1_vec);
                delta_v2(2*k) = norm(vc2_vec-v2_jie);
                deta_v(2*k) = delta_v1(2*k) + delta_v2(2*k);          % �����ٶ�����
                
                [v1_jie,v2_jie] = solve_lambert(r1_vec,r2_vec,tf,mu,lw,k,0);  %������֫ 
                delta_v1(2*k+1) = norm(v1_jie-vc1_vec);
                delta_v2(2*k+1) = norm(vc2_vec-v2_jie);
                deta_v(2*k+1) = delta_v1(2*k+1) + delta_v2(2*k+1);          % �����ٶ�����
            end
end 
min_v=min(deta_v);
end