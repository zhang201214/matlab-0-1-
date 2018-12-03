r1=1;
r2=1;
  theta=[20,20,30,-120,-40];
 %theta=[20,-20,35,150,-185,-15];
  theta=theta*pi/180;
%  theta=-51.4*pi/180;
miu=4*pi*pi;
miu_true=398600.436;
%对于GEO轨道而言
r1_true=42164.169637;
 k1=sqrt(miu_true/r1_true)/sqrt(miu/r1);
 k2=((2*pi*r1_true^(1.5))/sqrt(miu_true))/((2*pi*r1^(1.5))/sqrt(miu));
for j=1:size(theta,2)
d=sqrt(r1^2+r2^2-2*cos(theta(j))*r1*r2);
s=(r1+r2+d)/2;
am=s/2;
arfa1=2*asin(sqrt(s./(2*am)));
beta1=2*asin(sqrt((s-d)./(2*am)));
tm=sqrt(am.^3/miu).*((arfa1-sin(arfa1))-(beta1-sin(beta1)));
t=tm:0.01:13;
m=size(t,2);
for i=1:m
    deltav(i)=Cal_min(r1,r2,theta(j),t(i));
end
deltav1=deltav;
for i=1:m-1      %加入末端滑行的情况                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    if (deltav1(i)<deltav1(i+1))
        deltav1(i+1)=deltav1(i);
    end
end    
 plot(t,deltav);
 hold on;
  plot(t*k2/3600,deltav1*k1*1000);
%  plot(t,deltav1,'--');
% hold on;
%求a，b，c
cn=1;
cm=1;
for i=1:m-2
    if deltav1(i)>deltav(i+1) && deltav1(i+1)==deltav1(i+2)
        a(j,cn)=t(i+1);
        c(j,cn)=deltav1(i+1);
        cn=cn+1;
          plot(t(i+1)*k2/3600,deltav1(i+1)*k1*1000,'o');
%         plot(t(i+1),deltav1(i+1),'o');
    elseif deltav1(i)==deltav1(i+1) && deltav1(i+1)>deltav1(i+2)
        b(j,cm)=t(i+1);
          plot(b(j,cm)*k2/3600,deltav1(i+1)*k1*1000,'+');
%         plot(b(j,cm),deltav1(i+1),'+');
        cm=cm+1;   
    else
    end   
end 

  for p=1:size(b,1)          %用于tf与deltav的关系中，无视下降段，并将恒稳段连起来的情况
      a(p,1)=tm;
      for i=2:size(b,2)
          a(p,i)=b(p,i-1);
      end
  end
clear deltav;
clear deltav1;
end
  c=c;  
  a=a;   
  b=b;   