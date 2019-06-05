
%% 第一问
syms x1 y1 x2 y2  x3 y3 x4 y4
[x1,y1]=solve('y1^12-(1.5)^6*y1^6-(2.611104485/4)*x1^(-1)*(1.5)^6=0','y1^12-(1.8)^6*y1^6-(-1.061025201/4)*x1^(-1)*(1.8)^6=0');
[x2,y2]=solve('y2^12-(1.9)^6*y2^6-(-1.341402343/4)*x2^(-1)*(1.9)^6=0','y2^12-(2.3)^6*y2^6-(-1.481624726/4)*x2^(-1)*(2.3)^6=0');
[x3,y3]=solve('y3^12-(2.5)^6*y3^6-(-1.370806496/4)*x3^(-1)*(2.5)^6=0','y3^12-(3)^6*y3^6-(-0.912203381/4)*x3^(-1)*(3)^6=0');
[x4,y4]=solve('y4^12-(3.3)^6*y4^6-(-0.589300528/4)*x4^(-1)*(3.3)^6=0','y4^12-(3.8)^6*y4^6-(-0.191850231/4)*x4^(-1)*(3.8)^6=0');
x1=vpa(x1,8);
y1=vpa(y1,8);
x2=vpa(x2,8);
y2=vpa(y2,8);
x3=vpa(x3,8);
y3=vpa(y3,8);
x4=vpa(x4,8);
y4=vpa(y4,8);
A=[x1,x2,x3,x4;y1,y2,y3,y4];
mean(abs(A),2)
% x 代表势能井 y代表范德华半径

%% 第二问兰纳琼斯函数求x0
syms x1 x2
equ1=x1*exp(1.5)+x2*1.5.^(-7)-50.99339386==0;
equ2=x1*exp(1.6)+x2*1.6.^(-7)-29.45603784==0;
[x1,x2]=solve(equ1,equ2);
x1=vpa(x1,8)
x2=vpa(x2,8)
%%  第二问兰纳琼斯拟合
xdata= [1.5
1.6
1.7
1.8
1.9
2
2.1
2.2
2.3
2.4
2.5
2.6
2.7
2.8
2.9
3
3.1
3.2
3.3
3.4
3.5
3.6
3.7
3.8
3.9
4
4.1
]
ydata= [50.99339386
29.45603784
16.14900065
8.149689578
3.500426677
0.922380883
0.697720557
-0.077584849
-0.770026947
-1.103052932
-1.461661639
-1.673505428
-1.777123806
-1.897867831
-2.007891984
-2.132391571
-2.205109185
-2.159555125
-1.995174244
-1.790830878
-1.624578884
-1.494768094
-1.399044847
-1.302741441
-1.168752552
-0.985335632
-0.72218956
];
x0=[-1.4288937,980.68585];
[x,resnorm,residual,exitflag,output] = lsqcurvefit(@nihe,x0,xdata,ydata)
syms a b
equ1=a^(-1)*(6808.6)+48*b^12==0;
equ2=a^(-1)*(300.41)-24*b^6==0;
[a,b]=solve(equ1,equ2);
a=abs(vpa(a,8))
b=abs(vpa(b,8))
%a 是势能阱 b表示范德华
r=[1.5 1.8 1.9 2.3 2.5 3 3.3 3.8];
V1=arrayfun( @(r) 4*  1.104560997610864792051188487676*(  (  1.4987144912808439110738219226071/r)^12-(   1.4987144912808439110738219226071/r)^6) ,r );
V2=[2.611104485 -1.061025201 -1.341402343 -1.481624726 -1.370806496 -0.912203381 -0.589300528 -0.191850231];
subplot(1,2,1);
plot(xdata,nihe(x,xdata),'k','LineWidth',1);
hold on;
plot(xdata,ydata,'r-','LineWidth',1);
title('力与距离'); 
xlabel('距离(?)');

ylabel('力(eV/?)');
grid on;
legend('拟合图像','原始数据');
subplot(1,2,2);
plot(r,V1,'k','LineWidth',1);
hold on;
plot(r,V2,'r','LineWidth',1);
xlabel('距离(?)');
ylabel('势能(eV)');
grid on;
legend('拟合图像','原始数据');
title('势能与距离'); 
%% 第二问兰纳琼斯修正函数计算x0
syms x1 x2 x3
equ1=x1*(1.5)^(-13)+x2*(1.5)^(-7)+x3*(1.5)^(-5)-(50.99339386)==0;
equ2=x1*(1.6)^(-13)+x2*(1.6)^(-7)+x3*(1.6)^(-5)-(29.45603784)==0;
equ3=x1*(1.7)^(-13)+x2*(1.7)^(-7)+x3*(1.7)^(-5)-(16.14900065)==0;
[x1,x2,x3]=solve(equ1,equ2,equ3);
x1=vpa(x1,8)
x2=vpa(x2,8)
x3=vpa(x3,8)
%% 第二问兰纳琼斯修正函数优化
xdata= [1.5
1.6
1.7
1.8
1.9
2
2.1
2.2
2.3
2.4
2.5
2.6
2.7
2.8
2.9
3
3.1
3.2
3.3
3.4
3.5
3.6
3.7
3.8
3.9
4
4.1
]
ydata= [50.99339386
29.45603784
16.14900065
8.149689578
3.500426677
0.922380883
0.697720557
-0.077584849
-0.770026947
-1.103052932
-1.461661639
-1.673505428
-1.777123806
-1.897867831
-2.007891984
-2.132391571
-2.205109185
-2.159555125
-1.995174244
-1.790830878
-1.624578884
-1.494768094
-1.399044847
-1.302741441
-1.168752552
-0.985335632
-0.72218956
];
x0=[ -3574.4437, 2501.0917,-584.89598];
[x,resnorm,residual,exitflag,output] = lsqcurvefit(@nihebeta,x0,xdata,ydata)
syms a b c ;
equ1=-48*b^(12)+3435.7*a^(-1)==0;
equ2=24*b^(6)-2428.0*a^(-1)==0;
equ3=16*c+558.59*a^(-1)==0;
[a,b,c]=solve(equ1,equ2,equ3);
a=abs(vpa(a,8))
b=abs(vpa(b,8))
c=abs(vpa(c,8))
r=[1.5 1.8 1.9 2.3 2.5 3 3.3 3.8];
V1=arrayfun(@(r) 4* 142.98842545429849337779160123318*((    0.94396543650484804544479955467617/(r+0.8))^12- (  0.94396543650484804544479955467617/(r+0.8))^6-   0.24415874843770779545337745730649/((r+0.8)^4) ),r);
V2=[2.611104485 -1.061025201 -1.341402343 -1.481624726 -1.370806496 -0.912203381 -0.589300528 -0.191850231];
subplot(1,2,1)
plot(r,V1,'k','LineWidth',1);
hold on;
plot(r,V2,'r','LineWidth',1);
title('势能与距离');  
xlabel('距离(埃米)');
ylabel('势能(eV)');
grid on;
legend('拟合图像','原始数据');
subplot(1,2,2)
plot(xdata,nihebeta(x,xdata),'k','LineWidth',1);
hold on;
plot(xdata,ydata,'r','LineWidth',1);
xlabel('距离(埃米)');
ylabel('力(eV/埃米)');
grid on;
legend('拟合图像','原始数据');
title('力与距离'); 
%% 第二问  Buckingham函数拟合
xdata= [1.5
1.6
1.7
1.8
1.9
2
2.1
2.2
2.3
2.4
2.5
2.6
2.7
2.8
2.9
3
3.1
3.2
3.3
3.4
3.5
3.6
3.7
3.8
3.9
4
4.1
]
ydata= [50.99339386
29.45603784
16.14900065
8.149689578
3.500426677
0.922380883
0.697720557
-0.077584849
-0.770026947
-1.103052932
-1.461661639
-1.673505428
-1.777123806
-1.897867831
-2.007891984
-2.132391571
-2.205109185
-2.159555125
-1.995174244
-1.790830878
-1.624578884
-1.494768094
-1.399044847
-1.302741441
-1.168752552
-0.985335632
-0.72218956
];
% x0=[-1177.0728,-74.294143,1722.8966,-981.97255, 282.28737];

[x,resnorm,residual,exitflag,output] = lsqcurvefit(@nihe01,x0 , xdata,ydata)


  syms a b c ;%a 是A b是 phai c是c
  equ1=(119.238222017794)*b+a==0;
 equ2=(-322.353558882783)*b^2-a==0;
 equ3=-(1113.18395448538)+6*c==0;
 [a,b,c]=solve(equ1,equ2,equ3)
 a=abs(vpa(a,8))
 b=abs(vpa(b,8))
 c=abs(vpa(c,8))
 r=[1.5:0.1:4.1];
 y=@(r)  (44.106085377933482050138991326094*exp(-1*r/ 0.36989888503496387572866410664574)- 185.53065908089664048929989803582*r.^(-6))*(1)
m=[1.5 1.8 1.9 2.3 2.5 3 3.3 3.8];
n=[2.611104485 -1.061025201 -1.341402343 -1.481624726 -1.370806496 -0.912203381 -0.589300528 -0.191850231];
subplot(1,2,1);
plot(xdata,nihe01(x,xdata),'k','Linewidth',1);
hold on;
plot(xdata,ydata,'r','Linewidth',1);
title('Buckingham力与距离'); 
xlabel('距离(埃米)');
ylabel('力(eV/埃米)');
grid on;
legend('拟合图像','原始数据');
subplot(1,2,2);
plot(r,y(r),'k','Linewidth',1);
hold on;
plot(m,n,'r','Linewidth',1);
xlabel('距离(埃米)');
ylabel('势能(eV)');
grid on;
legend('拟合图像','原始数据');
title('Buckingham势能与距离'); 
 
%% 第三问
syms d  ;
f1=1*(8*hanshu(sqrt(3)*d/2));
f2=8*(hanshu(sqrt(3)*d/2)+3*hanshu(d)+3*hanshu(sqrt(2)*d)+hanshu(sqrt(3)*d));
f=f1+f2;
df1=diff(f,d);
subplot(2,2,1);
ezplot(df1);%ezplot用来plot未知变量
title('体心立方势能导数函数'); 
xlabel('晶胞参数d(?)');
ylabel('势能导数V(r）’');
grid on;
subplot(2,2,2);
ezplot(f);
title('体心立方势能函数'); 
xlabel('晶胞参数(?)');
ylabel('V(r）(eV/?)');
grid on;
f3=6*(8*hanshu(sqrt(3)*d/2)+4*hanshu(sqrt(5)*d/2)+1*hanshu(d)  );
f4=8*(3*hanshu(sqrt(2)*d/2)+3*hanshu(d)+1*hanshu(sqrt(2)*d)+3*hanshu(sqrt(3)*d)+3*hanshu(sqrt(6)*d/2));
m=f3+f4;
df2=diff(m,d);
subplot(2,2,3);
ezplot(df2);
title('面心立方势能导数函数'); 
xlabel('晶胞参数d(?)');
ylabel('势能导数V(r）’');
grid on;
subplot(2,2,4);
ezplot(m);
title('面心立方势能函数'); 
xlabel('晶胞参数(?)');
ylabel('V(r）(eV/?)');
grid on;
%% 第四问
syms x y z a b c ;
r0=1.17;
n=[0:10];%长/宽/高的铁原子个数为(n+1)个
x=n*2*r0/sqrt(3);
y=n*2*r0/sqrt(3);%设探针原子在铁块上平面的中心
z=4*n*r0/sqrt(3)+2;%探针原子与铁块上平面的中心的高度距离
sum1=0;
sum2=0;
sum3=0;
a=0;
b=0;
c=0;
d1=sqrt((x-4*a*r0/sqrt(3)).^2+(y-4*b*r0/sqrt(3)).^2+(z-4*c*r0/sqrt(3)).^2);%边角原子距离
d2=sqrt( (x-2*r0/sqrt(3)-4*a*r0/sqrt(3) ).^2  +  (y-2*r0/sqrt(3)-4*b*r0/sqrt(3)).^2  +  (z-2*r0/sqrt(3)-4*c*r0/sqrt(3) ).^2   );
F1=lihanshu(d1);%边角原子力
F2=lihanshu(d2);%中心原子力
thetax1=(x-4*a*r0/sqrt(3))/d1;%边角原子沿着x轴
thetay1=(y-4*b*r0/sqrt(3))/d1;%边角原子沿着y轴
thetaz1=(z-4*c*r0/sqrt(3))/d1;%边角原子沿着z轴
thetax2=(x-2*r0/sqrt(3)-4*a*r0/sqrt(3))/d2;%中心原子沿着x轴
thetay2=(y-2*r0/sqrt(3)-4*b*r0/sqrt(3))/d2;%中心原子沿着y轴
thetaz2=(z-2*r0/sqrt(3)-4*c*r0/sqrt(3))/d2;%中心原子沿着z轴
Fx=F1*thetax1+F2*thetax2;
Fy=F1*thetay1+F2*thetay2;
Fz=F1*thetaz1+F2*thetaz2;
for c=0:n
        for b=0:n
            for a=0:n
               sum1=sum1+Fx;
            end
        end
end
for c=0:n
        for b=0:n
            for a=0:n
               sum2=sum2+Fy;
            end
        end
end
for c=0:n
        for b=0:n
            for a=0:n
               sum3=sum3+Fz;
            end
        end
end
sum1;
sum2;
sum3;
sum=sqrt(sum1.^2+sum2.^2+sum3.^2);
subplot(2,3,1)
plot(n,sum1,'r','LineWidth',1);
xlabel('n晶胞数量');
ylabel('力(eV/?10^(-10)m)');
title('Fx-n的变化函数');
subplot(2,3,2)
plot(n,sum2,'b','LineWidth',1);
xlabel('n晶胞数量');
ylabel('力(eV/10^(-10)m)');
title('Fy-n的变化函数');
subplot(2,3,3)
plot(n,sum3,'k','LineWidth',1);
xlabel('n晶胞数量');
ylabel('力(eV/10^(-10)m)');
title('Fz-n的变化函数');
subplot(2,3,4)
plot(n,sum3,'k','LineWidth',1);
xlabel('n晶胞数量');
ylabel('力(eV/10^(-10)m)');
hold on;
plot(n,sum1,'r','LineWidth',1);
xlabel('n晶胞数量');
ylabel('力(eV/10^(-10)m)');
title('Fx/Fy-Fz的比较');
subplot(2,3,5)
plot(n,sum,'k','LineWidth',1);
xlabel('n晶胞数量');
ylabel('力(eV/10^(-10)m)');
title('合力的模的大小');

