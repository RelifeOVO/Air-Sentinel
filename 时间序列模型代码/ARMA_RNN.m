%平均移动法
clc;
clear all;
y=[0.35 0.33 0.29 0.19 0.23 0.24 0.39 0.37 0.21 0.21 0.21];
m=length(y);
n=[1,2]; % 自定义
for i=1:length(n)
   for j=1:m-n(i)+1
        yhat{i}(j)=sum(y(j:j+n(i)-1))/n(i);
   end
    y31(i)=yhat{i}(end);
    s(i)=sqrt(mean((y(n(i)+1:m)-yhat{i}(1:end-1)).^2));
end
y31,s


%加权平均
clc;
clear all;
y=[215 197 203 234 194 108 191 241 232 221 196 226 201 219 217 213 203 225 237 188 212 198 219 177 231 199 203];
w=[1/7;3/7;3/7];
m=length(y);n=3;
for i=1:m-n+1
   yhat(i)=y(i:i+n-1)*w;
    
end
yhat;
err=abs(y(n+1:m)-yhat(1:end-1))./y(n+1:m);
T_err=1-sum(yhat(1:end-1))/sum(y(n+1:m));
y1989=yhat(i)/(1-T_err);


%趋势移动平均法
clc;
clear all;
y=[216 199 222 218 217 259 206 230 255 221 214 212 219 224 210 205 186 249 214 228 211 226 219 238 217 205 206];
y=[676 825 774 716 940 1159 1384 1524 1668 1688 1958 2031 2234 2566 2820 3006 3093 3277 3514 3770 4107];
m1=length(y);
n=6;
for i=1:m1-n+1
   yhat1(i)=sum(y(i:i+n-1))/n;
end
yhat1;
m2=length(yhat1);
for i=1:m2-n+1
   yhat2(i)=sum(yhat1(i:i+n-1))/n;
end
yhat2;
plot(1:27,y,'*');
a21=2*yhat1(end)-yhat2(end);
b21=2*(yhat1(end)-yhat2(end))/(n-1);
y1986=a21+b21
y1987=a21+2*b21


%指数平滑法
clc,clear all;
yt=[216 199 222 218 217 259 206 230 255 221 214 212 219 224 210 205 186 249 214 228 211 226 219 238 217 205 206];
yt=[50 52 47 51 49 48 51 40 48 52 51 59];
n=length(yt); 
alpha=[0.2 0.5 0.8];m=length(alpha); 
yhat(1,1:m)=(yt(1)+yt(2))/2; 
for i=2:n 
 yhat(i,:)=alpha*yt(i-1)+(1-alpha).*yhat(i-1,:); 
end
yhat;
y1=yt';
err=sqrt(mean((repmat(y1,1,m)-yhat).^2))
xlswrite('yt',yhat) ;
yhat1988=alpha*yt(n)+(1-alpha).*yhat(n,:)


%三次指数平滑法
clc,clear;
yt=[20.04 20.06 25.72 34.61 51.77 55.92 80.65 131.11 148.58 162.67 232.26];
n=length(yt); 
alpha=0.3; st1_0=mean(yt(1:3)); st2_0=st1_0;st3_0=st1_0; 
st1(1)=alpha*yt(1)+(1-alpha)*st1_0; 
st2(1)=alpha*st1(1)+(1-alpha)*st2_0; 
st3(1)=alpha*st2(1)+(1-alpha)*st3_0; 
for i=2:n 
st1(i)=alpha*yt(i)+(1-alpha)*st1(i-1); 
st2(i)=alpha*st1(i)+(1-alpha)*st2(i-1); 
st3(i)=alpha*st2(i)+(1-alpha)*st3(i-1); 
end 
xlswrite('touzi.xls',[st1',st2',st3']) 
st1=[st1_0,st1];st2=[st2_0,st2];st3=[st3_0,st3]; 
a=3*st1-3*st2+st3; 
b=0.5*alpha/(1-alpha)^2*((6-5*alpha)*st1-2*(5-4*alpha)*st2+(4-3*alpha)*st3); 
c=0.5*alpha^2/(1-alpha)^2*(st1-2*st2+st3); 
yhat=a+b+c; 
xlswrite('touzi.xls',yhat','Sheet1','D1') 
plot(1:n,yt,'*',1:n,yhat(1:n),'O') 
legend('实际值','预测值') 
xishu=[c(n+1),b(n+1),a(n+1)]; 
yhat1990=polyval(xishu,2)


%自适应滤波法
clc,clear;
yt=[217 207.5 215 223 222 221.5 209 213 217 213 217 215];
m=length(yt);k=0.083;
N=12;Terr=10000;
w=ones(1,N)/N;
while abs(Terr)>0.00001
   Terr=[];
   for j=N+1;m-1
        yhat(j)=w*yt(j-1:-1:j-N)';
        err=yt(j)-yhat(j);
        Terr=[Terr,abs(err)];
        w=w+2*k*err*yt(j-1:-1:j-N);
    end
    Terr=max(Terr);
 end
 w,yhat
 
 
%趋势外推预测法——修正指数曲线法(以下三种类似，选S标准误差小的模型)
function chanliang
clc,clear;
global a b k
%yt=[217 207.5 215 223 222 221.5 209 213 217 213 217];
yt=[42.1 47.5 52.7 57.7 62.5 67.1 71.5 75.7 79.8 83.7 87.5 91.1 94.6 97.9 101.1];
n=length(yt);m=n/3;
%值得注意的是，并不是任何一组数据都可以用修正指数曲线拟合。采用前应对数据进行检验，检验方法是看给定数据的逐期增长量的比率是否接近某一常数b
cf=diff(yt);
for i=1:n-2
   bzh(i)=cf(i+1)/cf(i)
end
range=minmax(bzh)     %b的范围
s1=sum(yt(1:m)),s2=sum(yt(m+1:2*m)),s3=sum(yt(2*m+1:end))
b=((s3-s2)/(s2-s1))^(1/m)
a=(s2-s1)*(b-1)/(b*(b^m-1)^2)
k=(s1-a*b*(b^m-1)/(b-1))/m
y=yuce(1:18)
%定义预测函数
function y=yuce(t)
global a b k
y=k+a*b.^t;

Compertz 曲线:  初期增长缓慢，以后逐渐加快。当达到一定程度后，增长率又逐渐下降。
clc,clear
yuce=@(t,a,b,k)k*a.^(b.^t);
y=[42.1 47.5 52.7 57.7 62.5 67.1 71.5 75.7 79.8 83.7 87.5 91.1 94.6 97.9 101.1];
yt=log(y);n=length(yt);m=n/3;
s1=sum(yt(1:m)),s2=sum(yt(m+1:2*m)),s3=sum(yt(2*m+1:end))
b=((s3-s2)/(s2-s1))^(1/m)
a=(s2-s1)*(b-1)/(b*(b^m-1)^2)
k=(s1-a*b*(b^m-1)/(b-1))/m
a=exp(a)
k=exp(k)
y=yuce(1:18,a,b,k)


Logistic 曲线（生长曲线）
clc,clear
yuce=@(t,a,b,k) 1./(k+a*b.^t);
y=[42.1 47.5 52.7 57.7 62.5 67.1 71.5 75.7 79.8 83.7 87.5 91.1 94.6 97.9 101.1];
yt=1./y;n=length(yt);m=n/3;
s1=sum(yt(1:m)),s2=sum(yt(m+1:2*m)),s3=sum(yt(2*m+1:end))
b=((s3-s2)/(s2-s1))^(1/m)
a=(s2-s1)*(b-1)/(b*(b^m-1)^2)
k=(s1-a*b*(b^m-1)/(b-1))/m
y1=yuce(1:18,a,b,k)