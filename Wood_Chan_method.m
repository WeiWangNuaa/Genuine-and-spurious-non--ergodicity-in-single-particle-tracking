function FBM=Wood_Chan_method(N,H,h)


     
v=ceil(log(N)/log(2));
m=2*2^v;
% Ïà¹Øº¯Êý¾ØÕóc  gamma(t-s)=1/2*D^2*{(t-s-1)^2+(t-s+1)^2-2*(t-s)^2}
c=zeros(1,m);

for j=0:m-1
    
    if j<=m/2
     
     c(j+1)=1/2*(abs(j-1)^(2*H)+abs(j+1)^(2*H)-2*abs(j)^(2*H));
    else
        
     c(j+1)=1/2*(abs(m-j-1)^(2*H)+abs(m-j+1)^(2*H)-2*abs(m-j)^(2*H)); 
    end
    
end


lamda=fft(c);


U=normrnd(0,1,1,m/2);
V=normrnd(0,1,1,m/2);


for j=0:m/2
    if j==0
        a(j+1)=sqrt(lamda(j+1))*randn/sqrt(m);
    end
    if j==m/2
        a(j+1)=sqrt(lamda(j+1))*randn/sqrt(m);
    end
    if j>=1 && j<m/2
        a(j+1)=sqrt(lamda(j+1))/sqrt(2*m)*(complex(U(j+1),V(j+1)));
        a(m-j+1)=sqrt(lamda(j+1))/sqrt(2*m)*(complex(U(j+1),-V(j+1)));
    end
end



g=fft(a);

FGN=real(g(1:N));

FBM(1)=FGN(1);

for i=1:N
    if i>1
       FBM(i)=sum(real(FGN(1:i)));  
    end
   
end
FBM=(h).^(H)*FBM;


end
