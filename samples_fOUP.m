%% fractional Ornstein¨CUhlenbeck process

function [y]=samples_fOUP(N,dt,lambda,sigma,T,H,x0)

% y(1)=normrnd(0,sqrt(gamma(2*H+1)/2));
y(1)=x0;

dBH=fGN(T,H,dt);

for i=1:N
    

 y(i+1)=(1-lambda.*dt).*y(i)+sigma*dBH(i);
          
end

y=y(2:end);

end