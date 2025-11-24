function y=generate_levy_walk(T,dt,a,t0,num,x0)

beta=a;
N=T/dt;
u=unifrnd(0,1);
v(1)=1*(u<0.5)-1*(u>=0.5);
% v(1)=1;
n=1;


while n<=N+1

if num==2
u=unifrnd(0,1);
tau=t0*((u)^(-1/beta)-1);
n_tau=round(tau/dt);
v(n:n+n_tau)=v(n);
n=n+1+n_tau;
v(n)=-v(n-1);
end

if num==1
tau=exprnd(a);
n_tau=round(tau/dt);
v(n:n+n_tau)=v(n);
n=n+1+n_tau;
v(n)=-v(n-1);
end

end

y=x0+cumsum(v(1:N)*dt);





end