function X=generate_sample_RLMFBM_MEMORY_new(T,dt,alpha,V)


N=T/dt;
v=randn(1,N);

for j=1:length(V)

    
      i=[1:1:V(j)];    
      % r=alpha-1;  
      % w=((V(j)-i+1).^(r+1)-(V(j)-i).^(r+1)).*dt.^(r+1)./(r+1)./dt/gamma(alpha);
      s=2*alpha-1;
      w=(((V(j)-i+1).^s-(V(j)-i).^s)*dt^(s)/(s*dt)).^(1/2)/gamma(alpha);
    
      mfbm(j)=sum(w(1:V(j)).*v(1:V(j))*sqrt(dt));
      
end

X=mfbm;

end