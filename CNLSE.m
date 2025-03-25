function [Uout,Vout]=CNLSE(betai,beta3i,gammai,dbetai,zi,nzi,Ti,Uini,Vini)
A=2/3;
B=1/3;
Uin=Uini;
Vin=Vini;

beta=betai;
beta3=beta3i;
gamma=gammai;
dbeta=dbetai;


z=zi;
nz=nzi;
dz=z/nz;


T=Ti;
n=length(Uin);
dt=T/n;
t=((1:n)'-(n+1)/2)*dt;
w=2*pi*[(0:n/2-1),(-n/2:-1)]'/(dt*n);


uf=fft(Uin);
vf=fft(Vin);


uhs=(j*0.5*dbeta+j*0.5*beta*w.^2-(1/6)*j*beta3*w.^3);
vhs=(-j*0.5*dbeta+j*0.5*beta*w.^2-(1/6)*j*beta3*w.^3);


uhs=exp(uhs*0.5*dz);
vhs=exp(vhs*0.5*dz);


for xz=1:nz
    
if xz==1
    uh=ifft(uhs.*uf);
    vh=ifft(vhs.*vf);
else
    

uhs=(j*0.5*dbeta+j*0.5*beta*w.^2-(1/6)*j*beta3*w.^3);
vhs=(-j*0.5*dbeta+j*0.5*beta*w.^2-(1/6)*j*beta3*w.^3);


uhs=exp(uhs*0.5*dz);
vhs=exp(vhs*0.5*dz);
    


    uf=fft(ur);
    vf=fft(vr);
    uh=ifft(uhs.*uf);
    vh=ifft(vhs.*vf);
end;

    k1=j*gamma*(abs(uh).^2+A*abs(vh).^2).*uh+j*B*gamma*conj(uh).*((vh).^2);
    l1=j*gamma*(abs(vh).^2+A*abs(uh).^2).*vh+j*B*gamma*conj(vh).*((uh).^2);
    uh1=uh+dz*0.5*k1;
    vh1=vh+dz*0.5*l1;
    k2=j*gamma*(abs(uh1).^2+A*abs(vh1).^2).*uh1+j*B*gamma*conj(uh1).*((vh1).^2);
    l2=j*gamma*(abs(vh1).^2+A*abs(uh1).^2).*vh1+j*B*gamma*conj(vh1).*((uh1).^2);
    uh2=uh+dz*0.5*k2;
    vh2=vh+dz*0.5*l2;
    k3=j*gamma*(abs(uh2).^2+A*abs(vh2).^2).*uh2+j*B*gamma*conj(uh2).*((vh2).^2);
    l3=j*gamma*(abs(vh2).^2+A*abs(uh2).^2).*vh2+j*B*gamma*conj(vh2).*((uh2).^2);
    uh3=uh+dz*k3;
    vh3=vh+dz*l3;
    k4=j*gamma*(abs(uh3).^2+A*abs(vh3).^2).*uh3+j*B*gamma*conj(uh3).*((vh3).^2);
    l4=j*gamma*(abs(vh3).^2+A*abs(uh3).^2).*vh3+j*B*gamma*conj(vh3).*((uh3).^2);
    uhh=zeros(1,n);
    vhh=zeros(1,n);
    uhh=uh+(dz/6)*(k1+2*k2+2*k3+k4);
    vhh=vh+(dz/6)*(l1+2*l2+2*l3+l4);
    
    
    uf=fft(uhh);
    vf=fft(vhh);
    uh=ifft(uhs.*uf);
    vh=ifft(vhs.*vf);
    ur=uh;
    vr=vh;
end;


Uout=uh;
Vout=vh;

              