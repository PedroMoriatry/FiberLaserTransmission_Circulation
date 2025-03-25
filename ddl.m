function [Oout,Uout_Ge]=ddl(Oin,beta2,beta3,z,nz,T,ni)
Uout_Ge=zeros(ni,nz);
dz=z/nz;
nt = length(Oin);
dt=T/nt;
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt);
uf=fft(Oin);
hs =j*beta2*w.^2/2-j*beta3*w.^3/6;
hs = exp(hs*dz);
for xz=1:nz
    uh=ifft(hs.*uf);
    uf=fft(uh);   
    Uout_Ge(:,xz)=uh;
end;
Oout=uh;

    
    
    