function [S]=matrix(N,alpha, dt)
%% 2D BTCS matrix

format long

Nx=N;
Nz=N;
dx=(1)./(Nx-1);
dz=(1)./(Nz-1);
dx2=dx^2;
dz2=dz^2;

R = alpha*dt/dx2;

a=zeros(N*7);   % x's
b=zeros(N*7);   % y's
c=zeros(N*7);   % nz's
pos=1;
%% form S
for m = 1: Nx
    for j = 1: Nz       
        
        %% Interior
        if ((m>=1&&m<=Nx))&&(j>=1&&j<=(Nz))
            a(pos) = (m-1)*(Nz)+j;
            b(pos) = (m-1)*(Nz)+j;
            c(pos)=1 + 4*R;
            pos = pos+1;
        end
        
        %% NSEW in stencil
        if j<Nz
            a(pos) = (m-1)*(Nz)+j;
            b(pos) = (m-1)*(Nz)+j+1;
            c(pos)=-R; %% N
            pos = pos+1;
        end
        if j>1
            a(pos) = (m-1)*(Nz)+j;
            b(pos) = (m-1)*(Nz)+j-1;
            c(pos)=-R; %% S
            pos = pos+1;
        end
        if m<Nx
            a(pos) = (m-1)*(Nz)+j;
            b(pos) = m*(Nz)+j;
            c(pos)=-R; %% E
            pos = pos+1;
        end
        if m>1
            a(pos) = (m-1)*(Nz)+j;
            b(pos) = (m-2)*(Nz)+j;
            c(pos)=-R; %% W
            pos = pos+1;
        end

    end
a=a(1:pos-1);
b=b(1:pos-1);
c=c(1:pos-1);
S=sparse(a,b,c);

end
