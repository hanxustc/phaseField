clear all;

npoint=6;
Nx=2^7;
Ny=2^7;
dx=1.0;
dy=1.0;

a=1;b=2;c=3;d=4;
while (a~=b||b~=c||c~=d||d~=a)
    Voronoi(npoint,Nx,Ny);
    theta=Griding(Nx,Ny,dx,dy);
    a=theta(1);
    b=theta(Ny);
    c=theta((Nx-1)*Ny+1);
    d=theta(Nx*Ny);
end

fid1=fopen('poly_matrix.in','w');
ii=0;
for i=1:Nx;
    for j=1:Ny
        ii=ii+1;
        fprintf(fid1,'%10.6f ',theta(ii));
    end
    fprintf(fid1,'\n');
end