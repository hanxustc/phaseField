clear all;
format long;

load theta.out;
load theta12_num.out;

Nx=2^7;
Ny=2^7;
k=11;

theta0=theta(1:Nx*Ny);
theta_plot=theta((k-1)*Nx*Ny+1:k*Nx*Ny);
theta1=theta12_num((k-2)*Nx*Ny+1:(k-1)*Nx*Ny,1);
theta2=theta12_num((k-2)*Nx*Ny+1:(k-1)*Nx*Ny,2);
dtheta12=theta2-theta1;

%Initial%
fid1=fopen('theta0_matrix.out','w');
ii=0;
for i=1:Nx;
    for j=1:Ny
        ii=ii+1;
        theta0_matrix(i,j)=theta0(ii);
        fprintf(fid1,'%10.6f ',theta0_matrix(i,j));
    end
    fprintf(fid1,'\n');
end

%Errors%
fid2=fopen('errors.out','w');
for i=1:k-1
    ttheta1=0;
    ttheta2=0;
    for j=1:Nx*Ny
        ttheta1=ttheta1+abs(theta(Nx*Ny*(i-1)+j));
        ttheta2=ttheta2+abs(theta(Nx*Ny*i+j));
    end
    ttheta1=ttheta1/(Nx*Ny);
    ttheta2=ttheta2/(Nx*Ny);
    dtheta(i)=abs(ttheta2-ttheta1);
    fprintf(fid2,'%10.6f\n',dtheta(i));
end

%Contour plot%
fid3=fopen('theta_matrix.out','w');
ii=0;
for i=1:Nx;
    for j=1:Ny
        ii=ii+1;
        theta_matrix(i,j)=theta_plot(ii);
        fprintf(fid3,'%10.6f ',theta_matrix(i,j));
    end
    fprintf(fid3,'\n');
end

%theta12%
fid4=fopen('dtheta12_matrix.out','w');
ii=0;
for i=1:Nx;
    for j=1:Ny
        ii=ii+1;
        dtheta12_matrix(i,j)=dtheta12(ii);
        fprintf(fid4,'%10.6f ',dtheta12_matrix(i,j));
    end
    fprintf(fid4,'\n');
end