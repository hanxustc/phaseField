function theta0=Griding(Nx,Ny,dx,dy)
format long;

%Read file%
in=fopen('grains.inp','r');
epsilon=1.0e-1;
ndime=2;
npoint=fscanf(in,'%d',1);
nnode=fscanf(in,'%d',1);
nelem=fscanf(in,'%d',1);
grainNum=fscanf(in,'%d',1);

if (grainNum==5)
	theta_num=[-0.6 -0.5 -0.25 0.0 0.15];
end
if (grainNum==6)
	theta_num=[-0.6 -0.45 -0.3 -0.15 0.0 0.15];
end

for ipoint=1:npoint
    jpoint=fscanf(in,'%d',1);
    dummy=fscanf(in,'%lf %lf',[2,1]);
    for idime=1:ndime
        vnode(jpoint,idime)=dummy(idime);
    end
end

for ielem=1:nelem
    jelem=fscanf(in,'%d',1);
    dummy=fscanf(in,'%d',[nnode+1,1]);
    for inode=1:nnode+1
        velem(jelem,inode)=dummy(inode);
    end
end

for ielem=1:nelem
    jnode=0;
    for inode=1:nnode
        knode=velem(ielem,inode);
        if(knode~=0)
            jnode=jnode+1;
        end
    end
    nnode2(ielem)=jnode;
end

%Griding%
vgrid=zeros(Nx*Ny,2);
for i=1:Nx
    for j=1:Ny
        ii=(i-1)*Ny+j;
        vgrid(ii,1)=i*dx-dx/2;
        vgrid(ii,2)=j*dy-dy/2;
    end
end

%Distribute theta%
theta0=zeros(Nx*Ny,1);
for i=1:Nx
    for j=1:Ny
        ii=(i-1)*Ny+j;
        for ielem=1:nelem
            theta=0.0;
            igrain=velem(ielem,nnode+1);
            mnode=nnode2(ielem);
            for inode=1:mnode
                knode=velem(ielem,inode);
                xv1=vnode(knode,1);
                yv1=vnode(knode,2);
                jnode=velem(ielem,inode+1);
                if (inode==mnode)
                    jnode=velem(ielem,1);
                end
                xv2=vnode(jnode,1);
                yv2=vnode(jnode,2);
                p1x=xv1-vgrid(ii,1);
                p1y=yv1-vgrid(ii,2);
                p2x=xv2-vgrid(ii,1);
                p2y=yv2-vgrid(ii,2);
                x1=sqrt(p1x*p1x+p1y*p1y);
                x2=sqrt(p2x*p2x+p2y*p2y);
                if (x1*x2<=epsilon)
                    theta=2*pi;
                else
                    cos_theta=((p1x*p2x+p1y*p2y)/(x1*x2));
                    theta=theta+acos(cos_theta);
                end
            end
            if (abs(theta-2*pi)<=epsilon)
                theta0(ii)=theta_num(igrain);
            end
        end
    end
end

% %Generate plot%
% grid_out=fopen('grid.out','w');
% for ii=1:Nx*Ny
%     if (vgrid(ii,1)>Nx*dx)
%         vgrid(ii,1)=vgrid(ii,1)-Nx*dx;
%     end
%     if (vgrid(ii,2)>Ny*dy)
%         vgrid(ii,2)=vgrid(ii,2)-Ny*dy;
%     end
% 	fprintf(grid_out,'%14.6e %14.6e %14.6e\n',vgrid(ii,1),vgrid(ii,2),theta0(ii));
% end
% fclose(grid_out);
end