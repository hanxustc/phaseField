function []=Voronoi(npoint,Nx,Ny)
format long;

%Create files%
out=fopen('Voroni_vertices.out','w');
out1=fopen('grains.out','w');
out2=fopen('final_plot.p','w');
out3=fopen('original_points.out','w');
out4=fopen('cell.out','w');
out5=fopen('grains.inp','w');

%Create box%
xmax=Nx;
ymax=Ny;
x0=1;
y0=1;
extra=1.0;

%Generate random points%
x=xmax*rand(npoint,1);
y=ymax*rand(npoint,1);

%Periodic boundary%
for i=1:npoint
    j=npoint+i;
    x(j)=x(i);
    y(j)=y(i)-ymax;
end
for i=1:npoint
    j=2*npoint+i;
    x(j)=x(i)+xmax;
    y(j)=y(i)-ymax;
end
for i=1:npoint
    j=3*npoint+i;
    x(j)=x(i)+xmax;
    y(j)=y(i);
end
for i=1:npoint
    j=4*npoint+i;
    x(j)=x(i)+xmax;
    y(j)=y(i)+ymax;
end
for i=1:npoint
    j=5*npoint+i;
    x(j)=x(i);
    y(j)=y(i)+ymax;
end
for i=1:npoint
    j=6*npoint+i;
    x(j)=x(i)-xmax;
    y(j)=y(i)+ymax;
end
for i=1:npoint
    j=7*npoint+i;
    x(j)=x(i)-xmax;
    y(j)=y(i);
end
for i=1:npoint
    j=8*npoint+i;
    x(j)=x(i)-xmax;
    y(j)=y(i)-ymax;
end

%Generate voronoi diagram%
[c,f]=voronoin([x,y]);
nelem=size(f,1);
ncount=0;

for i=1:nelem
    flag=1;
    vnodes=f{i,:};
    nnode=size(vnodes,2);
    for j=1:nnode
        if(vnodes(j)==1)
            flag=0;
        end
    end
    if(flag==1)
        ncount=ncount+1;
        new_nnode(ncount)=nnode;
        for j=1:nnode
           new_vnode(ncount,j)=vnodes(j);
        end
    end
end

%Clip far outside voronoi elements%
new_nelem=0;
for i=1:ncount
	flag=0;
	for j=1:new_nnode(i)
        kk=new_vnode(i,j);
        if(kk~=0)
            if(c(kk,1)>=-extra&&c(kk,1)<=xmax+extra&&c(kk,2)>=-extra&&c(kk,2)<=ymax+extra)
                flag=1;
            end
        end
	end
    
	if(flag==1)
        new_nelem=new_nelem+1;
        jnode=0;
        for j=1:new_nnode(i)
            kk=new_vnode(i,j);
            if(kk~=0)
                jnode=jnode+1;
                new_vnode2(new_nelem,jnode)=new_vnode(i,j);
            end
        end
        new_nnode2(new_nelem)=jnode;
	end
end

%Assign grain numbers to voronoi elements%
epsilon=1.0e-4;
for isector=1:9
	for ipoint=1:npoint
        jpoint=(isector-1)*npoint+ipoint;
        for ielem=1:new_nelem
            theta=0.0;
            nnode=new_nnode2(ielem);
            for inode=1:nnode
                kk1=new_vnode2(ielem,inode);
                xv1=c(kk1,1);
                yv1=c(kk1,2);
                jnode=inode+1;
                if(inode==nnode)
                    jnode=1;
                end
                kk2=new_vnode2(ielem,jnode);
                xv2=c(kk2,1);
                yv2=c(kk2,2);
                p1x=(xv1-x(jpoint));
                p1y=(yv1-y(jpoint));
                p2x=(xv2-x(jpoint));
                p2y=(yv2-y(jpoint));
                x1=sqrt(p1x*p1x+p1y*p1y);
                x2=sqrt(p2x*p2x+p2y*p2y);
                cos_theta=((p1x*p2x+p1y*p2y)/(x1*x2));
                theta=theta+acos(cos_theta);
            end
            if(abs(theta-2*pi)<=epsilon)
                igrain(ielem)=ipoint;
            end
        end
	end
end

%Print origunal random points%
for i=1:npoint
    fprintf(out3,'%4.6e %14.6e\n',x(i),y(i));
end

%Simulation cell%
fprintf(out4,'%14.6e %14.6e\n',x0,y0);
fprintf(out4,'%14.6e %14.6e\n',xmax,y0);
fprintf(out4,'%14.6e %14.6e\n',xmax,ymax);
fprintf(out4,'%14.6e %14.6e\n',x0,ymax);
fprintf(out4,'%14.6e %14.6e\n',x0,y0);

%Print original voronoi diagram%
for i=1:ncount
	fprintf(out,'# i %d\n',i);
    for j=1:nnode
        kk=new_vnode(i,j);
        if(kk~=0)
            fprintf(out,'%14.6e %14.6e\n',c(kk,1),c(kk,2));
        end
    end
	kk=new_vnode(i,1);
	fprintf(out,'%14.6e %14.6e\n',c(kk,1),c(kk,2));
    fprintf(out,'\n');
end

%Final plot%
fprintf(out2,'set terminal pdfcairo size 50cm,50cm font "Times New Roman,60"\n');
fprintf(out2,'set output "plot.pdf"\n\n');
% fprintf(out2,'set xrange [0:');
% fprintf(out2,'%d]\n',xmax);
% fprintf(out2,'set yrange [0:');
% fprintf(out2,'%d]\n',ymax);
fprintf(out2,'set xtics font "Times New Roman,60"\n');
fprintf(out2,'set ytics font "Times New Roman,60"\n');
fprintf(out2,'set xlabel "x" font "Times New Roman,72"\n');
fprintf(out2,'set ylabel "y" font "Times New Roman,72"\n');
fprintf(out2,'set mxtics 5\n');
fprintf(out2,'set mytics 5\n\n');
fprintf(out2,'set key top right Left reverse font "Times New Roman,48"\n\n');

for i=1:new_nelem
    fprintf(out1,'# i %d %d\n',i,new_nnode2(i));
    nnode=size(new_vnode2,2);
	ncount=0;
	xcod=0.0;
	ycod=0.0;

	for j=1:nnode
        kk=new_vnode2(i,j);
        fprintf(out1,'# %5d %5d\n',j,kk);
        if(kk~=0)
            fprintf(out1,'%14.6e %14.6e\n',c(kk,1),c(kk,2));
            ncount=ncount+1;
            xcod=xcod+c(kk,1);
            ycod=ycod+c(kk,2);
        end
	end
	kk=new_vnode2(i,1);
	fprintf(out1,'%14.6e %14.6e\n',c(kk,1),c(kk,2));
	fprintf(out1,'\n');
	xcod=xcod/ncount;
	ycod=ycod/ncount;
    fprintf(out2,'set label "');
    fprintf(out2,'%d',i);
    fprintf(out2,'" at');
    fprintf(out2,'%14.6e,%14.6e\n',xcod,ycod);
    fprintf(out2,'set label "');
	fprintf(out2,'%d',igrain(i));
	fprintf(out2,'" at');
    fprintf(out2 ,'%14.6e,%14.6e',xcod+5.0,ycod+8.0);
	fprintf(out2,'\n');
end
fprintf(out2,'plot "grains.out" w line linetype 1 lw 6 lc 1, "cell.out" w line linetype 1 lw 2 lc 2, "original_points.out" w p pt 7 ps 3.0');

%Print grains input file%
[nn1,nn2]=size(c);
nnode=size(new_vnode2,2);
fprintf(out5,'%5d %5d %5d %5d\n',(nn1-1),nnode,new_nelem,npoint);
for i=2:nn1
	fprintf(out5,'%5d %14.6e %14.6e\n',i,c(i,1),c(i,2));
end
for i=1:new_nelem
	fprintf(out5,'%5d',i);
	for j=1:nnode
        fprintf(out5,'%5d',new_vnode2(i,j));
	end
    fprintf(out5,'%5d',igrain(i));
	fprintf(out5,'\n');
end
end