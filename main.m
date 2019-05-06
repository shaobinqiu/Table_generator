clear
cell=import_poscar('/home/qiusb/Documents/str/diamond_scell.vasp');
cell.coords=cell.coords-repmat(cell.coords(1,:),size(cell.coords,1),1);
coords=cell.coords*cell.lattice;
s=[];
for ii=1:size(cell.atomcount,1)
    s=[s repmat(ii,1,cell.atomcount(ii,1))];
end
%%%%%%%%%%%%%%%%%
r=0;
for ii=1:size(cell.lattice,1)
    r=max(r,sqrt(sum(cell.lattice(ii,:).^2)));
end
test_poits=[];
for ii=1:100
    for jj=1:100
        test_poits=[test_poits;r*sin(pi/100*ii)*cos(2*pi/100*jj) r*sin(pi/100*ii)*sin(2*pi/100*jj) r*cos(pi/100*ii)];
    end
end
%%%%%%%%%%%%%%%%%%
    lattice=cell.lattice;
 for ii=1:9
     lattice_test_points=[test_poits;lattice];
     k = convhull( lattice_test_points(:,1), lattice_test_points(:,2), lattice_test_points(:,3));
     K=unique(k);
     in=K(K<=size(test_poits,1));
     if size(in,1)>0
         for jj=-ii:ii
             for kk=-ii:ii
                 for ll=-ii:ii
                     lattice=[lattice;jj*cell.lattice(1,:)+kk*cell.lattice(2,:)+ll*cell.lattice(3,:)];
                 end
             end
         end
         lattice=unique(lattice,'rows');
     else
         break
     end
 end
%  close all
%  figure
%  plot3(test_poits(:,1),test_poits(:,2),test_poits(:,3))
%  hold on
%  plot3( lattice_test_points(K,1), lattice_test_points(K,2), lattice_test_points(K,3),'*')
 %%%%%%%%%%%%%%%%
 d_lattice=[];
 for ii=1:size(lattice,1)
     if sum(lattice(ii,:).^2)<r^2+0.5
         d_lattice=[d_lattice ii];
     end
 end
 lattice=lattice(d_lattice,:);
 %%%%%%%%%
 temp=nchoosek(1:size(lattice,1),3);
 lattice_r=[];
 for ii=1:size(temp,1)
     if abs(abs(det(lattice(temp(ii,:),:)))-det(cell.lattice))<0.1
         lattice_r=[lattice_r;lattice(temp(ii,:),:)];
     end
 end
 %%%%%%%%%%
 R_t=[];
 per=perms(1:3);
 for ii=1:size(lattice_r,1)/3
     for jj=1:size(per,1)
         lattice_ij=lattice_r(per(jj,:)+repmat(ii*3-3,1,3),:);
         R_ij=lattice_ij/cell.lattice;
         if sum(sum(abs(R_ij*R_ij'-eye(3))))<0.1 %&& det(R_ij)>0
             R_t=[R_t;lattice_ij/cell.lattice];
         end
     end
 end
 %%%%%%%%%%
 table=[1:size(cell.coords,1)];
R_rt=[];
T_rt=[];
 for ii=1:size(coords,1)
     t=coords(ii,:)-coords(1,:);
     coords_t=coords+repmat(-t,size(coords,1),1);
     for jj=1:size(R_t,1)/3
         coords_rt=coords_t*R_t(jj*3-2:jj*3,:);
         coords_rt_D=coords_rt/cell.lattice+repmat(cell.coords(ii,:),size(coords,1),1);
         n=0;
         table_temp=zeros(1,size(cell.coords,1));
         for kk=1:size(coords_rt_D,1)
             for ll=1:size(cell.coords,1)
                 if sum((coords_rt_D(kk,:)-cell.coords(ll,:)-round(coords_rt_D(kk,:)-cell.coords(ll,:))).^2)<0.1 && s(kk)==s(ll)
                     table_temp(1,ll)=kk;
                     n=n+1;
                 end
             end
         end
         if n==size(cell.coords,1)
             table=[table;table_temp];
             R_rt=[R_rt;R_t(jj*3-2:jj*3,:)];
             T_rt=[T_rt;t/cell.lattice];
         end
     end
 end
 table=unique(table,'rows');
 
 
     
 
 
 
 
 
 