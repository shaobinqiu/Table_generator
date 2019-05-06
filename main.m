clear
cell=import_poscar('/home/qiusb/Documents/str/diamond_pcell.vasp');
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
    lattice=cell.lattice;
 for ii=1:9
     ii
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
     else
         break
     end
 end
 close all
 figure
 plot3(test_poits(:,1),test_poits(:,2),test_poits(:,3))
 hold on
 plot3( lattice_test_points(K,1), lattice_test_points(K,2), lattice_test_points(K,3),'*')