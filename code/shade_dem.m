function shadows = shade_dem(dem, s)
% shade_dem calculates DEM shading according to Corripio (2003, IJGIS).
%
%   shadows = shade_dem(dem, s)
%
% This code was translated from the Fortran F90 code in the git
% repo https://github.com/cran/insol/blob/master/src/doshade.f90 and R
% package https://meteoexploration.com/R/insol/insol.pdf. Original
% reference: https://doi.org/10.1080/713811744

rows=size(dem.Z,1);
cols=size(dem.Z,2);

inversesunvector=-s/max(abs(s(1:2)));
normalsunvector(3)=sqrt(s(1)^2 + s(2)^2);
normalsunvector(1)=-s(1)*s(3)/normalsunvector(3);
normalsunvector(2)=-s(2)*s(3)/normalsunvector(3);

z=dem.Z;
l=dem.l;

if s(1)<0
    f_i=1;
else
    f_i=cols;
end

if s(2)<0
    f_j=1;
else
    f_j=rows;
end
shadows=ones(rows,cols);

j=f_j;
for i=1:cols
    n=0;
    zcompare=-inf;
    
    idx=round(i);
    idy=round(j);
    
    dx=0;
    dy=0;
    
    while idx>=1 && idx<=cols && idy>=1 && idy<=rows
       vectortoorigin(1)=dx*l;
       vectortoorigin(2)=dy*l;
       vectortoorigin(3)=z(idy,idx);
       zprojection=dot(vectortoorigin,normalsunvector);
       if zprojection<zcompare
           shadows(idy,idx)=0;
       else
           zcompare=zprojection;
       end
       n=n+1;
       
        idx=round(i+dx);
        idy=round(j+dy);
        dx=inversesunvector(1)*n;
        dy=inversesunvector(2)*n;
    end
end

i=f_i;
for j=1:rows
    n=0;
    zcompare=-inf;
    
    dx=0;
    dy=0;
    idx=i;
    idy=j;
    while idx>=1 && idx<=cols && idy>=1 && idy<=rows

    
       vectortoorigin(1)=dx*l;
       vectortoorigin(2)=dy*l;
       vectortoorigin(3)=z(idy,idx);
       
       zprojection=dot(vectortoorigin,normalsunvector);
       if zprojection<zcompare
           shadows(idy,idx)=0;
       else
           zcompare=zprojection;
       end
       n=n+1;
       dx=inversesunvector(1)*n;
        dy=inversesunvector(2)*n;
    
    idx=round(i+dx);
    idy=round(j+dy);
    end
end