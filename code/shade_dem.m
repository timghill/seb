function shadows=doshade(dem, s)
% function doshade computes DEM shading according to Corripio (2003,
% IJGIS). This code was translated from the Fortran F90 code in the git
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
% SUBROUTINE doshade(dem, sunvector, cols, rows, dl, sombra)
% IMPLICIT NONE
% INTEGER :: cols, rows, newshape(2)
% DOUBLE PRECISION :: z(cols,rows), sombra(cols, rows), dem(cols*rows),sunvector(3)
% INTEGER :: idx, jdy, n, i, j, f_i, f_j, casx, casy
% DOUBLE PRECISION :: inversesunvector(3), normalsunvector(3), vectortoorigin(3)
% DOUBLE PRECISION :: dl, dx, dy, zprojection, zcompare
% inversesunvector = -sunvector/Maxval(ABS(sunvector(1:2)))
% normalsunvector(3)=sqrt(sunvector(1)**2+sunvector(2)**2)
% normalsunvector(1)=-sunvector(1)*sunvector(3)/normalsunvector(3)
% normalsunvector(2)=-sunvector(2)*sunvector(3)/normalsunvector(3)
% newshape(1)=cols
% newshape(2)=rows
% z=reshape(dem,newshape)

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
% !*** casx is an integer, this makes the value large enough to compare effectively
% casx=NINT(1e6*sunvector(1))        
% casy=NINT(1e6*sunvector(2))
% SELECT CASE (casx)
% !******** case (:0) means sunvector(x) negative, 
% ! sun is on the West: beginning of grid cols
% CASE (:0)    
% f_i=1            !** fixed i_value
% CASE default
% f_i=cols
% END SELECT
% SELECT CASE (casy)
% !******** case (:0) sunvector(y) negative, 
% ! Sun is on the North: beginning of grid rows
% CASE (:0)
% f_j=1
% CASE default 
% f_j=rows
% END SELECT
% !******************* Grid scanning *******************************
% !*** the array sombra stores the shaded cells, it is set 
% !*** to 1 before the grid scanning.
% !*****************************************************************
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

        
%         idx;
        
%         idx
    
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

% sombra = 1
% j=f_j
% DO i=1, cols        
%     n = 0
%     zcompare = -HUGE(zcompare) !** initial value lower than any possible zprojection
%     DO 
%         dx=inversesunvector(1)*n
%         dy=inversesunvector(2)*n
%         idx = NINT(i+dx)
%         jdy = NINT(j+dy)
%         IF ((idx < 1) .OR. (idx > cols) .OR. (jdy < 1) .OR. (jdy > rows)) exit
%         vectortoorigin(1) = dx*dl
%         vectortoorigin(2) = dy*dl
%         VectortoOrigin(3) = z(idx,jdy)
%         zprojection = Dot_PRODUCT(vectortoorigin,normalsunvector)
%         IF (zprojection < zcompare) THEN 
%             sombra(idx,jdy) = 0 
%             ELSE
%             zcompare = zprojection
%         END IF  
%         n=n+1
%     END DO 
% END DO
% i=f_i

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
% DO j=1,rows
%     n = 0
%     zcompare = -HUGE(zcompare)  !** initial value lower than any possible zprojection
%     DO 
%         dx=inversesunvector(1)*n    
%         dy=inversesunvector(2)*n
%         idx = NINT(i+dx)
%         jdy = NINT(j+dy)
%         IF ((idx < 1) .OR. (idx > cols) .OR. (jdy < 1) .OR. (jdy > rows)) exit
%         vectortoorigin(1) = dx*dl
%         vectortoorigin(2) = dy*dl
%         VectortoOrigin(3) = z(idx,jdy)
%         zprojection = Dot_PRODUCT(vectortoorigin,normalsunvector)
%         IF (zprojection < zcompare) THEN 
%             sombra(idx,jdy) = 0 
%             ELSE
%             zcompare = zprojection
%         END IF  
%         n=n+1
%     END DO 
% END DO
% END SUBROUTINE doshade
%   