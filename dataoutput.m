% Kesson
if num==0

writeID = fopen( [ 'Res_strain.txt' ] , 'w' ) ;
a=Eps0{1};
b=Eps0{2};
c=Eps0{3};
x=0;
fprintf(writeID,'E=%e Nu=%e \r\n',E,nu);
fprintf(writeID,[num2str(length(voxelx)),'\r\n']);
for i=1:length(Eps0{1})
    for j=1:length(Eps0{1})
        x=x+1;
        fprintf(writeID,'%d %e %e %e \r\n',x,a(i,j),b(i,j),c(i,j));
    end
end
else

writeID = fopen( [ 'Eigenstrain' ] , 'w' ) ;
a=Eps0{1};
b=Eps0{2};
c=Eps0{3};
x=0;
fprintf(writeID,[num2str(length(voxelx)),'\r\n']);
for i=1:length(Eps0{1})
    for j=1:length(Eps0{1})
        x=x+1;
        fprintf(writeID,'%d %e %e %e \r\n',x,a(i,j),b(i,j),c(i,j));
    end
end
end

fclose('all');

writeID = fopen( 'run.bat' , 'w' ) ;
if num==0
 fprintf(writeID,[ 'call abaqus job=' inpFileName '   user=sigini']) ;
else
 fprintf(writeID,[ 'call abaqus job=' inpFileName '   user=uexpan']) ;
end
fclose('all')
system('run.bat')