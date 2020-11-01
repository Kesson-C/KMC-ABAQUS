% Kesson
inpFileName = ['PY_Example1_',sprintf('%03d',num)];
readID = fopen( [ inpFileName '.inp' ] , 'r' ) ;
writeID = fopen([ inpFileName 'x.inp' ] , 'w') ;
n=0;
Node={};
Element={};
node={};
element={};
k=0;
while ~feof(readID)
    str = fgetl(readID) ;  
    % reset every element to a single set
    findNodeSet = regexp(str, '*Node' ) ;
    findElementSet = regexp(str, '*Element, type=CPE4R' ) ;
    if k==1
    A=str2num(str);
    E=A(1)*10^3;
    nu=A(2);
    break
    end
    if length(str)>7
    if  str(1:8)=='*Elastic'
        k=1;
    end
    end
    if n==1 & str(1)~='*'
    Node{end+1} = str;
    elseif n==2 & str(1)~='*'
    Element{end+1} = str;
    elseif n==2 & str(1)=='*'
        n=3;
    end
if n<3
    if findNodeSet == 1 
        fprintf(writeID, str ) ;
        n=1;
    elseif findElementSet == 1 ;
        n=2;
            for i=1:length(Node)
                fprintf(writeID, Node{i}) ;
                node{end+1}=str2num(Node{i});
            end
    end
end
if n==3
          	for i=1:length(Element)
                fprintf(writeID, Element{i}) ;
                element{end+1}=str2num(Element{i});
                n=4;
            end
end
end
voxelx=[];
voxely=[];
for i=1:length(element)
    a=element{i};
    n1=node{a(2)};
    n2=node{a(3)};
    n3=node{a(4)};
    n4=node{a(5)};
    voxelx(end+1)=(n1(2)+n2(2)+n3(2)+n4(2))/4;
    voxely(end+1)=(n1(3)+n2(3)+n3(3)+n4(3))/4;
% writeID = fopen( 'run.bat' , 'w' ) ;
% fprintf(writeID,[ 'call abaqus job=' inpFileName '.inp']) ;
% fclose('all')
end
    for i=1:mesh(1)
        for j=1:mesh(2)
            xx(i,j) = voxelx(i+j);
            yy(i,j) = voxely(i+j);
        end
    end
            X0{1} = xx;
            X0{2} = yy;

fclose('all');