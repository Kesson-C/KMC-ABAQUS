% Kesson
inpFileName = ['PY_Example1_',sprintf('%03d',num)];
read = fopen([ inpFileName '.dat'] , 'r');
write = fopen( [inpFileName '.txt'], 'w');
n=0;
all={};
start=[];
ND={};
S11=[];
S12=[];
S22=[];
E11=[];
E12=[];
E22=[];
S={};
EE={};
while ~feof(read)
    n=n+1;
    str = char(fgets(read)) ;
    all{end+1}=str;
    if length(str)>11
        if strcmp(str(1:11),'    ELEMENT') | strcmp(str(1:11),'       NODE')
            start(end+1)=n;
        elseif strcmp(str(1:9),'*Instance')
            
        end 
    end
end

for i=start(1)+3:start(1)+length(element)+2
temp=all{i};
temp=str2num(temp(25:end));
S{end+1}=temp;
end


for i=1:length(S)
    ss = S{i};
    S11(end+1) = ss(5);
    S12(end+1) = ss(8);
    S22(end+1) = ss(6);
    E11(end+1) = ss(1);
    E12(end+1) = ss(4);
    E22(end+1) = ss(2);
end

    for i=1:mesh(1)
        for j=1:mesh(2)
            e11(i,j) = E11(i+j);
            e12(i,j) = E12(i+j);
            e22(i,j) = E22(i+j);
            s11(i,j) = S11(i+j);
            s12(i,j) = S12(i+j);
            s22(i,j) = S22(i+j);
        end
    end
            Eps0{1} = e11;
            Eps0{2} = e12;
            Eps0{3} = e22;
            Sigma{1} = s11;
            Sigma{2} = s12;
            Sigma{3} = s22;
fclose('all');