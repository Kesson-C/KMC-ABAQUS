% Kesson
inpFileName = ['PY_Example1_',sprintf('%03d',num)];
change = ['PY_Example1_',sprintf('%03d',num-1)];
readID = fopen( [ inpFileName '.inp' ] , 'r' ) ;
n=0;
str={};
while ~feof(readID)
    n=n+1;
    str{end+1} = char(fgets(readID)) ;
    if length(str{n})>11
        a=str{n};
        if strcmp(a(1:9),'*Instance')
            str{n} = ['*Instance, library=',change,', instance="BMG ribbon-1"'];
            ss=n;
        elseif strcmp(a(1:11),'Dummy_Right')
            str{n} = ['Dummy_Right, 1, 1, ',sprintf('%d',EpsAvg)]; %%% need modified
            sss=n;
        end
    end
end
fclose(readID);
readID = fopen( [ inpFileName '.inp' ] , 'w' ) ;
for i=1:n
    if i==ss
    fprintf(readID,'%s\n',str{i});
    elseif i==sss
    fprintf(readID,'%s\n',str{i});
    else
    fprintf(readID,'%s',str{i});
    end
end
fclose(readID);