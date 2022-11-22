% Make one file of KoriModel

clear;

fid=fopen('KoriModelAll.m','w');
fod=fopen('KoriModel.m');
while ~feof(fod)
    fprintf(fid,'\n%s',(fgetl(fod)));
end
fclose('all');

fid=fopen('KoriModelAll.m','a+');

list=dir('subroutines/*.m');
for i=1:length(list)
        fod=fopen(strcat('subroutines/',list(i).name));
        while ~feof(fod)
            fprintf(fid,'\n%s',(fgetl(fod)));
        end
end
fclose('all');

