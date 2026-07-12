function addfilelist(path)     %add all *.SAC.r files' names to the filelist

%cd(path)
d=dir([ path '/Data/*.SAC*']);
% i=size(d)
% filelist=""

% for j=1:1:size(d)
%     
%     filelist(j)=string(d(j).name);
% 
% end

fid = fopen([path '/Data/filelist'],'wt');
for i=1:size(d)
    fprintf(fid,'%s\n',d(i).name);
end
% filelist=filelist(2:end)
% fid = fopen('filelist','wt');
% fprintf(fid,'%s\n',filelist); 
fclose(fid);

% copyfile('filelist', 'Data');
end
