function initialize(root_address,workPath,project)

% if ~exist(workPath)
%     mkdir(workPath)
% end
% cd (workPath)
if ~exist(project)
    mkdir(project)
end
mkdir(project)
cd (project)
mkdir Input
mkdir Data
mkdir Fig
cd ..
% command = ['cp ',root_address,'funcLib/libBP/General_BP.m ./'];
% system(command);
end
