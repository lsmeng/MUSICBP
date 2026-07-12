function check_data(path, Band)
%this function is used to check for and remove NaN or Inf in the ret.xori
%it should be carried out before alignment
filename=strcat(path,'Input/data',num2str(Band-1))
load(filename);
siz=size(ret.xori);
index=[];
count=0;
for i=1:1:siz(1)
    if isnan(ret.xori(i,1)) | ret.xori(i,1)==Inf
        count=count+1;
        index(count)=i;      %record the row number of NaN or Inf
    end
end

if count~=0
    for i=1:1:count
        siz=size(ret.xori);
        I=[1:index(i)-1 index(i)+1:siz(1)];
        ret=orderAll(ret,0,I);
    end
    save(filename)
end

end