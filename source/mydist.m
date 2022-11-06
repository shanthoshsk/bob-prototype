function z = mydist(w,p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[S]=size(w);
for i=1:S
        z(i,:)=sqrt(sum((w(i,:)-p(1,:)).^2));
end
end

