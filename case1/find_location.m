function [j,i] = find_location(p,N)

temp = 0;
for i = 1:N
    temp = temp + i;
    if p <= temp
        j =p-(temp-i);
        break
    end
end
