function [cur_gccind,cur_disind] = ExtractGcc(B)

[S,C]=graphconncomp(B, 'WEAK', true);

maxind=-1;
maxsize=0;

size_v = zeros(0, S);

for k=1:S
    size_v(k)=size(find(C == k), 2);
end

[size_sort,I]=sort(size_v, 'descend');

cur_gccind = find(C == I(1));

cur_disind = zeros(0,0);

for k=2:S
    curind = find(C == I(k));
    cur_disind = [cur_disind curind];
end

% fprintf('\tgccsize\t%d\n', size(cur_gccind,2));

