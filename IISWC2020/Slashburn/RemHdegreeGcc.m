function [disind,gccind,topind] = RemHdegreeGcc(B,k,dir)


if nargin<3
    dir=1;
end

n = size(B,1);

if (dir == 1)
	D = sum(B,2);
	D = D + sum(B,1)';
else
	D=sum(B,2);
end
[Dsort,I]=sort(D);


topind = flipud(I(n-k+1:n));

B(topind, :) = 0;
B(:, topind) = 0;

[gccind,disind] = ExtractGcc(B);
topind = topind';

mask = ismember(disind, topind);
disind = disind(~mask);
