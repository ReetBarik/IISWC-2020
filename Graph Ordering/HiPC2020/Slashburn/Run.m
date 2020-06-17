function exitcode = Run(filename, undirected)
fid = fopen(filename);
C = textscan(fid,'%d %d\n');
fclose(fid);
if undirected == 1
    G = graph(C{1},C{2});
else
    G = digraph(C{1},C{2});
end

A = adjacency(G);
A = sparse(1.0 * ((A+A') > 0));
%%% Run SlashBurn
k=1;
dir=0;
[~,~,Ak] = SlashBurn(A, k, dir);

e = Edgelist(Ak);
fid = fopen(filename, 'w');
for row = 1:length(e)
    fprintf(fid, '%d %d\n', e(row), e(row,2));
end
fclose(fid);
exitcode = 0;
end