% 
% SlashBurn: shatter graph and reorder nodes to make a compact adjacency matrix.
%
% Parameter
%   AOrig : adjacency matrix of a graph. 
%   dir : 1 (for a directed graph) or 0 (undirected graph)
%   k : # of nodes to to shatter in each iteration
%
% Return values
%   niter : # of iteration
%   gccsize : a vector containing gcc sizes over iterations
%   Ak : newly reordered adjacency matrix
%
% Example:
% >> [niter,gccsize,Ak] = SlashBurn( oregon_orig, 128, 1);
%

function [niter,gccsize,Ak] = SlashBurn(AOrig, k, dir)

if nargin<3
    dir=0;
end

gccsize = zeros(0,0);
niter=0;
n = size(AOrig,2);
totalind = zeros(1,n);
cur_lpos = 1;
cur_rpos = n;
gccind = [1:n];
cur_gccsize = n;


while niter == 0 | gccsize > k
	niter = niter+1;
	fprintf('Iteration %d...\n', niter);

	A = AOrig(gccind,gccind);
	[disind,newgccind,topind] = RemHdegreeGcc(A,k,dir);

	topind_size = size(topind, 2);

	totalind(cur_lpos:cur_lpos + topind_size - 1) = gccind(topind);
	cur_lpos = cur_lpos + topind_size;
	totalind(cur_rpos - size(disind,2) + 1:cur_rpos) = gccind(disind);
	cur_rpos = cur_rpos - size(disind,2);
	newind = totalind;

	gccind = gccind(newgccind);

	newind(cur_lpos:cur_rpos) = gccind;
	
	Ak = AOrig(newind,newind);
	%spy(Ak);

	cur_gccsize = size(gccind, 2);

	gccsize = [
	gccsize
	size(gccind,2)
	];
end

fprintf('Total iterration = %d. \n', niter);

