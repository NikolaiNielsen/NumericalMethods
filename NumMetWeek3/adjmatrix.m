% function [A]=adjmatrix(n);
% This code generates a random adjacency matrix of programmable size n
% Generates the upper triangulat portion of the matrix by cycling through
% every upper triangular element and assign 1 or 0 with 50/50 probability
% then add the transpose to itself to get the lower triangular portion to
% match, row to column, the upper trangular portion-- this makes it adj...
% add weight w, w=.5 results in all zeros, w>>1 leads to mostly 1's
% w=1 leads to a 50/50 chance of a connection at each node
function [A]=adjmatrix(n,w)
A=zeros(n,n);
k=1;
for j=1:n-1 % Go through each row
    for k=j+1:n % Go through each column of the row chosen above
    A(j,k)=round(rand*w);
    if A(j,k)>1
       A(j,k)=1;
    end
    end
end
   % Now formulate the lower triangular portion to make it an adjacency
   % matrix
A=A+A';