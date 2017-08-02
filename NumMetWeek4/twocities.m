function a = twocities(n,c,p) % 100 2 0.4

x = random_modular_graph(n,c,p,1);

G = graph(x);

moduleOne = find(conncomp(G)==1);
moduleTwo = find(conncomp(G)==2);

for i = 1:4
    neighModuleOne = neighbors(G,moduleOne(i));
    neighModuleTwo = neighbors(G,moduleTwo(i));
    neigh = unique(cat(1,neighModuleOne,neighModuleTwo));
    newGuy = ones(length(neigh),1)*(size(G.Nodes,1)+1);
    G = addedge(G,newGuy,neigh,1);
    G = rmedge(G,moduleOne(i));
    G = rmedge(G,moduleTwo(i));
    clear neighModuleOne
    clear neighModuleTwo
    clear neigh
    clear newGuy
end

a = adjacency(G);