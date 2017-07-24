function [nodeInit,aInit,ls,lmax] = treatNetwork(name,randormost,verb)

    data = importdata(name);
    data = data(:,1:2);
    data = unique(sort(data,2),'rows');
    s = data(:,1);
    t = data(:,2);
    
    a = zeros(max(data(:)));
    for i = 1:length(s)
       a(s(i),t(i)) = 1; 
    end
    a(end,end) = 0;
    a = a+a';
    
    
    node = createnet(a,0);
    [node,ls{1}] = bfs(node);
    
    nodeInit = node;
    aInit = a;
    
    if randormost == 'rand'
        [ls,lmax] = removerand(node,a,verb);
    elseif randormost == 'most'
        [ls,lmax] = removemost(node,a,verb);
    else
        error('either rand or most')
    end

