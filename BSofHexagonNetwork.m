function P=BSofHexagonNetwork(n,r)  
% n : which tier to be constructed, n=1,2,3; 
% P(1,:): x positions of cells, 
% P(2,:): y postions of cells; 
% r: radius of hexagon cell

a=[3/2; sqrt(3)/2]*r;
b=[0;sqrt(3)]*r;
if(n==1)
    m=6;
    alfa1=[1 0 -1 -1 0 1];
    alfa2=[0 1 1 0 -1 -1];
    PP=repmat(a,1,m).*repmat(alfa1,2,1) + repmat(b,1,m).*repmat(alfa2,2,1); 
elseif(n==2)
    m=12;
    alfa1=[2 1 0 -1 -2 -2 -2 -1 0 1 2 2];
    alfa2=[0 1 2 2 2 1 0 -1 -2 -2 -2 -1];
    PP=repmat(a,1,m).*repmat(alfa1,2,1) + repmat(b,1,m).*repmat(alfa2,2,1); 
elseif(n==3)
    m=18;
    alfa1=[3 2 1 0 -1 -2 -3 -3 -3 -3 -2 -1 0 1 2 3 3 3];
    alfa2=[0 1 2 3 3 3 3 2 1 0 -1 -2 -3 -3 -3 -3 -2 -1];
    PP=repmat(a,1,m).*repmat(alfa1,2,1) + repmat(b,1,m).*repmat(alfa2,2,1); 
end
 P=PP(1,:) + 1i* PP(2,:);