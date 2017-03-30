function Distance=SystemPlot(nbrBS, nbrUE, Radius)
% nbrBS: Number of base stations in the area, consider 2 tiers
% nbrUE: Number of users per BS
% Radius: The radius of a cell

%Define the radius of a cell (m)
radius=Radius;

%Number of base stations in the area, consider two tiers
nbrBSs = nbrBS;

%Number of users per BS
K = nbrUE;

% %Pathloss exponent
% alpha = 3.76;

%Generate the BS locations in the area, Positions of BSs of the central cell and the two tiers aroun
%BSpositions = [0 BSofHexagonNetwork(1,radius) BSofHexagonNetwork(2,radius)];

% My code
BSpositions = [0 BSofHexagonNetwork(1,radius)];


%Compute alternative BS locations by using wrap around.
% wrapHorizontal = radius * [0 -7.5 -3 4.5 7.5 3 -4.5];
% wrapVertical = radius*sqrt(3)/2*[0 1 7.5 6.5 -1 -8 -7.5];
% wrapLocations = wrapHorizontal+ 1i*wrapVertical;  % Move the original nbrBSs BSs to the left, right, up, down, leftup corner, leftdown corner, rightup corner, rightdown corner.
% BSpositionsWrapped = repmat(BSpositions.',[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrBSs 1]);  % Each row of BSpositionsWrapped represents the wrapped positions of a BS, a nbrBss*7 matrix

% My code
% wrapHorizontal = radius*[0 0 0 3 3 -3 -3];
% wrapVertical = radius* [0 2*sqrt(3) -2*sqrt(3) sqrt(3) -sqrt(3) -sqrt(3) sqrt(3)];
% wrapLocations = wrapHorizontal+ 1i*wrapVertical;  % Move the original nbrBSs BSs to the left, right, up, down, leftup corner, leftdown corner, rightup corner, rightdown corner.
% BSpositionsWrapped = repmat(BSpositions.',[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrBSs 1]);  % Each row of BSpositionsWrapped represents the wrapped positions of a BS, a nbrBss*7 matrix

% My code 2 (only moving the central cell to all the surrounding positions)
% if you name the sourronding positions starting in the up central from 1
% to 12, the order is (0,1,7,3,5,11,9,4,10,2,6,8,12)
BSpositionsWrapped = ones(7,13);
centralBSWrappedVertical = radius*[0 2*sqrt(3) -2*sqrt(3) sqrt(3) -sqrt(3) sqrt(3) -sqrt(3) 0 0 (3/2)*sqrt(3) -(3/2)*sqrt(3) -(3/2)*sqrt(3) (3/2)*sqrt(3)];
centralBSWrappedHorizontal = radius*[0 0 0 3 3 -3 -3 3 -3 1.5 1.5 -1.5 -1.5];
BSpositionsWrapped(1,:) = centralBSWrappedHorizontal + 1i*centralBSWrappedVertical;
for i=2:nbrBS
    BSpositionsWrapped(i,:) = BSpositions(i);
end

%Prepare to put out UEs in the cells
UEpositions = zeros(K,nbrBSs);

distancesSquaredBSj = zeros(K*nbrBSs, nbrBSs);
%Go through all the cells
for l = 1:nbrBSs
    
%     disp([num2str(l) ' cells out of ' num2str(nbrBSs)]);
     
    %Put out K users in each cell
    vertices=exp(1j*pi/3*(1:6));
    UEpositions(:,l) = DropUEinHexagonCell(K,vertices,radius) + BSpositions(l);
 
    for j = 1:nbrBSs
        
        %Compute the distance from the users to BS j
        distancesSquaredBSj(1+(l-1)*K:l*K, j) = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);

    end
    
end
Distance=distancesSquaredBSj;

% %Plot the system
% figure;
% plot(real(BSpositionsWrapped),imag(BSpositionsWrapped),'b+');
% voronoi(real(BSpositionsWrapped(:)),imag(BSpositionsWrapped(:)),'b+'); %BSs and voroni regions
% axis equal;
% hold on;
% 
% howmanytoplot = K;
% for l = 1:nbrBSs
%     
%     plot(real(UEpositions(1:howmanytoplot,l)),imag(UEpositions(1:howmanytoplot,l)),'k*'); %User locations
%     
% end


