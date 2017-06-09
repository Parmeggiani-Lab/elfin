close all;
figure; 
xyrange = [0,400,0,400];
axis(xyrange); 
hold off;

grid on
box on
[x,y]=getpts();

plot(x,y,'-bx'); 
axis(xyrange)

P=[x,y];

D = (P-circshift(P, 1));
normaliseFactor = 40; % this should be set around mean pair-CoM distance
m = mean(sqrt(sum(D.^2, 2))) / normaliseFactor;

format long g
P=P./m;
P=[P, zeros(size(P,1), 1)]

cd('/Users/joy/src/elfin/bm/fun')

% save using dlmwrite('R3.csv', P*3, ' ')

% output figure styles
% linewidth 2, markersize 15
% manually set ranges

% xlabel('X'); ylabel('Y'); zlabel('Z')
% xticklabels([]); yticklabels([]); zticklabels([])
% set(gca, 'fontsize', 32)
% xticks(-400:25:400); yticks(-400:25:400); zticks(-400:25:400)
% lims = [(floor(min(X(:,1))/25)-1)*25 (ceil(max(X(:,1))/25)+1)*25
%     (floor(min(X(:,2))/25)-1)*25 (ceil(max(X(:,2))/25)+1)*25
% (floor(min(X(:,3))/25)-1)*25 (ceil(max(X(:,3))/25)+1)*25]';
% axis(lims(:))
% view([45, 45, 120]); axis equal; axis(lims(:))
% view([90,0,0]); axis equal; axis(lims(:))
% view([0,90,0]); axis equal; axis(lims(:))

% To plot modified CSV points:
% figure; hold all;
% plot3(X(:,1), X(:,2), X(:,3), '-b', 'linewidth', 2, 'markersize', 15)
% plot3(X(:,1), X(:,2), X(:,3), 'xr', 'linewidth', 2, 'markersize', 15)
% grid on; box on
% axis equal