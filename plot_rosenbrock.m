%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Post Processing for Rosenbrock problem --------------------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Objective Function (Rosenbrock) -----------------------------------
plotEnd = history.iterCount(end);
x_k = history.xVal;
f_k = history.fVal;

% axis
lx = -2;
ly = -1;
ux = 2;
uy = 4;

x = linspace(lx,ux,500);
y = linspace(ly,uy,500);
[X,Y] = meshgrid(x,y);


XX = [X(:) Y(:)]';
FF = zeros(1,size(XX,2));
for i = 1:size(XX,2)
   
    FF(i) = f(XX(:,i));
    
end
FF = reshape(FF,size(X));

clines = [0.01 0.1 0.5 1 10 25 50 100 150 250];


%% ---- Plot results ------------------------------------------------------
fig = figure(figcount);
figcount = figcount+1;
fig.Position = [40 80 900 700];

% Mesh and Frame Size
subplot(2,2,1)
semilogy(1:plotEnd,history.meshSize(1:plotEnd),'b.-','MarkerSize',5)
hold on
semilogy(1:plotEnd,history.frameSize(1:plotEnd),'r.-','MarkerSize',5)
semilogy(1:plotEnd,history.error.relGradError(1:plotEnd),'g.-','MarkerSize',5)
semilogy(1:plotEnd,history.error.relMeshGradError(1:plotEnd),'m.-','MarkerSize',5)
xlabel('Iteration')
hold off
% legend('\delta_k', '\Delta_k',...
%     'Location', 'southwest')
legend('\delta_k', '\Delta_k','relgrad Error','relgradmesh Error',...
    'Location', 'southwest')
title('Iteration Tolerances')
drawnow


% Iteration Type Count Histogram
vec = history.iterTypeCode;
% Determine which search steps were used
temp = 2;
check = [1 1];
if options.searchStep.quasiNewton
    % QN search used
    temp = temp + 2;
    check = [check 1 1];
else
    check = [check 0 0];
end
if options.searchStep.lineSearch
    % LS used
    temp = temp + 1;
    check = [check 1];
else
    check = [check 0];
end
if options.searchStep.samplingSearch
    % Sampling Search used
    temp = temp + 1;
    check = [check 1];
else
    check = [check 0];
end
check = logical(check);

% total counts
count = zeros(temp,1);
for i = 1:6
    count(i) = sum(vec == i);
end

typeCodes = [1:6];
% The possible search steps
typeNamesStr = ["Refine"; "Poll"; "QN"; "QN-LS"; "LS"; "Sample"];
% only plot search steps which were used
typeNames = categorical(typeNamesStr(check));
typeNames = reordercats(typeNames,typeNamesStr(check));
C = [[1 1 0]; [1 0 0]; [0 0 0]; [0.5 0.5 0.5]; [0 1 1]; [1 0 1]]; % iter type colors
C = C(check,:);

subplot(2,2,2)
b = bar(typeNames,count(check),'facecolor', 'flat');
b.CData = C;
xlabel('Step Type')
ylabel('Number of Steps')
title('Total Steps')


% Objective Function Value
subplot(2,2,3)
plot(1:plotEnd,history.fVal(1:plotEnd),'b.-','MarkerSize',5)
title('Objective Function f(x^k)')
xlabel('Iteration')
ylabel('f(x^k)')
drawnow


% Plot function solution at each iteration
subplot(2,2,4)
contourf(X,Y,FF,clines,'edgecolor','none')
colormap('winter')
axis([lx ux ly uy])
colorbar
set(gca,'ColorScale','log')
drawnow

C = {[1 1 0]; [1 0 0]; [0 0 0]; [0.5 0.5 0.5]; [0 1 1]; [1 0 1]}; % iter type colors {Refine, Poll, QN, LS, Sample}
hold on
plot(1,1,'kx','LineWidth',1,'MarkerSize',5)
plot([x0(1,:) history.xVal(1,1)],[x0(2,:) history.xVal(2,1)],'-',...
    'Color',C{history.iterTypeCode(1)},...
    'MarkerSize',5,...
    'MarkerEdgeColor',C{history.iterTypeCode(1)},...
    'MarkerFaceColor',C{history.iterTypeCode(1)})
plot(history.xVal(1,1),history.xVal(2,1),'*',...
    'MarkerSize',5,...
    'MarkerEdgeColor',C{history.iterTypeCode(1)},...
    'MarkerFaceColor',C{history.iterTypeCode(1)})
plot(x0(1),x0(2),'m*','MarkerSize',5)
plot(x0(1),x0(2),'mo','MarkerSize',7)
for k = 2:length(history.xVal)
    plot([history.xVal(1,k-1) history.xVal(1,k)],[history.xVal(2,k-1) history.xVal(2,k)],'-',...
        'Color',C{history.iterTypeCode(k)},...
        'MarkerSize',5,...
        'MarkerEdgeColor',C{history.iterTypeCode(k)},...
        'MarkerFaceColor',C{history.iterTypeCode(k)})
    plot(history.xVal(1,k),history.xVal(2,k),'*',...
        'MarkerSize',5,...
        'MarkerEdgeColor',C{history.iterTypeCode(k)},...
        'MarkerFaceColor',C{history.iterTypeCode(k)})
    drawnow
end
hold off
xlabel('x_1')
ylabel('x_2')
title('Rosenbrock Function')

