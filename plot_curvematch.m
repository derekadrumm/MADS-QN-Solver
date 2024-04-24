%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Post Processing for curvematching problem ------------------ %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Load solution -----------------------------------------------------
plotEnd = history.iterCount(end);
cpts_x = (x_sol(1:N))';
cpts_y = (x_sol(N+1:2*N))';

weights_sol = x_sol(2*N+1:end); % NURBS
% weights_sol = weights_init; % Bspline

cpts = [cpts_x; cpts_y];
% spline
[pnurbs,t] = pnurbspline(cpts,weights_sol,p,100,"none",[0,1,250]);

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


% plot initial spline
subplot(2,2,4)
plot(pnurbs_init(1,1),pnurbs_init(2,1),'go',pnurbs_init(1,end),pnurbs_init(2,end),'go')
hold on
plot(pnurbs_init(1,:),pnurbs_init(2,:),'g','Linewidth',1)
plot([cpts_x_init cpts_x_init(1)],[cpts_y_init cpts_y_init(1)],'g*--')

% plot curve to match
plot(curve(:,1),curve(:,2),'b','Linewidth',1)

% plot solution spline
plot(pnurbs(1,1),pnurbs(2,1),'ko',pnurbs(1,end),pnurbs(2,end),'ko')
plot(pnurbs(1,:),pnurbs(2,:),'k','Linewidth',1)
plot([cpts_x cpts_x(1)],[cpts_y cpts_y(1)],'r*--')

axis([lb ub lb ub])
title('Shape Plots')
drawnow
hold off



%% ---- Plot over iterations ----------------------------------------------
figure(figcount)
figcount = figcount + 1;

for i = 1:history.iterCount(end)
    
    cpts_x = (history.xVal(1:N,i))';
    cpts_y = (history.xVal(N+1:2*N,i))';
    
    weights_sol = history.xVal(2*N+1:end,i); % NURBS
%     weights_sol = weights_init; % Bspline
    
    cpts = [cpts_x; cpts_y];
    % spline
    [pnurbs,t,cpts_new] = pnurbspline(cpts,weights_sol,p,100,"none",[0,1,250]);
    
    % plot solution spline
    plot(pnurbs(1,1),pnurbs(2,1),'ko',pnurbs(1,end),pnurbs(2,end),'ko')
    hold on
    plot(pnurbs(1,:),pnurbs(2,:),'k','Linewidth',1)
    plot([cpts_new(1,:) cpts_new(1,1)],[cpts_new(2,:) cpts_new(2,1)],'r*--')
    text(cpts_new(1,:)+0.1, cpts_new(2,:), string(1:numel(cpts_new(1,:))));
    
    % plot curve to match
    plot(curve(:,1),curve(:,2),'b','Linewidth',1)
    
    axis([lb ub lb ub])
    title(['Iteration ', num2str(i)])
    drawnow
    hold off

end
