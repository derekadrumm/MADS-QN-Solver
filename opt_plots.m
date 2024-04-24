%%%% Plot MADS_SMF_CON iteration info
plotEnd = history.iterCount(end);

% Objective Function Value
figure(figcount)
figcount = figcount+1;
plot(1:plotEnd,history.fVal(1:plotEnd),'b.-','MarkerSize',5)
title('Objective Function f(x_k)')
xlabel('Iteration')
ylabel('f(x_k)')
drawnow


%% Individual Plots
% Mesh and Frame Size
figure(figcount)
figcount = figcount+1;
semilogy(1:plotEnd,history.meshSize(1:plotEnd),'b.-','MarkerSize',5)
hold on
semilogy(1:plotEnd,history.frameSize(1:plotEnd),'r.-','MarkerSize',5)
semilogy(1:plotEnd,history.error.relGradError(1:plotEnd),'g.-','MarkerSize',5)
semilogy(1:plotEnd,history.error.relMeshGradError(1:plotEnd),'m.-','MarkerSize',5)
% title('Mesh and Frame Size')
xlabel('Iteration')
hold off
legend('\delta^k', '\Delta^k','||relgrad||_{\infty}','||relgradmesh||_{\infty}',...
    'Location', 'southwest')
drawnow

% % Objective Function Value AND Surrogate Function Value
% figure(figcount)
% figcount = figcount+1;
% for k = 1:length(history.x_vals)
%     plot(1:k,history.fun_vals(1:k),'b.-','MarkerSize',5)
%     hold on
% %     plot(1:k,history.funsur_vals(1:k),'r.-','MarkerSize',5)
%     title('Objective Function f(x_k)')
%     xlabel('Iteration')
%     ylabel('f(x_k)')
%     drawnow
%     hold off
% end


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

figure(figcount)
figcount = figcount+1;
b = bar(typeNames,count(check),'facecolor', 'flat');
b.CData = C;
xlabel('Step Type')
ylabel('Number of Steps')


%% Subplots
fig = figure(figcount);
figcount = figcount+1;
fig.Position = [100 270 1320 500];

% Mesh and Frame Size
subplot(1,2,1)
semilogy(1:plotEnd,history.meshSize(1:plotEnd),'b.-','MarkerSize',5)
hold on
semilogy(1:plotEnd,history.frameSize(1:plotEnd),'r.-','MarkerSize',5)
semilogy(1:plotEnd,history.error.relGradError(1:plotEnd),'g.-','MarkerSize',5)
semilogy(1:plotEnd,history.error.relMeshGradError(1:plotEnd),'m.-','MarkerSize',5)
xlabel('Iteration')
hold off
% legend('\delta_k', '\Delta_k',...
%     'Location', 'southwest')
legend('\delta^k', '\Delta^k','||relgrad||_{\infty}','||relgradmesh||_{\infty}',...
    'Location', 'southwest')
drawnow

% % Objective Function Value AND Surrogate Function Value
% figure(figcount)
% figcount = figcount+1;
% for k = 1:length(history.x_vals)
%     plot(1:k,history.fun_vals(1:k),'b.-','MarkerSize',5)
%     hold on
% %     plot(1:k,history.funsur_vals(1:k),'r.-','MarkerSize',5)
%     title('Objective Function f(x_k)')
%     xlabel('Iteration')
%     ylabel('f(x_k)')
%     drawnow
%     hold off
% end


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

subplot(1,2,2)
b = bar(typeNames,count(check),'facecolor', 'flat');
b.CData = C;
xlabel('Step Type')
ylabel('Number of Steps')






