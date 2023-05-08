%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compiles results from LDES proposal simulations
% Kaden Plewe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

%% Load data and FF Classes

% Sim 1
FF_test1 = load('FF_test1.mat'); FF_test1 = FF_test1.FF_test;
FF_test1.thetaFolder = 'thetaTest1';

% Sim 2
FF_test2 = load('FF_test2.mat'); FF_test2 = FF_test2.FF_test;
FF_test2.thetaFolder = 'thetaTest2';

% Sim 3
FF_test3 = load('FF_test3.mat'); FF_test3 = FF_test3.FF_test;
FF_test3.thetaFolder = 'thetaTest3';

% Sim 4
FF_test4 = load('FF_test4.mat'); FF_test4 = FF_test4.FF_test;
FF_test4.thetaFolder = 'thetaTest4';

% Sim 5
FF_test5 = load('FF_test5.mat'); FF_test5 = FF_test5.FF_test;
FF_test5.thetaFolder = 'thetaTest5';

% Sim 6
FF_test6 = load('FF_test6.mat'); FF_test6 = FF_test6.FF_test;
FF_test6.thetaFolder = 'thetaTest6';

% Sim 7
FF_test7 = load('FF_test7.mat'); FF_test7 = FF_test7.FF_test;
FF_test7.thetaFolder = 'thetaTest7';

% Sim 8
FF_test8 = load('FF_test8.mat'); FF_test8 = FF_test8.FF_test;
FF_test8.thetaFolder = 'thetaTest8';

% Sim 9
FF_test9 = load('FF_test9.mat'); FF_test9 = FF_test9.FF_test;
FF_test9.thetaFolder = 'thetaTest9';

% Sim 10
FF_test10 = load('FF_test10.mat'); FF_test10 = FF_test10.FF_test;
FF_test10.thetaFolder = 'thetaTest10';

% Sim 11
FF_test11 = load('FF_test11.mat'); FF_test11 = FF_test11.FF_test;
FF_test11.thetaFolder = 'thetaTest11';

% Sim 12
FF_test12 = load('FF_test12.mat'); FF_test12 = FF_test12.FF_test;
FF_test12.thetaFolder = 'thetaTest12';

%% Compile data for plots
% energy loss (kJ)
% total energy in (J)

t1 = FF_test1.Fo2t(FF_test1.Fo, 1);
t2 = FF_test2.Fo2t(FF_test2.Fo, 1);
t3 = FF_test3.Fo2t(FF_test3.Fo, 1);
t4 = FF_test4.Fo2t(FF_test4.Fo, 1);
t5 = FF_test5.Fo2t(FF_test5.Fo, 1);
t6 = FF_test6.Fo2t(FF_test6.Fo, 1);
t7 = FF_test7.Fo2t(FF_test7.Fo, 1);
t8 = FF_test8.Fo2t(FF_test8.Fo, 1);
t9 = FF_test9.Fo2t(FF_test9.Fo, 1);
t10 = FF_test10.Fo2t(FF_test10.Fo, 1);
t11 = FF_test11.Fo2t(FF_test11.Fo, 1);
t12 = FF_test12.Fo2t(FF_test12.Fo, 1);

eTot = 0.95*FF_test1.rhopPack*FF_test1.cpp*FF_test1.Hp^2*pi*FF_test1.bp^2*(600-25);
qTot1 = zeros(length(FF_test1.Fo)-1, 1);
energyLoss1 = zeros(length(FF_test1.Fo)-1, 1);
if length(FF_test1.Fo) < FF_test1.ls, loadTheta(FF_test1, length(FF_test1.Fo));  
else, loadTheta(FF_test1, FF_test1.ls); end
for i = 2:length(FF_test1.Fo)                                      
    i_ = mod(i-1, FF_test1.ls) + 1;
    if i_ == 1 && i ~= length(FF_test1.Fo)
        loadTheta(FF_test1, i+FF_test1.ls-1);
    end
    qTot1(i-1) = FF_test1.theta{i_, 14}; 
    if i == 2
        energyLoss1(i-1) = 0;
    elseif i == 3
        energyLoss1(i-1) = trapz(t1(1:i-1), qTot1(1:i-1));
    else
        ft = simpsonIntegrator(FF_test1, t1(1:i-1));
        energyLoss1(i-1) = ft*qTot1(1:i-1);
    end
end

qTot2 = zeros(length(FF_test2.Fo)-1, 1);
energyLoss2 = zeros(length(FF_test2.Fo)-1, 1);
if length(FF_test2.Fo) < FF_test2.ls, loadTheta(FF_test2, length(FF_test2.Fo));  
else, loadTheta(FF_test2, FF_test2.ls); end
for i = 2:length(FF_test2.Fo)                                      
    i_ = mod(i-1, FF_test2.ls) + 1;
    if i_ == 1 && i ~= length(FF_test2.Fo)
        loadTheta(FF_test2, i+FF_test2.ls-1);
    end
    qTot2(i-1) = FF_test2.theta{i_, 14}; 
    if i == 2
        energyLoss2(i-1) = 0;
    elseif i == 3
        energyLoss2(i-1) = trapz(t2(1:i-1), qTot2(1:i-1));
    else
        ft = simpsonIntegrator(FF_test2, t2(1:i-1));
        energyLoss2(i-1) = ft*qTot2(1:i-1);
    end
end

qTot3 = zeros(length(FF_test3.Fo)-1, 1);
energyLoss3 = zeros(length(FF_test3.Fo)-1, 1);
if length(FF_test3.Fo) < FF_test3.ls, loadTheta(FF_test3, length(FF_test3.Fo));  
else, loadTheta(FF_test3, FF_test3.ls); end
for i = 2:length(FF_test3.Fo)                                      
    i_ = mod(i-1, FF_test3.ls) + 1;
    if i_ == 1 && i ~= length(FF_test3.Fo)
        loadTheta(FF_test3, i+FF_test3.ls-1);
    end
    qTot3(i-1) = FF_test3.theta{i_, 14}; 
    if i == 2
        energyLoss3(i-1) = 0;
    elseif i == 3
        energyLoss3(i-1) = trapz(t3(1:i-1), qTot3(1:i-1));
    else
        ft = simpsonIntegrator(FF_test3, t3(1:i-1));
        energyLoss3(i-1) = ft*qTot3(1:i-1);
    end
end

qTot4 = zeros(length(FF_test4.Fo)-1, 1);
energyLoss4 = zeros(length(FF_test4.Fo)-1, 1);
if length(FF_test4.Fo) < FF_test4.ls, loadTheta(FF_test4, length(FF_test4.Fo));  
else, loadTheta(FF_test4, FF_test4.ls); end
for i = 2:length(FF_test4.Fo)                                      
    i_ = mod(i-1, FF_test4.ls) + 1;
    if i_ == 1 && i ~= length(FF_test4.Fo)
        loadTheta(FF_test4, i+FF_test4.ls-1);
    end
    qTot4(i-1) = FF_test4.theta{i_, 14}; 
    if i == 2
        energyLoss4(i-1) = 0;
    elseif i == 3
        energyLoss4(i-1) = trapz(t4(1:i-1), qTot4(1:i-1));
    else
        ft = simpsonIntegrator(FF_test4, t4(1:i-1));
        energyLoss4(i-1) = ft*qTot4(1:i-1);
    end
end

qTot5 = zeros(length(FF_test5.Fo)-1, 1);
energyLoss5 = zeros(length(FF_test5.Fo)-1, 1);
if length(FF_test5.Fo) < FF_test5.ls, loadTheta(FF_test5, length(FF_test5.Fo));  
else, loadTheta(FF_test5, FF_test5.ls); end
for i = 2:length(FF_test5.Fo)                                      
    i_ = mod(i-1, FF_test5.ls) + 1;
    if i_ == 1 && i ~= length(FF_test5.Fo)
        loadTheta(FF_test5, i+FF_test5.ls-1);
    end
    qTot5(i-1) = FF_test5.theta{i_, 14}; 
    if i == 2
        energyLoss5(i-1) = 0;
    elseif i == 3
        energyLoss5(i-1) = trapz(t5(1:i-1), qTot5(1:i-1));
    else
        ft = simpsonIntegrator(FF_test5, t5(1:i-1));
        energyLoss5(i-1) = ft*qTot5(1:i-1);
    end
end

qTot6 = zeros(length(FF_test6.Fo)-1, 1);
energyLoss6 = zeros(length(FF_test6.Fo)-1, 1);
if length(FF_test6.Fo) < FF_test6.ls, loadTheta(FF_test6, length(FF_test6.Fo));  
else, loadTheta(FF_test6, FF_test6.ls); end
for i = 2:length(FF_test6.Fo)                                      
    i_ = mod(i-1, FF_test6.ls) + 1;
    if i_ == 1 && i ~= length(FF_test6.Fo)
        loadTheta(FF_test6, i+FF_test6.ls-1);
    end
    qTot6(i-1) = FF_test6.theta{i_, 14}; 
    if i == 2
        energyLoss6(i-1) = 0;
    elseif i == 3
        energyLoss6(i-1) = trapz(t6(1:i-1), qTot6(1:i-1));
    else
        ft = simpsonIntegrator(FF_test6, t6(1:i-1));
        energyLoss6(i-1) = ft*qTot6(1:i-1);
    end
end

qTot7 = zeros(length(FF_test7.Fo)-1, 1);
energyLoss7 = zeros(length(FF_test7.Fo)-1, 1);
if length(FF_test7.Fo) < FF_test7.ls, loadTheta(FF_test7, length(FF_test7.Fo));  
else, loadTheta(FF_test7, FF_test7.ls); end
for i = 2:length(FF_test7.Fo)                                      
    i_ = mod(i-1, FF_test7.ls) + 1;
    if i_ == 1 && i ~= length(FF_test7.Fo)
        loadTheta(FF_test7, i+FF_test7.ls-1);
    end
    qTot7(i-1) = FF_test7.theta{i_, 14}; 
    if i == 2
        energyLoss7(i-1) = 0;
    elseif i == 3
        energyLoss7(i-1) = trapz(t7(1:i-1), qTot7(1:i-1));
    else
        ft = simpsonIntegrator(FF_test7, t7(1:i-1));
        energyLoss7(i-1) = ft*qTot7(1:i-1);
    end
end

qTot8 = zeros(length(FF_test8.Fo)-1, 1);
energyLoss8 = zeros(length(FF_test8.Fo)-1, 1);
if length(FF_test8.Fo) < FF_test8.ls, loadTheta(FF_test8, length(FF_test8.Fo));  
else, loadTheta(FF_test8, FF_test8.ls); end
for i = 2:length(FF_test8.Fo)                                      
    i_ = mod(i-1, FF_test8.ls) + 1;
    if i_ == 1 && i ~= length(FF_test8.Fo)
        loadTheta(FF_test8, i+FF_test8.ls-1);
    end
    qTot8(i-1) = FF_test8.theta{i_, 14}; 
    if i == 2
        energyLoss8(i-1) = 0;
    elseif i == 3
        energyLoss8(i-1) = trapz(t8(1:i-1), qTot8(1:i-1));
    else
        ft = simpsonIntegrator(FF_test8, t8(1:i-1));
        energyLoss8(i-1) = ft*qTot8(1:i-1);
    end
end

qTot9 = zeros(length(FF_test9.Fo)-1, 1);
energyLoss9 = zeros(length(FF_test9.Fo)-1, 1);
if length(FF_test9.Fo) < FF_test9.ls, loadTheta(FF_test9, length(FF_test9.Fo));  
else, loadTheta(FF_test9, FF_test9.ls); end
for i = 2:length(FF_test9.Fo)                                      
    i_ = mod(i-1, FF_test9.ls) + 1;
    if i_ == 1 && i ~= length(FF_test9.Fo)
        loadTheta(FF_test9, i+FF_test9.ls-1);
    end
    qTot9(i-1) = FF_test9.theta{i_, 14}; 
    if i == 2
        energyLoss9(i-1) = 0;
    elseif i == 3
        energyLoss9(i-1) = trapz(t9(1:i-1), qTot9(1:i-1));
    else
        ft = simpsonIntegrator(FF_test9, t9(1:i-1));
        energyLoss9(i-1) = ft*qTot9(1:i-1);
    end
end

qTot10 = zeros(length(FF_test10.Fo)-1, 1);
energyLoss10 = zeros(length(FF_test10.Fo)-1, 1);
if length(FF_test10.Fo) < FF_test10.ls, loadTheta(FF_test10, length(FF_test10.Fo));  
else, loadTheta(FF_test10, FF_test10.ls); end
for i = 2:length(FF_test10.Fo)                                      
    i_ = mod(i-1, FF_test10.ls) + 1;
    if i_ == 1 && i ~= length(FF_test10.Fo)
        loadTheta(FF_test10, i+FF_test10.ls-1);
    end
    qTot10(i-1) = FF_test10.theta{i_, 14}; 
    if i == 2
        energyLoss10(i-1) = 0;
    elseif i == 3
        energyLoss10(i-1) = trapz(t10(1:i-1), qTot10(1:i-1));
    else
        ft = simpsonIntegrator(FF_test10, t10(1:i-1));
        energyLoss10(i-1) = ft*qTot10(1:i-1);
    end
end

qTot11 = zeros(length(FF_test11.Fo)-1, 1);
energyLoss11 = zeros(length(FF_test11.Fo)-1, 1);
if length(FF_test11.Fo) < FF_test11.ls, loadTheta(FF_test11, length(FF_test11.Fo));  
else, loadTheta(FF_test11, FF_test11.ls); end
for i = 2:length(FF_test11.Fo)                                      
    i_ = mod(i-1, FF_test11.ls) + 1;
    if i_ == 1 && i ~= length(FF_test11.Fo)
        loadTheta(FF_test11, i+FF_test11.ls-1);
    end
    qTot11(i-1) = FF_test11.theta{i_, 14}; 
    if i == 2
        energyLoss11(i-1) = 0;
    elseif i == 3
        energyLoss11(i-1) = trapz(t11(1:i-1), qTot11(1:i-1));
    else
        ft = simpsonIntegrator(FF_test11, t11(1:i-1));
        energyLoss11(i-1) = ft*qTot11(1:i-1);
    end
end

qTot12 = zeros(length(FF_test12.Fo)-1, 1);
energyLoss12 = zeros(length(FF_test12.Fo)-1, 1);
if length(FF_test12.Fo) < FF_test12.ls, loadTheta(FF_test12, length(FF_test12.Fo));  
else, loadTheta(FF_test12, FF_test12.ls); end
for i = 2:length(FF_test12.Fo)                                      
    i_ = mod(i-1, FF_test12.ls) + 1;
    if i_ == 1 && i ~= length(FF_test12.Fo)
        loadTheta(FF_test12, i+FF_test12.ls-1);
    end
    qTot12(i-1) = FF_test12.theta{i_, 14}; 
    if i == 2
        energyLoss12(i-1) = 0;
    elseif i == 3
        energyLoss12(i-1) = trapz(t12(1:i-1), qTot12(1:i-1));
    else
        ft = simpsonIntegrator(FF_test12, t12(1:i-1));
        energyLoss12(i-1) = ft*qTot12(1:i-1);
    end
end


%% energy loss histogram plot

eLoss = [energyLoss1(end); energyLoss2(end); energyLoss3(end); energyLoss4(end); ...
    energyLoss5(end); energyLoss6(end); energyLoss7(end); energyLoss8(end); ...
    energyLoss9(end); energyLoss10(end); energyLoss11(end); energyLoss12(end)]./3600;
eTot_ = eTot/1000/3600;
eLossPercent = 100*eLoss./eTot_;
names = categorical({'Sim 1', 'Sim 2', 'Sim 3', 'Sim 4', 'Sim 5', 'Sim 6', ...
    'Sim 7', 'Sim 8', 'Sim 9', 'Sim 10', 'Sim 11', 'Sim 12'});
names = reordercats(names, {'Sim 1', 'Sim 2', 'Sim 3', 'Sim 4', 'Sim 5', 'Sim 6', ...
    'Sim 7', 'Sim 8', 'Sim 9', 'Sim 10', 'Sim 11', 'Sim 12'});
f = figure('Units', 'normalized', 'Position', [0 0 0.5 0.25], 'Visible', 'on');
set(f, 'defaultAxesColorOrder', [[0, 0, 0]; [0, 0, 0]]);
b = bar(names, eLoss, 'BarWidth', 0.7); hold on;
ylim([0, 50]);
b.FaceColor = 'k';
ylabel('Total Energy Loss ($kWh_{t}$)', 'interpreter', 'latex', 'FontSize', 14);

yyaxis right
b = bar(names, eLossPercent, 'BarWidth', 0.7);
ylabel('Total Energy Loss (\%)', 'interpreter', 'latex', 'FontSize', 14);
ylim([0, 50*100/eTot_]);
b.FaceColor = 'k';

set(gca, 'Color', [1 1 1], 'TickLabelInterpreter', 'latex', 'FontSize', 14);

%% Temporal energy loss plots

f = figure('Units', 'normalized', 'Position', [0 0 0.2 0.2], 'Visible', 'on');

plot(t1/3600, [0; energyLoss1]./3600, '.k'); hold on;
plot(t2/3600, [0; energyLoss2]./3600, '.k');
plot(t3/3600, [0; energyLoss3]./3600, '.k');
plot(t4/3600, [0; energyLoss4]./3600, '.k');
plot(t5/3600, [0; energyLoss5]./3600, '.k');
plot(t6/3600, [0; energyLoss6]./3600, '.k');
plot(t7/3600, [0; energyLoss7]./3600, '.k');
plot(t8/3600, [0; energyLoss8]./3600, '.k');
plot(t9/3600, [0; energyLoss9]./3600, '.k');
plot(t10/3600, [0; energyLoss10]./3600, '.k');
plot(t11/3600, [0; energyLoss11]./3600, '.k');
plot(t12/3600, [0; energyLoss12]./3600, '.k');

ylabel('Total Energy Loss ($kWh_{t}$)', 'interpreter', 'latex', 'FontSize', 14);
xlabel('$t$ (h)', 'interpreter', 'latex', 'FontSize', 14);
set(gca, 'Color', [1 1 1], 'TickLabelInterpreter', 'latex', 'FontSize', 14);
xlim([0, 165]);
ylim([0, 35]);

%% Sample contour plot 
kCharge = 31;
kHold = 160;
kDischarge = 220;
tCharge = FF_test1.Fo2t(FF_test1.Fo(kCharge), 1);
tHold = FF_test1.Fo2t(FF_test1.Fo(kHold), 1);
tDischarge = FF_test1.Fo2t(FF_test1.Fo(kDischarge), 1);
climits = [400, 600];

% Charge
loadTheta(FF_test1, FF_test1.ls);
knCharge = mod(kCharge-1, FF_test1.ls) + 1;
r_Charge = [-fliplr(FF_test1.theta{knCharge, 3}), FF_test1.theta{knCharge, 3}];
z_Charge = [FF_test1.theta{knCharge, 2}, max(FF_test1.theta{knCharge, 2}), ...
    1 + 1.000001*FF_test1.h];
thetaA_Charge = FF_test1.theta{knCharge, 12};
theta_Charge  = [fliplr(FF_test1.theta{knCharge, 1}), FF_test1.theta{knCharge, 1}; ...
                    thetaA_Charge*ones(2, length(r_Charge))];

% Hold
loadTheta(FF_test1, kHold);
knHold = mod(kHold-1, FF_test1.ls) + 1;
r_Hold = [-fliplr(FF_test1.theta{knHold, 3}), FF_test1.theta{knHold, 3}];
z_Hold = [FF_test1.theta{knHold, 2}, max(FF_test1.theta{knHold, 2}), ...
    1 + 1.000001*FF_test1.h];
thetaA_Hold = FF_test1.theta{knHold, 12};
theta_Hold  = [fliplr(FF_test1.theta{knHold, 1}), FF_test1.theta{knHold, 1}; ...
                    thetaA_Hold*ones(2, length(r_Hold))];

% Discharge
loadTheta(FF_test1, kDischarge);
knDischarge = mod(kDischarge-1, FF_test1.ls) + 1;
r_Discharge = [-fliplr(FF_test1.theta{knDischarge, 3}), FF_test1.theta{knDischarge, 3}];
z_Discharge = [FF_test1.theta{knDischarge, 2}, max(FF_test1.theta{knDischarge, 2}), ...
    1 + 1.000001*FF_test1.h];
thetaA_Discharge = FF_test1.theta{knDischarge, 12};
theta_Discharge  = [fliplr(FF_test1.theta{knDischarge, 1}), FF_test1.theta{knDischarge, 1}; ...
                    thetaA_Discharge*ones(2, length(r_Discharge))];

fzri = figure('Units', 'normalized', 'Position', 0.2*[0 0 2*2*FF_test1.b 1]);
format compact;

% plot charging contour
pos1 = [0.1300-0.01 0.1100 0.2134 0.8150];
subplot(1, 3, 1, 'Position', pos1);
[R, Z] = meshgrid(r_Charge, z_Charge); 
temp = FF_test1.theta2T(theta_Charge);

pzri1 = surf(R, Z, temp, 'EdgeColor', 'interp', 'FaceColor', 'interp');
caxis(climits);
xlim([min(R(:)), max(R(:))])
ylim([min(Z(:)), max(Z(:))])
pbaspect([2*FF_test1.b, 1, 1]);  % figure sized proportional to aspect ratio
view(0, 90);
colormap(jet);
set(pzri1, 'linestyle', 'none');
title(sprintf('$t$ = %1.2f h', tCharge/3600), 'interpreter', 'latex', ...
                            'FontSize', 14);
set(gca, 'Color', [1 1 1], 'TickLabelInterpreter', 'latex', 'FontSize', 14);
ylabel('$z/H$', 'FontWeight', 'bold', 'interpreter', 'latex', 'FontSize', 14);

% plot holding contour
pos2 = [0.4108-0.05 0.1100 0.2134 0.8150];
subplot(1, 3, 2, 'Position', pos2);
[R, Z] = meshgrid(r_Hold, z_Hold); 
temp = FF_test1.theta2T(theta_Hold);

pzri1 = surf(R, Z, temp, 'EdgeColor', 'interp', 'FaceColor', 'interp');
caxis(climits);
xlim([min(R(:)), max(R(:))])
ylim([min(Z(:)), max(Z(:))])
pbaspect([2*FF_test1.b, 1, 1]);  % figure sized proportional to aspect ratio
view(0, 90);
colormap(jet);
set(pzri1, 'linestyle', 'none');
title(sprintf('$t$ = %1.2f h', tHold/3600), 'interpreter', 'latex', ...
                            'FontSize', 14);
set(gca, 'Color', [1 1 1], 'TickLabelInterpreter', 'latex', 'FontSize', 14, ...
    'YTickLabel',[]);
xlabel('$r/H$', 'FontWeight', 'bold', 'interpreter', 'latex', 'FontSize', 14);


% plot discharging contour
pos3 = [0.6916-0.09 0.1100 0.2134 0.8150];
subplot(1, 3, 3, 'Position', pos3);
[R, Z] = meshgrid(r_Discharge, z_Discharge); 
temp = FF_test1.theta2T(theta_Discharge);

pzri1 = surf(R, Z, temp, 'EdgeColor', 'interp', 'FaceColor', 'interp');
caxis(climits);
xlim([min(R(:)), max(R(:))])
ylim([min(Z(:)), max(Z(:))])
pbaspect([2*FF_test1.b, 1, 1]);  % figure sized proportional to aspect ratio
view(0, 90);
set(pzri1, 'linestyle', 'none');
title(sprintf('$t$ = %1.2f h', tDischarge/3600), 'interpreter', 'latex', ...
                            'FontSize', 14); 
set(gca, 'Color', [1 1 1], 'TickLabelInterpreter', 'latex', 'FontSize', 14, ...
    'YTickLabel',[]);

h = axes(fzri, 'visible', 'off');
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';
% ylabel(h, '$z/H$', 'FontWeight', 'bold', 'interpreter', 'latex', 'FontSize', 14);
% xlabel(h, '$r/H$', 'FontWeight', 'bold', 'interpreter', 'latex', 'FontSize', 14);
c = colorbar(h, 'Position', [0.93-0.1 0.168 0.022 0.7]);
c.Label.String = '$T$ ($^\circ$C)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 14;
c.TickLabelInterpreter = 'latex'; 
colormap(h, 'jet');
caxis(h, climits);
