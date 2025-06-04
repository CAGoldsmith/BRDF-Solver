close all
clear all

% addr = 'C:\Users\Clarissa G\Documents\MATLAB\Design Drosophibot\Scale Drosophila Male\k_spring = 3.8e-08\T_swing = 0.1\allVars.mat';
addr = 'C:\Users\Clarissa G\Documents\MATLAB\Design Drosophibot\Scale Drosophila Male\Flat\k_spring = 3.8e-08\T_swing = 0.1\phi_I = 0.5\phi_C = 0.5\allVars.mat';
load(addr)
typX = 10e-6*ones(3,1);
% opts = optimoptions(@fsolve,'StepTolerance',1e-10,'OptimalityTolerance',1e-10,'TypicalX',typX);
opts = optimoptions(@fsolve,'StepTolerance',1e-20,'OptimalityTolerance',1e-20,'FunctionTolerance',1e-20,'Display','off');

jointNames = {'Thorax-Coxa','Coxa-Trochanter','Trochanter-Femur','Femur-Tibia','Tibia-Tarsus'};
colors = {'#0072BD','#D95319','#EDB120'};
timeVectorS = linspace(0,swingDuration/stanceDuty,n);
saveAddr = 'C:\Users\Clarissa G\Documents\MATLAB\Design Drosophibot\Scale Drosophila Male\Flat\k_spring = 3.8e-08\T_swing = 0.1\phi_I = 0.5\phi_C = 0.5\BrianData\Combined Forces';
if ~exist(saveAddr)
    mkdir(saveAddr)
end

% saveAddr = 'G:\Other computers\NeuroMINT Lab Computer\MATLAB\Design Drosophibot\Scale Drosophila Male\k_spring = 3.8e-08\T_swing = 0.1\BrianData\Combined Forces\';
ranges = {['A2:D' num2str(2+n)], ['F2:I' num2str(2+n)], ['K2:N' num2str(2+n)], ['P2:S' num2str(2+n)], ['U2:X' num2str(2+n)]};
titleRanges = {'A1:D1', 'F1:I1', 'K1:N1', 'P1:S1', 'U1:X1'};

for l = 1:6
    fprintf(['Beginning ' legNames{l} '\n'])
    for s = 1:numSegs(l)
        fprintf(['Segment ' num2str(s) '\n'])
        for t = 1:n
            rVec = jointPoints{t}{l}(:,s+1) - jointPoints{t}{l}(:,s);
            solverEq = @(x) tauForceEq(Tjoint{l}{s}(:,t),rVec,x);
            x0 = [0;0;0];
            tauForce{l}{s}(:,t) = (fsolve(solverEq,x0,opts))*(10^-6);
        end
        totForce{l}{s} = tauForce{l}{s} + Fjoint{l}{s};
    end
    for ss = 1:5
        if ss == 1
            totForceComb{l}{ss} = totForce{l}{1} + totForce{l}{2} + totForce{l}{3};
        else
            totForceComb{l}{ss} = totForce{l}{ss+2};
        end
        for t=1:n
            totForceCombMags{l}{ss}(t) = norm(totForceComb{l}{ss}(:,t));
        end
        totForceCombMax(l,ss) = max(totForceCombMags{l}{ss});

        figure
        tl{l,ss} = tiledlayout(3,1,'TileSpacing','compact','Padding','loose');
        title(tl{l,ss},[legNames{l} ' ' jointNames{ss}],'Interpreter','tex','FontWeight', 'bold')
        subtitle(tl{l,ss},{[springKPlotTitle ' ' '^{Nm}/_{rad}' ',   ' 'T_{step}= ' num2str(swingDuration/stanceDuty) ' s'], ['Max. Force Magnitude: ' num2str(totForceCombMax(l,ss)) ' N']},'Interpreter','tex')
        nexttile
        plot(timeVectorS,totForceComb{l}{ss}(1,:),'LineWidth',1.5,'Color',colors{1})
        grid on
        xlim([0 timeVectorS(end)])
        xlabel('Time (s)')
        ylabel({'X Component (N)','Post <-- --> Ant'})

        nexttile
        plot(timeVectorS,totForceComb{l}{ss}(2,:),'LineWidth',1.5,'Color',colors{2})
        grid on
        xlim([0 timeVectorS(end)])
        xlabel('Time (s)')
        ylabel({'Y (Lateral) Component (N)'})

        nexttile
        plot(timeVectorS,totForceComb{l}{ss}(3,:),'LineWidth',1.5,'Color',colors{3})
        grid on
        xlim([0 timeVectorS(end)])
        xlabel('Time (s)')
        ylabel({'Z Component (N)','Ventral <-- --> Dorsal'})
        saveas(tl{l,ss},[saveAddr '\' legNames{l} '_' jointNames{ss}])
        close gcf
        % tab{l,ss} = array2table([timeVectorS', totForceComb{l}{ss}'],'VariableNames',{'Time (s)','X Component (N)','Y Component (N)','Z Component (N)'});
        % writetable(tab{l,ss},[saveAddr '\forceData.xls'],'Sheet',l,'Range',ranges{ss})
        % writematrix(jointNames{ss},[saveAddr '\forceData.xls'],'Sheet',l,'Range',titleRanges{ss})
    end
end




