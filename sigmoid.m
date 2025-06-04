close all
clear all

thetaMax = 1;
thetaMin = -1;
thetaRest = 0;

steep = 250;
R = 20e-3;

%Always equal to thetaRest when Uext = Uflx (symmetrical only)
% amp = thetaMax-thetaMin;
% offset = thetaRest;
% thetaExt = @(U) amp./(1 + exp(steep*(20e-3-U)))+offset;
% thetaFlx = @(U) -amp./(1 + exp(steep*(20e-3-U)))+offset;

%theta = thetaMax when Uext = R, Uflx = 0, and vice-versa.
ampExt = 2*(thetaMax-thetaRest);
ampFlx = 2*(thetaMin-thetaRest);
offset = thetaRest/2;
thetaExt = @(U,S) ampExt./(1 + exp(S*(20e-3-U)))+offset;
thetaFlx = @(U,S) ampFlx./(1 + exp(S*(20e-3-U)))+offset;

Utest = 0:.1e-3:R;

figure
plot(Utest,thetaExt(Utest,steep));
hold on
plot(Utest,thetaFlx(Utest,steep));
xlabel('U')
ylabel('\theta')
legend('ext','flx')


[Uext,Uflx] = meshgrid(Utest,Utest);
UextV = Uext(:);
UflxV = Uflx(:);

figure
surf(Uext,Uflx,thetaExt(Uext,steep)+thetaFlx(Uflx,steep),'edgealpha',0)
xlabel('U_{ext}')
ylabel('U_{flx}')
view(0,90)
axis equal

negDiff = (thetaRest - thetaMin);
posDiff = (thetaMax - thetaRest);
negLevs = thetaMin:negDiff/5:thetaRest;
posLevs = thetaRest:posDiff/5:thetaMax;
levs = union(negLevs,posLevs);

figure
contour(Uext,Uflx,thetaExt(Uext,steep)+thetaFlx(Uflx,steep),levs)
hold on
% plot(
xlabel('U_{ext}')
ylabel('U_{flx}')
axis equal

f = @(S) thetaExt(Utest,steep) + thetaFlx(Utest,S) - thetaRest;
g = @(S) f(S)*f(S)';

sTry = 100:25:1000;
gTry = zeros(size(sTry));
for i=1:length(sTry)
    gTry(i) = g(sTry(i));
end

figure
plot(sTry,gTry)
xlabel('steepness')
ylabel('err')

% steep = sTry(gTry == min(gTry));
steepNew = fminunc(g,steep);
hold on
plot(steepNew,g(steepNew),'rx')
title(sprintf('steepness = %3.3f',steepNew))


figure
contour(Uext,Uflx,thetaExt(Uext,steep)+thetaFlx(Uflx,steepNew),levs)
hold on
% plot(
xlabel('U_{ext}')
ylabel('U_{flx}')
axis equal

thExt = thetaExt(Utest,steep) + offset;
figure
plot(Utest,thExt);
hold on

thTest = linspace(min(thExt),max(thExt));
UextCom = R-1./steep.*log(ampExt./(thTest - thetaRest) - 1);

hold on
plot(thTest,UextCom)