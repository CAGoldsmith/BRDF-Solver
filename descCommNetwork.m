function [scalarErr,matrixErr,gTransMax,gRotMax,PEPposRest,PEPnegRest,PEP,PEPposMat,PEPnegMat] = descCommNetwork(kSynTrans,kSynRot,PEPposBias,PEPnegBias,invertTrans,invertRot,Erest,delEdep,R,vTrans,vRot,pep,toPlot)
    
    delEhyp = -100;
    gHyp = -R/delEhyp;
    vNegTrans = (Erest + R + min(max( (vTrans - Erest)/R, 0), 1)*gHyp*delEhyp)./(1 + min(max( (vTrans - Erest)/R, 0), 1)*gHyp);
    vNegRot = (Erest + R + min(max( (vRot - Erest)/R, 0), 1)*gHyp*delEhyp)./(1 + min(max( (vRot - Erest)/R, 0), 1)*gHyp);
    
    gTransMax = abs(kSynTrans)*R/(delEdep - abs(kSynTrans)*R);
    if invertTrans
        gTransPos = gTransMax*min(max( (vNegTrans - Erest)/R, 0), 1);
        gTransNeg = gTransMax*min(max( (vTrans - Erest)/R, 0), 1);
    else
        gTransPos = gTransMax*min(max( (vTrans - Erest)/R, 0), 1);
        gTransNeg = gTransMax*min(max( (vNegTrans - Erest)/R, 0), 1);
    end
    
    gRotMax = abs(kSynRot)*R/(delEdep - abs(kSynRot)*R);
    if invertRot
        gRotPos = gRotMax*min(max( (vNegRot - Erest)/R, 0), 1);
        gRotNeg = gRotMax*min(max( (vRot - Erest)/R, 0), 1);
    else
        gRotPos = gRotMax*min(max( (vRot - Erest)/R, 0), 1);
        gRotNeg = gRotMax*min(max( (vNegRot - Erest)/R, 0), 1);
    end
    
    bmZero = zeros(size(vTrans));
    bmZero(ceil(length(vTrans(:))/2)) = 1;
    
    PEPpos = @(I)(I + gTransPos*delEdep + gRotPos*delEdep)./(1 + gTransPos + gRotPos);
    PEPneg = @(I)(I + gTransNeg*delEdep + gRotNeg*delEdep)./(1 + gTransNeg + gRotNeg);
    
    biasPos = @(I) sum(PEPpos(I).*bmZero,'all');
    biasNeg = @(I) sum(PEPneg(I).*bmZero,'all');
    
    try PEPposRest = fzero(biasPos,PEPposBias);
    catch
        keyboard
    end
    PEPnegRest = fzero(biasNeg,PEPnegBias);
    
    gDepCopy = R/(delEdep - R);
    PEPposMat = PEPpos(PEPposRest);
    PEPnegMat = PEPneg(PEPnegRest);
    gPos = gDepCopy*min(max( (PEPpos(PEPposRest) - Erest)/R, 0), 1);
    gNeg = gDepCopy*min(max( (PEPneg(PEPnegRest) - Erest)/R, 0), 1);
    
    PEP = (gPos*delEdep + gNeg*delEdep)./(1 + gPos + gNeg);
    
    if toPlot
        figure
        surf(vTrans,vRot,PEPpos(PEPposRest))
        hold on
        surf(vTrans,vRot,PEPneg(PEPnegRest))
        drawnow

        figure
        surf(vTrans,vRot,pep)
        title('desired pep data')

        figure
        surf(vTrans,vRot,PEP)
        title('network''s pep data')
        
        drawnow
    end
    
    matrixErr = pep - PEP;
    scalarErr = 1/2*sum(matrixErr.^2,'all');
end