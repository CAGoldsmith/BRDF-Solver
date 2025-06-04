function F = quartCircEqs(x,prevX,prevY,prevZ,m,chordL,circR,circX,circY,circZ)
F(1) = sqrt((x(1)-prevX)^2 + (x(2)-prevY)^2 + (prevZ + m*(x(1)-prevX)-prevZ)^2)-chordL;
F(2) = (x(1)-circX)^2 + (x(2)-circY)^2 + (prevZ + m*(x(1)-prevX)-circZ)^2 - circR^2;
end