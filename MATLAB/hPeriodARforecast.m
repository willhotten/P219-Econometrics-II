function yf = hPeriodARforecast(y,d,Phi,p,h)
    [T,n] = size(y);
    F_companion = [Phi;[eye(p-1) zeros(p-1 ,1)]];
    c_companion = [d; zeros(p-1 ,1)];
    F_Power = zeros(size(F_companion));
    for j = 1:h
        F_Power = F_Power + F_companion^(j-1);
    end
    A = F_Power*c_companion ;
    B = F_companion^h ;
    yWithLags = [y mlag2(y,p-1)];
    yf = repmat (A',T,1) + yWithLags*B';
    yf = yf (:,1);
end