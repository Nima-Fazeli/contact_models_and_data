function [xn,xp,y] = drawE( m11, m12, m22, p, v2)
    
    qu = linspace(-1.05*v2,1.05*v2,200);
    
    y = zeros(length(qu),1);
    xn = y;
    xp = y;
    
    for i=1:length(qu)
        y(i) = qu(i);
        rs = roots([m11,2*m12*qu(i),m22*qu(i)^2-p]);
        xn(i) = rs(1);
        xp(i) = rs(2);
    end
end



