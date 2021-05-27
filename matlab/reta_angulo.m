function [conj_x,conj_y] = reta_angulo(x,y,L,alpha)
%desenha reta a partir do angulo e do comprimento
x2=x+(L*cosd(alpha));
y2=y+(L*sind(alpha));
conj_x = [x x2];
conj_y = [y y2];
end

