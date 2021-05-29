function [x,y] = reta_angulo(x_i,y_i,lineLength,alpha)
%desenha reta a partir do angulo e do comprimento
x(1) = x_i;
y(1) = y_i;
x(2) = x(1) + lineLength * cosd(alpha);
y(2) = y(1) + lineLength * sind(alpha);
end

