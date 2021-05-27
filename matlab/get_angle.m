function [ang] = get_angle(x,y)
%pega o angula de um ponto pra origem
ang = rad2deg(angle(x+y*1i));
end