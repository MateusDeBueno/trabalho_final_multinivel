function [g,h] = reta_angulo_gh(g_i,h_i,L,angul)
%dado um angulo e um comprimento de uma linha no plano alpha beta, retorna
%o ponto g e h para essa reta
g(1) = g_i;
h(1) = h_i;
a = g_i + L * cosd(angul);
b = h_i + L * sind(angul);
[Va,Vb,Vc] = transformada_clarke_inv(a,b,0);
[g(2), h(2)] = transformada_fast(Va,Vb,Vc);
end

