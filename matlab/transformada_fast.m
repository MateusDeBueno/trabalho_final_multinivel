function [g,h] = transformada_fast(Va,Vb,Vc)
Vab = Va-Vb;
Vbc = Vb-Vc;
Vca = Vc-Va;

g=(2.0*Vab - 1.0*Vbc - 1.0*Vca)/3.0;
h=(-1.0*Vab + 2.0*Vbc - 1.0*Vca)/3.0;
end