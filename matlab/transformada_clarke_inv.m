function [Va,Vb,Vc] = transformada_clarke_inv(a,b,g)
    Va = a + g;
    Vb = -0.5*a + sqrt(3)*0.5*b + g;
    Vc = -0.5*a - sqrt(3)*0.5*b + g;
end

