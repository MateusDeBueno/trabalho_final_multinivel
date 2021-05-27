function [alpha,beta,gama] = transformada_clarke(Va,Vb,Vc)
alpha = 2/3*(Va-0.5*Vb-0.5*Vc);
beta = 2/3*(sqrt(3)/2*Vb-sqrt(3)/2*Vc);
gama = 2/3*(Va/sqrt(2)+Vb/sqrt(2)+Vc/sqrt(2));
end