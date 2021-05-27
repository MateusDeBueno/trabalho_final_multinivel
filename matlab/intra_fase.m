clear all
close all
clc

N = 7; %numero de niveis
Vdc = 1;
V = (0:1:N-1)*Vdc - floor(N/2); % opcoes de tensoes de fase

%cria as tens es de fase
for j = 1:N
    Va(j) = V(j);
    Vb(j) = V(j);
    Vc(j) = V(j);
end

n_vetores = 0;
vetor_alpha = zeros(1,N^3,1); %cria vetor para alphas
vetor_beta = zeros(1,N^3,1); %cria vetor para betas
vetor_gama = zeros(1,N^3,1); %cria vetor para gamas
vetor = string(1:N^3); %vetor de string para salvar as tensoes de linha que geram a cada vetor

%preenche o vetor alpha e o vetor beta com seus respectivos valores
for k = 1:N %varre a tens o Va
    for j = 1:N %varre a tens o Vb
        for i = 1:N %varre a tens o Vc
            n_vetores = n_vetores+1; %quantos vetores existem
            [vetor_alpha(n_vetores), vetor_beta(n_vetores), vetor_gama(n_vetores)] = transformada_clarke(Va(k), Vb(j), Vc(i));
            char_Va = int2str(Va(k)/Vdc);
            char_Vb = int2str(Vb(j)/Vdc);
            char_Vc = int2str(Vc(i)/Vdc);
            vetor(n_vetores) = append(char_Va, char_Vb, char_Vc); %salvando qual conjuento de tensoes de fase gera qual vetor
        end
    end
end


vetor_alpha = round(vetor_alpha,10); %arrendonda com 10 casa decimais, estava dando problema com a funcao unique
vetor_beta = round(vetor_beta,10);

num_redundancias = zeros(1,n_vetores);
matrix_vector = [vetor_alpha', vetor_beta', num_redundancias']; %primeira coluna alpha, segunda coluna beta, terceira coluna num redundancias
matrix_uniq = unique(matrix_vector, 'rows'); %pega somente os unicos
dados = num2cell(matrix_uniq,1); %transforma a matrix em celula, para poder salvar string junto com numeros
n_vet_unic = 3*N*(N-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix

for z = 1:n_vetores
    alpha_talvez_unico = vetor_alpha(z);
    beta_talvez_unico = vetor_beta(z);
    for j = 1:n_vet_unic
        if (alpha_talvez_unico == dados{1}(j) && beta_talvez_unico == dados{2}(j))
            dados{3}(j) = dados{3}(j)+1; %contando o numero de redundancias por vetor
            dados{4}(j,dados{3}(j)) = vetor(z); %salva a combinacao de tensao de linha que leva a esse vetor
        end
    end
end

%plot do numero de redundancias
figure
scatter(dados{1},dados{2},25,dados{3},'filled')
stg = blanks(1);
for z = 1:n_vet_unic
    for i=1:dados{3}(z,1)
        stg = append(stg, newline,dados{4}(z,i));
    end
    text(dados{1}(z),dados{2}(z),stg,'FontSize',7)
    stg = erase(stg,stg);
end
grid on
xlabel('Alpha')
ylabel('Beta')
save_figure("joint_phase")


figure
plot3(vetor_alpha, vetor_beta, vetor_gama, 'o')
% hold on
% [X,Y] = meshgrid(-4:0.5:4,-4:0.5:4);
% Z = X*0+3;
% surf(X,Y,Z)
% hold off
xlabel('Alpha')
ylabel('Beta')
zlabel('Gama')






figure
gscatter(dados{1},dados{2},dados{3})
stg = blanks(1);
for z = 1:n_vet_unic
    numero_redundancias = int2str(dados{3}(z,1));
    text(dados{1}(z),dados{2}(z),numero_redundancias,'FontSize',11)
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Numero de redundancias para cada posição no mapa Alpha Beta')



