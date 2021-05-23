clear all
close all
clc

% N = 7; %numero de niveis
% Vdc = 1;
% 
% Vbig = [-2 0 2];
% Vsmall = [-1 0 1];
% size_b = length(Vbig); %tamanho do grande
% size_s = length(Vsmall); %tamanho do pequeno
% size_t = size_b*size_s; %tamanho total
% vetor = string(1:size_t); %vetor de string para salvar as tensoes de linha que levam a cada vetor
% V = string(1:size_t);
% n = 1;
% %vetores de uma fase
% for k = 1:size_b %varre a tensao Vbig
%     for j = 1:size_s %varre a tensao Vsmall
%         V(n) = Vbig(k) + Vsmall(j);
%         char_Vbig = int2str(Vbig(k));
%         char_Vsmall = int2str(Vsmall(j));
%         vetor(n) = append(char_Vbig, char_Vsmall);
%         n=n+1;
%     end
% end

Vdc_big = 2;
Nbig = 3;
Vbig = (0:1:Nbig-1)*Vdc_big;  % opcoes de tensoes de fase


%cria os vetores das tensoes de fase
Va_big = zeros(1,Nbig);
Vb_big = zeros(1,Nbig);
Vc_big = zeros(1,Nbig);
for j = 1:Nbig %atribui os valores para as tensoes de fase
    Va_big(j) = Vbig(j);
    Vb_big(j) = Vbig(j);
    Vc_big(j) = Vbig(j);
end


n_vetores_big = 0; %conta o numero de vetores totais, deve ser igual a nˆ3
vetor_alpha_big = zeros(1,Nbig^3,1); %cria vetor para alphas
vetor_beta_big = zeros(1,Nbig^3,1); %cria vetor para betas
vetor_gama_big = zeros(1,Nbig^3,1); %cria vetor para betas
vetor_big = string(1:Nbig^3);
%preenche o vetor alpha e o vetor beta com seus respectivos valores
for k = 1:Nbig %varre a tensao Va
    for j = 1:Nbig %varre a tensao Vb
        for i = 1:Nbig %varre a tensao Vc
        n_vetores_big = n_vetores_big+1; %quantos vetores existem
        [vetor_alpha_big(n_vetores_big), vetor_beta_big(n_vetores_big), vetor_gama_big(n_vetores_big)] = transformada_clarke(Va_big(k), Vb_big(j), Vc_big(i));

        %salva o char somente para saber quais tensoes estavam
        %aplicadas
        char_Va_big = int2str(Va_big(k));
        char_Vb_big = int2str(Vb_big(j));
        char_Vc_big = int2str(Vc_big(i));
        vetor_big(n_vetores_big) = append(char_Va_big, char_Vb_big, char_Vc_big); %salvando qual conjuento de tensoes de fase gera qual vetor
        end
    end
end



num_redundancias_big = zeros(1,n_vetores_big);
matrix_vector_big = [vetor_alpha_big', vetor_beta_big', num_redundancias_big']; %primeira coluna alpha, segunda coluna beta
matrix_uniq_big = unique(matrix_vector_big, 'rows');
dados_big = num2cell(matrix_uniq_big,1);

n_vet_unic_big = 3*Nbig*(Nbig-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix


for z = 1:n_vetores_big
    alpha_talvez_unico_big = vetor_alpha_big(z);
    beta_talvez_unico_big = vetor_beta_big(z);
    for j = 1:n_vet_unic_big
        if (alpha_talvez_unico_big == dados_big{1}(j) && beta_talvez_unico_big == dados_big{2}(j))
            dados_big{3}(j) = dados_big{3}(j)+1; %contando o numero de redundancias por vetor
            dados_big{4}(j,dados_big{3}(j)) = vetor_big(z);
        end
    end
end




%plot do numero de redundancias
figure
scatter(dados_big{1},dados_big{2},25,dados_big{3},'filled')
stg = blanks(1);
for z = 1:n_vet_unic_big
    for i=1:dados_big{3}(z,1)
        stg = append(stg, newline,dados_big{4}(z,i));
    end
    text(dados_big{1}(z),dados_big{2}(z),stg,'FontSize',7)
    stg = erase(stg,stg);
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Gráfico alpha beta da celular maior')


figure
plot3(vetor_alpha_big, vetor_beta_big, vetor_gama_big, 'o')
xlabel('Alpha')
ylabel('Beta')
zlabel('Gama')
title('Gráfico alpha beta gama da celular maior')



%% small
Vdc_small = 1;
Nsmall = 3;
Vsmall = (0:1:Nsmall-1)*Vdc_small;  % opcoes de tensoes de fase


%cria os vetores das tensoes de fase
Va_small = zeros(1,Nsmall);
Vb_small = zeros(1,Nsmall);
Vc_small = zeros(1,Nsmall);
for j = 1:Nsmall %atribui os valores para as tensoes de fase
    Va_small(j) = Vsmall(j);
    Vb_small(j) = Vsmall(j);
    Vc_small(j) = Vsmall(j);
end


n_vetores_small = 0; %conta o numero de vetores totais, deve ser igual a nˆ3
vetor_alpha_small = zeros(1,Nsmall^3,1); %cria vetor para alphas
vetor_beta_small = zeros(1,Nsmall^3,1); %cria vetor para betas
vetor_gama_small = zeros(1,Nsmall^3,1); %cria vetor para betas
vetor_small = string(1:Nsmall^3);
comprimento = zeros(1,Nsmall^3,1);

%preenche o vetor alpha e o vetor beta com seus respectivos valores
for k = 1:Nsmall %varre a tensao Va
    for j = 1:Nsmall %varre a tensao Vb
        for i = 1:Nsmall %varre a tensao Vc
        n_vetores_small = n_vetores_small+1; %quantos vetores existem
        [vetor_alpha_small(n_vetores_small), vetor_beta_small(n_vetores_small), vetor_gama_small(n_vetores_small)] = transformada_clarke(Va_small(k), Vb_small(j), Vc_small(i));
        comprimento(n_vetores_small) = sqrt(vetor_alpha_small(n_vetores_small)^2+vetor_beta_small(n_vetores_small)^2);

        %salva o char somente para saber quais tensoes estavam
        %aplicadas
        char_Va_small = int2str(Va_small(k));
        char_Vb_small = int2str(Vb_small(j));
        char_Vc_small = int2str(Vc_small(i));
        vetor_small(n_vetores_small) = append(char_Va_small, char_Vb_small, char_Vc_small); %salvando qual conjuento de tensoes de fase gera qual vetor
        end
    end
end


num_redundancias_small = zeros(1,n_vetores_small);
matrix_vector_small = [vetor_alpha_small', vetor_beta_small', num_redundancias_small', comprimento'];
%coluna 1 = alpha
%coluna 2 = beta
%coluna 3 = num redundancias
%coluna 4 = comprimento
matrix_uniq_small = unique(matrix_vector_small, 'rows');
dados_small = num2cell(matrix_uniq_small,1);


n_vet_unic_small = 3*Nsmall*(Nsmall-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix


for z = 1:n_vetores_small
    alpha_talvez_unico_small = vetor_alpha_small(z);
    beta_talvez_unico_small = vetor_beta_small(z);
    for j = 1:n_vet_unic_small
        if (alpha_talvez_unico_small == dados_small{1}(j) && beta_talvez_unico_small == dados_small{2}(j))
            dados_small{3}(j) = dados_small{3}(j)+1; %contando o numero de redundancias por vetor
            dados_small{5}(j,dados_small{3}(j)) = vetor_small(z);
        end
    end
end


%plot do numero de redundancias
figure
scatter(dados_small{1},dados_small{2},25,dados_small{3},'filled')
stg = blanks(1);
for z = 1:n_vet_unic_small
    for i=1:dados_small{3}(z,1)
        stg = append(stg, newline,dados_small{5}(z,i));
    end
    text(dados_small{1}(z),dados_small{2}(z),stg,'FontSize',7)
    stg = erase(stg,stg);
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Gráfico alpha beta da celular menor')


figure
plot3(vetor_alpha_small, vetor_beta_small, vetor_gama_small, 'o')
xlabel('Alpha')
ylabel('Beta')
zlabel('Gama')
title('Gráfico alpha beta gama da celular menor')




figure
for k=1:size(dados_small{1})-1
    hold on
    plot([dados_small{1}(k),dados_small{2}(k)],[dados_small{1}(k+1),dados_small{2}(k+1)])
    hold off
end

cnt=0;
for k=1:size(dados_small{1})
    if (max(dados_small{4})==dados_small{4}(k))
        %salva o ponto
        cnt=cnt+1;
        alpha_externo(cnt)=dados_small{1}(k);
        beta_externo(cnt)=dados_small{2}(k);
    end
end

figure
patch(alpha_externo,beta_externo,'red')

figure
scatter(beta_externo,alpha_externo,'red')
