clear all
close all
clc

%numero de celulas por fase
num_celula = 2;
nivel = 2^(num_celula+1)-1;
%essa calculo so serve para verificar numero de niveis por fase

% %% ANALISE DO GRANDE
% %CRIANDO TENSOES DE CADA CELULA
% %CELULA GRANDE
% Vdc_big = 2;
% Nbig = 3;
% Vbig = (0:1:Nbig-1)*Vdc_big-2;  % opcoes de tensoes de fase
% 
% %cria os vetores das tensoes de fase
% Va_big = zeros(1,Nbig);
% Vb_big = zeros(1,Nbig);
% Vc_big = zeros(1,Nbig);
% for j = 1:Nbig %atribui os valores para as tensoes de fase
%     Va_big(j) = Vbig(j);
%     Vb_big(j) = Vbig(j);
%     Vc_big(j) = Vbig(j);
% end
% 
% %cria vetores necessarios
% n_vetores_big = 0;
% vetor_alpha_big = zeros(1,Nbig^3,1); %cria vetor para alphas
% vetor_beta_big = zeros(1,Nbig^3,1); %cria vetor para betas
% vetor_gama_big = zeros(1,Nbig^3,1); %cria vetor para gamas
% vetor_big = string(1:Nbig^3); %vetor de string para salvar as tensoes de linha que geram a cada vetor
% comprimento = zeros(1,Nbig^3,1); %cria vetor para gamas
% 
% %preenche o vetor alpha e o vetor beta com seus respectivos valores
% for k = 1:Nbig %varre a tens o Va
%     for j = 1:Nbig %varre a tens o Vb
%         for i = 1:Nbig %varre a tens o Vc
%             n_vetores_big = n_vetores_big+1; %quantos vetores existem
%             [vetor_alpha_big(n_vetores_big), vetor_beta_big(n_vetores_big), vetor_gama_big(n_vetores_big)] = transformada_clarke(Va_big(k), Vb_big(j), Vc_big(i));
%             char_Va_big = int2str(Va_big(k));
%             char_Vb_big = int2str(Vb_big(j));
%             char_Vc_big = int2str(Vc_big(i));
%             vetor_big(n_vetores_big) = append(char_Va_big, char_Vb_big, char_Vc_big); %salvando qual conjuento de tensoes de fase gera qual vetor
%             comprimento(n_vetores_big) = ((vetor_alpha_big(n_vetores_big))^2 + (vetor_beta_big(n_vetores_big))^2)^0.5;
%         end
%     end
% end
% 
% 
% vetor_alpha_big = round(vetor_alpha_big,10); %arrendonda com 10 casa decimais, estava dando problema com a funcao unique
% vetor_beta_big = round(vetor_beta_big,10);
% comprimento = round(comprimento,10);
% 
% num_redundancias_big = zeros(1,n_vetores_big);
% matrix_vector_big = [vetor_alpha_big', vetor_beta_big', num_redundancias_big', comprimento']; %primeira coluna alpha, segunda coluna beta, terceira coluna num redundancias
% matrix_uniq_big = unique(matrix_vector_big, 'rows'); %pega somente os unicos
% dados_big = num2cell(matrix_uniq_big,1); %transforma a matrix em celula, para poder salvar string junto com numeros
% n_vet_unic_big = 3*Nbig*(Nbig-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix
% %MATRIX COMPLETA
% %coluna 1 = alpha
% %coluna 2 = beta
% %coluna 3 = num redundancias
% %coluna 4 = comprimento
% 
% 
% for z = 1:n_vetores_big
%     alpha_talvez_unico_big = vetor_alpha_big(z);
%     beta_talvez_unico_big = vetor_beta_big(z);
%     for j = 1:n_vet_unic_big
%         if (alpha_talvez_unico_big == dados_big{1}(j) && beta_talvez_unico_big == dados_big{2}(j))
%             dados_big{3}(j) = dados_big{3}(j)+1; %contando o numero de redundancias por vetor
%             dados_big{5}(j,dados_big{3}(j)) = vetor_big(z); %salva a combinacao de tensao de linha que leva a esse vetor
%         end
%     end
% end
% 
% 
% %plot do numero de redundancias
% figure
% hold on
% 
% scatter(dados_big{1},dados_big{2},25,dados_big{3},'filled')
% stg = blanks(1);
% for z = 1:n_vet_unic_big
%     for i=1:dados_big{3}(z,1)
%         stg = append(stg, newline,dados_big{5}(z,i));
%     end
%     text(dados_big{1}(z),dados_big{2}(z),stg,'FontSize',7)
%     stg = erase(stg,stg);
% end
% grid on
% xlabel('Alpha')
% ylabel('Beta')
% % for k=1:12
% %     angle = k*30;
% %     [conj_x,conj_y] = reta_angulo(0,0,2,angle);
% %     plot(conj_x,conj_y)
% % end
% hold off

%CELULA GRANDE
Vdc_big = 2;
Nbig = 3;
Vbig = (0:1:Nbig-1)*Vdc_big-1;  % opcoes de tensoes de fase

%cria os vetores das tensoes de fase
Va_big = zeros(1,Nbig);
Vb_big = zeros(1,Nbig);
Vc_big = zeros(1,Nbig);
for j = 1:Nbig %atribui os valores para as tensoes de fase
    Va_big(j) = Vbig(j);
    Vb_big(j) = Vbig(j);
    Vc_big(j) = Vbig(j);
end

%cria vetores necessarios
n_vetores_big = 0;
vetor_alpha_big = zeros(1,Nbig^3,1); %cria vetor para alphas
vetor_beta_big = zeros(1,Nbig^3,1); %cria vetor para betas
vetor_gama_big = zeros(1,Nbig^3,1); %cria vetor para gamas
vetor_big = string(1:Nbig^3); %vetor de string para salvar as tensoes de linha que geram a cada vetor
comprimento_big = zeros(1,Nbig^3,1); %cria vetor para comprimento de cada vetor
angulo_big = zeros(1,Nbig^3,1); %cria vetor para angulo de cada vector

%preenche o vetor alpha e o vetor beta com seus respectivos valores
for k = 1:Nbig %varre a tens o Va
    for j = 1:Nbig %varre a tens o Vb
        for i = 1:Nbig %varre a tens o Vc
            n_vetores_big = n_vetores_big+1; %quantos vetores existem
            [vetor_alpha_big(n_vetores_big), vetor_beta_big(n_vetores_big), vetor_gama_big(n_vetores_big)] = transformada_clarke(Va_big(k), Vb_big(j), Vc_big(i));
            char_Va_big = int2str(Va_big(k));
            char_Vb_big = int2str(Vb_big(j));
            char_Vc_big = int2str(Vc_big(i));
            vetor_big(n_vetores_big) = append(char_Va_big, char_Vb_big, char_Vc_big); %salvando qual conjuento de tensoes de fase gera qual vetor
            comprimento_big(n_vetores_big) = ((vetor_alpha_big(n_vetores_big))^2 + (vetor_beta_big(n_vetores_big))^2)^0.5;
            angulo_big(n_vetores_big) = get_angle(vetor_alpha_big(n_vetores_big),vetor_beta_big(n_vetores_big));
        end
    end
end


vetor_alpha_big = round(vetor_alpha_big,10); %arrendonda com 10 casa decimais, estava dando problema com a funcao unique
vetor_beta_big = round(vetor_beta_big,10);
comprimento_big = round(comprimento_big,10);
angulo_big = round(angulo_big,10);

num_redundancias_big = zeros(1,n_vetores_big);
matrix_vector_big = [vetor_alpha_big', vetor_beta_big', num_redundancias_big', comprimento_big', angulo_big']; %primeira coluna alpha, segunda coluna beta, terceira coluna num redundancias
matrix_uniq_big = unique(matrix_vector_big, 'rows'); %pega somente os unicos
dados_big = num2cell(matrix_uniq_big,1); %transforma a matrix em celula, para poder salvar string junto com numeros
n_vet_unic_big = 3*Nbig*(Nbig-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix
%MATRIX COMPLETA
%coluna 1 = alpha
%coluna 2 = beta
%coluna 3 = num redundancias
%coluna 4 = comprimento
%coluna 5 = angulo

for z = 1:n_vetores_big
    alpha_talvez_unico_big = vetor_alpha_big(z);
    beta_talvez_unico_big = vetor_beta_big(z);
    for j = 1:n_vet_unic_big
        if (alpha_talvez_unico_big == dados_big{1}(j) && beta_talvez_unico_big == dados_big{2}(j))
            dados_big{3}(j) = dados_big{3}(j)+1; %contando o numero de redundancias por vetor
            dados_big{6}(j,dados_big{3}(j)) = vetor_big(z); %salva a combinacao de tensao de linha que leva a esse vetor
        end
    end
end

figure
hold on

scatter(dados_big{1},dados_big{2},25,dados_big{3},'filled')
stg = blanks(1);
for z = 1:n_vet_unic_big
    for i=1:dados_big{3}(z,1)
        stg = append(stg, newline,dados_big{6}(z,i));
    end
    text(dados_big{1}(z),dados_big{2}(z),stg,'FontSize',7)
    stg = erase(stg,stg);
end
grid on
xlabel('Alpha')
ylabel('Beta')
% for k=1:12
%     angle = k*30;
%     [conj_x,conj_y] = reta_angulo(0,0,2,angle);
%     plot(conj_x,conj_y)
% end
hold off


figure
scatter(dados_big{1},dados_big{2})
stg = blanks(1);
for z = 1:n_vet_unic_big
    angulo = int2str(dados_big{5}(z,1));
    comprimento = num2str(round(dados_big{4}(z,1),2));
    text(dados_big{1}(z),dados_big{2}(z),angulo,'FontSize',11)
    text(dados_big{1}(z),dados_big{2}(z)-0.2,comprimento,'FontSize',11)
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Angulo de cada vetor')

%pegar os valores longos de interesse
comprimeto_big = dados_big{4};
n_vetores_longos_big = ceil(Nbig/2);
max_big = zeros(1,n_vetores_longos_big);
for k=1:n_vetores_longos_big
    max_big(k) = max(comprimeto_big);
    comprimeto_big = comprimeto_big(comprimeto_big~=max_big(k));
end

%coloca na ordem crescente a partir do angulo (coluna 5)
matrix_uniq_angle_big = sortrows(matrix_uniq_big,5);
figure
hold on
cont_big=0;
scatter(dados_big{1},dados_big{2})
for k=1:n_vet_unic_big
    for j=1:n_vetores_longos_big
        if (matrix_uniq_angle_big(k,4)==max_big(j))
            scatter(matrix_uniq_angle_big(k,1),matrix_uniq_angle_big(k,2),'filled','red')
            cont_big=cont_big+1;
            text(matrix_uniq_angle_big(k,1),matrix_uniq_angle_big(k,2), int2str(cont_big))
            alphas_externos_big(cont_big) = matrix_uniq_angle_big(k,1);
            betas_externos_big(cont_big) = matrix_uniq_angle_big(k,2);
        end
    end
end
hold off

centering_a = 1;
centering_b = 2;
figure
hold on
h = fill(alphas_externos_big,betas_externos_big,'r');
set(h,'facealpha',.5)
h = fill(centering_a+alphas_externos_big,centering_b+betas_externos_big,'r');
set(h,'facealpha',.5)
hold off
%% ANALISE DO PEQUENO

%CELULA PEQUENA
Vdc_small = 1;
Nsmall = 3;
Vsmall = (0:1:Nsmall-1)*Vdc_small-1;  % opcoes de tensoes de fase

%cria os vetores das tensoes de fase
Va_small = zeros(1,Nsmall);
Vb_small = zeros(1,Nsmall);
Vc_small = zeros(1,Nsmall);
for j = 1:Nsmall %atribui os valores para as tensoes de fase
    Va_small(j) = Vsmall(j);
    Vb_small(j) = Vsmall(j);
    Vc_small(j) = Vsmall(j);
end

%cria vetores necessarios
n_vetores_small = 0;
vetor_alpha_small = zeros(1,Nsmall^3,1); %cria vetor para alphas
vetor_beta_small = zeros(1,Nsmall^3,1); %cria vetor para betas
vetor_gama_small = zeros(1,Nsmall^3,1); %cria vetor para gamas
vetor_small = string(1:Nsmall^3); %vetor de string para salvar as tensoes de linha que geram a cada vetor
comprimento_small = zeros(1,Nsmall^3,1); %cria vetor para comprimento de cada vetor
angulo_small = zeros(1,Nsmall^3,1); %cria vetor para angulo de cada vector

%preenche o vetor alpha e o vetor beta com seus respectivos valores
for k = 1:Nsmall %varre a tens o Va
    for j = 1:Nsmall %varre a tens o Vb
        for i = 1:Nsmall %varre a tens o Vc
            n_vetores_small = n_vetores_small+1; %quantos vetores existem
            [vetor_alpha_small(n_vetores_small), vetor_beta_small(n_vetores_small), vetor_gama_small(n_vetores_small)] = transformada_clarke(Va_small(k), Vb_small(j), Vc_small(i));
            char_Va_small = int2str(Va_small(k));
            char_Vb_small = int2str(Vb_small(j));
            char_Vc_small = int2str(Vc_small(i));
            vetor_small(n_vetores_small) = append(char_Va_small, char_Vb_small, char_Vc_small); %salvando qual conjuento de tensoes de fase gera qual vetor
            comprimento_small(n_vetores_small) = ((vetor_alpha_small(n_vetores_small))^2 + (vetor_beta_small(n_vetores_small))^2)^0.5;
            angulo_small(n_vetores_small) = get_angle(vetor_alpha_small(n_vetores_small),vetor_beta_small(n_vetores_small));
        end
    end
end


vetor_alpha_small = round(vetor_alpha_small,10); %arrendonda com 10 casa decimais, estava dando problema com a funcao unique
vetor_beta_small = round(vetor_beta_small,10);
comprimento_small = round(comprimento_small,10);
angulo_small = round(angulo_small,10);

num_redundancias_small = zeros(1,n_vetores_small);
matrix_vector_small = [vetor_alpha_small', vetor_beta_small', num_redundancias_small', comprimento_small', angulo_small']; %primeira coluna alpha, segunda coluna beta, terceira coluna num redundancias
matrix_uniq_small = unique(matrix_vector_small, 'rows'); %pega somente os unicos
dados_small = num2cell(matrix_uniq_small,1); %transforma a matrix em celula, para poder salvar string junto com numeros
n_vet_unic_small = 3*Nsmall*(Nsmall-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix
%MATRIX COMPLETA
%coluna 1 = alpha
%coluna 2 = beta
%coluna 3 = num redundancias
%coluna 4 = comprimento
%coluna 5 = angulo

for z = 1:n_vetores_small
    alpha_talvez_unico_small = vetor_alpha_small(z);
    beta_talvez_unico_small = vetor_beta_small(z);
    for j = 1:n_vet_unic_small
        if (alpha_talvez_unico_small == dados_small{1}(j) && beta_talvez_unico_small == dados_small{2}(j))
            dados_small{3}(j) = dados_small{3}(j)+1; %contando o numero de redundancias por vetor
            dados_small{6}(j,dados_small{3}(j)) = vetor_small(z); %salva a combinacao de tensao de linha que leva a esse vetor
        end
    end
end

figure
hold on

scatter(dados_small{1},dados_small{2},25,dados_small{3},'filled')
stg = blanks(1);
for z = 1:n_vet_unic_small
    for i=1:dados_small{3}(z,1)
        stg = append(stg, newline,dados_small{6}(z,i));
    end
    text(dados_small{1}(z),dados_small{2}(z),stg,'FontSize',7)
    stg = erase(stg,stg);
end
grid on
xlabel('Alpha')
ylabel('Beta')
% for k=1:12
%     angle = k*30;
%     [conj_x,conj_y] = reta_angulo(0,0,2,angle);
%     plot(conj_x,conj_y)
% end
hold off


figure
scatter(dados_small{1},dados_small{2})
stg = blanks(1);
for z = 1:n_vet_unic_small
    angulo = int2str(dados_small{5}(z,1));
    comprimento = num2str(round(dados_small{4}(z,1),2));
    text(dados_small{1}(z),dados_small{2}(z),angulo,'FontSize',11)
    text(dados_small{1}(z),dados_small{2}(z)-0.2,comprimento,'FontSize',11)
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Angulo de cada vetor')

%pegar os valores longos de interesse
comprimeto_small = dados_small{4};
n_vetores_longos_small = ceil(Nsmall/2);
max_small = zeros(1,n_vetores_longos_small);
for k=1:n_vetores_longos_small
    max_small(k) = max(comprimeto_small);
    comprimeto_small = comprimeto_small(comprimeto_small~=max_small(k));
end

%coloca na ordem crescente a partir do angulo (coluna 5)
matrix_uniq_angle_small = sortrows(matrix_uniq_small,5);
figure
hold on
cont_small=0;
scatter(dados_small{1},dados_small{2})
for k=1:n_vet_unic_small
    for j=1:n_vetores_longos_small
        if (matrix_uniq_angle_small(k,4)==max_small(j))
            scatter(matrix_uniq_angle_small(k,1),matrix_uniq_angle_small(k,2),'filled','red')
            cont_small=cont_small+1;
            text(matrix_uniq_angle_small(k,1),matrix_uniq_angle_small(k,2), int2str(cont_small))
            alphas_externos_small(cont_small) = matrix_uniq_angle_small(k,1);
            betas_externos_small(cont_small) = matrix_uniq_angle_small(k,2);
        end
    end
end
hold off

centering_a = 1;
centering_b = 2;
figure
hold on
h = fill(alphas_externos_small,betas_externos_small,'r');
set(h,'facealpha',.5)
h = fill(centering_a+alphas_externos_small,centering_b+betas_externos_small,'r');
set(h,'facealpha',.5)
hold off


%% sistema completo


