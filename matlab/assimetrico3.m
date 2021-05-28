clear all
close all
clc

%numero de celulas por fase
num_celula = 2;
nivel = 2^(num_celula+1)-1; %isso para V e 2V
%essa calculo so serve para verificar numero de niveis por fase

%% ANALISE DO GRANDE
%CRIANDO TENSOES DE CADA CELULA
%CELULA GRANDE
Vdc_big = 2;
Nbig = 3;
Vbig = (0:1:Nbig-1)*Vdc_big - ceil(Nbig/2);  % opcoes de tensoes de fase

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
title('Possiveis estados para a celula grande')


% ESSE PLOT EH DOS ANGULOS E COMPRIMENTOS DE CADA CELULA GRANDE
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
title('Angulo e comprimento de cada vetor da celula grande')


%pegar os valores longos de interesse
comprimeto_big = dados_big{4};
n_vetores_longos_big = ceil(Nbig/2);
max_big = zeros(1,n_vetores_longos_big);
for k=1:n_vetores_longos_big
    max_big(k) = max(comprimeto_big);
    comprimeto_big = comprimeto_big(comprimeto_big~=max_big(k));
end

%coloca na ordem crescente a partir do angulo (coluna 5), dessa maneira eh
%possivel printar as regioes em torno de um ponto, as regioes coloridas
%possuem capacidade plena de PWM
matrix_uniq_angle_big = sortrows(matrix_uniq_big,5);
cont_big=0;
for k=1:n_vet_unic_big
    for j=1:n_vetores_longos_big
        if (matrix_uniq_angle_big(k,4)==max_big(j))
            cont_big=cont_big+1;
            alphas_externos_big(cont_big) = matrix_uniq_angle_big(k,1);
            betas_externos_big(cont_big) = matrix_uniq_angle_big(k,2);
        end
    end
end


%pegar os valores longos de interesse
comprimeto_big = dados_big{4};
n_vetores_longos_big = ceil(Nbig/2);
max_big = zeros(1,n_vetores_longos_big);
for k=1:n_vetores_longos_big
    max_big(k) = max(comprimeto_big);
    comprimeto_big = comprimeto_big(comprimeto_big~=max_big(k));
end

%coloca na ordem crescente a partir do angulo (coluna 5), dessa maneira eh
%possivel printar as regioes em torno de um ponto
matrix_uniq_angle_big = sortrows(matrix_uniq_big,5);
cont_big=0;
for k=1:n_vet_unic_big
    for j=1:n_vetores_longos_big
        if (matrix_uniq_angle_big(k,4)==max_big(j))
            cont_big=cont_big+1;
            alphas_externos_big(cont_big) = matrix_uniq_angle_big(k,1);
            betas_externos_big(cont_big) = matrix_uniq_angle_big(k,2);
        end
    end
end






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
title('Possiveis estados para uma celula pequena')

% % Angulo e comprimento de cada estado da celula pequena
% figure
% scatter(dados_small{1},dados_small{2})
% stg = blanks(1);
% for z = 1:n_vet_unic_small
%     angulo = int2str(dados_small{5}(z,1));
%     comprimento = num2str(round(dados_small{4}(z,1),2));
%     text(dados_small{1}(z),dados_small{2}(z),angulo,'FontSize',11)
%     text(dados_small{1}(z),dados_small{2}(z)-0.2,comprimento,'FontSize',11)
% end
% grid on
% xlabel('Alpha')
% ylabel('Beta')
% title('Angulo e comprimento de cada vetor da celula pequena')

%pegar os valores longos de interesse
comprimeto_small = dados_small{4};
n_vetores_longos_small = ceil(Nsmall/2);
max_small = zeros(1,n_vetores_longos_small);
for k=1:n_vetores_longos_small
    max_small(k) = max(comprimeto_small);
    comprimeto_small = comprimeto_small(comprimeto_small~=max_small(k));
end

%coloca na ordem crescente a partir do angulo (coluna 5), dessa maneira eh
%possivel printar as regioes em torno de um ponto, as regioes coloridas
%possuem capacidade plena de PWM
matrix_uniq_angle_small = sortrows(matrix_uniq_small,5);
cont_small=0;
for k=1:n_vet_unic_small
    for j=1:n_vetores_longos_small
        if (matrix_uniq_angle_small(k,4)==max_small(j))
            cont_small=cont_small+1;
            alphas_externos_small(cont_small) = matrix_uniq_angle_small(k,1);
            betas_externos_small(cont_small) = matrix_uniq_angle_small(k,2);
        end
    end
end

%% SISTEMA COMPLETO

%variavel para contar o numero de vetores grandes
n_vetores_Big = 0;

%variavel para contar o numero de vetores totais (grande e pequenos)
n_vetores=0;

%numero de vetores totais (grande e pequenos)
n = Nbig^3*Nsmall^3;

%CRIANDO VARIAS VARIAVEIS PARA SEREM USADAS NO LOOP QUE VARRE TENSOES
vetor_alpha = zeros(1,n,1); %cria vetor para alphas
vetor_beta = zeros(1,n,1); %cria vetor para betas
vetor_gama = zeros(1,n,1); %cria vetor para betas
Va = zeros(1,n,1); %cria vetor para Va
Vb = zeros(1,n,1); %cria vetor para Vb
Vc = zeros(1,n,1); %cria vetor para Vc
conjunto = zeros(1,n,1); %cria vetor para contar os estados da celula grande
vetor_vabc = string(1:n); %vetor de string para salvar as tensoes de fase que geram a cada vetor
vetor_vbig = string(1:n); %vetor de string para salvar as tensoes da celula grande
vetor_vsmall = string(1:n); %vetor de string para salvar as tensoes da celula pequena
vetor_g = zeros(1,n,1); %cria vetor para g
vetor_h = zeros(1,n,1); %cria vetor para h
comp = zeros(1,n,1); %cria vetor para salvar o comprimento de cada vetor
ang = zeros(1,n,1); %cria vetor para salvar o angulo de cada vetor

%VARRENDO AS TENSOES
for a = 1:Nbig %varre a tensao Va_big
    for b = 1:Nbig %varre a tensao Vb_big
        for c = 1:Nbig %varre a tensao Vc_big
            n_vetores_Big = n_vetores_Big+1; %quantos vetores big existem
            for d = 1:Nsmall %varre a tensao Va_small
                for e = 1:Nsmall %varre a tensao Vb_small
                    for f = 1:Nsmall %varre a tensao Vc_small
                        n_vetores = n_vetores+1; %quantos vetores existem
                        
                        %somandos as tensoes do big com as tensoes do small
                        Va(n_vetores) = Va_big(a) + Va_small(d);
                        Vb(n_vetores) = Vb_big(b) + Vb_small(e);
                        Vc(n_vetores) = Vc_big(c) + Vc_small(f);
                        
                        %conjunto, eh difinido pela posicao do vetor big
                        conjunto(n_vetores) = n_vetores_big;
                        
                        %calculo do vetor alpha beta
                        [vetor_alpha(n_vetores), vetor_beta(n_vetores), vetor_gama(n_vetores)] = transformada_clarke(Va(n_vetores), Vb(n_vetores), Vc(n_vetores));
                        
                        %calcula angulo e comprimento
                        comp(n_vetores) = ((vetor_alpha(n_vetores))^2 + (vetor_beta(n_vetores))^2)^0.5;
                        ang(n_vetores) = get_angle(vetor_alpha(n_vetores),vetor_beta(n_vetores));
                        
                        %salvar posicao
                        char_Va = int2str(Va(n_vetores));
                        char_Vb = int2str(Vb(n_vetores));
                        char_Vc = int2str(Vc(n_vetores));
                        vetor_vabc(n_vetores) = append(char_Va, char_Vb, char_Vc); %salvando qual conjuento de tensoes de fase gera qual vetor
                        
                        %salva posicao Vbig
                        char_Va_big = int2str(Va_big(a));
                        char_Vb_big = int2str(Vb_big(b));
                        char_Vc_big = int2str(Vc_big(c));
                        vetor_vbig(n_vetores) = append(char_Va_big, char_Vb_big, char_Vc_big); %salvando conjuento de tensoes da celula grande

                        %salva posicao Vsmall
                        char_Va_small = int2str(Va_small(d));
                        char_Vb_small = int2str(Vb_small(e));
                        char_Vc_small = int2str(Vc_small(f));
                        vetor_vsmall(n_vetores) = append(char_Va_small, char_Vb_small, char_Vc_small); %salvando conjuento de tensoes da celula pequena
                        
                        %transformada fast
                        [vetor_g(n_vetores), vetor_h(n_vetores)] = transformada_fast(Va(n_vetores), Vb(n_vetores), Vc(n_vetores));
                    end
                end
            end
        end
    end
end

%aproximar para evitar erros, 10 casas decimais
vetor_alpha = round(vetor_alpha,10);
vetor_beta = round(vetor_beta,10);
vetor_gama = round(vetor_gama,10);
comp = round(comp,10);
ang = round(ang,10);

%vetor para numero de redundancias no plano alpha beta
num_redundancias = zeros(1,n_vetores);

matrix_vetor = [vetor_alpha', vetor_beta', vetor_gama', num_redundancias', comp', ang'];
%MATRIX QUE SERA REDUZIDA
%coluna 1 = alpha
%coluna 2 = beta
%coluna 3 = gama
%coluna 4 = num redundancias
%coluna 5 = comprimento
%coluna 6 = angulo
%coluna 7 = tensao Vbig+Vsmall
%coluna 8 = tensao Vbig
%coluna 9 = tensao Vsmall
%coluna 10 = g
%coluna 11 = h

matrix_uniq = unique(matrix_vetor, 'rows');
dados = num2cell(matrix_uniq,1);
n_vet_unic = length(matrix_uniq);

cnt = 0;
%conta o numero de redundancias
for z = 1:n_vetores
    alpha_talvez_unico = vetor_alpha(z);
    beta_talvez_unico = vetor_beta(z);
    for j = 1:n_vet_unic
        if (alpha_talvez_unico == dados{1}(j) && beta_talvez_unico == dados{2}(j))
            dados{4}(j) = dados{4}(j)+1; %contando o numero de redundancias por vetor
            dados{7}(j,dados{4}(j)) = vetor_vabc(z);
            dados{8}(j,dados{4}(j)) = vetor_vbig(z);
            dados{9}(j,dados{4}(j)) = vetor_vsmall(z);
            dados{10}(j,dados{4}(j)) = vetor_g(z);
            dados{11}(j,dados{4}(j)) = vetor_h(z);
            cnt=cnt+1;
        end
    end
end


figure
gscatter(dados{1},dados{2},dados{4})
stg = blanks(1);
for z = 1:n_vet_unic
    numero_redundancias = int2str(dados{4}(z,1));
    text(dados{1}(z),dados{2}(z),numero_redundancias,'FontSize',11)
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Numero de redundancias para cada posição no mapa Alpha Beta')
% save_figure('numero_redundancias')


% % PLOT DO ANGULO DE CADA VETOR DO SISTEMA COMPLETO
% figure
% scatter(dados{1},dados{2})
% stg = blanks(1);
% for z = 1:n_vet_unic
%     angulo = int2str(dados{6}(z,1));
%     comprimento = num2str(round(dados{5}(z,1),2));
%     text(dados{1}(z),dados{2}(z),angulo,'FontSize',11)
% %     text(dados{1}(z),dados{2}(z)-0.2,comprimento,'FontSize',11)
% end
% grid on
% xlabel('Alpha')
% ylabel('Beta')
% title('Angulo de cada vetor')



%%
%A partir dos dados tirados ate entao nota-se que o sistema pode fazer PWM
%durante todo o sistema
figure
angulo_inicial = 15;
hold on
for k=1:n_vet_unic_big
    h = fill(dados_big{1}(k)+alphas_externos_small,dados_big{2}(k)+betas_externos_small,'r');
    set(h,'facealpha',.4)
end
stg = blanks(1);
for z = 1:n_vet_unic_big
    for i=1:dados_big{3}(z,1)
        stg = append(stg, newline,dados_big{6}(z,i));
    end
    text(dados_big{1}(z),dados_big{2}(z),stg,'FontSize',12)
    stg = erase(stg,stg);
end

for k=-180:30:180
    [t1, t2]= reta_angulo(0,0,4,k+angulo_inicial);
    plot(t1, t2, 'b')
end
r = 3.4;
x=0;
y=0;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit, 'g');
scatter(dados_big{1},dados_big{2}, 50, 'filled', 'yellow')
scatter(dados{1},dados{2}, 15,'filled','black')
hold off
xlabel('Alpha')
ylabel('Beta')


figure
scatter(dados{1},dados{2},25,dados{4},'filled')
stg = blanks(1);
k=0;
for z = 1:n_vet_unic
    k = k+1;
    if (dados{5}(k)>=2.8)
        for i=1:dados{4}(z,1)
            stg = append(stg, newline,dados{7}(z,i));
        end
        text(dados{1}(z),dados{2}(z),stg,'FontSize',8)
        stg = erase(stg,stg);
    end
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Tensao total aplicada em cada fase')


figure
hold on
for k=-180:30:180
    [t1, t2]= reta_angulo(0,0,4,k+angulo_inicial);
    plot(t1, t2, 'b')
end
scatter(dados{1},dados{2},'filled')
stg = blanks(1);
k=0;
for z = 1:n_vet_unic
    k = k+1;
    if (dados{5}(k)>=2.8)
        for i=1:dados{4}(z,1)
            stg = append(stg, newline,dados{8}(z,i));
        end
        text(dados{1}(z),dados{2}(z),stg,'FontSize',8)
        stg = erase(stg,stg);
    end
end
scatter(alphas_externos_big,betas_externos_big, 'filled')
hold off
grid on
xlabel('Alpha')
ylabel('Beta')
title('Tensao aplicada em cada fase pela celula grande')


%% VALIDACAO NO TEMPO

f_rede = 60; %freq do motor
Vref = 3.4; %a tensao de referencia maxima eh 3.4607
fs = 720*2; %frequencia de amostragem (7 segmentos)
n_pontos_t = fs/f_rede;
n_pontos = n_pontos_t+1;
time = linspace(1/(2*fs),(1/f_rede+1/(2*fs)),n_pontos); %vetor de tempo


%criando diversos vetores para preencher
VA = zeros(1,n_pontos);
VB = zeros(1,n_pontos);
VC = zeros(1,n_pontos);
g_ref = zeros(1,n_pontos);
h_ref = zeros(1,n_pontos);
Vul_g = zeros(1,n_pontos);
Vul_h = zeros(1,n_pontos);
Vlu_g = zeros(1,n_pontos);
Vlu_h = zeros(1,n_pontos);
Vll_g = zeros(1,n_pontos);
Vll_h = zeros(1,n_pontos);
Vuu_g = zeros(1,n_pontos);
Vuu_h = zeros(1,n_pontos);
V1_g = zeros(1,n_pontos);
V1_h = zeros(1,n_pontos);
V2_g = zeros(1,n_pontos);
V2_h = zeros(1,n_pontos);
V3_g = zeros(1,n_pontos);
V3_h = zeros(1,n_pontos);
delta_V1 = zeros(1,n_pontos);
delta_V2 = zeros(1,n_pontos);
delta_V3 = zeros(1,n_pontos);
cont_1 = zeros(1,n_pontos);
cont_2 = zeros(1,n_pontos);
cont_3 = zeros(1,n_pontos);




%amostrando a cada 1/fs
for k=1:length(time)
    VA(k) = Vref*sin(time(k)*2*pi*f_rede);
    VB(k) = Vref*sin(time(k)*2*pi*f_rede+2*pi/3);
    VC(k) = Vref*sin(time(k)*2*pi*f_rede+4*pi/3);
    
    [g_ref(k), h_ref(k)] = transformada_fast(VA(k), VB(k), VC(k));
    
    [Vul_g(k),Vul_h(k),Vlu_g(k),Vlu_h(k),Vll_g(k),Vll_h(k),Vuu_g(k),Vuu_h(k),V1_g(k),V1_h(k),V2_g(k),V2_h(k),V3_g(k),V3_h(k),delta_V1(k),delta_V2(k),delta_V3(k)] = aproximacoes_fast(g_ref(k),h_ref(k),nivel);
    
    for p=1:nivel %zero ate o nuumero de niveis
        pretenso_Va_1 = p-1;
        pretenso_Vb_1 = p-1-V1_g(k);
        pretenso_Vc_1 = p-1-V1_g(k)-V1_h(k);
        if (pretenso_Va_1>=0 && pretenso_Va_1<7 && pretenso_Vb_1>=0 && pretenso_Vb_1<7 && pretenso_Vc_1>=0 && pretenso_Vc_1<7)
            cont_1(k) = cont_1(k)+1;
            Va_1(k,cont_1(k))=pretenso_Va_1;
            Vb_1(k,cont_1(k))=pretenso_Vb_1;
            Vc_1(k,cont_1(k))=pretenso_Vc_1;
        end
        pretenso_Va_2 = p-1;
        pretenso_Vb_2 = p-1-V2_g(k);
        pretenso_Vc_2 = p-1-V2_g(k)-V2_h(k);  
        if (pretenso_Va_2>=0 && pretenso_Va_2<7 && pretenso_Vb_2>=0 && pretenso_Vb_2<7 && pretenso_Vc_2>=0 && pretenso_Vc_2<7)
            cont_2(k) = cont_2(k)+1;
            Va_2(k,cont_2(k))=pretenso_Va_2;
            Vb_2(k,cont_2(k))=pretenso_Vb_2;
            Vc_2(k,cont_2(k))=pretenso_Vc_2;
        end
        pretenso_Va_3 = p-1;
        pretenso_Vb_3 = p-1-V3_g(k);
        pretenso_Vc_3 = p-1-V3_g(k)-V3_h(k);     
        if (pretenso_Va_3>=0 && pretenso_Va_3<7 && pretenso_Vb_3>=0 && pretenso_Vb_3<7 && pretenso_Vc_3>=0 && pretenso_Vc_3<7)
            cont_3(k) = cont_3(k)+1;
            Va_3(k,cont_3(k))=pretenso_Va_3;
            Vb_3(k,cont_3(k))=pretenso_Vb_3;
            Vc_3(k,cont_3(k))=pretenso_Vc_3;
        end
    end
    
    %pegando os alphas, betas e gamas de cada ponto a partir do gh
    [vetor_alpha_ref(k), vetor_beta_ref(k), vetor_gama_ref(k)] = transformada_clarke(VA(k), VB(k), VC(k));
    [vetor_alpha1(k), vetor_beta1(k), vetor_gama1(k)] = transformada_clarke(Va_1(k,1), Vb_1(k,1), Vc_1(k,1));
    [vetor_alpha2(k), vetor_beta2(k), vetor_gama2(k)] = transformada_clarke(Va_2(k,1), Vb_2(k,1), Vc_2(k,1));
    [vetor_alpha3(k), vetor_beta3(k), vetor_gama3(k)] = transformada_clarke(Va_3(k,1), Vb_3(k,1), Vc_3(k,1));
end

%grafico dos pontos utilizados no plano alpha beta
figure
hold on
scatter(dados{1},dados{2})
scatter(vetor_alpha_ref, vetor_beta_ref, 'filled', 'r')
plot(vetor_alpha_ref, vetor_beta_ref)
scatter(vetor_alpha1, vetor_beta1, 'filled', 'b')
scatter(vetor_alpha2, vetor_beta2, 'filled', 'b')
scatter(vetor_alpha3, vetor_beta3, 'filled', 'b')
% for k=1:length(g_ref)-1 %escreve as razoes ciclicas
%     text(vetor_alpha1(k), vetor_beta1(k),num2str(round(delta_V1(k),2)),'FontSize',8)
%     text(vetor_alpha2(k), vetor_beta2(k),num2str(round(delta_V2(k),2)),'FontSize',8)
%     text(vetor_alpha3(k), vetor_beta3(k),num2str(round(delta_V2(k),2)),'FontSize',8)
% end
% matriz_repeticao = [];
% vet_alpha_plot = unique([vetor_alpha1, vetor_alpha2, vetor_alpha3]);
% vet_beta_plot = unique([vetor_beta1, vetor_beta2, vetor_beta3]);
for k=1:length(vetor_alpha_ref)-1 %enumera os pontos de amostragem
    text(vetor_alpha_ref(k), vetor_beta_ref(k), int2str(k))
    text(vetor_alpha1(k)-.05, vetor_beta1(k)-.05,int2str(k),'FontSize',8)
    text(vetor_alpha2(k)+.05, vetor_beta2(k)+.05,int2str(k),'FontSize',8)
    text(vetor_alpha3(k), vetor_beta3(k),int2str(k),'FontSize',8)
end
xlabel('Alpha')
ylabel('Beta')
for k=-180:30:180
    [t1, t2]= reta_angulo(0,0,4,k+angulo_inicial);
    plot(t1, t2, 'b')
end
hold off

%grafico dos pontos utilizados no plano G H
figure
hold on
scatter(dados{10}(:,1),dados{11}(:,1))
scatter(g_ref, h_ref,'filled','red')
plot(g_ref, h_ref)
scatter(V1_g, V1_h,'filled','blue')
scatter(V2_g, V2_h,'filled','blue')
scatter(V3_g, V3_h,'filled','blue')
% for k=1:length(g_ref)-1 %escreve as razoes ciclicas
%     text(V1_g(k), V1_h(k),num2str(round(delta_V1(k),2)),'FontSize',8)
%     text(V2_g(k), V2_h(k),num2str(round(delta_V2(k),2)),'FontSize',8)
%     text(V3_g(k), V3_h(k),num2str(round(delta_V3(k),2)),'FontSize',8)
% end
for k=1:length(g_ref)-1 %enumera os pontos de amostragem
    text(g_ref(k), h_ref(k),int2str(k))
    text(V1_g(k), V1_h(k),int2str(k),'FontSize',8)
    text(V2_g(k)-.05, V2_h(k)-.05,int2str(k),'FontSize',8)
    text(V3_g(k)+.05, V3_h(k)+.05,int2str(k),'FontSize',8)
end
hold off
xlabel('G')
ylabel('H')


figure
b = bar([cont_1(1:end-1)', cont_2(1:end-1)', cont_3(1:end-1)']);
% b = bar([cont_1(1:end-1)']);
legend('V1', 'V2', 'V3')
title('Opcoes de vetores para cada ponto de amostragem')




%%
%precisao de cada intervalo de amostragem
precision = 10000;
mid = round(precision/2);
Va_emula_total = [];
Vb_emula_total = [];
Vc_emula_total = [];

%criando vetores para analise em um periodo da rede
segmento = zeros(1, precision);
Va_emula = zeros(1, precision);
Vb_emula = zeros(1, precision);
Vc_emula = zeros(1, precision);

%variando em um periodo da rede, cada periodo de amostragem possue precisao
%igual a variavel precision
for j=1:length(time)-1
    %converte as razoes ciclicas
    del_V1 = delta_V1(j)*precision;
    del_V2 = delta_V2(j)*precision;
    del_V3 = delta_V3(j)*precision;

    %pegando as primeiras opcoes de tensao de todas as fases, essa etapa q deve
    %ser repensada

    %tensoes para o vetor V1
    Va1 = Va_1(j,1);
    Vb1 = Vb_1(j,1);
    Vc1 = Vc_1(j,1);

    %tensoes para o vetor V2
    Va2 = Va_2(j,1);
    Vb2 = Vb_2(j,1);
    Vc2 = Vc_2(j,1);

    %tensoes para o vetor V2
    Va3 = Va_3(j,1);
    Vb3 = Vb_3(j,1);
    Vc3 = Vc_3(j,1);

    %7 segmentos
    %ordem V1 V2 V3 V1 V3 V2 V1
    del_seg1 = del_V1/4;
    del_seg2 = del_V2/2;
    del_seg3 = del_V3/2;
    del_seg4 = del_V1/2;
    del_seg5 = del_V3/2;
    del_seg6 = del_V2/2;
    del_seg7 = del_V1/4;

    for k=1:precision
        if (del_seg1>=k)
            segmento(k) = 1;
        elseif ((del_seg1+del_seg2)>=k && k<=mid)
            segmento(k) = 2;
        elseif ((del_seg1+del_seg2+del_seg3)>=k && k<=mid)
            segmento(k) = 3;
        elseif ((del_seg1+del_seg2+del_seg3)<=k && (del_seg1+del_seg2+del_seg3+del_seg4)>=k)
            segmento(k) = 4;
        elseif ((del_seg1+del_seg2+del_seg3+del_seg4+del_seg5)>=k && k>=mid)
            segmento(k) = 5;
        elseif ((del_seg1+del_seg2+del_seg3+del_seg4+del_seg5+del_seg6)>=k && k>=mid)
            segmento(k) = 6;
        else
            segmento(k) = 7;
        end
        %7 segmentos
        %ordem V1 V2 V3 V1 V3 V2 V1
        if (segmento(k)==1 || segmento(k) ==4 || segmento(k)==7)
            Va_emula(k) = Va1;
            Vb_emula(k) = Vb1;
            Vc_emula(k) = Vc1;
        elseif (segmento(k)==2 || segmento(k) ==6)
            Va_emula(k) = Va2;
            Vb_emula(k) = Vb2;
            Vc_emula(k) = Vc2;
        else
            Va_emula(k) = Va3;
            Vb_emula(k) = Vb3;
            Vc_emula(k) = Vc3;
        end
    end
    Va_emula_total = [Va_emula_total, Va_emula];
    Vb_emula_total = [Vb_emula_total, Vb_emula];
    Vc_emula_total = [Vc_emula_total, Vc_emula];
end

figure
plot([Va_emula_total, Va_emula_total]-3)
title('Dois periodos da tensao da fase A emulada')
xlim([1 length([Va_emula_total, Va_emula_total])])
% figure
% subplot(3,1,1);
% plot(Va_emula_total)
% xlim([1 length(Va_emula_total)])
% title('Va')
% subplot(3,1,2); 
% plot(Vb_emula_total)
% xlim([1 length(Va_emula_total)])
% title('Vb')
% subplot(3,1,3); 
% plot(Vc_emula_total)
% xlim([1 length(Va_emula_total)])
% title('Vc')
% 
% 
% figure
% subplot(3,1,1);
% plot(Va_emula_total-Vb_emula_total)
% xlim([1 length(Va_emula_total)])
% title('Vab')
% subplot(3,1,2); 
% plot(Vb_emula_total-Vc_emula_total)
% xlim([1 length(Va_emula_total)])
% title('Vbc')
% subplot(3,1,3); 
% plot(Vc_emula_total-Va_emula_total)
% xlim([1 length(Va_emula_total)])
% title('Vca')