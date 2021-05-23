clear all;
close all;
clc;

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


n_vetores_big = 0; %conta o numero de vetores big
% vetor_alpha_big = zeros(1,Nbig^3,1); %cria vetor para alphas
% vetor_beta_big = zeros(1,Nbig^3,1); %cria vetor para betas
% vetor_gama_big = zeros(1,Nbig^3,1); %cria vetor para gama
% vetor_big = string(1:Nbig^3);

n_vetores=0;
n = Nbig^3*Nsmall^3;

vetor_alpha = zeros(1,n,1); %cria vetor para alphas
vetor_beta = zeros(1,n,1); %cria vetor para betas
vetor_gama = zeros(1,n,1); %cria vetor para betas
vetor = string(1:n); %vetor de string para salvar as tensoes de linha que levam a cada vetor
Va = zeros(1,n,1); %cria vetor para Va
Vb = zeros(1,n,1); %cria vetor para Va
Vc = zeros(1,n,1); %cria vetor para Va
conjunto = zeros(1,n,1); %cria vetor para Va

for a = 1:Nbig %varre a tensao Va_big
    for b = 1:Nbig %varre a tensao Vb_big
        for c = 1:Nbig %varre a tensao Vc_big
            n_vetores_big = n_vetores_big+1; %quantos vetores big existem
            for d = 1:Nsmall %varre a tensao Va_small
                for e = 1:Nsmall %varre a tensao Vb_small
                    for f = 1:Nsmall %varre a tensao Vc_small
                        n_vetores = n_vetores+1; %quantos vetores existem
                        Va(n_vetores) = Va_big(a) + Va_small(d);
                        Vb(n_vetores) = Vb_big(b) + Vb_small(e);
                        Vc(n_vetores) = Vc_big(c) + Vc_small(f);
                        conjunto(n_vetores) = n_vetores_big;
                        [vetor_alpha(n_vetores), vetor_beta(n_vetores), vetor_gama(n_vetores)] = transformada_clarke(Va(n_vetores), Vb(n_vetores), Vc(n_vetores));
                    end
                end
            end
        end
    end
end

%aproximar para evitar erros
vetor_alpha = round(vetor_alpha,10);
vetor_beta = round(vetor_beta,10);
vetor_gama = round(vetor_gama,10);

num_redundancias = zeros(1,n_vetores);

matrix_vector = [vetor_alpha', vetor_beta', vetor_gama', num_redundancias', conjunto'];
%coluna 1 = alpha
%coluna 2 = beta
%coluna 3 = gama
%coluna 4 = num redundancias
%coluna 5 = conjunto (estado do vetor big) //isso so esta na matrix_vector
matrix_uniq = unique(matrix_vector, 'rows');

%%
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
            %dados{4}(j,dados{3}(j)) = vetor(z);
            cnt=cnt+1;
        end
    end
end

% %conta o numero de estados da celula grande para cada ponto
% for z = 1:n_vetores
%     alpha_talvez_unico = vetor_alpha(z);
%     beta_talvez_unico = vetor_beta(z);
%     conjunto_talvez = conjunto(z);
%     for j = 1:n_vet_unic
%         if (alpha_talvez_unico == dados{1}(j) && beta_talvez_unico == dados{2}(j) && conjunto_talvez~=dados{5}(j))
%             dados{6}(j) = dados{6}(j)+1; %contando o numero de redundancias por vetor
%         end
%     end
% end


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
% 
figure
scatter3(dados{1},dados{2},dados{5})



figure
scatter(dados{1},dados{2},(dados{5}+1)*10)

% figure
% area(dados{1}, dados{2})




% figure
% subplot(3,1,1);
% plot(Va)
% subplot(3,1,2); 
% plot(Vb)
% subplot(3,1,3); 
% plot(Vc)
