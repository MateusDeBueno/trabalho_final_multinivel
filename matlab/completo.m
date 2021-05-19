clear all
close all
clc

N = 7; %numero de niveis
Vdc = 1;

% num_coluna = N-1; %numero de colunas de chaves
% estados_chaves= dec2bin(0:(2^num_coluna-1)) - '0'; %Binary matrix table
% tensoes_de_linha = zeros(1,2^num_coluna); %tensao gerada por cada combinacao do conversor
% for k=1:2^num_coluna %k coluna
%     for j=1:num_coluna %j linha
%         tensoes_de_linha(k)= tensoes_de_linha(k) + estados_chaves(k,j); %a tensao gerada eh a soma de chaves acionadas
%     end
% end
% 
% %adiciona o nome correto para cada coluna da tabela
% for k=1:num_coluna
%     column_name{k} = ['Sp' int2str(k)];
% end
% column_name{num_coluna+1} = 'Tensao';
% table = array2table([estados_chaves, tensoes_de_linha'], 'VariableNames', column_name); %crio uma tabela com as combinacoes nas chaves Sp e a tensao de saida

tensoes_de_linha = [0 1 1 2 2 3 3 4 4 5 5 6];

V = sort(tensoes_de_linha); %todas as combinacoes de chaves geram as seguintes tensoes

n = length(V); %opcoes de chaveamento intra fase

%cria os vetores das tensoes de fase
Va = zeros(1,n);
Vb = zeros(1,n);
Vc = zeros(1,n);
for j = 1:n %atribui os valores para as tensoes de fase
    Va(j) = V(j);
    Vb(j) = V(j);
    Vc(j) = V(j);
end

n_vetores = 0; %conta o numero de vetores totais, deve ser igual a nˆ3
vetor_alpha = zeros(1,n^3,1); %cria vetor para alphas
vetor_beta = zeros(1,n^3,1); %cria vetor para betas
vetor = string(1:n^3); %vetor de string para salvar as tensoes de linha que levam a cada vetor

%preenche o vetor alpha e o vetor beta com seus respectivos valores
for k = 1:n %varre a tens o Va
    for j = 1:n %varre a tens o Vb
        for i = 1:n %varre a tens o Vc
        n_vetores = n_vetores+1; %quantos vetores existem
        [vetor_alpha(n_vetores), vetor_beta(n_vetores)] = transformada_clarke(Va(k), Vb(j), Vc(i)); %calcula alpha e beta

        %salva o char somente para saber quais tensoes estavam
        %aplicadas
        char_Va = int2str(Va(k)/Vdc);
        char_Vb = int2str(Vb(j)/Vdc);
        char_Vc = int2str(Vc(i)/Vdc);
        vetor(n_vetores) = append(char_Va, char_Vb, char_Vc); %salvando qual conjuento de tensoes de fase gera qual vetor
        end
    end
end

vetor_alpha = round(vetor_alpha,10);
vetor_beta = round(vetor_beta,10);

num_redundancias = zeros(1,n_vetores);
matrix_vector = [vetor_alpha', vetor_beta', num_redundancias']; %primeira coluna alpha, segunda coluna beta
matrix_uniq = unique(matrix_vector, 'rows');
dados = num2cell(matrix_uniq,1);

n_vet_unic = 3*N*(N-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix
%%


for z = 1:n_vetores
    alpha_talvez_unico = vetor_alpha(z);
    beta_talvez_unico = vetor_beta(z);
    for j = 1:n_vet_unic
        if (alpha_talvez_unico == dados{1}(j) && beta_talvez_unico == dados{2}(j))
            dados{3}(j) = dados{3}(j)+1; %contando o numero de redundancias por vetor
            dados{4}(j,dados{3}(j)) = vetor(z);
        end
    end
end


figure
gscatter(dados{1},dados{2},dados{3})
stg = blanks(1);
for z = 1:n_vet_unic
    numero_redundancias = int2str(dados{3}(z,1));
    text(dados{1}(z),dados{2}(z),numero_redundancias,'FontSize',11)
end
grid on
%title("Diagrama fase dos vetores espaciais para inversor de " + N + " n veis")
xlabel('Alpha')
ylabel('Beta')
%save_figure("intra_fase" + N)
