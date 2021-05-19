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
vetor_g = zeros(1,N^3,1); %cria vetor para alphas
vetor_h = zeros(1,N^3,1); %cria vetor para betas
vetor = string(1:N^3); %vetor de string para salvar as tensoes de linha que geram a cada vetor

%preenche o vetor alpha e o vetor beta com seus respectivos valores
for k = 1:N %varre a tens o Va
    for j = 1:N %varre a tens o Vb
        for i = 1:N %varre a tens o Vc
            n_vetores = n_vetores+1; %quantos vetores existem
            [vetor_g(n_vetores), vetor_h(n_vetores)] = transformada_fast(Va(k), Vb(j), Vc(i));
            char_Va = int2str(Va(k)/Vdc);
            char_Vb = int2str(Vb(j)/Vdc);
            char_Vc = int2str(Vc(i)/Vdc);
            Vab(n_vetores) = Va(k) - Vb(j);
            Vbc(n_vetores) = Vb(j) - Vc(i);
            Vca(n_vetores) = Vc(i) - Va(k);
            vetor(n_vetores) = append(char_Va, char_Vb, char_Vc); %salvando qual conjuento de tensoes de fase gera qual vetor
        end
    end
end




num_redundancias = zeros(1,n_vetores);
posicoes = zeros(1,n_vetores);
matrix_vector = [vetor_g', vetor_h', num_redundancias', Vab', Vbc', Vca'];
%coluna 1 = g
%coluna 2 = g
%coluna 3 = numero de redundancias
%coluna 4 = tensao Vab
%coluna 5 = tensao Vbc
%coluna 6 = tensao Vca
%coluna 7 = estados do conversor
matrix_uniq = unique(matrix_vector, 'rows');
dados = num2cell(matrix_uniq,1);

n_vet_unic = 3*N*(N-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix


for z = 1:n_vetores
    g_talvez_unico = vetor_g(z);
    h_talvez_unico = vetor_h(z);
    for j = 1:n_vet_unic
        if (g_talvez_unico == dados{1}(j) && h_talvez_unico == dados{2}(j))
            dados{3}(j) = dados{3}(j)+1; %contando o numero de redundancias por vetor
            dados{7}(j,dados{3}(j)) = vetor(z);
        end
    end
end


figure
scatter3(dados{4},dados{5},dados{6})
stg = blanks(1);
for z = 1:n_vet_unic
    numero_redundancias = int2str(dados{3}(z,1));
    text(dados{4}(z),dados{5}(z), dados{6}(z), numero_redundancias,'FontSize',11)
end
xlabel('Vab')
ylabel('Vbc')
zlabel('Vca')





figure
gscatter(dados{1},dados{2},dados{3})
stg = blanks(1);
for z = 1:n_vet_unic
    numero_redundancias = int2str(dados{3}(z,1));
    text(dados{1}(z),dados{2}(z),numero_redundancias,'FontSize',11)
end
grid on
%title("Diagrama fase dos vetores espaciais para inversor de " + N + " n veis")
xlabel('g')
ylabel('h')
%save_figure("intra_fase" + N)





%plot do numero de redundancias
figure
scatter(dados{1},dados{2},25,dados{3},'filled')
stg = blanks(1);
for z = 1:n_vet_unic
    for i=1:dados{3}(z,1)
        stg = append(stg, newline,dados{7}(z,i));
    end
    text(dados{1}(z),dados{2}(z),stg,'FontSize',7)
    stg = erase(stg,stg);
end
grid on



%% VALIDACAO NO TEMPO


f_rede = 50; %freq da rede
m = 3.47; %indice de modulacao
n_pontos = 1000;
time = linspace(0,1/f_rede,n_pontos); %vetor de tempo


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
delta_ul = zeros(1,n_pontos);
delta_lu = zeros(1,n_pontos);
delta_ll = zeros(1,n_pontos);
delta_uu = zeros(1,n_pontos);

for k=1:length(time)
    VA(k) = m*sin(time(k)*2*pi*f_rede);
    VB(k) = m*sin(time(k)*2*pi*f_rede-2*pi*f_rede/3);
    VC(k) = m*sin(time(k)*2*pi*f_rede-4*pi*f_rede/3);
    
    [g_ref(k), h_ref(k)]= transformada_fast(VA(k), VB(k), VC(k));
    
    Vul_g(k) = ceil(g_ref(k));
    Vul_h(k) = floor(h_ref(k));
    Vlu_g(k) = floor(g_ref(k));
    Vlu_h(k) = ceil(h_ref(k));
    Vll_g(k) = floor(g_ref(k));
    Vll_h(k) = floor(h_ref(k));
    Vuu_g(k) = ceil(g_ref(k));
    Vuu_h(k) = ceil(h_ref(k));

    
    %ïdentificar quais os 3 pontos usar, P1 e P2 sao sempre usados
    %pore, eh utilizado ou P3 ou P4
    V1_g(k) = Vul_g(k);
    V1_h(k) = Vul_h(k);
    V2_g(k) = Vlu_g(k);
    V2_h(k) = Vlu_h(k);    
    if ((g_ref(k)+h_ref(k))-(Vul_g(k)+Vul_h(k)))>0 %escolhe Vuu
        V3_g(k) = Vuu_g(k);
        V3_h(k) = Vuu_h(k);
        delta_ul(k)=g_ref(k)-Vll_g(k);
        delta_lu(k)=h_ref(k)-Vll_h(k);
        delta_ll(k)=1-delta_ul(k)-delta_lu(k);  
    else                                         %escolhe Vll
        V3_g(k) = Vll_g(k);
        V3_h(k) = Vll_h(k);
        delta_ul(k)=-(h_ref(k)-Vuu_h(k));
        delta_lu(k)=-(g_ref(k)-Vuu_g(k));
        delta_uu(k)=1-delta_ul(k)-delta_lu(k);
    end
end

% figure
% plot(VA)
% hold on
% plot(VB)
% plot(VC)
% hold off
% xlim([0 n_pontos])
% legend('Va', 'Vb', 'Vc')
% title('Tensoes de linha')

figure
plot(g_ref)
hold on
plot(h_ref)
hold off
xlim([0 n_pontos])
legend('g', 'h')

% figure
% plot(V1_g)
% hold on
% plot(V1_h)
% plot(V2_g)
% plot(V2_h)
% plot(V3_g)
% plot(V3_h)
% hold off
% xlim([0 n_pontos])
% legend('P1_g', 'P1_h', 'P2_g', 'P2_h', 'P3_g', 'P3_h')


figure
scatter(dados{1},dados{2})
hold on
plot(g_ref,h_ref)
hold off
grid on



figure
scatter(V1_g, V1_h, 50)
hold on
scatter(V2_g, V2_h, 30)
scatter(V3_g, V3_h, 10)
plot(g_ref,h_ref)
hold off
legend('V1', 'V2', 'V3', 'Vref')

