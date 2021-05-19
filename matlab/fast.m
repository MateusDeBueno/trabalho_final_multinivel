clear all
close all
clc

N = 7; %numero de niveis
Vdc = 1;
V = (0:1:N-1)*Vdc - 3; % opcoes de tensoes de fase

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
matrix_vector = [vetor_g', vetor_h', num_redundancias', Vab', Vbc', Vca']; %primeira coluna alpha, segunda coluna beta
matrix_uniq = unique(matrix_vector, 'rows');
dados = num2cell(matrix_uniq,1);

n_vet_unic = 3*N*(N-1)+1; %numero de vetores unicos, eh o mesmo comprimendo da matrix


for z = 1:n_vetores
    g_talvez_unico = vetor_g(z);
    h_talvez_unico = vetor_h(z);
    for j = 1:n_vet_unic
        if (g_talvez_unico == dados{1}(j) && h_talvez_unico == dados{2}(j))
            dados{3}(j) = dados{3}(j)+1; %contando o numero de redundancias por vetor
            dados{4}(j,dados{3}(j)) = vetor(z);
        end
    end
end

dados{5} = Vab;
dados{6} = Vbc;
dados{7} = Vca;

figure
scatter3(dados{5},dados{6},dados{7})
stg = blanks(1);
for z = 1:n_vetores
    numero_redundancias = int2str(dados{3}(z,1));
    text(dados{5}(z),dados{6}(z), dados{7}(z), numero_redundancias,'FontSize',11)
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
        stg = append(stg, newline,dados{4}(z,i));
    end
    text(dados{1}(z),dados{2}(z),stg,'FontSize',7)
    stg = erase(stg,stg);
end
grid on



%% VALIDACAO NO TEMPO


f_rede = 50; %freq da rede
m = 0.5; %indice de modulacao
n_pontos = 1000;
time = linspace(0,1/f_rede,n_pontos);
VA = zeros(1,n_pontos);
VB = zeros(1,n_pontos);
VC = zeros(1,n_pontos);


for k=1:length(time)
    VA(k) = 3*sin(time(k)*2*pi*f_rede);
    VB(k) = 3*sin(time(k)*2*pi*f_rede-2*pi*f_rede/3);
    VC(k) = 3*sin(time(k)*2*pi*f_rede-4*pi*f_rede/3);
    
    [g, h]= transformada_fast(VA, VB, VC);
    
    P1_g = ceil(g);
    P1_h = floor(h);
    P2_g = floor(g);
    P2_h = ceil(h);
    P3_g = floor(g);
    P3_h = floor(h);
    P4_g = ceil(g);
    P4_h = ceil(h);
end

figure
plot(VA)
hold on
plot(VB)
plot(VC)
hold off
xlim([0 n_pontos])
legend('Va', 'Vb', 'Vc')
title('Tensoes de linha')

figure
plot(g)
hold on
plot(h)
hold off
xlim([0 n_pontos])
legend('g', 'h')

figure
plot(P1_g)
hold on
plot(P1_h)
plot(P2_g)
plot(P2_h)
plot(P3_g)
plot(P3_h)
plot(P4_g)
plot(P4_h)
hold off
xlim([0 n_pontos])
legend('P1_g', 'P1_h', 'P2_g', 'P2_h', 'P3_g', 'P3_h', 'P4_g', 'P4_h')

figure
scatter(dados{1},dados{2})
hold on
plot(g,h)
hold off
grid on