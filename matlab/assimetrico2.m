clear all;
close all;
clc;

%numero de celulas por fase
num_celula = 2;
nivel = 2^(num_celula+1)-1;
%essa calculo so serve para verificar numero de niveis por fase

%CRIANDO TENSOES DE CADA CELULA
%CELULA GRANDE
Vdc_big = 2;
Nbig = 3;
Vbig = (0:1:Nbig-1)*Vdc_big-2;  % opcoes de tensoes de fase

%cria os vetores das tensoes de fase
Va_big = zeros(1,Nbig);
Vb_big = zeros(1,Nbig);
Vc_big = zeros(1,Nbig);
for j = 1:Nbig %atribui os valores para as tensoes de fase
    Va_big(j) = Vbig(j);
    Vb_big(j) = Vbig(j);
    Vc_big(j) = Vbig(j);
end

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

%variavel para contar o numero de vetores grandes
n_vetores_big = 0;

%variavel para contar o numero de vetores totais (grande e pequenos)
n_vetores=0;

%numero de vetores totais (grande e pequenos)
n = Nbig^3*Nsmall^3;

%CRIANDO VARIAS VARIAVEIS PARA SEREM USADAS NO LOOP QUE VARRE TENSOES
vetor_alpha = zeros(1,n,1); %cria vetor para alphas
vetor_beta = zeros(1,n,1); %cria vetor para betas
vetor_gama = zeros(1,n,1); %cria vetor para betas
Va = zeros(1,n,1); %cria vetor para Va
Vb = zeros(1,n,1); %cria vetor para Va
Vc = zeros(1,n,1); %cria vetor para Va
conjunto = zeros(1,n,1); %cria vetor para Va
vetor_vabc = string(1:n); %vetor de string para salvar as tensoes de fase que geram a cada vetor
vetor_vbig = string(1:n); %vetor de string para salvar as tensoes da celula grande
vetor_vsmall = string(1:n); %vetor de string para salvar as tensoes da celula pequena
vetor_g = zeros(1,n,1); %cria vetor para g
vetor_h = zeros(1,n,1); %cria vetor para h

%VARRENDO AS TENSOES
for a = 1:Nbig %varre a tensao Va_big
    for b = 1:Nbig %varre a tensao Vb_big
        for c = 1:Nbig %varre a tensao Vc_big
            n_vetores_big = n_vetores_big+1; %quantos vetores big existem
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

%o numero de redundancias eh dado pelo tamanho da matrix unique
num_redundancias = zeros(1,n_vetores);

matrix_vector = [vetor_alpha', vetor_beta', vetor_gama', num_redundancias', conjunto'];
%MATRIX COMPLETA
%coluna 1 = alpha
%coluna 2 = beta
%coluna 3 = gama
%coluna 4 = num redundancias
%coluna 5 = conjunto (estado do vetor big) //isso so esta na matrix_vector

matrix_quase_uniq = [vetor_alpha', vetor_beta', vetor_gama', num_redundancias'];
%MATRIX QUE SERA REDUZIDA
%coluna 1 = alpha
%coluna 2 = beta
%coluna 3 = gama
%coluna 4 = num redundancias
%coluna 5 = tensao Vab
%coluna 6 = tensao Vbig
%coluna 7 = tensao Vsmall
%coluna 8 = g
%coluna 9 = h

matrix_uniq = unique(matrix_quase_uniq, 'rows');

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
            dados{5}(j,dados{4}(j)) = vetor_vabc(z);
            dados{6}(j,dados{4}(j)) = vetor_vbig(z);
            dados{7}(j,dados{4}(j)) = vetor_vsmall(z);
            dados{8}(j,dados{4}(j)) = vetor_g(z);
            dados{9}(j,dados{4}(j)) = vetor_h(z);
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

figure
scatter(dados{1},dados{2},25,dados{4},'filled')
stg = blanks(1);
for z = 1:n_vet_unic
    for i=1:dados{4}(z,1)
        stg = append(stg, newline,dados{5}(z,i));
    end
    text(dados{1}(z),dados{2}(z),stg,'FontSize',5)
    stg = erase(stg,stg);
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Tensao aplicada em cada fase no mapa Alpha Beta')
% save_figure('tensao_total')

figure
scatter(dados{1},dados{2},25,dados{4},'filled')
stg = blanks(1);
for z = 1:n_vet_unic
    for i=1:dados{4}(z,1)
        stg = append(stg, newline,dados{6}(z,i));
    end
    text(dados{1}(z),dados{2}(z),stg,'FontSize',7)
    stg = erase(stg,stg);
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Tensao aplicada em cada fase na celula grande no mapa Alpha Beta')
% save_figure('tensao_grande')

figure
scatter(dados{1},dados{2},25,dados{4},'filled')
stg = blanks(1);
for z = 1:n_vet_unic
    for i=1:dados{4}(z,1)
        stg = append(stg, newline,dados{7}(z,i));
    end
    text(dados{1}(z),dados{2}(z),stg,'FontSize',5)
    stg = erase(stg,stg);
end
grid on
xlabel('Alpha')
ylabel('Beta')
title('Tensao aplicada em cada fase na celula pequena no mapa Alpha Beta')
% save_figure('tensao_pequena')



figure
gscatter(dados{8},dados{9},dados{4})
stg = blanks(1);
for z = 1:n_vet_unic
    numero_redundancias = int2str(dados{4}(z,1));
    text(dados{8}(z),dados{9}(z),numero_redundancias,'FontSize',11)
end
grid on
xlabel('g')
ylabel('h')


%% VALIDACAO NO TEMPO


f_rede = 50; %freq da rede
Vref = 3.4; %a tensao de referencia maxima eh 3.4607
fs = 1200; %frequencia de chaveamento (7 segmentos)
n_pontos_t = fs/f_rede;
n_pontos = n_pontos_t+1;
time = linspace(1/(2*fs),(1/f_rede+1/(2*fs)),n_pontos); %vetor de tempo
% time = 0:0.1:1;

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

for k=1:length(time)
    VA(k) = Vref*sin(time(k)*2*pi*f_rede);
    VB(k) = Vref*sin(time(k)*2*pi*f_rede-2*pi*f_rede/3);
    VC(k) = Vref*sin(time(k)*2*pi*f_rede-4*pi*f_rede/3);
    
    [g_ref(k), h_ref(k)]= transformada_fast(VA(k), VB(k), VC(k));
    
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
end





% %cria o gif
% fig = figure;
% filename = 'vetores_usados.gif';
% segundos = 20; %quantos segundos pra um ciclo da animacao
% precisao_color = 100;
% matrix_color = [zeros(1,precisao_color)' zeros(1,precisao_color)' zeros(1,precisao_color)'];
% for k=1:precisao_color
%     matrix_color(k,1)=k/abs(precisao_color);
%     k=k+1;
% end
% for N = 1:n_pontos
%       scatter(V1_g(N), V1_h(N),[], [delta_V1(N) 0 0],'filled')
%       colorbar;
%       hold on
%       scatter(V2_g(N), V2_h(N),[], [delta_V2(N) 0 0],'filled')
%       scatter(V3_g(N), V3_h(N),[], [delta_V3(N) 0 0],'filled')
%       scatter(g_ref(N), h_ref(N),'filled','blue')
%       %quiver(0, 0, g_ref(n), h_ref(n))
%       plot(g_ref,h_ref)
%       scatter(vetor_g,vetor_h,'black')
%       hold off
%       xlabel('G');
%       ylabel('H');
% %       xlim([-(N-1) N-1]);
% %       ylim([-(N-1) N-1]);
%       colormap(matrix_color)
%       caxis([0 1])
%       colorbar;
%       color = colorbar;
%       ylabel(color, 'Razão cíclica')
%       %drawnow
%       frame = getframe(fig);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if n == 1
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',segundos/n_pontos);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',segundos/n_pontos);
%       end
% end




figure
scatter(dados{8}(:,1),dados{9}(:,1))
hold on
scatter(g_ref, h_ref,'filled','red')
plot(g_ref, h_ref)
scatter(V1_g, V1_h,'filled','blue')
scatter(V2_g, V2_h,'filled','blue')
scatter(V3_g, V3_h,'filled','blue')
hold off
save_figure('aqui2')

% figure
% subplot(2,1,1);
% hold on
% plot(Va_1(:,1))
% plot(Va_2(:,1))
% plot(Va_3(:,1))
% hold off
% subplot(2,1,2); 
% plot(VA)
% title('VA')
% 
% figure
% subplot(2,1,1);
% hold on
% plot(Vb_1(:,1))
% plot(Vb_2(:,1))
% plot(Vb_3(:,1))
% hold off
% subplot(2,1,2); 
% plot(VB)
% title('VB')
% 
% figure
% subplot(2,1,1);
% hold on
% plot(Vc_1(:,1))
% plot(Vc_2(:,1))
% plot(Vc_3(:,1))
% hold off
% subplot(2,1,2); 
% plot(VC)
% title('VC')


figure
b = bar([cont_1(1:end-1)', cont_2(1:end-1)', cont_3(1:end-1)']);
% b = bar([cont_3(1:end-1)']);
legend('V1', 'V2', 'V3')
title('Opcoes de vetores para cada ponto de amostragem')

%%

initial_pont = 3;
m = initial_pont;
m2 = initial_pont+n_pontos_t/2;
%analisando ponto a ponto
figure
scatter(dados{8}(:,1),dados{9}(:,1))
hold on
scatter(g_ref(m), h_ref(m),'filled','red')
plot(g_ref(m), h_ref(m))
text(g_ref(m),h_ref(m),'Ref1','FontSize',9)
scatter(V1_g(m), V1_h(m),'filled','blue')
text(V1_g(m),V1_h(m),'V1','FontSize',9)
scatter(V2_g(m), V2_h(m),'filled','blue')
text(V2_g(m),V2_h(m),'V2','FontSize',9)
scatter(V3_g(m), V3_h(m),'filled','blue')
text(V3_g(m),V3_h(m),'V3','FontSize',9)
scatter(g_ref(m2), h_ref(m2),'filled','red')
plot(g_ref(m2), h_ref(m2))
text(g_ref(m2),h_ref(m2),'Ref2','FontSize',9)
scatter(V1_g(m2), V1_h(m2),'filled','blue')
text(V1_g(m2),V1_h(m2),'V1','FontSize',9)
scatter(V2_g(m2), V2_h(m2),'filled','blue')
text(V2_g(m2),V2_h(m2),'V2','FontSize',9)
scatter(V3_g(m2), V3_h(m2),'filled','blue')
text(V3_g(m2),V3_h(m2),'V3','FontSize',9)
hold off
contador=0;
text(-5.5,-3-contador,"Ref1 m="+m,'FontSize',9)
for k=1:cont_1(m)
    contador=contador+0.5;
    junta = append(int2str(Va_1(m,k)-3),int2str(Vb_1(m,k)-3),int2str(Vc_1(m,k)-3));
    text(-5.5,-3-contador,"V1: ("+junta+")",'FontSize',9)
end
text(-3.8,-3-contador,"D1="+round(delta_V1(m),2),'FontSize',9)
for k=1:cont_2(m)
    contador=contador+0.5;
    junta = append(int2str(Va_2(m,k)-3),int2str(Vb_2(m,k)-3),int2str(Vc_2(m,k)-3));
    text(-5.5,-3-contador,"V2: ("+junta+")",'FontSize',9)
end
text(-3.8,-3-contador,"D2="+round(delta_V2(m),2),'FontSize',9)
for k=1:cont_3(m)
    contador=contador+0.5;
    junta = append(int2str(Va_3(m,k)-3),int2str(Vb_3(m,k)-3),int2str(Vc_3(m,k)-3));
    text(-5.5,-3-contador,"V3: ("+junta+")",'FontSize',9)
end
text(-3.8,-3-contador,"D3="+round(delta_V3(m),2),'FontSize',9)

contador2=-0.3;
text(4.2,5.5-contador2,"Ref2 m="+m2,'FontSize',9)
for k=1:cont_1(m2)
    contador2=contador2+0.5;
    junta = append(int2str(Va_1(m2,k)-3),int2str(Vb_1(m2,k)-3),int2str(Vc_1(m2,k)-3));
    text(4.2,5.5-contador2,"V1: ("+junta+")",'FontSize',9)
end
text(2.8,5.5-contador2,"D1="+round(delta_V1(m2),2),'FontSize',9)
for k=1:cont_2(m2)
    contador2=contador2+0.5;
    junta = append(int2str(Va_2(m2,k)-3),int2str(Vb_2(m2,k)-3),int2str(Vc_2(m2,k)-3));
    text(4.2,5.5-contador2,"V2: ("+junta+")",'FontSize',9)
end
text(2.8,5.5-contador2,"D2="+round(delta_V2(m2),2),'FontSize',9)
for k=1:cont_3(m2)
    contador2=contador2+0.5;
    junta = append(int2str(Va_3(m2,k)-3),int2str(Vb_3(m2,k)-3),int2str(Vc_3(m2,k)-3));
    text(4.2,5.5-contador2,"V3: ("+junta+")",'FontSize',9)
end
text(2.8,5.5-contador2,"D3="+round(delta_V3(m2),2),'FontSize',9)
% save_figure("m_igual"+m)


%%
kg = 2.3;
kh = 3.4041;
% kg = 3.7453;
% kh = 3.7453;
figure
hold on
%PRIMEIRO QUADRANTE
[Vul_g_pt,Vul_h_pt,Vlu_g_pt,Vlu_h_pt,Vll_g_pt,Vll_h_pt,Vuu_g_pt,Vuu_h_pt,V1_g_pt,V1_h_pt,V2_g_pt,V2_h_pt,V3_g_pt,V3_h_pt,delta_V1_pt,delta_V2_pt,delta_V3_pt,cont_1_pt,Va_1_pt,Vb_1_pt,Vc_1_pt,cont_2_pt,Va_2_pt,Vb_2_pt,Vc_2_pt,cont_3_pt,Va_3_pt,Vb_3_pt,Vc_3_pt] = aproximacoes_fast(kg,kh,nivel);
scatter(kg,kh,'filled','green')
scatter(Vul_g_pt,Vul_h_pt,'filled','blue')
text(Vul_g_pt,Vul_h_pt,"Vul",'FontSize',9)
scatter(Vlu_g_pt,Vlu_h_pt,'filled','blue')
text(Vlu_g_pt,Vlu_h_pt,"Vlu",'FontSize',9)
scatter(Vll_g_pt,Vll_h_pt,'filled','blue')
text(Vll_g_pt,Vll_h_pt,"Vll",'FontSize',9)
scatter(Vuu_g_pt,Vuu_h_pt,'filled','blue')
text(Vuu_g_pt,Vuu_h_pt,"Vuu",'FontSize',9)
text((Vll_g_pt + Vul_g_pt)/2,Vll_h_pt-1,"1quad",'FontSize',9)
scatter(V1_g_pt,V1_h_pt, 10,'filled','red')
scatter(V2_g_pt,V2_h_pt, 10,'filled','red')
scatter(V3_g_pt,V3_h_pt, 10,'filled','red')
text(V1_g_pt,V1_h_pt-0.2,"D=" + num2str(delta_V1_pt),'FontSize',9)
text(V2_g_pt,V2_h_pt-0.2,"D=" + num2str(delta_V2_pt),'FontSize',9)
text(V3_g_pt,V3_h_pt-0.2,"D=" + num2str(delta_V3_pt),'FontSize',9)
% text(V1_g_pt,V1_h_pt-0.4,num2str(cont_1_pt),'FontSize',9)
% text(V2_g_pt,V2_h_pt-0.4,num2str(cont_2_pt),'FontSize',9)
% text(V3_g_pt,V3_h_pt-0.4,num2str(cont_3_pt),'FontSize',9)
contador_pt = -0.4;
for k=1:cont_1_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_1_pt(k)-3),int2str(Vb_1_pt(k)-3),int2str(Vc_1_pt(k)-3));
    text(V1_g_pt,V1_h_pt-0.6-contador_pt,"V1: ("+junta_pt+")",'FontSize',9)
end
contador_pt = -0.4;
for k=1:cont_2_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_2_pt(k)-3),int2str(Vb_2_pt(k)-3),int2str(Vc_2_pt(k)-3));
    text(V2_g_pt,V2_h_pt-0.6-contador_pt,"V2: ("+junta_pt+")",'FontSize',9)
end
contador_pt = -0.4;
for k=1:cont_3_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_3_pt(k)-3),int2str(Vb_3_pt(k)-3),int2str(Vc_3_pt(k)-3));
    text(V3_g_pt,V3_h_pt-0.6-contador_pt,"V3: ("+junta_pt+")",'FontSize',9)
end


%SEGUNDO QUADRANTE
[Vul_g_pt,Vul_h_pt,Vlu_g_pt,Vlu_h_pt,Vll_g_pt,Vll_h_pt,Vuu_g_pt,Vuu_h_pt,V1_g_pt,V1_h_pt,V2_g_pt,V2_h_pt,V3_g_pt,V3_h_pt,delta_V1_pt,delta_V2_pt,delta_V3_pt,cont_1_pt,Va_1_pt,Vb_1_pt,Vc_1_pt,cont_2_pt,Va_2_pt,Vb_2_pt,Vc_2_pt,cont_3_pt,Va_3_pt,Vb_3_pt,Vc_3_pt] = aproximacoes_fast(-kg,kh,nivel);
scatter(-kg,kh,'filled','green')
scatter(Vul_g_pt,Vul_h_pt,'filled','blue')
text(Vul_g_pt,Vul_h_pt,"Vul",'FontSize',9)
scatter(Vlu_g_pt,Vlu_h_pt,'filled','blue')
text(Vlu_g_pt,Vlu_h_pt,"Vlu",'FontSize',9)
scatter(Vll_g_pt,Vll_h_pt,'filled','blue')
text(Vll_g_pt,Vll_h_pt,"Vll",'FontSize',9)
scatter(Vuu_g_pt,Vuu_h_pt,'filled','blue')
text(Vuu_g_pt,Vuu_h_pt,"Vuu",'FontSize',9)
text((Vul_g_pt + Vll_g_pt)/2,Vll_h_pt-1,"2quad",'FontSize',9)
scatter(V1_g_pt,V1_h_pt, 10,'filled','red')
scatter(V2_g_pt,V2_h_pt, 10,'filled','red')
scatter(V3_g_pt,V3_h_pt, 10,'filled','red')
text(V1_g_pt,V1_h_pt-0.2,"D=" + num2str(delta_V1_pt),'FontSize',9)
text(V2_g_pt,V2_h_pt-0.2,"D=" + num2str(delta_V2_pt),'FontSize',9)
text(V3_g_pt,V3_h_pt-0.2,"D=" + num2str(delta_V3_pt),'FontSize',9)
% text(V1_g_pt,V1_h_pt-0.4,num2str(cont_1_pt),'FontSize',9)
% text(V2_g_pt,V2_h_pt-0.4,num2str(cont_2_pt),'FontSize',9)
% text(V3_g_pt,V3_h_pt-0.4,num2str(cont_3_pt),'FontSize',9)
contador_pt = -0.4;
for k=1:cont_1_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_1_pt(k)-3),int2str(Vb_1_pt(k)-3),int2str(Vc_1_pt(k)-3));
    text(V1_g_pt,V1_h_pt-0.6-contador_pt,"V1: ("+junta_pt+")",'FontSize',9)
end
contador_pt = -0.4;
for k=1:cont_2_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_2_pt(k)-3),int2str(Vb_2_pt(k)-3),int2str(Vc_2_pt(k)-3));
    text(V2_g_pt,V2_h_pt-0.6-contador_pt,"V2: ("+junta_pt+")",'FontSize',9)
end
contador_pt = -0.4;
for k=1:cont_3_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_3_pt(k)-3),int2str(Vb_3_pt(k)-3),int2str(Vc_3_pt(k)-3));
    text(V3_g_pt,V3_h_pt-0.6-contador_pt,"V3: ("+junta_pt+")",'FontSize',9)
end

%TERCEIRO QUADRANTE
[Vul_g_pt,Vul_h_pt,Vlu_g_pt,Vlu_h_pt,Vll_g_pt,Vll_h_pt,Vuu_g_pt,Vuu_h_pt,V1_g_pt,V1_h_pt,V2_g_pt,V2_h_pt,V3_g_pt,V3_h_pt,delta_V1_pt,delta_V2_pt,delta_V3_pt,cont_1_pt,Va_1_pt,Vb_1_pt,Vc_1_pt,cont_2_pt,Va_2_pt,Vb_2_pt,Vc_2_pt,cont_3_pt,Va_3_pt,Vb_3_pt,Vc_3_pt] = aproximacoes_fast(-kg,-kh,nivel);
scatter(-kg,-kh,'filled','green')
scatter(Vul_g_pt,Vul_h_pt,'filled','blue')
text(Vul_g_pt,Vul_h_pt,"Vul",'FontSize',9)
scatter(Vlu_g_pt,Vlu_h_pt,'filled','blue')
text(Vlu_g_pt,Vlu_h_pt,"Vlu",'FontSize',9)
scatter(Vll_g_pt,Vll_h_pt,'filled','blue')
text(Vll_g_pt,Vll_h_pt,"Vll",'FontSize',9)
scatter(Vuu_g_pt,Vuu_h_pt,'filled','blue')
text(Vuu_g_pt,Vuu_h_pt,"Vuu",'FontSize',9)
text((Vul_g_pt + Vll_g_pt)/2,Vll_h_pt+1,"3quad",'FontSize',9)
scatter(V1_g_pt,V1_h_pt, 10,'filled','red')
scatter(V2_g_pt,V2_h_pt, 10,'filled','red')
scatter(V3_g_pt,V3_h_pt, 10,'filled','red')
text(V1_g_pt,V1_h_pt-0.2,"D=" + num2str(delta_V1_pt),'FontSize',9)
text(V2_g_pt,V2_h_pt-0.2,"D=" + num2str(delta_V2_pt),'FontSize',9)
text(V3_g_pt,V3_h_pt-0.2,"D=" + num2str(delta_V3_pt),'FontSize',9)
% text(V1_g_pt,V1_h_pt-0.4,num2str(cont_1_pt),'FontSize',9)
% text(V2_g_pt,V2_h_pt-0.4,num2str(cont_2_pt),'FontSize',9)
% text(V3_g_pt,V3_h_pt-0.4,num2str(cont_3_pt),'FontSize',9)
contador_pt = -0.4;
for k=1:cont_1_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_1_pt(k)-3),int2str(Vb_1_pt(k)-3),int2str(Vc_1_pt(k)-3));
    text(V1_g_pt,V1_h_pt-0.6-contador_pt,"V1: ("+junta_pt+")",'FontSize',9)
end
contador_pt = -0.4;
for k=1:cont_2_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_2_pt(k)-3),int2str(Vb_2_pt(k)-3),int2str(Vc_2_pt(k)-3));
    text(V2_g_pt,V2_h_pt-0.6-contador_pt,"V2: ("+junta_pt+")",'FontSize',9)
end
contador_pt = -0.4;
for k=1:cont_3_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_3_pt(k)-3),int2str(Vb_3_pt(k)-3),int2str(Vc_3_pt(k)-3));
    text(V3_g_pt,V3_h_pt-0.6-contador_pt,"V3: ("+junta_pt+")",'FontSize',9)
end

%QUARTO QUADRANTE
[Vul_g_pt,Vul_h_pt,Vlu_g_pt,Vlu_h_pt,Vll_g_pt,Vll_h_pt,Vuu_g_pt,Vuu_h_pt,V1_g_pt,V1_h_pt,V2_g_pt,V2_h_pt,V3_g_pt,V3_h_pt,delta_V1_pt,delta_V2_pt,delta_V3_pt,cont_1_pt,Va_1_pt,Vb_1_pt,Vc_1_pt,cont_2_pt,Va_2_pt,Vb_2_pt,Vc_2_pt,cont_3_pt,Va_3_pt,Vb_3_pt,Vc_3_pt] = aproximacoes_fast(kg,-kh,nivel);
scatter(kg,-kh,'filled','green')
scatter(Vul_g_pt,Vul_h_pt,'filled','blue')
text(Vul_g_pt,Vul_h_pt,"Vul",'FontSize',9)
scatter(Vlu_g_pt,Vlu_h_pt,'filled','blue')
text(Vlu_g_pt,Vlu_h_pt,"Vlu",'FontSize',9)
scatter(Vll_g_pt,Vll_h_pt,'filled','blue')
text(Vll_g_pt,Vll_h_pt,"Vll",'FontSize',9)
scatter(Vuu_g_pt,Vuu_h_pt,'filled','blue')
text(Vuu_g_pt,Vuu_h_pt,"Vuu",'FontSize',9)
text((Vul_g_pt + Vll_g_pt)/2,Vll_h_pt+1,"4quad",'FontSize',9)
scatter(V1_g_pt,V1_h_pt, 10,'filled','red')
scatter(V2_g_pt,V2_h_pt, 10,'filled','red')
scatter(V3_g_pt,V3_h_pt, 10,'filled','red')
text(V1_g_pt,V1_h_pt-0.2,"D=" + num2str(delta_V1_pt),'FontSize',9)
text(V2_g_pt,V2_h_pt-0.2,"D=" + num2str(delta_V2_pt),'FontSize',9)
text(V3_g_pt,V3_h_pt-0.2,"D=" + num2str(delta_V3_pt),'FontSize',9)
% text(V1_g_pt,V1_h_pt-0.4,num2str(cont_1_pt),'FontSize',9)
% text(V2_g_pt,V2_h_pt-0.4,num2str(cont_2_pt),'FontSize',9)
% text(V3_g_pt,V3_h_pt-0.4,num2str(cont_3_pt),'FontSize',9)
contador_pt = -0.4;
for k=1:cont_1_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_1_pt(k)-3),int2str(Vb_1_pt(k)-3),int2str(Vc_1_pt(k)-3));
    text(V1_g_pt,V1_h_pt-0.6-contador_pt,"V1: ("+junta_pt+")",'FontSize',9)
end
contador_pt = -0.4;
for k=1:cont_2_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_2_pt(k)-3),int2str(Vb_2_pt(k)-3),int2str(Vc_2_pt(k)-3));
    text(V2_g_pt,V2_h_pt-0.6-contador_pt,"V2: ("+junta_pt+")",'FontSize',9)
end
contador_pt = -0.4;
for k=1:cont_3_pt
    contador_pt=contador_pt+0.2;
    junta_pt = append(int2str(Va_3_pt(k)-3),int2str(Vb_3_pt(k)-3),int2str(Vc_3_pt(k)-3));
    text(V3_g_pt,V3_h_pt-0.6-contador_pt,"V3: ("+junta_pt+")",'FontSize',9)
end
scatter(dados{8}(:,1),dados{9}(:,1))
% xlim([-3.5, 3.5])
% ylim([-3.5, 3.5])
hold off
grid;





