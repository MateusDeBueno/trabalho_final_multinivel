clear all
close all
clc

N = 7; %numero de niveis
Vdc = 1;
V = (0:1:N-1)*Vdc - floor(N/2);  % opcoes de tensoes de fase
%- floor(N/2);

%cria as tensoes de fase
Va = zeros(1,N);
Vb = zeros(1,N);
Vc = zeros(1,N);
for j = 1:N
    Va(j) = V(j);
    Vb(j) = V(j);
    Vc(j) = V(j);
end

n_vetores = 0;
vetor_g = zeros(1,N^3,1); %cria vetor para g
vetor_h = zeros(1,N^3,1); %cria vetor para h
vetor = string(1:N^3); %vetor de string para salvar as tensoes de linha que geram a cada vetor
Vab = zeros(1,N^3,1); 
Vbc = zeros(1,N^3,1); 
Vca = zeros(1,N^3,1); 

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
%coluna 2 = h
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


% figure
% scatter3(dados{4},dados{5},dados{6})
% stg = blanks(1);
% for z = 1:n_vet_unic
%     numero_redundancias = int2str(dados{3}(z,1));
%     text(dados{4}(z),dados{5}(z), dados{6}(z), numero_redundancias,'FontSize',11)
% end
% xlabel('Vab')
% ylabel('Vbc')
% zlabel('Vca')
% % view(65,25)



%cria o gif
fig = figure;
filename = 'tensoes_linha.gif';
frames_linha = 180; %quantos frames neste gif
for n = 1:frames_linha
      scatter3(dados{4},dados{5},dados{6})
%       stg = blanks(1);
%       for z = 1:n_vet_unic
%           numero_redundancias = int2str(dados{3}(z,1));
%           text(dados{4}(z),dados{5}(z), dados{6}(z), numero_redundancias,'FontSize',11)
%       end
      xlabel('Vab')
      ylabel('Vbc')
      zlabel('Vca')
      view(n,25)
      %drawnow
      frame = getframe(fig);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',10/frames_linha);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',10/frames_linha);
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
save_figure('aqui')


%% VALIDACAO NO TEMPO


f_rede = 50; %freq da rede
m = 3.4; %indice de modulacao ??????????????????
n_pontos = 100;
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
delta_V1 = zeros(1,n_pontos);
delta_V2 = zeros(1,n_pontos);
delta_V3 = zeros(1,n_pontos);

for k=1:length(time)
    VA(k) = m*sin(time(k)*2*pi*f_rede);
    VB(k) = m*sin(time(k)*2*pi*f_rede-2*pi*f_rede/3);
    VC(k) = m*sin(time(k)*2*pi*f_rede-4*pi*f_rede/3);
    
    [g_ref(k), h_ref(k)]= transformada_fast(VA(k), VB(k), VC(k));
    
%     Vul_g(k) = ceil(g_ref(k));
%     Vul_h(k) = floor(h_ref(k));
%     Vlu_g(k) = floor(g_ref(k));
%     Vlu_h(k) = ceil(h_ref(k));
%     Vll_g(k) = floor(g_ref(k));
%     Vll_h(k) = floor(h_ref(k));
%     Vuu_g(k) = ceil(g_ref(k));
%     Vuu_h(k) = ceil(h_ref(k));
% 
%     
%     %REESCREVER
%     V1_g(k) = Vul_g(k);
%     V1_h(k) = Vul_h(k);
%     V2_g(k) = Vlu_g(k);
%     V2_h(k) = Vlu_h(k);    
%     if ((g_ref(k)+h_ref(k))-(Vul_g(k)+Vul_h(k)))>0 %escolhe Vuu
%         V3_g(k) = Vuu_g(k);
%         V3_h(k) = Vuu_h(k);
%         delta_ul(k)=-(h_ref(k)-Vuu_h(k));
%         delta_lu(k)=-(g_ref(k)-Vuu_g(k));
%         delta_uu(k)=1-delta_ul(k)-delta_lu(k);
%         delta_V1(k) = delta_ul(k);
%         delta_V2(k) = delta_lu(k);
%         delta_V3(k) = delta_uu(k);
%     else                                           %escolhe Vll
%         V3_g(k) = Vll_g(k);
%         V3_h(k) = Vll_h(k);
%         delta_ul(k)=g_ref(k)-Vll_g(k);
%         delta_lu(k)=h_ref(k)-Vll_h(k);
%         delta_ll(k)=1-delta_ul(k)-delta_lu(k); 
%         delta_V1(k) = delta_ul(k);
%         delta_V2(k) = delta_lu(k);
%         delta_V3(k) = delta_ll(k);
%     end
    [Vul_g(k),Vul_h(k),Vlu_g(k),Vlu_h(k),Vll_g(k),Vll_h(k),Vuu_g(k),Vuu_h(k),V1_g(k),V1_h(k),V2_g(k),V2_h(k),V3_g(k),V3_h(k),delta_V1(k),delta_V2(k),delta_V3(k)] = aproximacoes_fast(g_ref(k),h_ref(k),N);

end


% DELETAR
V1 = [V1_g; V1_h];
V2 = [V2_g; V2_h];
V3 = [V3_g; V3_h];


% scatter(V1')


% figure
% plot(VA)
% hold on
% plot(VB)
% plot(VC)
% hold off
% xlim([0 n_pontos])
% legend('Va', 'Vb', 'Vc')
% title('Tensoes de linha')



% figure
% plot(g_ref)
% hold on
% plot(h_ref)
% hold off
% xlim([0 n_pontos])
% legend('g', 'h')



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



% figure
% scatter(dados{1},dados{2})
% hold on
% scatter(g_ref,h_ref)
% hold off
% grid on


% teste = 50;
% figure
% scatter(dados{1},dados{2})
% hold on
% scatter(V1_g(teste), V1_h(teste),[], [delta_V1(teste) 0 0],'filled')
% scatter(V2_g(teste), V2_h(teste),[], [delta_V2(teste) 0 0],'filled')
% scatter(V3_g(teste), V3_h(teste),[], [delta_V3(teste) 0 0],'filled')
% scatter(g_ref(teste),h_ref(teste),'filled','blue')
% precisao_color = 100;
% matrix_color = [zeros(1,precisao_color)' zeros(1,precisao_color)' zeros(1,precisao_color)'];
% for k=1:precisao_color
%     matrix_color(k,1)=k/abs(precisao_color);
%     k=k+1;
% end
% colormap(matrix_color)
% caxis([0 1])
% colorbar;
% hold off







%cria o gif
fig = figure;
filename = 'vetores_usados3.gif';
segundos = 20; %quantos segundos pra um ciclo da animacao
precisao_color = 100;
matrix_color = [zeros(1,precisao_color)' zeros(1,precisao_color)' zeros(1,precisao_color)'];
for k=1:precisao_color
    matrix_color(k,1)=k/abs(precisao_color);
    k=k+1;
end
for n = 1:n_pontos
      scatter(V1_g(n), V1_h(n),[], [delta_V1(n) 0 0],'filled')
      colorbar;
      hold on
      scatter(V2_g(n), V2_h(n),[], [delta_V2(n) 0 0],'filled')
      scatter(V3_g(n), V3_h(n),[], [delta_V3(n) 0 0],'filled')
      scatter(g_ref(n), h_ref(n),'filled','blue')
      %quiver(0, 0, g_ref(n), h_ref(n))
      plot(g_ref,h_ref)
      scatter(dados{1},dados{2},'black')
      hold off
      xlabel('G');
      ylabel('H');
      xlim([-(N-1) N-1]);
      ylim([-(N-1) N-1]);
      colormap(matrix_color)
      caxis([0 1])
      colorbar;
      color = colorbar;
      ylabel(color, 'Razão cíclica')
      %drawnow
      frame = getframe(fig);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',segundos/n_pontos);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',segundos/n_pontos);
      end
end

