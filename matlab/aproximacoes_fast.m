function [Vul_g,Vul_h,Vlu_g,Vlu_h,Vll_g,Vll_h,Vuu_g,Vuu_h,V1_g,V1_h,V2_g,V2_h,V3_g,V3_h,delta_V1,delta_V2,delta_V3] = aproximacoes_fast(g_ref,h_ref, nivel)
    if (g_ref>=0 && h_ref>=0)%G e H positivos, PRIMEIRO quadrante
        Vul_g = ceil(g_ref);
        Vul_h = floor(h_ref);
        Vlu_g = floor(g_ref);
        Vlu_h = ceil(h_ref);
        Vll_g = floor(g_ref);
        Vll_h = floor(h_ref);
        Vuu_g = ceil(g_ref);
        Vuu_h = ceil(h_ref);
    end
    if (g_ref<=0 && h_ref>=0)%G negativo e H positivo, SEGUNDO quadrante
        Vul_g = floor(g_ref);
        Vul_h = floor(h_ref);
        Vlu_g = ceil(g_ref);
        Vlu_h = ceil(h_ref);
        Vll_g = ceil(g_ref);
        Vll_h = floor(h_ref);
        Vuu_g = floor(g_ref);
        Vuu_h = ceil(h_ref);
    end
    if (g_ref<=0 && h_ref<=0)%G e H negativos, TERCEIRO quadrante
        Vul_g = floor(g_ref);
        Vul_h = ceil(h_ref);
        Vlu_g = ceil(g_ref);
        Vlu_h = floor(h_ref);
        Vll_g = ceil(g_ref);
        Vll_h = ceil(h_ref);
        Vuu_g = floor(g_ref);
        Vuu_h = floor(h_ref);
    end
    if (g_ref>=0 && h_ref<=0)%G positivo e H negativo, QUARTO quadrante
        Vul_g = ceil(g_ref);
        Vul_h = ceil(h_ref);
        Vlu_g = floor(g_ref);
        Vlu_h = floor(h_ref);
        Vll_g = floor(g_ref);
        Vll_h = ceil(h_ref);
        Vuu_g = ceil(g_ref);
        Vuu_h = floor(h_ref);
    end
    
    %REESCREVER
    V1_g = Vul_g;
    V1_h = Vul_h;
    V2_g = Vlu_g;
    V2_h = Vlu_h;    
    if ((abs(g_ref)+abs(h_ref))-(abs(Vul_g)+abs(Vul_h)))>0%escolhe Vuu
        V3_g = Vuu_g;
        V3_h = Vuu_h;
        delta_V1 = -(abs(h_ref)-abs(Vuu_h));
        delta_V2 = -(abs(g_ref)-abs(Vuu_g));
        delta_V3 = 1-delta_V1-delta_V2;
    else%escolhe Vll
        V3_g = Vll_g;
        V3_h = Vll_h;
        delta_V1 = abs(g_ref)-abs(Vll_g);
        delta_V2 = abs(h_ref)-abs(Vll_h);
        delta_V3 = 1-delta_V1-delta_V2; 
    end
end

