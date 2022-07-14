clear all
clc
close all
rng(27);
%% mainGeral
% % Declaração de variaveis globais
global MR NL VOLT 
% Global linprog
global DBAR NBAR PLOAD ANGLE Bref
% linha
global DLIN NLIN SB EB r x G B Bor FLIM 
global SIG NCV NCVCIR NCIRCV 
global fij fji
% gerador
global DGER NGER BARPG PGMIN PGMAX CPG
% constantes gerais
global PB NITER PAREI DIFMAX TOL
% PL
global CX AX BX Aiq Biq Vub Vlb NVAR NCAeq NLAeq NLBeq

global PLOAD_b QLOAD_b NBUS QLOAD PGWMAX DGW NGW NOGW

% GLOBAL DMR SOLVER
global PGEN R_LIN R_GER R_BAR FR TO FRref TOref

% GLOBAL FMINCON
global OBJF DADOSDMR DTEN DVIO MONLT QGEN NUMGERFIC BARGERFIC DADOSFMINCON CPG Bsh NOGEN

global linhasMonitoradas

%% ----------Carregando os Dados----------
 dir_atual = cd;   % Carregando diretorio atual do programa
 cd Dados;         % Acessando a pasta de casos
% Carregando Caso
IEEE118_v
load('Series_Merrick_sep_pu.mat', 'wind');
load('SAM_WindFarmsOutput_PU.mat');
cd (dir_atual);   % Retornando ao diretorio do programa
%% leitura de Dados
% Escolher função objetivo aqui
OBJF = 2;
inpdatlinprog;
inpdatdmr;
OBJF = 1;
INPDAT;
% Limpando dados
clear DBAR;
clear DGER;
clear DLIN;
clear DTEN;
%% Declaração de variaveis  
Nc = 1000; %Numero de casos a serem rodados, usar um valor alto Nc=2000;
dia = 0; %Dia do ano no arquivo de dados
hora = 18; %Hora do dia, para pegar apenas carga pesada
ictg = 1;

Max_FLUX = zeros(NLIN,NLIN); % armazenar o maior fluxo em cada linha;

Sum_FLUX = zeros(NLIN,NLIN); % armazenar a soma dos fluxos que ultrapassaram o limite da linha; 
TOL_BETA = 0.01; % tolerância de 1%

Min_VOLT = 2*ones(NLIN,NBUS);

% Rampas de geração
MVu = zeros(NLIN,NGER);
MVd= zeros(NLIN,NGER);

% ResultLinprog = zeros();
FOBS = zeros(NLIN,Nc);
% RESULTFMINCON = zeros(Nc,(NBUS+NGER)*2);
ResultLinprog = zeros(Nc,NGER+NBAR);

LambdaLinhas = zeros(Nc,NLIN);
fluxoGeral = zeros(Nc,NLIN);

PgwGeral = zeros(Nc,3);
cargaGeral = zeros(Nc, NBAR);

eolicasEUA = [pout_amazon_wind pout_Atla_VIII pout_block_island pout_Diamond_Vista pout_fenner pout_High_Sheldon ...
    pout_madison pout_maple pout_Timber_road pout_Top_Crop pout_waymart];

%% loop for
for ic = 1:1:Nc
    a = 0.95; 
    b = 1.05;
    variacao = a + (b-a)*rand();
    
    rA = normrnd(1,0.01,1,NBUS); %Agora usando a distribuição normal
    PLOAD = PLOAD_b.*rA*variacao;
    QLOAD = QLOAD_b.*rA*variacao;
    
    cargaGeral(ic,:) = PLOAD;
    PGW = PGWMAX.*wind(24*dia+hora,:); %Usando dados reais
% PGW = PGWMAX.*eolicasEUA(ic*dia+hora,:);
%     PGW = PGWMAX.*[pout_fenner(ic*dia+hora) pout_High_Sheldon(ic*dia+hora) pout_Atla_VIII(ic*dia+hora)];
    PgwGeral(ic,:) = PGW; 
%     PGW = PGWMAX.*wblrnd(0.8,5,1,3); %Usando a distribuição de Weibull
    hora = hora + 1;
    if hora>20
        hora = 18;
        dia = dia + 1;
    end
    
    for i= 1:NGW
        k = NOGW(i);
        PLOAD(k) = PLOAD(k) - PGW(i);
    end
    
    % Rodar linprog
    % salvar o melhor resultado do linprog
    if ic>1
        PG_ant = MR(1:NGER)'; % Geração anterior
    end
    OPF_linprog;
    LambdaLinhas(ic,:) = Vd.ineqlin;
    FOBS(ictg,ic) = FOB;
    MR = Vp;
%   Salvar o resultado do linprog
    ResultLinprog(ic,:) = MR';
    ANGLE = MR(NGER+1:NGER+NBAR)';
    PGEN = MR(1:NGER)';
    
    PG = MR(1:NGER)';
    if ic > 1 % somente a partir da 2a. rodada
        for i=1:NGER
            Difer(i) = PG(i) - PG_ant(i); % diferença entre o atual e anterior
            if (Difer(i) < 0) % rampa down
                Difer(i) = - Difer(i);
                if Difer(i) > MVd(ictg,i)
                    MVd(ictg,i) = Difer(i);
                end
            end
            if (Difer(i) > 0) % rampa up
                if Difer(i) > MVu(ictg,i)
                    MVu(ictg,i) = Difer(i);
                end
            end
        end
    end
    
    FR   = FRref;                % de
    TO   = TOref;                % para
% % Verifica fluxo nas linhas
    %--------- Cálculo dos fluxos potência ativa – Depois da linprog convergida com perdas sem limite de fluxo
    FLIM2 = FLIM+10000;
for i=1:max(size(DVIO))
    FLIM2(DVIO(i,1)) = DVIO(i,2)/PB;
end
    for I=1:NL
        NUM_L(I) = I;
        k = FR(I);
        m = TO(I); 
        Tk  = ANGLE(k); Tm = ANGLE(m); Gkm = G(I); Bkm = Bor(I);
        FPKM     = Bkm*(Tk-Tm) + 0.5*Gkm*(Tk-Tm)^2;
        FPMK     = Bkm*(Tm-Tk) + 0.5*Gkm*(Tm-Tk)^2;
        FPIJS = FPKM * PB;
        FPJIS = FPMK * PB;
        
        if FPIJS >= FPJIS
            MAIOR = FPIJS;
        else
            MAIOR = FPJIS;
        end
        
        fluxoGeral(ic,I) = MAIOR;
        
         if MAIOR > Max_FLUX(ictg,I)
            Max_FLUX(ictg,I) = MAIOR;
        end
        
        
       if MAIOR > FLIM2(I)*PB
            Difer = MAIOR - FLIM2(I)*PB;
            Sum_FLUX(ictg,I) = Sum_FLUX(ictg,I) + Difer;
       end
       
       if I==8
          Ft(:,ic) =  MAIOR/PB;
       end
       
    end
%     
%     if ic > 10 % ic é o número de iterações
%         % fluxo em uma LT em pu 
%         % a função teste Ft recebe o fluxo em uma LT, mas pode substituir pela função desejada
%         BETA = std(Ft)/mean(Ft);  % Desvio padrão sobre a média, 
%     end
%     
%     % Rodar DMR Solver
%     Main_DMR_EEFC
%     for i=1:NBUS
%         if (VOLT(i)<Min_VOLT(ictg,i))
%             Min_VOLT(ictg,i)    =  VOLT(i);
%         end
%     end
% 
%     
%     if ic > 10 % RODADAS é o número de iterações  
%       if BETA <= TOL_BETA 
%           break
%       end
%     end

    % verificar se tem violação nas linhas de transmissão
    % Finalizar com o FMINCON
%     if verificaViolacao()
%         disp("ocorreu violação")
%         OBJF = 2;
%    OPF_FMINCON_ALU
%    RESULTFMINCON(ic,:) = DADOSFMINCON;
%     end
end

%% Corte de carga
corteCarga = zeros(Nc,max(size(NUMGERFIC)));
for i=1:max(size(NUMGERFIC))
%     corteCarga(:,i) = RESULTFMINCON(:,2*NBUS+NUMGERFIC(i));
    corteCarga(:,i) = ResultLinprog(:,NUMGERFIC(i));
end

sumCorteCarga = sum(corteCarga);

[sumCorteCargaOrdenado,numGerFicCorte] = sort(-sumCorteCarga);
sumCorteCargaOrdenado = -sumCorteCargaOrdenado;
%% FOB médio
mediaFobs = zeros(NLIN,2);
for i=1:NLIN
    mediaFobs(i,1) = mean(FOBS(i,:));
    mediaFobs(i,2) = i;
end

[mediaFobsOrdenado,linhaFob] = sort(-mediaFobs(:,1));
mediaFobsOrdenado = -mediaFobsOrdenado;

%% Afundamento de Tensão
minimoVolt = 2;
minVoltGeral = zeros(NBAR,2);
for i=1:NBAR
    if Min_VOLT(1,i)<=minimoVolt
        minimoVolt = Min_VOLT(1,i);
        minVoltGeral(i,1) = Min_VOLT(1,i);
        minVoltGeral(i,2) = i;
    end
    minimoVolt = 2;
end

%  [Max_FLUX2, I_lin] = sort(Max_FLUX);%%Ranking. Returns a sort index I_lin which specifies how the elements of MAX_FLUX were rearranged to obtain the sorted output.
% [Sun_FLUX2, I_lin] = sort(Sum_FLUX);%%Ranking.
% 
% [ORDE_VOLT, I_bus] = sort(Min_VOLT);%%Ranking. Returns a sort index I_bus which specifies how the elements of ORDE _VOLT were rearranged to obtain the sorted output.
% [PIOR_VOLT, I_bus] = min(Min_VOLT); %%Returns the indices into operating dimension corresponding to the minimum values.
%% violação de fluxo, deixando mais facil de ver
violacaoPercentualFluxo = zeros(NLIN,3);

violacaoPercentualFluxo(:,:) = violacao(Max_FLUX(1,:), PB*FLIM2');
ssssss=1;

% sumViolacaoFlux = zeros(NLIN,2);
% for i = 1:NLIN
%     sumViolacaoFlux(i,1) = i;
%     sumViolacaoFlux(i,2) =  sum(violacaoPercentualFluxo(:,3,i));
% end
% 
% [sumViolacaoFluxOrdenado,numLinhaFlux] = sort(-sumViolacaoFlux(:,2));
% sumViolacaoFluxOrdenado = -sumViolacaoFluxOrdenado;
%% MVU//MVD
MVdGeral = zeros(NGER,2);
maior = 0;
for i=1:NGER
    if MVd(1,i)>=maior
        maior = MVd(1,i);
        MVdGeral(i,1) = MVd(1,i);
        MVdGeral(i,2) = BARPG(i);
    end
    maior = 0;
end

MVuGeral = zeros(NGER,2);
for i=1:NGER
    if MVu(1,i)>=maior
        maior = MVu(1,i);
        MVuGeral(i,1) = MVu(1,i);
        MVuGeral(i,2) = BARPG(i);
    end
    maior = 0;
end

sumLambdaLinhas = sum(LambdaLinhas);
[sumLambdaLinhasOrdenado,numLinhaLambda] = sort(-sumLambdaLinhas);
sumLambdaLinhasOrdenado = -sumLambdaLinhasOrdenado;
%% Graficos fluxos de linhas selecionadas
de = 100;
ate = 200;
cor = ["r" "g" "b" "k" "m"];
for j=0:6
    figure
    legenda(1:size(linhasMonitoradas,2)/5) = "Vazio";
    x1 = linspace(de,ate,ate-de +1);
    for i=1:5
        plot(x1,fluxoGeral(de:ate,linhasMonitoradas(i+j*5)),'lineWidth',2,'color',cor(i))
        hold on
        legenda(2*(i-1)+1) = "Fluxo na Linha: " + linhasMonitoradas(i+j*5);
        x2 = linspace(FLIM(linhasMonitoradas(i+j*5))*PB,FLIM(linhasMonitoradas(i+j*5))*PB,ate-de+1);
        plot(x1,x2, cor(i));
        legenda(2*(i-1)+2) = "Limite da Linha: " + linhasMonitoradas(i+j*5);
    end
    legend(legenda)
    title('Fluxo por linha')
    xlabel('Caso')
    ylabel('Fluxo(MW)')
    clear legenda
end
%% Geradores convencionais 
somaGeracao = sum(ResultLinprog(1:end,1:NGER),2);
tempo = linspace(de,ate,ate-de+1);
figure
j=0;
for i=1:NGER
    if PB*max(ResultLinprog(:,i))>1
        j = j+1;
        plot(tempo,PB*ResultLinprog(de:ate,i))
        legenda(j) = "Gerador Barra: " + NOGEN(i);
        if any(i==NUMGERFIC())
            legenda(j) = "Gerador Fictício " + NOGEN(i);
        end
        hold on
    end
end
hold off
% plot(linspace(1,Nc,Nc),PB*ResultLinprog(:,1),linspace(1,Nc,Nc),PB*ResultLinprog(:,2),linspace(1,Nc,Nc),...
% PB*ResultLinprog(:,3), linspace(1,Nc,Nc),PB*ResultLinprog(:,4), linspace(1,Nc,Nc),PB*ResultLinprog(:,5),...
% linspace(1,Nc,Nc),PB*somaGeracao(:))
% 
legend(legenda)
title('Geradores convencionais')
xlabel('Caso')
ylabel('Geração(MW)')
clear legenda
%% Geração eolica
somaEolicas = sum(PgwGeral(),2);
figure
legenda(1:size(PgwGeral,2)+1)="Vazio";
for i=1:size(PgwGeral,2)
plot(tempo,PB*PgwGeral(de:ate,i))
legenda(i) = "Eólica " + i;
hold on
end
legenda(i+1) = "Geração Eólica Total";
plot(tempo,PB*somaEolicas(de:ate))
legend(legenda)
title('Geradores Eólicos')
xlabel('Caso')
ylabel('Geração(MW)')
hold off 
%% Grafico de carga
figure
vetorCarga = PB*sum(cargaGeral,2);
plot(linspace(de,ate,ate-de+1),vetorCarga(de:ate))
title('Carga do Sistema')
xlabel('Caso')
ylabel('Carga(MW)')

ssss=1;


