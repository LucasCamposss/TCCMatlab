%-------------------------------------------------------------
%------------ Subrotina para entrada de dados-----------------
%-------------------------------------------------------------


function inpdat

% barra
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
global CX Aeq Beq Aiq Biq Vub Vlb NVAR NCAeq NLAeq NLBeq NUMGERFIC BARGERFIC DVIO linhasMonitoradas

global numCorteVento barCorteVento


%----------Dados das barras----------

[NBAR,AUX]=size(DBAR);
PB = 100;

for i=1:NBAR
   PLOAD(i)=1.0*DBAR(i,10)/100; %%%IMPORTANTE(EU REDUZI A CARGA PRA TENTAR RODAR O LINPROG DO 118) 
end


%----------Dados das linhas----------

[NLIN,AUX]=size(DLIN);

for i=1:NLIN
   SB(i)  = DLIN(i,1); %Vetor da barra de partida (STARTBUS)
   EB(i)  = DLIN(i,2); %Vetor da barra de chegada (ENDBUS);
   r(i)   = DLIN(i,3)/100; % resistência série
   x(i)   = DLIN(i,4)/100; % reatância série
   G(i)   = r(i)/(r(i)^2+x(i)^2); % condutância série
   B(i)   = x(i)/(r(i)^2+x(i)^2); % susceptância
   
   Bor(i) = 1/x(i); % susceptância série (gamma)
   FLIM(i)=DLIN(i,10)/100; %Fluxo limite em pu-MW
%   Mudanças pra conseguir rodar o IEEE118
   if FLIM(i)==0
    FLIM(i) = 50;
   end
end
for i=1:max(size(DVIO))
    FLIM(DVIO(i,1)) = DVIO(i,2)/100;
    linhasMonitoradas(i) = DVIO(i,1);
end

%----------Dados dos geradores----------

[NGER,AUX]=size(DGER);
NUMEROGERFIC = 0;
numCorte = 0;
for i=1:NGER
   BARPG(i) = DGER(i,1); % Apontador barra/gerador. BARPG(2)=4 significa BAR do GER 2 é 4
   PGMIN(i) = DGER(i,2)/PB;
   PGMAX(i) = DGER(i,3)/PB;
   CPG(i)   = DGER(i,4); % Custo operacional do gerador
   %%%%%%%%%%% Acrescentar: Déficit e perdas
    if CPG(i) < 400 && CPG(i) > 0
        CPG(i) = 1 + rand();
    elseif CPG(i) < 0
        numCorte = numCorte+1;
        numCorteVento(numCorte) = i;
        barCorteVento(numCorte) = BARPG(i);
    else
        NUMEROGERFIC = NUMEROGERFIC+1;
        NUMGERFIC(NUMEROGERFIC) = i;
        BARGERFIC(NUMEROGERFIC) = BARPG(i);
    end
    
end



%Inicializando variáveis

NVAR = NGER + NBAR; %Número de variáveis

