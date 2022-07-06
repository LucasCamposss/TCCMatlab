function inpdat

% Subrotina para entrada de dados

global DBAR DTEN NBUS BREF VOLT  VTMIN VTMAX ANGLE PLOAD QLOAD LBDP LBDQ QBAR
global DLIN NL FR TO G B Bsh TAP CAPOR DVIO MONLT
global DGER NGER NOGEN NGERD PGEN PGMIN PGMAX QGEN QGMIN QGMAX COSTPG
global PB X LAMBDA CPU_TIME OBJF R_BAR R_GER R_LIN
global PLOAD_b QLOAD_b
global DGW NGW NOGW NGERW PGW PGWMIN PGWMAX DGW


% Edimar: 22/02/2013

%-------------- Conjunto de dados de Barra (c�digo DBAR) --------------
PB = 100;

[NBUS,NCOL] = size(DBAR);            % n�mero de barras

for i = 1:NBUS
    
    TIPO = DBAR(i,2);                % tipo de barra
    
    if(TIPO == 2)
        BREF = i;                    % barra de refer�ncia
    end
    
    GBTEN(i) = DBAR(i,3);              % Grupo base de tens�o
    
    VOLT(i)  = DBAR(i,4);            % tens�o das barras
    ANGLE(i) = DBAR(i,5)*pi/180;     % �ngulo das barras em GRAUS
    PLOAD(i) = 1.0*DBAR(i,10)/PB;        % carga ativa
    QLOAD(i) = 1.0*DBAR(i,11)/PB;        % carga reativa
    QBAR(i)  = 0;%DBAR(i,12)/PB;        % inje��o fixa de reativo 
    % Capacitor (+)
    % Indutor   (-)
    QLOAD(i) = QLOAD(i) - QBAR(i);      % carga reativa total
    
    
end
PLOAD_b = PLOAD;
QLOAD_b = QLOAD;

%------------- Conjunto de dados de Linhas(c�digo DLIN) ---------------------

[NL,NCOL] = size(DLIN);

for i = 1:NL
    
    FR(i)   = DLIN(i,1);                % de
    TO(i)   = DLIN(i,2);                % para
    NCIR(i) = i;                        % n�mero do circuito
    ROR     = DLIN(i,3)/100;            % resist�ncia s�rie
    XOR     = DLIN(i,4)/100;            % reat�ncia s�rie
    
    DENO    = (ROR^2+XOR^2);
    G(i)    = +ROR/DENO;                % condut�ncia s�rie
    B(i)    = -XOR/DENO;                % suscept�ncia s�rie
    Bsh(i)  = 0.5*DLIN(i,5)/PB;         % suscept�ncia shunt
    TAP(i) = DLIN(i,6);                 % tap do trafo
    
    if(TAP(i) == 0)
        TAP(i) = 1;
    end
    
    TAP(i)    = 1/TAP(i);
    TAPMIN(i) = DLIN(i,7);             % tap m�nimo do trafo
    TAPMAX(i) = DLIN(i,8);             % tap m�ximo do trafo
    CAPOR(i)  = DLIN(i,10)/PB;         % tap m�ximo do trafo  
    MONLT(i) = 0;       % Se 0: Limite da LT n�o � monitorado
end

%------------ Monitoramento de Circuito Violado ------------------

[NVIO,NCOL] = size(DVIO);

for i=1:NVIO
    K        = DVIO(i,1);
    MONLT(K) = 1;       % Se 1: Limite da LT � monitorado
    CAPOR(K) = DVIO(i,2)/PB;
end


%------------ Conjunto de dados de tens�o (c�digo DTEN) -----------------------------

[NGBTEN,NCOL] = size(DTEN);         % n�mero de grupos de tens�o

for i=1:NGBTEN
    GBT(i)   = DTEN(i,1);
    TENMIN(i)= DTEN(i,2);
    TENMAX(i)= DTEN(i,3);
end

for i=1:NBUS
    for j=1:NGBTEN
        if(GBTEN(i)==GBT(j))
            VTMIN(i)=TENMIN(j);
            VTMAX(i)=TENMAX(j);        
        end
    end
end

%----------- Conjunto de dados de Gera��o Ativa (c�digo DGER) ----------

[NGER,NCOL] = size(DGER);           % n�mero de geradores

for i = 1:NGER
    
    NOGEN(i)        = DGER(i,1);             % barras geradoras
    NGERD(NOGEN(i)) = i;                     % associa a barra ao respectivo gerador
    VGEN(i)         = VOLT(NOGEN(i));        % tensao especificada para o gerador
    PGMIN(i)        = DGER(i,2)/PB;          % Limite m�nimo de gera��o ativa
    PGMAX(i)        = DGER(i,3)/PB;          % Limite m�ximo de gera��o ativa
  
    if OBJF == 1             % Minimo Custo de Gera��o
     COSTPG(i) = DGER(i,4);  % custo de gera��o US$/MW
    end
    
    if OBJF == 2             % Minimas Perdas
     COSTPG(i) = 100; 
    end
    
    
    PGEN(i)= 0.1; %DBAR(NOGEN(i),6)/PB;     % gera��o Ativa
    QGEN(i)=DBAR(NOGEN(i),7)/PB;     % gera��o reativa
    
    QGMIN(i)=DBAR(NOGEN(i),8)/PB;    % limite m�nimo de gera��o de reativo
    QGMAX(i)=DBAR(NOGEN(i),9)/PB;    % limite m�ximo de gera��o de reativo

%     -----Conjunto de dados de Gera��o E�lica (c�digo DGW)------
    [NGW,NCOL] = size(DGW);
    
    for j = 1:NGW
        NOGW(j)         = DGW(j,1);
        NGERW(NOGW(j))  = j;
        
        PGWMIN(j)        = DGW(j,2)/PB;
        PGWMAX(j)        = DGW(j,3)/PB;
    end
end





