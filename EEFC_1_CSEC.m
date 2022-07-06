function f = EEFC_1_CSEC(X)

% Subrotina para calculo do Jacobiano e dos resíduos - Load Flow Convencional

global DBAR DTEN NBUS BREF VOLT VTMIN VTMAX ANGLE PLOAD QLOAD QBAR TIPO
global DLIN NL FR TO G B Bsh TAP CAPOR
global DGER NGER NOGEN NGERD PGEN PGMIN PGMAX QGEN QGMIN QGMAX NPG
global PB CPU_TIME NEQ NVAL NTOT R_BAR R_GER R_LIN S_FP S_FQ  
global VREF ANGREF NBR
% Associadas as barra de controle

global NBC VCON BCON


%% ============== Atribuição das variáveis

%-----------  Atribuição das variáveis de tensão e ângulo

V    = X(1:NBUS);
TETA = X(NBUS+1:2*NBUS);

%-----------  Atribuição da variável PGEN/QGEN da barra de referência

k       = NGERD(BREF); % Número da barra onde está o gerador de referencia
PGEN(k) = X(2*NBUS+1);
QGEN(k) = X(2*NBUS+2);

%----------- Cálculo dos PGEN's: Despacho proporcional
if NGER > 1
    for i = 1:NGER     
       
       if TIPO(NOGEN(i)) == 1;
       BG      = NGERD(BREF); 
       prop    = PGMAX(i)/PGMAX(BG);
       PGEN(i) = prop * PGEN(BG); %Despacho proporcional
        
       end
    end    
end

%-----------  Atribuição das variáveis das barras de controle

AUX2 = 0;
if NBC > 0 % tem barra de controle 
    
    for i = 1:NBUS
        
       if TIPO(i) == 1 
           
       AUX2     = AUX2 + 1;
       BC       = NGERD(i);       
       QGEN(BC) = X(2*NBUS+NBR+AUX2);  
       
       end    
    end
end
PG   = PGEN;
QG   = QGEN;

% ========================================================================
%% Restrição Não Linear de Igualdade - EEFC

%---------------------- Contribuição dos fluxos de potência

S_FP(1:NBUS,1) = 0;
S_FQ(1:NBUS,1) = 0;

for I=1:NL
    k  = FR(I);
    m  = TO(I);
    Vk = V(k);
    Vm = V(m);
    Tk = TETA(k);
    Tm = TETA(m);
    
    Gkm   = G(I);
    Bkm   = B(I);
    Bshkm = Bsh(I);
 
%==========================================================================
%
%                 Equações de Fluxos 
%
%==========================================================================

%--------- Fluxo Ativo de k-m
    FPKM = Gkm*Vk^2 - Vk*Vm*(Gkm*cos(Tk - Tm) + Bkm*sin(Tk - Tm));
    
%--------- Fluxo Ativo de m-k 
    FPMK = Gkm*Vm^2 - Vk*Vm*(Gkm*cos(Tk - Tm) - Bkm*sin(Tk - Tm));
    
%--------- Fluxo Reativo de k-m 
    FQKM = -(Bkm+Bshkm)*Vk^2 + Vk*Vm* (Bkm*cos(Tk - Tm) - Gkm*sin(Tk - Tm));

%--------- Fluxo Reativo de m-k  
    FQMK = Vk*Vm*(Bkm*cos(Tk - Tm) + Gkm*sin(Tk - Tm)) - Vm^2*(Bkm + Bshkm);
    
%----------- Cálculo da Soma de Fluxo Ativo em cada barra
 
    S_FP(k,1)  = S_FP(k,1) + FPKM;
    S_FP(m,1)  = S_FP(m,1) + FPMK;    
                    
%----------- Cálculo da Soma de Fluxo Reativo em cada barra
 
    S_FQ(k,1)  = S_FQ(k,1) + FQKM;
    S_FQ(m,1)  = S_FQ(m,1) + FQMK;    
end
 
%% ================== Contribuições para as equações


%---------------------- Contribuição das Cargas Ativa e Reativa
    
for k = 1:NBUS   
    hp(k,1) = - PLOAD(k) - S_FP(k,1);
    hq(k,1) = - QLOAD(k) - S_FQ(k,1);
end

%---------------------- Contribuição dos Geradores
    
for i = 1:NGER
    
    k = NOGEN(i); % Número da barra onde está o gerador 
    
    hp(k,1) = hp(k,1) + PG(i);
    hq(k,1) = hq(k,1) + QG(i);
  
end

%---------------------- Contribuição da barra de referência
    
    hr(1,1) = V(BREF)    - VREF;
    hr(2,1) = TETA(BREF) - ANGREF;

        
%--------------------- Contribuição das barras de controle
 

if NBC > 0  
    
  for i = 1:NBC  
      
     BC      = BCON(i);
     hc(i,1) = V(BC) - VCON(i); 
     
  end
  
  f = [hp; hq; hr; hc];
    
end



end

