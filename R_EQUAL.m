function [cin,ceq]= R_EQUAL(x)

global NBUS BREF PLOAD QLOAD 
global NL FR TO G B Bsh TAP CAPOR DVIO MONLT
global NGER NOGEN NGERD 
global PB OBJF
% Flags
global FLG_LIM

%=============================================

V    = x(1:NBUS);
TETA = x(NBUS+1:2*NBUS);
PG   = x(2*NBUS+1:2*NBUS+NGER);
QG   = x(2*NBUS+NGER+1:2*NBUS+2*NGER);

if OBJF==4 % referente a Qinj
Qinj(1:NBUS,1) = x(2*NBUS+2*NGER+1:2*NBUS+2*NGER+NBUS); % Referente a Qinj
end

V = V';
TETA = TETA';
PG = PG';
QG = QG';

% ========================================================================
% Restrição Não Linear de Desigualdade - Fluxo Violado
cin    = [];
NCV    = 0;
NCVCIR = [];
NCIRCV = [];

% ========================================================================
% Restrição Não Linear de Desigualdade - EEFC

%---------------------- Contribuição dos fluxos de potência

Sum_fluxP(1:NBUS,1) = 0;
Sum_fluxQ(1:NBUS,1) = 0;

for I=1:NL
    k = FR(I);
    m = TO(I);
    Vk = V(k);
    Vm = V(m);
    Tk = TETA(k);
    Tm = TETA(m);
    
    Gkm = G(I);
    Bkm = B(I);
    
    Bshkm = Bsh(I);
       
 % Cálculo da Soma de Fluxo Ativo em cada barra
    
    COSkm = cos(Tk - Tm);
    SENkm = sin(Tk - Tm);
    
    FluxP_km = + Gkm*Vk^2 - Vk*Vm * (Gkm*COSkm + Bkm*SENkm);
    FluxP_mk = + Gkm*Vm^2 - Vm*Vk * (Gkm*COSkm - Bkm*SENkm);
 
    Sum_fluxP(k,1)  = Sum_fluxP(k,1) + FluxP_km;
    Sum_fluxP(m,1)  = Sum_fluxP(m,1) + FluxP_mk;
                    
 % Cálculo da Soma de Fluxo Reativo em cada barra
 
    FluxQ_km = -(Bkm+Bshkm)*Vk^2 + Vk*Vm* (Bkm*COSkm - Gkm*SENkm);
    FluxQ_mk = -(Bkm+Bshkm)*Vm^2 + Vm*Vk* (Bkm*COSkm + Gkm*SENkm);
    
    Sum_fluxQ(k,1)  = Sum_fluxQ(k,1) + FluxQ_km;
    Sum_fluxQ(m,1)  = Sum_fluxQ(m,1) + FluxQ_mk;
% ======================================================================== 
%                 Tratamento de Fluxo Violado 
% ========================================================================
  if (MONLT(I) == 1) % Entra se a LT I estiver na lista de monitoradas

    if (FluxP_km >= FluxP_mk)
        
        NCV=NCV+1;
        NCVCIR(I) = NCV;             % Apontadores
        NCIRCV(NCV) = I;
        
        cin(NCV) = FluxP_km - CAPOR(I);        
    end
    if (FluxP_km < FluxP_mk)
        
        NCV=NCV+1;
        NCVCIR(I) = NCV;      % Apontadores
        NCIRCV(NCV) = I;
        
        cin(NCV) = FluxP_mk - CAPOR(I);        
    end

  end
end

%---------------------- Contribuição das Cargas Ativa e Reativa
    
for k = 1:NBUS   
    hp(k,1) = - PLOAD(k) - Sum_fluxP(k,1);
    hq(k,1) = - QLOAD(k) - Sum_fluxQ(k,1);
    if OBJF==4      % referente a Qinj
    hq(k,1) = hq(k,1) + Qinj(k,1);
    end
end


%---------------------- Contribuição dos Geradores
    
for i = 1:NGER
    
    k = NOGEN(i); % Número da barra onde está o gerador 
    
    hp(k,1) = hp(k,1) + PG(i,1);
    hq(k,1) = hq(k,1) + QG(i,1);
end

ceq = [hp; hq];

if (FLG_LIM == 0)
cin=[];
end
end
