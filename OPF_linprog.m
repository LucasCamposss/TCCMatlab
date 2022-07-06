%----------Constantes do Problema----------

NITER     = 0;             % n�mero de itera��es
NITER_MAX = 10;

PAREI  = 10;            % convergencia 
PB     = 100;           % pot�ncia base do sistema
DIFMAX = 10;
TOL    = 0.0001;        % toler�ncia

%----Montagem da AX, BX, CX e outros vetores-----

%        NLAX=> Num. de linhas de AX = N�mero de Barras
%        NCAX=> Num. de colunas de AX = |NGER + NBAR-1 |
%        NLBX=> Num. de limhas de BX
%
NCAX = NGER + NBAR;
NLAX = NBAR + 1; % +1 devido  areferencia
NLBX = NLAX;
AX   = zeros(NLAX,NCAX);
BX   = zeros(NLAX,1);

for i=1:NCAX;
   CX(i)  = 0;
   Vlb(i) = 0;		
   Vub(i) = 0;
end

%---------Contribui�ao dos Geradores em AX, CX, VLB e VUB-----

for i=1:NGER;
    bar = BARPG(i);
    AX(bar,i) = 1;  % Representa os geradores
    CX(i)     = CPG(i); %custo de gera��o (US$/MWh)
    Vub(i)    = PGMAX(i);
    Vlb(i)    = 0;       % Referente a minima gera�ao
end
%-----------Contribui�ao da Rede em AX, Vub e Vlb-------------
Bbus = zeros(NBAR,NBAR);
for i = 1:NLIN
    NUMLIN = i;
    Gama  = Bor(NUMLIN);        % Representa as susceptancias das linhas
    FR    = SB(NUMLIN);         % de
    TO    = EB(NUMLIN);         % para
    Bbus(FR,FR) = Bbus(FR,FR) + Gama; 
    Bbus(TO,TO) = Bbus(TO,TO) + Gama;
    Bbus(FR,TO) = Bbus(FR,TO) - Gama;
    Bbus(TO,FR) = Bbus(TO,FR) - Gama;
end

AX(1:NBAR,NGER+1:NGER+NBAR) = - Bbus;

for i = NGER+1:NGER+NBAR;  
      Vlb(i) = - pi;
      Vub(i) = + pi;
end
%----------------Contribui�ao das Cargas em BX----------------

BX(1:NBAR)=PLOAD(1:NBAR);

%---------------- Contribui��o da Barra de Refer�ncia -------

AX(NLAX,NCAX) = 1;

X0=[];  % Chute inicial!!

NEI = NBAR + 1; % N�mero de equa��es de igualdade.

%-----------------------chamada do PL-----------------------

Ai=[];
Bi=[];
tic;
[Vp,fval,exitflag,output,LAMBDA] = linprog(CX,Ai,Bi,AX,BX,Vlb,Vub);

%--------------- Solu��o encontrada ---------------------
for I=1:NGER            % Geradores
   DESLOC = I;
   PG(I)  = PB * Vp(DESLOC);
end
for I=1:NBAR            % �ngulos
   DESLOC = NGER+I;
   ANGLE(I)  = Vp(DESLOC);
   ANG_G(I)  = ANGLE(I)*180/pi;
end
for i=1:NLIN
    rl = r(i);
    g  = G(i);
    y  = Bor(i);
    FR = SB(i);
    TO = EB(i);
    
    ang    = ANGLE(FR) - ANGLE(TO);
    fij(i) = y * (ANGLE(FR)-ANGLE(TO));
    fji(i) = y * (ANGLE(TO)-ANGLE(FR));
    
   
end

%============== Incluindo perdas e limite de fluxo ===========
ERRO = 0;
ITER = 0;

%----------------------loop do programa-----------------------

while (NITER<=NITER_MAX) & (DIFMAX>=TOL)
   
   NITER = NITER+1;

   %--------------- Inclus�o do Fluxo Violado e perdas ---------------
   PINJ = zeros(1,NBAR);
   for i=1:NLIN
    rl = r(i);
    g  = G(i);
    y  = Bor(i);
    FR = SB(i);
    TO = EB(i);
    
    ang       = ANGLE(FR)-ANGLE(TO);
    ang2      = ang*ang;
    perdas(i) = g*ang2;
    fij(i)    = y*(ANGLE(FR)-ANGLE(TO)) + g*ang2/2;
    fji(i)    = y*(ANGLE(TO)-ANGLE(FR)) + g*ang2/2;
    Bi(i)     = FLIM(i) - perdas(i);  % Termo independente 
    
%==============================================================          
    if ( fij(i) >= fji(i) )  % lado i-j positivo
%==============================================================
                                      
%-Contribui�ao dos Circuitos Violados em Ai, Bi

        Ai(i,NCAX)   = 0;            % Criando a linha com zero 
        Ai(i,NGER+FR)= + Bor(i);     % Coeficiente de teta-i
        Ai(i,NGER+TO)= - Bor(i);     % Coeficiente de teta j
   
    end
%==============================================================
    if ( fji(i) >= fij(i) )  % lado j-i positivo
%==============================================================                   
               
%-Contribui�ao dos Circuitos Violados em Ai, Bi

        Ai(i,NCAX)   = 0;            % Criando a linha com zero 
        Ai(i,NGER+FR)= - Bor(i);     % Coeficiente de teta-i
        Ai(i,NGER+TO)= + Bor(i);     % Coeficiente de teta j
    end
       
   %-----------modifica��es em BX devido as perdas-------------
    PINJ(FR)=PINJ(FR)+perdas(i)/2;
    PINJ(TO)=PINJ(TO)+perdas(i)/2;
    end
    BX(1:NBAR) = PLOAD(1:NBAR) + PINJ(1:NBAR);    % vetor de perdas em BX
   
   %-----------------------chamada do PL-----------------------

   [Vp,FOB, exitflag,output,Vd] = linprog(CX,Ai,Bi,AX,BX,Vlb,Vub);
%    if exitflag~=1
%         Ai=[];
%         Bi=[];
%         [Vp,FOB, exitflag,output,Vd] = linprog(CX,Ai,Bi,AX,BX,Vlb,Vub);
%         falhouLinprog = 1;
%         break
%    end

   %-----------------valor para teste de converg�ncia----------

   DIFMAX=0;
for i=1:NBAR
    DIFANG = abs(ANGLE(i)-Vp(NGER+i));
    if(DIFANG > DIFMAX)
        DIFMAX = DIFANG;
    end
end
ANGLE(1:NBAR) = Vp(NGER+1:NGER+NBAR);
ERRO(NITER) = DIFMAX;
ITER(NITER) = NITER;

end   % end do while

toc;

lbd_BUS = - Vd.eqlin(1:NBAR);
lbd_CV  = - Vd.ineqlin;
pilo    = - Vd.lower;
piup    = - Vd.upper;

% %----------------------Sa�da de dados-------------------------
% 
% %---------------Calculo dos valores das vari�veis-------------
% for I=1:NGER            % Geradores
%    DESLOC = I;
%    PG(I)  = PB * Vp(DESLOC);
% end
% for I=1:NBAR            % �ngulos
%    DESLOC = NGER+I;
%    ANGLE(I)  = Vp(DESLOC);
%    ANG_G(I)  = ANGLE(I)*180/pi;
% end
% 
% for i=1:NLIN  % fluxo nas linhas------------------
%     rl=r(i);
%     denom=(r(i)^2+x(i)^2);
%     g=rl/denom;
%     FR=SB(i);                    % de
%     TO=EB(i);                    % para
%     ang=ANGLE(FR)-ANGLE(TO);
%     ang2=ang*ang;  
%     y = Bor(i);
%     fij(i)= y*(ANGLE(FR)-ANGLE(TO)) + g*ang2/2;
%     fji(i)= y*(ANGLE(TO)-ANGLE(FR)) + g*ang2/2;
% end
% PERTOT = PB*sum(perdas);
% for i=1:NLIN
%    LINHA(i)=i;
% end
% for i=1:NBAR
%    BARRA(i)=i;
% end
% FOB = FOB*PB;
% disp('  ');
% disp('Fun��o objetivo: FOB = ');disp([FOB]);
% 
% disp('  ');
% disp('Itera��es = ');disp([NITER]);
% 
% disp('Perda total de pot�ncia do sistema:');disp([PERTOT]);
% disp('  ');
%    
% disp('----------DADOS DOS GERADORES-----------');
% disp('  ');
% disp('    Barra  Capacidade  Gera��o       cost  ');
% disp([BARPG',PGMAX'*PB,PG',CPG']);
% disp('  ');
% disp('------------DADOS DAS LINHAS------------');
% disp('  ');
% disp('     Linha       De      Para     fij  fji    Limite');
% disp([LINHA',SB',EB',fij'*PB,fji'*PB,FLIM'*100]);
% disp('  ');
% 
% disp('Dados das barras');
% disp('    Barra     Teta     Lambda')
% disp([BARRA',ANG_G',lbd_BUS(1:NBAR)]);
% disp('Linha(s) com viola��o de fluxo: ');
% disp('   Linviol    Lambda');
% disp([LINHA',lbd_CV]);













