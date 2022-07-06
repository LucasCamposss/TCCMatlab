% =======================================================================
% Objective Function
% =======================================================================
function f = F_OBJ(x)

global NBUS NGER COSTPG OBJF PGDC DADOSDMR CPG

%==========================================

V    = x(1:NBUS);
TETA = x(NBUS+1:2*NBUS);
PG(1,1:NGER)   = x(2*NBUS+1:2*NBUS+NGER);
QG   = x(2*NBUS+NGER+1:2*NBUS+2*NGER);


if OBJF==1 
    f = sum(COSTPG.*PG); % minimal Cost
end


if OBJF==2
%     f = sum(100*((DADOSDMR(2*NBUS+1:2*NBUS+NGER)-PG).^2));
%     f = sum(COSTPG.*PG)/10000; % minimal losses
%     f = f + 100000*(aux);  % minimal generation deviation sum(CPG.*PG)+
    f = sum(CPG.*PG); % minimal Cost
end
if OBJF==3
    f = 100*PG*PG';  % minimal generation deviation
end
if OBJF==4
    Qinj(1:NBUS,1) = x(2*NBUS+2*NGER+1:2*NBUS+2*NGER+NBUS); % Referente a Qinj
    f = 10*Qinj'*Qinj;  % minimal reactive power suply
    aux = PGDC-PG;
    f = f + 100*(aux*aux');  % minimal generation deviation
end
if OBJF==6
    aux = (DADOSDMR(2*NBUS+1:2*NBUS+NGER)-PG).^2;
    f = f + 100*(sum(aux));  % minimal generation deviation
end


