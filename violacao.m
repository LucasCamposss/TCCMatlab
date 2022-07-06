% MAIORES VIOLAÇÕES DE FLUXO
function f = violacao(maiorFluxo, limiteLinha)
    violacaoPercentual = zeros(max(size(maiorFluxo)),3);
    for i=1:max(size(maiorFluxo))
        violacaoPercentual(i,2) = i;
        if maiorFluxo(i)-limiteLinha(i)>0
            violacaoPercentual(i,1) = (maiorFluxo(i)-limiteLinha(i))/limiteLinha(i);
            violacaoPercentual(i,3) = (maiorFluxo(i)-limiteLinha(i));
        end
    end
    f = violacaoPercentual;
end