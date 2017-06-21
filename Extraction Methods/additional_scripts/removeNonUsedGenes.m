function modelNew=removeNonUsedGenes(model)
    % Remove genes from a model that do not longer correspond to any 
    % reaction. Also updates GPR rules.
    
    genes=unique(model.genes);
    if length(model.genes) ~= length(genes)
        disp('Some genes have identical IDs')
    end
    model.rules=[];
    model.rxnGeneMat=[];
    model.genes=[];
    warning off all
    for m=1:numel(model.rxns)
        model = changeGeneAssociation(model,model.rxns{m,1},model.grRules{m,1},genes,genes);
    end
    warning on all
    modelNew=model;
end