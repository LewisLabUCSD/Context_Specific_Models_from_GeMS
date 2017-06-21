function modelNew=removeIsoGenes(model)
    % Gets rid of the Entrez suffices ".x" in the gene identifiers in the
    % entire model, also updates GPR rules.
    % Input: model with genes in the form "entrezID.x"
    % Output: model where all genes are in the form "enrezID"
    
    genes={};
    for i = 1:numel(model.genes)
        str = model.genes{i};
        genes{i} = str(1:end-2);
    end
    genes = unique(genes);
        
    model.rules=[];
    model.rxnGeneMat=[];
    model.genes=[];
    warning off all
    for m=1:numel(model.rxns)
        grRule = model.grRules{m,1};
        grRule = regexprep(grRule,'[.]\d*','');
        model = changeGeneAssociation(model,model.rxns{m,1},grRule,genes,genes);
    end
    warning on all
    modelNew=model;
end