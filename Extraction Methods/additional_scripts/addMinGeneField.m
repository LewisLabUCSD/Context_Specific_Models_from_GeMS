function newModel = addMinGeneField(model, gene_expr, gene_id, th, rem_unk)
    %Find the genes that corresponded to "low exprssion" class and were
    %attempted to be removed. Genes < th and genes with no data if rem_unk = 1
    %were attempted to be removed
    %Returns model with field of minimized genes
    lowGenes = gene_id(gene_expr < th);
    minGenes = intersect(model.genes, lowGenes);
    if rem_unk == 1
        unkGenes = setdiff(model.genes,gene_id);
        minGenes = [minGenes;unkGenes];
    end
    newModel = model;
    newModel.minGenes = minGenes;
end