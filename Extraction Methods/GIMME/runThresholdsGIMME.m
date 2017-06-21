function runThresholdsGIMME(figName, modName)
    initCobraToolbox
    load(['parsedGPR_u_',modName,'.mat'])
    load(['model_u_',modName,'.mat'])
    load(['gene_expr_u_',modName,'.mat'])
    load(['gene_id_u_',modName,'.mat'])
    load(['corrRxn_u_',modName,'.mat'])
    load(['parsedGPR_c_',modName,'.mat'])
    load(['model_c_',modName,'.mat'])
    load(['gene_expr_c_',modName,'.mat'])
    load(['gene_id_c_',modName,'.mat'])
    load(['corrRxn_c_',modName,'.mat'])
    load(['ths_',modName,'.mat'])
    load(['model_s_',modName,'.mat'])
    load(['gene_expr_s_',modName,'.mat'])
    load(['gene_id_s_',modName,'.mat'])
    load(['corrRxn_s_',modName,'.mat'])
    load(['parsedGPR_s_',modName,'.mat'])
    load(['growthRate_',modName,'.mat'])

    
    if strcmp(figName,'U')
        tol = 1e-6;
        figName = [figName,'B'];
        disp('UNCONSTRAINED MODEL')
        expressionCol = mapGeneToRxn(model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u);
        singleRun(model_u, gene_expr_u, gene_id_u, expressionCol, ths.p10, 0.9, figName, 1, modName, tol)
        singleRun(model_u, gene_expr_u, gene_id_u, expressionCol, ths.mean, 0.9, figName, 2, modName, tol)
        singleRun(model_u, gene_expr_u, gene_id_u, expressionCol, ths.p25, 0.9, figName, 3, modName, tol)
        singleRun(model_u, gene_expr_u, gene_id_u, expressionCol, ths.p50, 0.9, figName, 4, modName, tol)
    end
    if strcmp(figName,'C')
        tol = 1e-8;
        disp('CONSTRAINED MODEL')
        figName = [figName,'B'];
        expressionCol = mapGeneToRxn(model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c);
        singleRun(model_c, gene_expr_c, gene_id_c, expressionCol, ths.p10, 0.9, figName, 1, modName, tol)
        singleRun(model_c, gene_expr_c, gene_id_c, expressionCol, ths.mean, 0.9, figName, 2, modName, tol)
        singleRun(model_c, gene_expr_c, gene_id_c, expressionCol, ths.p25, 0.9, figName, 3, modName, tol)
        singleRun(model_c, gene_expr_c, gene_id_c, expressionCol, ths.p50, 0.9, figName, 4, modName, tol)
    end
    if strcmp(figName,'S')
        tol = 1e-6;
        disp('SEMI-CONSTRAINED MODEL')
        figName = [figName,'B'];
        expressionCol = mapGeneToRxn(model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s);
        singleRun(model_s, gene_expr_s, gene_id_s, expressionCol, ths.p10, 0.9, figName, 1, modName, tol)
        singleRun(model_s, gene_expr_s, gene_id_s, expressionCol, ths.mean, 0.9, figName, 2, modName, tol)
        singleRun(model_s, gene_expr_s, gene_id_s, expressionCol, ths.p25, 0.9, figName, 3, modName, tol)
        singleRun(model_s, gene_expr_s, gene_id_s, expressionCol, ths.p50, 0.9, figName, 4, modName, tol)
    end
end

function singleRun(model, gene_exp, gene_names, expressionCol, ut, obj_frac, figName, id, modName, tol)
    tName = ['GIMME_',figName,num2str(id),'_',modName];
    disp(tName)
    disp('RUNNING GIMME...')
    cMod= call_GIMME(model, expressionCol, ut, obj_frac, tol);
    cMod.name = tName;
    eval([tName,'= addMinGeneField(cMod,gene_exp,gene_names,ut,0)']);
    save([tName,'.mat'],tName)
end