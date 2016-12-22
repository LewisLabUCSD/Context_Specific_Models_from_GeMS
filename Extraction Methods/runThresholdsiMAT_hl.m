function runThresholdsiMAT_hl(figName, bb, modName, hl)

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

if strcmp(hl,'H')
    blb = blb*1.5;
elseif strcmp(hl,'L')
    blb = blb*0.5;
else
    quit
end

if strcmp(figName,'U')
    %UNCONSTRAINED
    tol = 1e-6;
    core = {};
    runtime = 3600;
    if strcmp(bb,'B')
        model_u = changeRxnBounds(model_u, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model_u = addBiomassSinks(model_u);
        figName = [figName,'F'];
    end
    expressionCol = mapGeneToRxn(model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u);
    epsil = 1;
    figName = [figName,'_',hl,'_'];
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p10, ths.p10, 1, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.mean, ths.mean, 2, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p25, ths.p25, 3, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p50, ths.p50, 4, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.mean, ths.p10, 5, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p25, ths.p10, 6, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p50, ths.p10, 7, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p25, ths.mean, 8, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p50, ths.mean, 9, modName, tol, runtime)
    run_iMat(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p50, ths.p25, 10, modName, tol, runtime)
end

if strcmp(figName,'C')
    tol = 1e-8;
    %CONSTRAINED
    core = {};
    runtime = 7200;
    if strcmp(bb,'B')
        model_c = changeRxnBounds(model_c, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model_c = addBiomassSinks(model_c);
        figName = [figName,'F'];
    end
    if strcmp(bb,'H')
        model_c = changeRxnBounds(model_c, 'Biomass_reaction', 1e-3, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'H'];
    end
    expressionCol = mapGeneToRxn(model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c);
    epsil = 1e-6;
    figName = [figName,'_',hl,'_'];
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p10, ths.p10, 1, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.mean, ths.mean, 2, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p25, ths.p25, 3, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p50, ths.p50, 4, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.mean, ths.p10, 5, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p25, ths.p10, 6, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p50, ths.p10, 7, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p25, ths.mean, 8, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p50, ths.mean, 9, modName, tol, runtime)
    run_iMat(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p50, ths.p25, 10, modName, tol, runtime)
end

if strcmp(figName,'S')
    %SEMI-CONSTRAINED
    tol = 1e-6;
    core = {};
    runtime = 3600;
    if strcmp(bb,'B')
        model_s = changeRxnBounds(model_s, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model_s = addBiomassSinks(model_s);
        figName = [figName,'F'];
    end
    expressionCol = mapGeneToRxn(model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s);
    epsil = 1;
    figName = [figName,'_',hl,'_'];
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p10, ths.p10, 1, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.mean, ths.mean, 2, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p25, ths.p25, 3, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p50, ths.p50, 4, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.mean, ths.p10, 5, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p25, ths.p10, 6, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p50, ths.p10, 7, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p25, ths.mean, 8, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p50, ths.mean, 9, modName, tol, runtime)
    run_iMat(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p50, ths.p25, 10, modName, tol, runtime)
end
exit;
end

function run_iMat(core, model, gene_exp, gene_names, expressionCol, figName, epsil, lt, ht, id, modName, tol, runtime)
    tName = ['iMAT_',figName, num2str(id),'_',modName];
    disp(tName)
    cMod = call_iMAT(model, core, epsil, expressionCol, lt, ht, tol, [tName,'.txt'], runtime);
    cMod.name = tName;
    eval([tName,'= addMinGeneField(cMod,gene_exp,gene_names,lt,0);']);
    save([tName,'.mat'],tName)
end