function runThresholdsMBA(figName, bb, modName, ext)
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

ext = num2str(ext);

if strcmp(figName,'U')
    %UNCONSTRAINED
    tol = 1e-6;
    core = {};
    if strcmp(bb, 'B')
        model_u = changeRxnBounds(model_u, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_u = addBiomassSinks(model_u);
        figName = [figName,'F'];
    end
    expressionCol = mapGeneToRxn(model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u);
    epsil = 0.5;
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p10, ths.p10, 1, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.mean, ths.mean, 2, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p25, ths.p25, 3, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p50, ths.p50, 4, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.mean, ths.p10, 5, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p25, ths.p10, 6, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p50, ths.p10, 7, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p25, ths.mean, 8, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p50, ths.mean, 9, modName, ext, tol)
    run_MBA(core, model_u, gene_expr_u, gene_id_u, expressionCol, figName, epsil, ths.p50, ths.p25, 10, modName, ext, tol)
end

if strcmp(figName,'C')
    %CONSTRAINED
    tol = 1e-8;
    core = {};
    if strcmp(bb, 'B')
        model_c = changeRxnBounds(model_c, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_c = addBiomassSinks(model_c);
        figName = [figName,'F'];
    end
    if strcmp(bb, 'H')
        model_c = changeRxnBounds(model_c, 'Biomass_reaction', 1e-3, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'H'];
    end
    epsil = 0.5;
    expressionCol = mapGeneToRxn(model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c);
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p10, ths.p10, 1, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.mean, ths.mean, 2, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p25, ths.p25, 3, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p50, ths.p50, 4, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.mean, ths.p10, 5, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p25, ths.p10, 6, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p50, ths.p10, 7, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p25, ths.mean, 8, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p50, ths.mean, 9, modName, ext, tol)
    run_MBA(core, model_c, gene_expr_c, gene_id_c, expressionCol, figName, epsil, ths.p50, ths.p25, 10, modName, ext, tol)
end
if strcmp(figName,'S')
    %SEMI-CONSTRAINED
    tol = 1e-6;
    core = {};
    if strcmp(bb, 'B')
        model_s = changeRxnBounds(model_s, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        core = {'Biomass_reaction','DM_atp(c)'};
        figName = [figName,'B'];
    end
    if strcmp(bb, 'F')
        model_s = addBiomassSinks(model_s);
        figName = [figName,'F'];
    end
    epsil = 0.5;
    expressionCol = mapGeneToRxn(model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s);
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p10, ths.p10, 1, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.mean, ths.mean, 2, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p25, ths.p25, 3, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p50, ths.p50, 4, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.mean, ths.p10, 5, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p25, ths.p10, 6, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p50, ths.p10, 7, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p25, ths.mean, 8, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p50, ths.mean, 9, modName, ext, tol)
    run_MBA(core, model_s, gene_expr_s, gene_id_s, expressionCol, figName, epsil, ths.p50, ths.p25, 10, modName, ext, tol)
end
exit;
end

function run_MBA(core, model, gene_exp, gene_names, expressionCol, figName, epsil, mt, ut, id, modName, ext, tol)
    tName = ['MBA',ext,'_',figName, num2str(id),'_',modName];
    disp(tName)
    %Make sure output model is consistent, if not, run again
    incon = true;
    nrun = 0;
    while incon && nrun < 10
        cMod = call_MBA(model, expressionCol, core, mt, ut, epsil, tol);
        nrun = nrun + 1;
        inactiveRxns = CheckModelConsistency(cMod, tol);
        if isempty(inactiveRxns)
            incon = false;
        end
    end
    if incon
        errmsg = 'Output model is inconsistent';
        save(['INC_',tName,'.mat'],'errmsg')
    end
    disp(['Number of rxns: ',num2str(numel(cMod.rxns))])
    cMod.name = tName;
    eval([tName,'= addMinGeneField(cMod, gene_exp, gene_names, mt, 1);']);
    save([tName,'.mat'],tName)
    disp(' ')
end