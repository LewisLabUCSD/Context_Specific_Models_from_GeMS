function runThresholdsFastCore_hl(figName, bb, modName, hl)
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
    end
    
    if strcmp(figName,'U')
        %UNCONSTRAINED
        core = {};
        if strcmp(bb,'B')
            model_u = changeRxnBounds(model_u, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_u = addBiomassSinks(model_u);
            figName = [figName,'F'];
        end
        epsil = 1e-6;
        scaling = 1e3;
        figName = [figName,'_',hl,'_'];
        expressionCol = mapGeneToRxn(model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u);
        singleRun(core, gene_expr_u, gene_id_u, expressionCol, figName, model_u, epsil, scaling, ths.p10, 1, modName)
        singleRun(core, gene_expr_u, gene_id_u, expressionCol, figName, model_u, epsil, scaling, ths.mean, 2, modName)
        singleRun(core, gene_expr_u, gene_id_u, expressionCol, figName, model_u, epsil, scaling, ths.p25, 3, modName)
        singleRun(core, gene_expr_u, gene_id_u, expressionCol, figName, model_u, epsil, scaling, ths.p50, 4, modName)
    end

    if strcmp(figName,'C')
        %CONSTRAINED
        core = {};
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
        epsil = 1e-8;
        scaling = 1e3;
        figName = [figName,'_',hl,'_'];
        expressionCol = mapGeneToRxn(model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c);
        singleRun(core, gene_expr_c, gene_id_c, expressionCol, figName, model_c, epsil, scaling, ths.p10, 1, modName)
        singleRun(core, gene_expr_c, gene_id_c, expressionCol, figName, model_c, epsil, scaling, ths.mean, 2, modName)
        singleRun(core, gene_expr_c, gene_id_c, expressionCol, figName, model_c, epsil, scaling, ths.p25, 3, modName)
        singleRun(core, gene_expr_c, gene_id_c, expressionCol, figName, model_c, epsil, scaling, ths.p50, 4, modName)
    end
    
    if strcmp(figName,'S')
        %CONSTRAINED
        core = {};
        if strcmp(bb,'B')
            model_s = changeRxnBounds(model_s, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_s = addBiomassSinks(model_s);
            figName = [figName,'F'];
        end
        epsil = 1e-6;
        scaling = 1e3;
        figName = [figName,'_',hl,'_'];
        expressionCol = mapGeneToRxn(model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s);
        singleRun(core, gene_expr_s, gene_id_s, expressionCol, figName, model_s, epsil, scaling, ths.p10, 1, modName)
        singleRun(core, gene_expr_s, gene_id_s, expressionCol, figName, model_s, epsil, scaling, ths.mean, 2, modName)
        singleRun(core, gene_expr_s, gene_id_s, expressionCol, figName, model_s, epsil, scaling, ths.p25, 3, modName)
        singleRun(core, gene_expr_s, gene_id_s, expressionCol, figName, model_s, epsil, scaling, ths.p50, 4, modName)
    end

end

function singleRun(core, gene_exp, gene_names, expressionCol, figName, model, epsil, scaling, th, id, modName)
    tName = ['FastCore_',figName, num2str(id),'_',modName];
    disp(tName)
    cMod = call_fastcore(model, expressionCol, core, th, epsil, scaling);
    cMod.name = tName;
    eval([tName,'= addMinGeneField(cMod, gene_exp, gene_names, th, 1);']);
    save([tName,'.mat'],tName)
end