function runThresholdsmCADRE_hl(figName, bb, modName, hl) 
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
        tol = 1e-6;
        core = {};
        if strcmp(bb,'B')
            model_u = changeRxnBounds(model_u, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_u = addBiomassSinks(model_u);
            figName = [figName,'F'];
            bmi = find(ismember(model_u.rxns,'Biomass_reaction'));
            model_u.rxns{bmi} = 'Biomass_reaction_rename';
        end
        figName = [figName,'_',hl,'_'];
        disp('UNCONSTRAINED MODEL')
        singleRun(core, model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u, ths.p10, 1, 1/3, figName, 1, modName, tol)
        singleRun(core, model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u, ths.mean, 1, 1/3, figName, 2, modName, tol)
        singleRun(core, model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u, ths.p25, 1, 1/3, figName, 3, modName, tol)
        singleRun(core, model_u, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u, ths.p50, 1, 1/3, figName, 4, modName, tol)
    end
    if strcmp(figName,'C')
        tol = 1e-8;
        disp('CONSTRAINED MODEL')
        core = {};
        if strcmp(bb,'B')
            model_c = changeRxnBounds(model_c, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_c = addBiomassSinks(model_c);
            figName = [figName,'F'];
            bmi = find(ismember(model_c.rxns,'Biomass_reaction'));
            model_c.rxns{bmi} = 'Biomass_reaction_rename';
        end
        if strcmp(bb,'H')
            model_c = changeRxnBounds(model_c, 'Biomass_reaction', 1e-3, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'H'];
        end
        figName = [figName,'_',hl,'_'];
        singleRun(core, model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c, ths.p10, 1, 1/3, figName, 1, modName, tol)
        singleRun(core, model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c, ths.mean, 1, 1/3, figName, 2, modName, tol)
        singleRun(core, model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c, ths.p25, 1, 1/3, figName, 3, modName, tol)
        singleRun(core, model_c, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c, ths.p50, 1, 1/2, figName, 4, modName, tol)
    end
    if strcmp(figName,'S')
        tol = 1e-6;
        core = {};
        if strcmp(bb,'B')
            model_s = changeRxnBounds(model_s, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
            core = {'Biomass_reaction','DM_atp(c)'};
            figName = [figName,'B'];
        end
        if strcmp(bb,'F')
            model_s = addBiomassSinks(model_s);
            figName = [figName,'F'];
            bmi = find(ismember(model_s.rxns,'Biomass_reaction'));
            model_s.rxns{bmi} = 'Biomass_reaction_rename';
        end
        disp('SEMI-CONSTRAINED MODEL')
        figName = [figName,'_',hl,'_'];
        singleRun(core, model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s, ths.p10, 1, 1/3, figName, 1, modName, tol)
        singleRun(core, model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s, ths.mean, 1, 1/3, figName, 2, modName, tol)
        singleRun(core, model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s, ths.p25, 1, 1/3, figName, 3, modName, tol)
        singleRun(core, model_s, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s, ths.p50, 1, 1/3, figName, 4, modName, tol)
    end
    exit;
end

function singleRun(core, model, gene_names, gene_exp, parsedGPR, corrRxn, ht, mcheck, eta, figName, id, modName, tol)
    tName = ['mCADRE_',figName,num2str(id),'_',modName];
    disp(tName)
    disp('RUNNING mCADRE...')
    cMod=call_mCADRE(model, gene_names, gene_exp, parsedGPR, corrRxn, core, ht, mcheck, eta, tol);
    inactiveRxns = CheckModelConsistency(cMod, tol);
    %Make sure output model is consistent
    if ~isempty(inactiveRxns)
        errmsg = 'Output model is inconsistent';
        save(['INC_',tName,'.mat'],'errmsg')
    end
    cMod.name = tName;
    eval([tName,'= addMinGeneField(cMod, gene_exp, gene_names, ht, 1);']);
    save([tName,'.mat'],tName)
end