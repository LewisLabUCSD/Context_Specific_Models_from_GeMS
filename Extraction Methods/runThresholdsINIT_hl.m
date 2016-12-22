function runThresholdsINIT_hl(figName, bb, modName, hl)
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
    runtime = 3600;
    if strcmp(bb,'B')
        model = changeRxnBounds(model_u, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_u);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    expressionCol = mapGeneToRxn(model, gene_id_u, gene_expr_u, parsedGPR_u, corrRxn_u);
end

if strcmp(figName,'C')
    tol = 1e-8;
    runtime = 7200;
    if strcmp(bb,'B')
        model = changeRxnBounds(model_c, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_c);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    if strcmp(bb,'H')
        model = changeRxnBounds(model_c, 'Biomass_reaction', 1e-3, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'H'];
    end
    expressionCol = mapGeneToRxn(model, gene_id_c, gene_expr_c, parsedGPR_c, corrRxn_c);
end

if strcmp(figName,'S')
    tol = 1e-6;
    runtime = 3600;
    if strcmp(bb,'B')
        model = changeRxnBounds(model_s, 'Biomass_reaction', blb, 'l'); %Force biomass and ATP demand to be active
        figName = [figName,'B'];
    end
    if strcmp(bb,'F')
        model = addBiomassSinks(model_s);
        bmsInd=~cellfun(@isempty,strfind(model.rxns,'BMS_'));
        figName = [figName,'F'];
    end
    expressionCol = mapGeneToRxn(model, gene_id_s, gene_expr_s, parsedGPR_s, corrRxn_s);
end


%w1
w1 = zeros(length(expressionCol),1);
w1(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p10);
w1(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w1(w1 < -max(w1)) = -max(w1);
if strcmp(bb,'B') || strcmp(bb,'H')
    w1(ismember(model.rxns, {'Biomass_reaction','DM_atp(c)'})) = max(w1); %Biomass and ATP demand get high weight
else
    w1(bmsInd) = 0;
end

%w2
w2 = zeros(length(expressionCol),1);
w2(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.mean);
w2(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w2(w2 < -max(w2)) = -max(w2);
if strcmp(bb,'B') || strcmp(bb,'H')
    w2(ismember(model.rxns, {'Biomass_reaction','DM_atp(c)'})) = max(w2); %Biomass and ATP demand get high weight
end

%w3
w3 = zeros(length(expressionCol),1);
w3(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p25);
w3(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w3(w3 < -max(w3)) = -max(w3);
if strcmp(bb,'B') || strcmp(bb,'H')
    w3(ismember(model.rxns, {'Biomass_reaction','DM_atp(c)'})) = max(w3); %Biomass and ATP demand get high weight
else
    w3(bmsInd) = 0;
end

%w4
w4 = zeros(length(expressionCol),1);
w4(expressionCol >= 0) = 5*log(expressionCol(expressionCol>=0)/ths.p50);
w4(expressionCol < 0) = -2; % "unknown" entries get a weight of -2
w4(w4 < -max(w4)) = -max(w4);
if strcmp(bb,'B') || strcmp(bb,'H')
    w4(ismember(model.rxns, {'Biomass_reaction','DM_atp(c)'})) = max(w4); %Biomass and ATP demand get high weight
else
    w4(bmsInd) = 0;
end
    
if strcmp(figName,'UB') || strcmp(figName,'UF')
    %UNCONSTRAINED
    figName = [figName,'_',hl,'_'];
    epsil = 1;
    run_INIT(model, gene_expr_u, gene_id_u, ths.p25, figName, epsil, w1, 1, modName, tol, runtime);
    run_INIT(model, gene_expr_u, gene_id_u, ths.p50, figName, epsil, w2, 2, modName, tol, runtime);
    run_INIT(model, gene_expr_u, gene_id_u, ths.p10, figName, epsil, w3, 3, modName, tol, runtime);
    run_INIT(model, gene_expr_u, gene_id_u, ths.mean, figName, epsil, w4, 4, modName, tol, runtime);
end

if strcmp(figName,'CB') || strcmp(figName,'CF') || strcmp(figName,'CH')
    %CONSTRAINED
    figName = [figName,'_',hl,'_'];
    epsil = 1e-6;
    run_INIT(model, gene_expr_c, gene_id_c, ths.p25, figName, epsil, w1, 1, modName, tol, runtime);
    run_INIT(model, gene_expr_c, gene_id_c, ths.p50, figName, epsil, w2, 2, modName, tol, runtime);
    run_INIT(model, gene_expr_c, gene_id_c, ths.p10, figName, epsil, w3, 3, modName, tol, runtime);
    run_INIT(model, gene_expr_c, gene_id_c, ths.mean, figName, epsil, w4, 4, modName, tol, runtime);
end
if strcmp(figName,'SB') || strcmp(figName,'SF')
    %SEMI-CONSTRAINED
    figName = [figName,'_',hl,'_'];
    epsil = 1;
    run_INIT(model, gene_expr_s, gene_id_s, ths.p25, figName, epsil, w1, 1, modName, tol, runtime);
    run_INIT(model, gene_expr_s, gene_id_s, ths.p50, figName, epsil, w2, 2, modName, tol, runtime);
    run_INIT(model, gene_expr_s, gene_id_s, ths.p10, figName, epsil, w3, 3, modName, tol, runtime);
    run_INIT(model, gene_expr_s, gene_id_s, ths.mean, figName, epsil, w4, 4, modName, tol, runtime);

end
exit;
end

function run_INIT(model, gene_exp, gene_names, th, figName, epsil, w, id, modName, tol, runtime)
    tName = ['INIT_',figName, num2str(id),'_',modName];
    disp(tName)
    cMod = call_INIT(model, epsil, w, tol, [tName,'.txt'], runtime);
    cMod.name = tName;
    eval([tName,'= addMinGeneField(cMod, gene_exp, gene_names, th, 1);']);
    save([tName,'.mat'],tName)
end