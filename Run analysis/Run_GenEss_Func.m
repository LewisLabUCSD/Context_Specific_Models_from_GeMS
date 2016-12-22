
function [GenEss_score,Func_score] = Run_GenEss_Func(modelName,CellLine)
%Run_GenEss_Func Performs analysis for Gene essentiality and Metabolic
%capabilities
    %linearMOMA
    %
    % [GenEss_score,Func_score] = Run_GenEss_Func(model,CellLine,savefile)
    %
    %INPUT
    % modelName        COBRA model structure of the extracted context-specific
    %                  model
    % CellLine         Cell Line of the extracted context-specific model
    %                  (A375, HL60, KBM7 or K562)
    %
    %OUTPUTS
    % GenEss_score     Computed Gene Essentiality score
    % Func_score       Functionality score (pourcentage of capabilities
    %                  recovered by the context-specific model)


global  depl_p 
    
    load(['model_u_',CellLine,'.mat']);
    load(['gene_expr_u_',CellLine,'.mat']);
    load(['gene_id_u_',CellLine,'.mat']);
    depl_p = findModelDepletionGenes_loc(model_u,['depletion_ratios_',CellLine,'.xls']);


        %% Load the context-specific model to analyze
        load([modelName,'.mat']);
        tmB = eval(modelName);
        eval(['clear ',modelName]);
  
        %% Check gene-essentiality
        maxgr=0.01;
        tmBcRed = tmB;
        tmBcRed = changeRxnBounds_loc(tmBcRed,'Biomass_reaction',0,'l');
        grRatios = singleGeneDeletion_loc(tmBcRed,'FBA',depl_p.genes);
    
        mEssP = depl_p.values(grRatios <= maxgr);
        mNonP = depl_p.values(grRatios > maxgr);
        
        GenEss_score = ranksum(mEssP,mNonP,'tail','left');
        [totTest, numTests] = testFunctionality(tmB);
        Func_score=totTest/numTests*100;
        
end

function [totTest, numTests] = testFunctionality(model)
    [~,bmMets,~] = xlsread('biomass_rxn.xls');
    bmMets = setdiff(bmMets,{'atp[c]','adp[c]','pi[c]','h2o[c]','h[c]'});
    numTests = numel(bmMets) + 1;
    comp = {'[c]','[e]','[g]','[l]','[m]','[n]','[r]','[x]'};
    totTest = 0;
    for i = 1:numel(bmMets)
        cMet = bmMets{i}(1:end-3);
        cTest = 0;
        for j = 1:numel(comp)
            met = [cMet,comp{j}];
            if ~isempty(intersect(model.mets,met))
                nm = addReaction(model,['BMS_',met],{met},-0.5,false,0,1000,0,'Biomass Metabolite Sink');
                nm = changeObjective(nm, ['BMS_',met]);
                sol = optimizeCbModel(nm);
                if sol.f > 1e-8;
                    cTest = 1;
                    break;
                end
            end
        end
        totTest = cTest + totTest;
    end
    
    warning off all
    nm = addReaction(model,'BMS_atp(c)',{'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'},[-0.5,-0.5,0.5,0.5,0.5],false,0,1000,0,'Biomass Metabolite Sink'); %Make a copy of ATP demand
    nm = changeRxnBounds(nm, 'DM_atp_c_',0,'l');
    warning on all
    nm = changeObjective(nm, 'BMS_atp(c)');
    sol = optimizeCbModel(nm);
    if sol.f > 1e-8;
        totTest = totTest + 1;
    end
    
    disp(['Functionality test: ',num2str(totTest),'/',num2str(numel(bmMets)+1)])
end

%%
function [grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion_loc(model,method,geneList,verbFlag)
%singleGeneDeletion Performs single gene deletion analysis using FBA, MOMA or
    %linearMOMA
    %
    % [grRatio,grRateKO,grRateWT,delRxns,hasEffect] = singleGeneDeletion(model,method,geneList,verbFlag)
    %
    %INPUT
    % model         COBRA model structure including gene-reaction associations
    %
    %OPTIONAL INPUT
    % method        Either 'FBA', 'MOMA', or 'lMOMA' (Default = 'FBA')
    % geneList      List of genes to be deleted (default = all genes)
    % verbFlag      Verbose output (Default false)
    %
    %OUTPUTS
    % grRatio       Computed growth rate ratio between deletion strain and wild type
    % grRateKO      Deletion strain growth rates (1/h)
    % grRateWT      Wild type growth rate (1/h)
    % hasEffect     Does a gene deletion affect anything (i.e. are any reactions
    %               removed from the model)
    % delRxns       List of deleted reactions for each gene KO
    % fluxSolution  FBA/MOMA/lMOMA fluxes for KO strains
    %
    % Markus Herrgard 8/7/06

    if (nargin < 2)
        method = 'FBA';
    end
    if (nargin < 3)
        geneList = model.genes;
    else
        if (isempty(geneList))
            geneList = model.genes;
        end
    end
    if (nargin < 4)
        verbFlag = false;
    end

    nGenes = length(model.genes);
    nDelGenes = length(geneList);
    solWT = optimizeCbModel(model); 
    grRateWT = solWT.f;

    if grRateWT <= 1e-2 
        warning('Very small GR')
    end

    grRateKO = ones(nDelGenes,1)*grRateWT;
    grRatio = ones(nDelGenes,1); 
    hasEffect = true(nDelGenes,1);
    fluxSolution = zeros(length(model.rxns),nDelGenes);
    delRxns = cell(nDelGenes,1);
    if (verbFlag)  
        fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n','No','Perc','Name','Growth rate','Rel. GR');
    end

    for i = 1:nDelGenes
        if any(ismember(model.genes,geneList{i}))
            [modelDel,hasEffect(i),constrRxnNames] = deleteModelGenes(model,geneList{i}); 
            delRxns{i} = constrRxnNames;
        else
            delRxns{i} = {};
            hasEffect(i) = false;
        end
        if (hasEffect(i))
            switch method
                case 'lMOMA'
                    solKO = linearMOMA(model,modelDel,'max');
                case 'MOMA'
                    solKO = MOMA_loc(model,modelDel,'max',false,true);
                otherwise
                    solKO = optimizeCbModel(modelDel,'max');
            end
            if (solKO.stat == 1)
                grRateKO(i) = solKO.f;
                fluxSolution(:,i) = solKO.x;
            else
                grRateKO(i) = NaN;
            end
        end
        if (verbFlag)
            fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n',i,100*i/nDelGenes,geneList{i},grRateKO(i),grRateKO(i)/grRateWT*100);
        end
    end

    grRatio = grRateKO/grRateWT;
    grRatio(isnan(grRatio)) = 0; 

end

%%
function model = changeRxnBounds_loc(model,rxnNameList,value,boundType)
%changeRxnBounds Change upper or lower bounds of a reaction or a set of
%reactions
%
% model = changeRxnBounds(model,rxnNameList,value,boundType)
%
%INPUTS
% model         COBRA model structure
% rxnNameList   List of reactions (cell array or string)
% value         Bound values
%               Can either be a vector or a single scalar value if the same
%               bound value is to be assinged to all reactions
%
%OPTIONAL INPUT
% boundType     'u' - upper, 'l' - lower, 'b' - both (Default = 'b')
%               Bound type can either be a cell array of strings or a 
%               string with as many letters as there are reactions in 
%               rxnNameList
%
%OUTPUT
% model         COBRA model structure with modified reaction bounds
%
% Markus Herrgard 4/21/06

if (nargin < 4)
    boundType = 'b';
end

if ((length(value) ~= length(rxnNameList) && length(value) > 1) || (length(boundType) ~= length(rxnNameList) && length(boundType) > 1))
   error('Inconsistent lenghts of arguments: rxnNameList, value & boundType'); 
end

rxnID = findRxnIDs(model,rxnNameList);

% Remove reactions that are not in the model
if (iscell(rxnNameList))
    missingRxns = rxnNameList(rxnID == 0);
    for i = 1:length(missingRxns)
        fprintf('Reaction %s not in model\n',missingRxns{i}); 
    end
    if (length(boundType) > 1)
        boundType = boundType(rxnID ~= 0);
    end
    if (length(value) > 1)
        value = value(rxnID ~= 0);
    end
    rxnID = rxnID(rxnID ~= 0);    
end

if (isempty(rxnID) || sum(rxnID) == 0)
    display('No such reaction in model');
else
    nRxns = length(rxnID);
    if (length(boundType) > 1)
        if (length(value) == 1)
            value = repmat(value,nRxns,1);
        end
        for i = 1:nRxns
            switch lower(boundType{i})
                case 'u'
                    model.ub(rxnID(i)) = value(i);
                case 'l'
                    model.lb(rxnID) = value(i);
                case 'b'
                    model.lb(rxnID) = value(i);
                    model.ub(rxnID) = value(i);
            end
        end
    else
        switch lower(boundType)
            case 'u'
                model.ub(rxnID) = value;
            case 'l'
                model.lb(rxnID) = value;
            case 'b'
                model.lb(rxnID) = value;
                model.ub(rxnID) = value;
        end
    end
end
end

function genes_ratios = findModelDepletionGenes_loc(model,filename)

    % Load depletion (CRISPR) and ID data
    [gGeck,gEntr] = textread('gecko_id_to_entrez.txt','%s %s');        
    [~,~,raw] = xlsread(filename);
    gDepl = raw(:,1);
    for i = 1:length(gDepl)
        if isnumeric(gDepl{i})
            gDepl{i} = num2str(gDepl{i});
        end
    end
   
    dRatio = raw(:,2);    
    % Only get the genes which can be translated from gecko to entrez
    [~, entrInds] = ismember(gDepl,gGeck);
    gEntr = gEntr(entrInds(entrInds~=0));
    dRatio = cell2mat(dRatio(find(entrInds))); 
    
    % Only get the genes which we have depletion data for and are in the
    % model
    [members, geneInds] = ismember(model.genes, gEntr);
    mRatio = dRatio(geneInds(geneInds~=0));
    mGenes = model.genes(members);

    genes_ratios.values = mRatio;
    genes_ratios.genes = mGenes;
    
end