function new_model = addBiomassSinks(model)
    [~,bmMets,~] = xlsread('biomass_rxn.xls');
    bmMets = setdiff(bmMets,{'atp[c]','adp[c]','pi[c]','h2o[c]','h[c]'});
    
    new_model = addReaction(model,'BMS_atp(c)',{'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'},[-0.5,-0.5,0.5,0.5,0.5],false,0,1000,0,'Biomass Metabolite Sink'); %Make a copy of ATP demand
    for i = 1:numel(bmMets)
        new_model = addReaction(new_model,['BMS_',bmMets{i}],{bmMets{i}},-0.5,false,0,1000,0,'Biomass Metabolite Sink');
    end
end
