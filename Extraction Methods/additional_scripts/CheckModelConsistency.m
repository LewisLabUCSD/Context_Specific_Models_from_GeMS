function inactiveRxns = CheckModelConsistency(model, tol)
    % Consistency check using fastcc
    % Input:
    %   model - COBRA model struct
    %   tol - tolerance by which reactions are defined inactive after model extraction
    %         (recommended lowest value 1e-8 since solver tolerance is 1e-9)
    % Output:
    %   inactiveRxns - cell array with the names of inconsistent reactions
    
    is_active = fastcccheck(model, tol);
    inactiveRxns = setdiff(model.rxns, model.rxns(is_active));
end
