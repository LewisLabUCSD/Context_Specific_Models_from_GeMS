function tissueModel = call_MBA(model, expressionCol, core, mt, ut, epsil, tol)
    %Run MBA with specified parameters
    %Input:
    %  model - consistent COBRA model struct (input model)
    %  expressionCol - reaction expression level. Must be same length as model.rxns
    %  core - cell with reaction names (strings) that are manually put in
    %         the high confidence core
    %  mt - reactions with expression above this threshold are medium
    %       confidence (expression threshold)
    %  ut - reactions with expression above this threshold are high confidence
    %       (expression threshold)
    %  epsil - trade-off between removing medium confidence and
    %          low confidence reactions (usually 0.5)
    %  tol - tolerance by which reactions are defined inactive after model extraction
    %        (recommended lowest value 1e-8 since solver tolerance is 1e-9)


    % Get High expression core and medium expression core
    indH = find(expressionCol > ut);
    indM = find(expressionCol >= mt & expressionCol <= ut);
    CH = union(model.rxns(indH),core); %#ok<FNDSB>
    CM = model.rxns(indM); %#ok<FNDSB>
    
    NC = setdiff(model.rxns,union(CH,CM));
    
    %Biomass metabolite sinks do not have to be pruned
    nonBmsInd=cellfun(@isempty,strfind(NC, 'BMS_'));
    NC = NC(nonBmsInd);
    
    %MBA
    PM = model;
    removed = {};
    while ~isempty(NC)
        ri = randi([1,numel(NC)],1);
        r = NC{ri};
        inactive = CheckConsistency(PM, r, tol);
        eH = intersect(inactive, CH);
        eM = intersect(inactive, CM);
        eX = setdiff(inactive,union(CH,CM));
        if numel(eH)==0 && numel(eM) < epsil*numel(eX)
            PM = removeRxns(PM, inactive);
            NC = setdiff(NC,inactive);
            removed = union(removed,inactive);
        else
            NC = setdiff(NC,r);
        end
    end


    tissueModel = removeNonUsedGenes(PM);
    
    is_active = fastcc(tissueModel, tol);
    inactiveRxns = setdiff(tissueModel.rxns, tissueModel.rxns(is_active));
    if ~isempty(inactiveRxns)
        warning('Extracted model is not consistent, this might be caused by (numerical) issues in fastcc consistency checks')
    end
end

function inactiveRxns = CheckConsistency(model, r, epsilon)
    model = removeRxns(model, r);
    sol = optimizeCbModel(model);
    if sol.stat ~= 1
        inactiveRxns = model.rxns;
    else
        is_active = fastcc(model, epsilon);
        inactiveRxns = setdiff(model.rxns, model.rxns(is_active));
    end  
    inactiveRxns = [inactiveRxns;r];
end

function A = fastcc( model, epsilon ) 
    % A = fastcc( model, epsilon )
    %
    % The FASTCC algorithm for testing the consistency of an input model
    % Output A is the consistent part of the model

    % (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
    %     LCSB / LSRU, University of Luxembourg
    %tic

    N = (1:numel(model.rxns));
    I = find(model.rev==0);

    A = [];

    % start with I
    J = intersect( N, I ); %fprintf('|J|=%d  ', numel(J));
    V = LP7( J, model, epsilon ); 
    Supp = find( abs(V) >= 0.99*epsilon );  
    A = Supp;  %fprintf('|A|=%d\n', numel(A));
    incI = setdiff( J, A );    
    if ~isempty( incI )
        %fprintf('\n(inconsistent subset of I detected)\n');
    end
    J = setdiff( setdiff( N, A ), incI);  %fprintf('|J|=%d  ', numel(J));

    % reversible reactions
    flipped = false;
    singleton = false;        
    while ~isempty( J )
        if singleton
            Ji = J(1);
            V = LP3( Ji, model ) ; 
        else
            Ji = J;
            V = LP7( Ji, model, epsilon ) ; 
        end    
        Supp = find( abs(V) >= 0.99*epsilon );  
        A = union( A, Supp);  %fprintf('|A|=%d\n', numel(A)); 
        if ~isempty( intersect( J, A ))
            J = setdiff( J, A );     %fprintf('|J|=%d  ', numel(J));
            flipped = false;
        else
            JiRev = setdiff( Ji, I );
            if flipped || isempty( JiRev )
                flipped = false;
                if singleton
                    J = setdiff( J, Ji );  
                    %fprintf('\n(inconsistent reversible reaction detected)\n');
                else
                    singleton = true;
                end
            else
                model.S(:,JiRev) = -model.S(:,JiRev);
                tmp = model.ub(JiRev);
                model.ub(JiRev) = -model.lb(JiRev);
                model.lb(JiRev) = -tmp;
                flipped = true;  %fprintf('(flip)  ');
            end
        end
    end

    if numel(A) == numel(N)
        %fprintf('\nThe input model is consistent.\n'); 
    end
    %toc
end

function V = LP3( J, model )
    % V = LP3( J, model )
    %
    % CPLEX implementation of LP-3 for input set J (see FASTCORE paper)

    % (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
    %     LCSB / LSRU, University of Luxembourg


    [m,n] = size(model.S);

    % objective
    f = zeros(1,n);
    f(J) = -1;

    % equalities
    Aeq = model.S;
    beq = zeros(m,1);

    % bounds
    lb = model.lb;
    ub = model.ub;
    
    % Set up problem
    LPproblem.A = Aeq;
    LPproblem.b = beq;
    LPproblem.c = f;
    LPproblem.lb = lb;
    LPproblem.ub = ub;
    LPproblem.osense = 1;
    LPproblem.csense(1:m,1) = 'E';

    %V = cplexlp(f,[],[],Aeq,beq,lb,ub);
    sol = solveCobraLP(LPproblem);
    V = sol.full;
end

function V = LP7( J, model, epsilon )
    % V = LP7( J, model, epsilon )
    %
    % CPLEX implementation of LP-7 for input set J (see FASTCORE paper)

    % (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
    %     LCSB / LSRU, University of Luxembourg

    nj = numel(J);
    [m,n] = size(model.S);

    % x = [v;z]

    % objective
    f = -[zeros(1,n), ones(1,nj)];

    % equalities
    Aeq = [model.S, sparse(m,nj)];
    beq = zeros(m,1);

    % inequalities
    Ij = sparse(nj,n); 
    Ij(sub2ind(size(Ij),(1:nj)',J(:))) = -1;
    Aineq = sparse([Ij, speye(nj)]);
    bineq = zeros(nj,1);

    % bounds
    lb = [model.lb; zeros(nj,1)];
    ub = [model.ub; ones(nj,1)*epsilon];

    % Set up problem
    LPproblem.A = [Aeq;Aineq];
    LPproblem.b = [beq;bineq];
    LPproblem.c = f;
    LPproblem.lb = lb;
    LPproblem.ub = ub;
    LPproblem.osense = 1;
    LPproblem.csense(1:m,1) = 'E';
    LPproblem.csense(m+1:length(bineq)+m,1) = 'L';
    
    sol = solveCobraLP(LPproblem);
    if sol.stat == 1
        x = sol.full;    
        %x = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub);
        V = x(1:n);
    else
        V = zeros(n,1);
    end
end