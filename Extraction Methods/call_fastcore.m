function tissueModel = call_fastcore(model, expressionCol, core, th, tol, scaling) 
    %Run FASTCORE with specified parameters
    %  model - consistent COBRA model struct (input model)
    %  expressionCol - reaction expression level. Must be same length as model.rxns
    %  core - cell with reaction names (strings) that are manually put in
    %         the core
    %  th - reactions with expression above this threshold put in the core
    %       (expression threshold)
    %  epsil - trade-off between removing medium confidence and
    %          low confidence reactions (usually 0.5)
    %  tol - tolerance by which reactions are defined inactive after model extraction
    %        (recommended lowest value 1e-8 since solver tolerance is 1e-9)
    %  scaling - scaling constant

    %This script is a modified version of the original, changes are
    %indicated with "FIX" as comment
    
    % (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
    %     LCSB / LSRU, University of Luxembourg
    model_orig = model;
    
    C = find(expressionCol >= th);
    C = union(C, find(ismember(model.rxns, core)));

    tic

    N = 1:numel(model.rxns);
    I = find(model.rev==0);

    A = [];
    flipped = false;
    singleton = false;  

    % start with I
    J = intersect( C, I ); %fprintf('|J|=%d  ', length(J));
    P = setdiff( N, C);
    Supp = findSparseMode( J, P, singleton, model, tol, scaling );
    if ~isempty( setdiff( J, Supp ) ) 
      fprintf ('Error: Inconsistent irreversible core reactions.\n');
      return;
    end
    A = Supp;  %fprintf('|A|=%d\n', length(A));
    J = setdiff( C, A ); %fprintf('|J|=%d  ', length(J));

    % main loop     
    while ~isempty( J )
        P = setdiff( P, A);
        Supp = findSparseMode( J, P, singleton, model, tol, scaling );
        A = union( A, Supp );   %fprintf('|A|=%d\n', length(A)); 
        if ~isempty( intersect( J, A ))
            J = setdiff( J, A );     %fprintf('|J|=%d  ', length(J));
            flipped = false;
        else
            if singleton
                JiRev = setdiff(J(1),I);
            else
                JiRev = setdiff(J,I);
            end
            if flipped || isempty( JiRev )
                if singleton
                    fprintf('\nError: Global network is not consistent.\n');
                    return
                else
                  flipped = false;
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
    %fprintf('|A|=%d\n', length(A));
    
    toRem = setdiff(model.rxns,model.rxns(A));
    tissueModel = removeRxns(model_orig, toRem);
    tissueModel = removeNonUsedGenes(tissueModel);
    toc
end

function Supp = findSparseMode( J, P, singleton, model, tol, scaling )
    %
    % Supp = findSparseMode( J, P, singleton, model, tol )
    %
    % Finds a mode that contains as many reactions from J and as few from P
    % Returns its support, or [] if no reaction from J can get flux above tol

    % (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
    %     LCSB / LSRU, University of Luxembourg

    Supp = [];
    if isempty( J ) 
      return;
    end

    if singleton
      V = LP7( J(1), model, tol );
    else
      V = LP7( J, model, tol );
    end

    K = intersect( J, find(V >= 0.99*tol) );   

    if isempty( K ) 
      return;
    end

    V = LP9( K, P, model, tol, scaling );

    Supp = find( abs(V) >= 0.99*tol );
end

function V = LP7( J, model, tol )
    % V = LP7( J, model, tol )
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
    ub = [model.ub; ones(nj,1)*tol];

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

function V = LP9( K, P, model, tol, scaling )
    % V = LP9( K, P, model, tol )
    %
    % CPLEX implementation of LP-9 for input sets K, P (see FASTCORE paper)

    % (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
    %     LCSB / LSRU, University of Luxembourg


    scalingfactor = scaling; %FIX: used to be 1e5, but since tol is smaller, different value

    V = [];
    if isempty(P) || isempty(K)
        return;
    end

    np = numel(P);
    nk = numel(K);
    [m,n] = size(model.S);

    % x = [v;z]

    % objective
    f = [zeros(1,n), ones(1,np)];

    % equalities
    Aeq = [model.S, sparse(m,np)];
    beq = zeros(m,1);

    % inequalities
    Ip = sparse(np,n); Ip(sub2ind(size(Ip),(1:np)',P(:))) = 1;
    Ik = sparse(nk,n); Ik(sub2ind(size(Ik),(1:nk)',K(:))) = 1;
    Aineq = sparse([[Ip, -speye(np)]; ...
                    [-Ip, -speye(np)]; ...
                    [-Ik, sparse(nk,np)]]);
    bineq = [zeros(2*np,1); -ones(nk,1)*tol*scalingfactor];

    % bounds
    lb = [model.lb; zeros(np,1)] * scalingfactor;
    ub = [model.ub; max(abs(model.ub(P)),abs(model.lb(P)))] * scalingfactor;

    % Set up problem
    LPproblem.A = [Aeq;Aineq];
    LPproblem.b = [beq;bineq];
    LPproblem.c = f;
    LPproblem.lb = lb;
    LPproblem.ub = ub;
    LPproblem.osense = 1;
    LPproblem.csense(1:m,1) = 'E';
    LPproblem.csense(m+1:length(bineq)+m,1) = 'L';
    
    %FIX: use gurobi as solver
    sol = solveCobraLP(LPproblem);
    if sol.stat == 1
        x = sol.full;    
        %x = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub);
        V = x(1:n);
    else
        V = zeros(n,1);
    end
end


