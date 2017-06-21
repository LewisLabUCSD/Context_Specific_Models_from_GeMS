function A = fastcccheck( model, tol ) 
    % Fast consistency check, using Gurobi5 as solver
    % Input:
    %   model - COBRA model struct
    %   tol - tolerance by which reactions are defined inactive after model extraction
    %         (recommended lowest value 1e-8 since solver tolerance is 1e-9)
    % Output:
    %   A - vector of indices of consistent reactions in the model

    % Implementation from:
    % (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
    %     LCSB / LSRU, University of Luxembourg
    
    N = (1:numel(model.rxns));
    I = find(model.rev==0);

    A = [];

    % start with I
    J = intersect( N, I ); %fprintf('|J|=%d  ', numel(J));
    V = LP7( J, model, tol ); 
    Supp = find( abs(V) >= 0.99*tol );  
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
            V = LP7( Ji, model, tol ) ; 
        end    
        Supp = find( abs(V) >= 0.99*tol );  
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