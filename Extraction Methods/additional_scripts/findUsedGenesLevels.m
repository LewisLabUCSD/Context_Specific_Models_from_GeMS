function [gene_id, gene_expr] = findUsedGenesLevels(model,num)
    % Find gene expression levels from the measured data (mRNA-Seq)
    % Input:
    %   model - COBRA model struct (trimmed gene suffices, strings)
    %   num - Nx2 matrix where the first column are the Entrez gene IDs
    %         (integers) and the second column are the corresponding 
    %          expression levels (from the mRNA-seq experiment).

    genes_ID = zeros(length(model.genes),1);
    for i = 1:length(model.genes)
        genes_ID(i) = str2num(model.genes{i});
    end
    cnts = -1*ones(length(genes_ID),2);
    cnts(:,1) = genes_ID;

    for i = 1:length(genes_ID)
        cur_ID = genes_ID(i);
        flag = 0;
        for ii = 1:length(num)
            if num(ii,1) == cur_ID
                % In case multiple expression levels for genes are found in
                % the data, add values (happens only rarely)
                if flag == 1
                    cnts(i,2) = cnts(i,2) + num(ii,2);
                end
                if flag == 0
                    flag = 1;
                    cnts(i,2) = num(ii,2);
                end
            end
        end
    end

    data_inds = find(cnts(:,2)~= -1);
    gene_expr = cnts(data_inds,2);
    gene_id = model.genes(data_inds);
end

