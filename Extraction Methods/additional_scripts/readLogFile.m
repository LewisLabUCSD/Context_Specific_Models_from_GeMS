function [relGap, time] = readLogFile(fname)
    % Read a log file from gurobi MILP to show convergeance
    % Input:
    %   fname - name of log file (string)
    % Output:
    %   relGap - relative gap (objevtive and incumbent solution)
    %   time - time corresponding to relGap
    
    fid = fopen(fname);
    tline = fgetl(fid);
    cont = true;
    started = false;
    nempty = 0;
    relGap = [];
    time = [];
    while ischar(tline) && cont
        tline = fgetl(fid);
        if started && isempty(tline)
            nempty = nempty + 1;
        end
        if nempty > 1
            cont = false;
        end
        if started && cont && ~isempty(tline)
            ls = strsplit(tline, ' ');
            numPerc = ls{end-2};
            numTime = ls{end};
            if numel (numPerc) > 1
                perc = str2num(numPerc(1:end-1));
                tim = str2num(numTime(1:end-1));
                relGap(end+1) = perc;
                time(end+1) = tim;
            end
        end
        if ~started && ~isempty(regexp(tline,'Expl Unexpl', 'once'))
            started = true;
        end
    end
    fclose(fid);
end
