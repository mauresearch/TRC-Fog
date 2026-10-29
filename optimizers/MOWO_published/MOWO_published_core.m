%% Multi‑objective Walrus optimizer (MOWO)
function [x,f] = MOWO_published_core(params,MultiObj)
flage=1;
name = MultiObj.name;                 % Problem name
numOfObj = MultiObj.numOfObj;         % Set the number of objective function
evaluate_objective = MultiObj.fun;    % Evaluate the fitness value of the objective function
D = MultiObj.nVar;                    % Set the dimensional of problem at hand
LB = MultiObj.var_min;                % Lower boundary constraint
UB = MultiObj.var_max;                % Upper boundary constraint
Max_iteration = params.maxgen;        % Set the maximum number of generation (GEN)
SearchAgents_no = params.Np;          % Set the population size (Search Agent)
ishow = 1;
if isfield(params,'showPlot'), showPlot = params.showPlot; else, showPlot = true; end
if isfield(params,'verbose'),  verbose  = params.verbose;  else, verbose  = true; end
Nr = params.Nr;
if flage == 1
    chromosome = initialize_variables(SearchAgents_no, numOfObj, D, LB, UB,evaluate_objective);
end
intermediate_chromosome = non_domination_sort_mod(chromosome, numOfObj, D);
Pop = replace_chromosome(intermediate_chromosome, numOfObj,D,Nr);
M=numOfObj;
%% Algorithm Variables
K = D+M;
Walrus_pos = Pop(:,1:K+1);
Walrus_pos_ad = zeros(SearchAgents_no,K);

Iteration = 1;% Loop counter
while Iteration <= Max_iteration
    for i = 1:SearchAgents_no
        j = floor(rand()* SearchAgents_no) + 1;
        while j == i
            j = floor(rand()* SearchAgents_no) + 1;
        end
        %% 挑选出Best\Second个体
        Walrus_Best = Walrus_pos(1,:);
        Walrus_Second = Walrus_pos(2,:);
        %% WO
        pp = 0.4;                                     % Proportion of females
        F_number = round(SearchAgents_no*pp);         % Number of females
        M_number = F_number;                          % Number of males
        C_number = SearchAgents_no-F_number-M_number; % Number of children

        Alpha=1-Iteration/Max_iteration;
        Beta=1-1/(1+exp((1/2*Max_iteration-Iteration)/Max_iteration*10));
        A=2*Alpha; % A decreases linearly fron 2 to 0
        r1=rand();
        R=2*r1-1;
        Danger_signal=A*R;
        r2=rand();
        Satey_signal=r2;

        if abs(Danger_signal)>=1 % Migration (exploration)
            r3=rand();
            Rs=size(Walrus_pos,1);
            Migration_step=(Beta*r3^2)*(Walrus_pos(randperm(Rs),1:D)-Walrus_pos(randperm(Rs),1:D));
            Walrus_posNew=Walrus_pos(i,1:D)+Migration_step(i,1:D);

        elseif abs(Danger_signal)<1 % Reproduction (exploitation)
            if Satey_signal>=0.5 % Roosting behavior
                % [speed-up, identical result] The original loops below overwrite
                % Walrus_posNew in every pass, so only their LAST pass matters:
                %   for ii = 1:M_number ... for jj = M_number+1:M_number+F_number ... end, end
                %   for kk = SearchAgents_no-C_number+1:SearchAgents_no
                %       P=rand; o=Walrus_Best(:,1:D)+Walrus_pos(kk,1:D).*levyFlight(D);
                %       Walrus_posNew = P*(o-Walrus_pos(kk,1:D)); end
                % Here only the last pass is computed. The random numbers of the
                % earlier juvenile passes are still drawn in the same order, so the
                % random stream - and therefore every result - is exactly unchanged.
                if M_number >= 1
                    if F_number >= 1
                        ii = M_number; jj = M_number + F_number;
                        Walrus_posNew = Walrus_pos(jj,1:D)+Alpha.*(Walrus_pos(ii,1:D)-Walrus_pos(jj,1:D))+(1-Alpha).*(Walrus_Best(:,1:D)-Walrus_pos(jj,1:D));
                    else
                        Walrus_posNew = LB+hal(M_number,7).*(UB-LB);
                    end
                end
                kks = SearchAgents_no-C_number+1:SearchAgents_no;
                for q = 1:numel(kks)-1
                    rand; levyFlight(D);          % draws of the skipped passes
                end
                if ~isempty(kks)
                    kk = kks(end);
                    P=rand;
                    o=Walrus_Best(:,1:D)+Walrus_pos(kk,1:D).*levyFlight(D);
                    Walrus_posNew = P*(o-Walrus_pos(kk,1:D));
                end
            end
            % Foraging behavior
            if Satey_signal<0.5 && abs(Danger_signal)>=0.5 % Fleeing behavior
                r4=rand;
                Walrus_posNew=Walrus_pos(i,1:D)*R-abs(Walrus_Best(:,1:D)-Walrus_pos(i,1:D))*r4^2;
            end

            if Satey_signal<0.5 && abs(Danger_signal)<0.5 % Gathering behavior
                for jj=1:D
                    theta1=rand();
                    a1=Beta*rand()-Beta;
                    b1=tan(theta1.*pi);
                    X1=Walrus_Best(jj)-a1*b1*abs(Walrus_Best(jj)-Walrus_pos(i,jj));

                    theta2=rand();
                    a2=Beta*rand()-Beta;
                    b2=tan(theta2.*pi);
                    X2=Walrus_Second(jj)-a2*b2*abs(Walrus_Second(jj)-Walrus_pos(i,jj));

                    Walrus_posNew(:,jj)=(X1+X2)/2;
                end
            end
        end
        % Handling constraints
        Walrus_posNew = bound(Walrus_posNew(:,1:D),UB,LB);
        % Evaluate function at Walrus_posNew
        Walrus_posNew(:,D + 1: K) = evaluate_objective(Walrus_posNew(:,1:D));

        % For the first trail Walrus_posNew
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (Walrus_posNew(:,D+k)<Walrus_pos(i,D+k))
                dom_less = dom_less + 1;
            elseif (Walrus_posNew(:,D+k)== Walrus_pos(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end
        if dom_more == 0 && dom_equal ~= M % If trial vector Walrus_posNew dominates
            % target vector Xi. Replace Xi by Walrus_posNew in current Walrus_pos system and
            % add Xi to advanced population 2
            Walrus_pos_ad(i,1:K) = Walrus_pos(i,1:K); % add Xi to advanced Walrus_pos
            Walrus_pos(i,1:K) = Walrus_posNew(:,1:K); % replace Xi by Walrus_posNew
        else % else Add Xi (trial vector) to advanced Walrus_pos
            Walrus_pos_ad(i,1:K) = Walrus_posNew;
        end
    end % end for i

    Walrus_pos_com = [Walrus_pos(:,1:K) ; Walrus_pos_ad];
    intermediate_Walrus_pos = non_domination_sort_mod(Walrus_pos_com, M, D);
    Pop  = replace_chromosome(intermediate_Walrus_pos, M,D,Nr);
    Walrus_pos = Pop;

    DOMINATED = checkDomination(Pop(:,D+1:D+M));
    Pop  = Pop(~DOMINATED,:);
    Walrus_pos1 = Pop(:,1:K+1);

    if verbose && rem(Iteration, ishow) == 0
        disp(['MOWO Iteration ' num2str(Iteration) ': Number of solutions in the archive = ' num2str(size(Walrus_pos1,1))]);
    end
    Iteration = Iteration+1;
    % Plot Pareto Front
    if ~showPlot, continue; end
    h_fig = figure(1);

    if(numOfObj == 2)
        pl_data = Walrus_pos1(:,D+1:D+M); % extract data to plot
        POS_fit = sortrows(pl_data,2);
        figure(h_fig);
        try delete(h_par); end
        h_par = plot(POS_fit(:,1),POS_fit(:,2),'ob'); hold on;
        if(isfield(MultiObj,'truePF'))
            try delete(h_pf); end
            h_pf = plot(MultiObj.truePF(:,1),MultiObj.truePF(:,2),'.k'); hold on;
        end
        title(name);
        grid on; xlabel('f1'); ylabel('f2');
        set(gca,'FontSize',20,'linewidth',2);
        drawnow;
    end
    if(numOfObj == 3)
        pl_data = Walrus_pos1(:,D+1:D+M); % extract data to plot
        POS_fit = sortrows(pl_data,3);
        figure(h_fig);
        try delete(h_par); end

        h_par = plot3(POS_fit(:,1),POS_fit(:,2),POS_fit(:,3),'ob'); hold on;

        if(isfield(MultiObj,'truePF'))
            try delete(h_pf); end
            h_pf = plot3(MultiObj.truePF(:,1),MultiObj.truePF(:,2),MultiObj.truePF(:,3),'.k'); hold on;
        end
        title(name);
        grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
        set(gca,'FontSize',20,'linewidth',2);
        drawnow;
    end
end
f = Walrus_pos1(:,D+1:D+M);
x = Walrus_pos1(:,1:D);
if showPlot
    hold off;
    if(isfield(MultiObj,'truePF'))
        legend('MOWO','truePF');
    else
        legend('MOWO');
    end
end
end
%% Check the boundary limit
function a=bound(a,ub,lb)
a(a>ub)=ub(a>ub); a(a<lb)=lb(a<lb);
end
function d = dominates(x,y)
d = all(x<=y,2) & any(x<y,2);
end
function dom_vector = checkDomination(fitness)
Np = size(fitness,1);
dom_vector = zeros(Np,1);
all_perm = nchoosek(1:Np,2);    % Possible permutations
all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];

d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
dominated_particles = unique(all_perm(d==1,2));
dom_vector(dominated_particles) = 1;
end