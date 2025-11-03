
%Step-by-step to create the HaeconGEMv0.2



% Load template model
H_model = load('model/Human-GEM.mat')
H_model = H_model.ihuman
W_model = load('model/Worm-GEM.mat')
W_model = W_model.wormGEM



%Load orthologue pairs table
orthopair = readcell('ultimate_orthologues_supportByGeneDF_3_80_justPairs.csv','Delimiter','\t')


%Create new model
haecon_model = getModelFromOrthology(H_model,orthopair)
haecon_model.id = 'Haecon-GEM'
save('model/Haecon_GEMv03.mat','haecon_model')



%Starting filling gaps
setRavenSolver('gurobi');


haecon_model = setParam(haecon_model,'obj', 'MAR04413', 1);
haecon_model = setParam(haecon_model,'lb', 'MAR04413', 0.01);





%Penalizing donnors
penalize = @(M) ( ...
(sum(M.S~=0,1)==1) | ...
reshape(startsWith(lower(string(M.rxns(:))),'ex_'),1,[]) | ...
( isfield(M,'rxnNames')   & reshape(contains(lower(string(M.rxnNames(:))),   'exchange'), 1,[]) ) | ...
( isfield(M,'subSystems') & reshape(arrayfun(@(i) ...
( i<=numel(M.subSystems) && ...
( ( ischar(M.subSystems{i})    && contains(lower(string(M.subSystems{i})),    'transport') ) || ...
( iscell(M.subSystems{i})    && any(contains(lower(string([M.subSystems{i}{:}])), 'transport')) ) ) ), ...
1:numel(M.rxns)), 1,[]) ) ...
);



% Indices to penalize in each donor
ix1 = find(penalize(W_model));
ix2 = find(penalize(H_model));




% Initialize and apply reaction scores
rxnScores = cell(2,1);
rxnScores{1} = -1   * ones(numel(W_model.rxns),1);   % worm preferred
rxnScores{2} = -0.2 * ones(numel(H_model.rxns),1);   % human allowed
rxnScores{1}(ix1) = rxnScores{1}(ix1) - 50;          % penalize exchanges
rxnScores{2}(ix2) = rxnScores{2}(ix2) - 50;


% Biomass bounds fix and objective and detect (if needed), fix bounds, and set objective
bioID = 'MAR04413'; %This is tho default obective at HumanGEM
% set biomass as objective
haecon_model = setParam(haecon_model,'obj',bioID,1);
i = find(strcmp(haecon_model.rxns,bioID));


% Fix: ensure UB >= LB and positive UB
lbTarget = 1e-3;                                    % your 0.001
haecon_model = setParam(haecon_model,'lb',i,lbTarget);
if ~isfinite(haecon_model.ub(i)) || haecon_model.ub(i) < lbTarget
    haecon_model = setParam(haecon_model,'ub',i,1000);  % allow production
end
haecon_model.rev(i) = false;                         % biomass forward-only

      
% targeted gap-filling (enforce biomass >= 0.001)
useModelConstraints = true;


% MILP params (the version you used)
params = struct('timeLimit',1800,'printLevel',0,'solver','gurobi','mipGap',1e-4);

%Run filling gaps
[newConnected, cannotConnect, addedRxns, newModel, exitFlag] = fillGaps(haecon_model, {W_model, H_model}, false, useModelConstraints, false, rxnScores, params);


fprintf('exitFlag=%d, added %d reactions\n', exitFlag, numel(addedRxns));



% Filling gaps using the objective in WormGEM
bioID = 'MAR00021';                             % your biomass rxn
model2 = setParam(newModel,'obj',bioID,1);              % set objective
model2 = setParam(model2,'lb', bioID, 0.001);         % enforce growth target
model2 = setParam(model2,'ub', bioID, 1000);
model2.rev(bioID) = false;

% 2nd filling gaps
[newConnected2, cannotConnect2, addedRxns2, Model3, exitFlag2] = fillGaps(model2, {W_model, H_model},false, useModelConstraints, false, rxnScores, params);
exitFlag2
fprintf('exitFlag=%d, added %d reactions\n', exitFlag2, numel(addedRxns2));

% find the reactions currently set as objective in Model3
objIdx = find(Model3.c ~= 0);
table(Model3.rxns(objIdx), Model3.rxnNames(objIdx), Model3.c(objIdx), 'VariableNames',{'rxn','name','coef'})



save('model/Haecon_GEMv04.mat','newModel')
save('model/Haecon_GEMv05.mat','Model3')