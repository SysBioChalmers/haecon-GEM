% CHRR flux sampling for the HaeconGEM stage/sex models.
% Run from the repo root so 'model/contextualized_models' resolves.
% Requires the COBRA Toolbox and a working LP solver (Gurobi used here).

MODEL_DIR = 'model/contextualized_models';
MODELS = {
    'egg_c.xml'
    'L1_c.xml'
    'L2_c.xml'
    'L3_c.xml'
    'L4f_c.xml'
    'L4m_c.xml'
    'Am_c.xml'
    'Af_c.xml'
    'L4mf_union_c.xml'
};

NUM_CHAINS = 4;
SAMPLES_PER_CHAIN = 1000;
GROWTH_FLOOR_FRACTION = 0.1;   % lower bound on biomass = 10% of the model's own FBA optimum
SEED_BASE = 20260806;

changeCobraSolver('gurobi', 'LP');

for i = 1:numel(MODELS)
    name = erase(MODELS{i}, '.xml');
    model = ravenCobraWrapper(importModel(fullfile(MODEL_DIR, MODELS{i})));

    fba = optimizeCbModel(model, 'max');
    if fba.stat ~= 1
        fprintf('%s: infeasible, skipping\n', name);
        continue
    end

    objRxn = model.rxns{find(model.c ~= 0, 1)};
    model = changeRxnBounds(model, objRxn, GROWTH_FLOOR_FRACTION * fba.f, 'l');
    model = changeObjective(model, objRxn, 0);

    blocked = findBlockedReaction(model, 'FVA');
    model = removeRxns(model, blocked);

    rng(SEED_BASE + i);
    [samples1, roundedPolytope] = chrrSampler(model, [], SAMPLES_PER_CHAIN, 1, [], false, 100);

    chainMedians = median(samples1, 2);
    for chain = 2:NUM_CHAINS
        rng(SEED_BASE + i*1000 + chain);
        samplesK = chrrSampler(model, [], SAMPLES_PER_CHAIN, 1, roundedPolytope, false, 100);
        chainMedians = [chainMedians, median(samplesK, 2)];
    end

    writetable(table(model.rxns, mean(chainMedians, 2), 'VariableNames', {'rxnID', 'medianFlux'}), ...
        sprintf('%s_median.csv', name));
    writetable(array2table(samples1', 'VariableNames', model.rxns), ...
        sprintf('%s_samples.csv', name));

    fprintf('%s: done, %d reactions, %d blocked removed\n', name, numel(model.rxns), numel(blocked));
end
