

function ravenModel = convertCobraToRaven(cobraModel)

    ravenModel = struct();

    % Basic identifiers
    ravenModel.id          = cobraModel.description;
    ravenModel.name        = cobraModel.description;
    ravenModel.description = '';

    % Core fields
    ravenModel.rxns        = cobraModel.rxns;
    ravenModel.rxnNames    = cobraModel.rxnNames;
    ravenModel.mets        = cobraModel.mets;
    ravenModel.metNames    = cobraModel.metNames;

    ravenModel.S           = cobraModel.S;
    ravenModel.lb          = cobraModel.lb;
    ravenModel.ub          = cobraModel.ub;
    ravenModel.c           = cobraModel.c;
    ravenModel.b           = cobraModel.b;

    % rev sometimes int64 → convert to double logical
    ravenModel.rev         = double(cobraModel.rev);

    % GPR and genes
    ravenModel.grRules     = cobraModel.grRules;
    ravenModel.genes       = cobraModel.genes;
    ravenModel.rxnGeneMat  = cobraModel.rxnGeneMat;

    % Compartments
    ravenModel.comps       = cobraModel.comps;
    ravenModel.compNames   = cobraModel.compNames;

    % Metabolite properties
    ravenModel.metComps    = double(cobraModel.metComps);
    ravenModel.metFormulas = cobraModel.metFormulas;
    ravenModel.metCharges  = cobraModel.metCharges;

    % Reaction metadata
    ravenModel.subSystems         = cobraModel.subSystems;
    ravenModel.rxnNotes           = cobraModel.rxnNotes;
    ravenModel.rxnReferences      = cobraModel.rxnReferences;

    % EC numbers mapping (COBRA → RAVEN)
    if isfield(cobraModel, 'rxnECNumbers')
        ravenModel.eccodes = cobraModel.rxnECNumbers;
    else
        ravenModel.eccodes = cell(size(cobraModel.rxns));
    end

    % Confidence scores (convert if needed)
    if iscell(cobraModel.rxnConfidenceScores)
        try
            ravenModel.rxnConfidenceScores = cellfun(@str2double, cobraModel.rxnConfidenceScores);
        catch
            ravenModel.rxnConfidenceScores = zeros(numel(cobraModel.rxns),1);
        end
    else
        ravenModel.rxnConfidenceScores = cobraModel.rxnConfidenceScores;
    end

    % Optional / missing fields (create empty placeholders)
    nRxns = numel(cobraModel.rxns);
    nMets = numel(cobraModel.mets);
    nGenes = numel(cobraModel.genes);

    ravenModel.geneShortNames = cell(nGenes,1);
    ravenModel.inchis         = cell(nMets,1);
    ravenModel.version        = '';
    ravenModel.annotation     = struct();

end
