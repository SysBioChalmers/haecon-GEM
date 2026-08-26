

W_model = load('model/Worm-GEM.mat')
W_model = W_model.wormGEM


%Load orthologue pairs table
T = readtable('data/orthologues_supported/suppOrth_table_bygene_3_80_elegans.csv', ...
              'Delimiter','\t', ...
              'FileType','text', ...
              'TextType','string')

orthopair = cellfun(@char, table2cell(T(:,1:2)), 'UniformOutput', false)
haecon_model = getModelFromOrthology(W_model,orthopair)
haecon_model.id = 'Haecon-GEM'
save('model/Haecon_GEMv03_fromWorm.mat','haecon_model')