%% Main_Ottimizzato – Orchestratore modulare per il confronto di metodi di segmentazione
%
% Il presente script esegue il batch processing sull'intero dataset di volumi
% MRI (Task01_BrainTumour) e coordina la valutazione comparativa di
% Otsu, K-Means e Watershed.
%
% ARCHITETTURA AGGIORNATA (Grid Search Per-Method):
%   - La grid search trova il miglior (filtro, config) PER OGNI METODO
%     separatamente: Otsu, K-Means e Watershed vengono ognuno ottimizzato
%     in modo indipendente sul validation set (80%).
%   - La valutazione finale confronta "il miglior Otsu possibile" vs
%     "il miglior K-Means possibile" vs "il miglior Watershed possibile",
%     ciascuno con il proprio filtro e configurazione ottimali.
%   - Questo risponde alla domanda: "Quale algoritmo, quando ottimizzato,
%     raggiunge la performance più alta sul test set (20%)?"
%   - Il confronto è equo perché ogni metodo ha avuto le stesse opportunità
%     di ottimizzazione (stessi filtri candidati, stesse configurazioni candidate).
%
% Se ESEGUI_GRID_SEARCH = false, viene usato un filtro e una config di default
% per ogni metodo (Median 3x3, configurazione intermedia) sull'intero dataset.

clear; clc; close all;

% =========================================================================
% ACQUISIZIONE DEL DATASET E VERIFICA DI INTEGRITÀ
% =========================================================================
cartella_immagini = "/Users/federicoraimondi/Library/CloudStorage/OneDrive-PolitecnicodiBari/Università/Magistrale/1 ° Anno/1° Semestre/Guerriero/Progetto_Image_Processing/data/Task01_BrainTumour/immaginiMRI/";
cartella_labels   = "/Users/federicoraimondi/Library/CloudStorage/OneDrive-PolitecnicodiBari/Università/Magistrale/1 ° Anno/1° Semestre/Guerriero/Progetto_Image_Processing/data/Task01_BrainTumour/GroundTruth/";

lista_immagini = dir(fullfile(cartella_immagini, '*.nii.gz'));
lista_labels   = dir(fullfile(cartella_labels,   '*.nii.gz'));

num_pazienti = length(lista_immagini);
fprintf('Pazienti individuati nel subset: %d\n\n', num_pazienti);

% Controllo di consistenza 1: stesso numero di immagini e label
if num_pazienti ~= length(lista_labels)
    error('Il numero di immagini (%d) non coincide con il numero di label (%d)!', ...
        num_pazienti, length(lista_labels));
end

% Controllo di consistenza 2: corrispondenza 1:1 dei nomi file
for k = 1:num_pazienti
    [~, nome_base_img, ext_img] = fileparts(lista_immagini(k).name);
    if strcmp(ext_img, '.gz'), [~, nome_base_img, ~] = fileparts(nome_base_img); end
    [~, nome_base_lbl, ext_lbl] = fileparts(lista_labels(k).name);
    if strcmp(ext_lbl, '.gz'), [~, nome_base_lbl, ~] = fileparts(nome_base_lbl); end
    if ~strcmp(nome_base_img, nome_base_lbl)
        error('Disallineamento: %s vs %s', lista_immagini(k).name, lista_labels(k).name);
    end
end
fprintf('Integrità del subset verificata: ogni immagine ha la propria label.\n\n');

% =========================================================================
% SPLIT VALIDATION / TEST (80% / 20%)
% =========================================================================
% VALIDATION SET (80%): usato ESCLUSIVAMENTE dalla grid search per trovare
%   il miglior (filtro, config) per ogni metodo. Non compare nella valutazione finale.
% TEST SET (20%): usato ESCLUSIVAMENTE per la valutazione finale.
%   Non viene mai usato per la selezione di filtro o config.
% Questa separazione garantisce l'assenza di data leakage.

rng('default');                         % per riproducibilità
indici_casuali = randperm(num_pazienti);
num_valid      = round(num_pazienti * 0.8);
idx_valid      = indici_casuali(1:num_valid);
idx_test       = indici_casuali(num_valid+1:end);

fprintf('Pazienti di validazione: %d\n', num_valid);
fprintf('Pazienti di test:        %d\n\n', num_pazienti - num_valid);

% =========================================================================
% DEFINIZIONE DEI METODI DI SEGMENTAZIONE DA CONFRONTARE
% =========================================================================
metodi_da_testare = {'Otsu', 'K-Means', 'Watershed'};

% =========================================================================
% GRID SEARCH PER-METHOD (OPZIONALE) — SOLO SUL VALIDATION SET
% =========================================================================
ESEGUI_GRID_SEARCH = true;   % false per usare i default senza grid search

if ESEGUI_GRID_SEARCH
    % Definizione dei filtri candidati da esplorare
    filtri_candidati = [
        struct('nome', 'Median 3x3',          'tipo', 'median',   'param', [3 3]);
        struct('nome', 'Gaussian sigma=0.5',  'tipo', 'gaussian', 'param', 0.5);
        struct('nome', 'Wiener 3x3',          'tipo', 'wiener',   'param', [3 3]);
        ];

    % Grid search per-method sul validation set.
    % Restituisce best_models: struct array con il miglior (filtro, config) per ogni metodo.
    best_models = grid_search_denoising(metodi_da_testare, ...
        cartella_immagini, cartella_labels, ...
        lista_immagini(idx_valid), lista_labels(idx_valid), ...
        num_valid, filtri_candidati);

    % La valutazione finale avviene SOLO sul test set (mai visto durante la grid search)
    lista_valutazione_img = lista_immagini(idx_test);
    lista_valutazione_lbl = lista_labels(idx_test);
    num_valutazione       = length(idx_test);

else
    % Default: filtro Median e configurazione intermedia per ogni metodo.
    % Nessuna ottimizzazione eseguita → potenza computazionale minima.
    default_filtro = struct('nome', 'Median 3x3', 'tipo', 'median', 'param', [3 3]);
    num_metodi     = numel(metodi_da_testare);
    best_models    = struct();

    for m = 1:num_metodi
        best_models(m).nome    = metodi_da_testare{m};
        best_models(m).filtro  = default_filtro;
        best_models(m).Accuracy_best = NaN;   % non calcolato (nessuna grid search)
        switch metodi_da_testare{m}
            case 'Otsu',      best_models(m).config = 3;    % 3 soglie (default)
            case 'K-Means',   best_models(m).config = 3;    % k=3 (default)
            case 'Watershed', best_models(m).config = 0.95; % p=0.95 (default)
        end
    end

    % Senza grid search si valuta sull'intero dataset (no leakage perché non si è
    % effettuata nessuna selezione su alcun sottinsieme)
    lista_valutazione_img = lista_immagini;
    lista_valutazione_lbl = lista_labels;
    num_valutazione       = num_pazienti;
end

% =========================================================================
% VALUTAZIONE FINALE — OGNI METODO CON IL PROPRIO BEST (FILTRO, CONFIG)
% =========================================================================
fprintf('\n============================================================\n');
if ESEGUI_GRID_SEARCH
    fprintf('   VALUTAZIONE FINALE SUL TEST SET (%d pazienti)\n', num_valutazione);
else
    fprintf('   VALUTAZIONE FINALE SULL''INTERO DATASET (%d pazienti)\n', num_valutazione);
end
fprintf('============================================================\n\n');

Risultati = cell(1, numel(metodi_da_testare));

for m = 1:numel(metodi_da_testare)
    nome_metodo = best_models(m).nome;
    filtro_best = best_models(m).filtro;
    config_best = best_models(m).config;

    fprintf('\n======================================================\n');
    fprintf('   METODO:  %s\n', nome_metodo);
    fprintf('   Filtro:  %s\n', filtro_best.nome);
    fprintf('   Config:  %s\n', format_config_main(nome_metodo, config_best));
    if ESEGUI_GRID_SEARCH && ~isnan(best_models(m).F1_best)
        fprintf('   F1 val:  %.4f\n', best_models(m).F1_best);
    end
    fprintf('======================================================\n');

    [Dataset_Acc, Dataset_Sens, Dataset_Spec, Dataset_Jacc, Dataset_F1] = ...
        esegui_pipeline_metodo(nome_metodo, cartella_immagini, cartella_labels, ...
        lista_valutazione_img, lista_valutazione_lbl, num_valutazione, ...
        filtro_best, config_best);

    Risultati{m}.nome        = nome_metodo;
    Risultati{m}.filtro      = filtro_best.nome;
    Risultati{m}.config      = config_best;
    Risultati{m}.Accuracy    = Dataset_Acc;
    Risultati{m}.Sensitivity = Dataset_Sens;
    Risultati{m}.Specificity = Dataset_Spec;
    Risultati{m}.Jaccard     = Dataset_Jacc;
    Risultati{m}.F1          = Dataset_F1;
end

% =========================================================================
% CONFRONTO FINALE — TABELLA RIEPILOGATIVA
% =========================================================================
fprintf('\n\n******************************************************\n');
fprintf('   ELABORAZIONE COMPLETATA - RISULTATI A CONFRONTO\n');
fprintf('******************************************************\n\n');

fprintf('%-12s  %-22s  %-14s  %6s  %6s  %6s  %6s  %6s\n', ...
    'Metodo', 'Filtro', 'Config', 'Acc', 'Sens', 'Spec', 'Jacc', 'F1');
fprintf('%s\n', repmat('-', 1, 84));

for m = 1:numel(metodi_da_testare)
    Media_Acc  = mean(Risultati{m}.Accuracy);
    Media_Sens = mean(Risultati{m}.Sensitivity);
    Media_Spec = mean(Risultati{m}.Specificity);
    Media_Jacc = mean(Risultati{m}.Jaccard);
    Media_F1   = mean(Risultati{m}.F1);
    cfg_label  = format_config_main(Risultati{m}.nome, Risultati{m}.config);
    fprintf('%-12s  %-22s  %-14s  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        Risultati{m}.nome, Risultati{m}.filtro, cfg_label, ...
        Media_Acc, Media_Sens, Media_Spec, Media_Jacc, Media_F1);
end

fprintf('\n');

% =========================================================================
% FUNZIONE DI SUPPORTO — Formatta la config come stringa leggibile
% =========================================================================
function s = format_config_main(metodo, cfg)
    switch metodo
        case 'Otsu',      s = sprintf('%d soglie', cfg);
        case 'K-Means',   s = sprintf('k=%d', cfg);
        case 'Watershed', s = sprintf('p=%.2f', cfg);
        otherwise,        s = sprintf('%g', cfg);
    end
end