%% Main_Ottimizzato – Orchestratore modulare per il confronto di metodi di segmentazione
%
% BRANCH: feature/tumor-slice-filter
%
% Questo branch implementa l'Oracle Sanity Check: la valutazione finale
% viene condotta SOLO sulle slice che contengono tumore nella Ground Truth.
% La Grid Search usa F1-Score come criterio di ottimizzazione (non Accuracy).
%
% Differenze rispetto al branch main:
%   - Grid Search: max(F1) invece di max(Accuracy)  [condiviso con main]
%   - Valutazione Test Set: solo_slice_positive = true  [solo questo branch]
%   - Output: unica tabella Oracle  [solo questo branch]
%
% Vedere §6.4 della relazione per la discussione metodologica.

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

if num_pazienti ~= length(lista_labels)
    error('Il numero di immagini (%d) non coincide con il numero di label (%d)!', ...
        num_pazienti, length(lista_labels));
end

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
rng('default');
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
% Criterio: max(F1-Score medio) — immune all'accuracy paradox
% =========================================================================
ESEGUI_GRID_SEARCH = true;

if ESEGUI_GRID_SEARCH
    filtri_candidati = [
        struct('nome', 'Median 3x3',          'tipo', 'median',   'param', [3 3]);
        struct('nome', 'Gaussian sigma=0.5',  'tipo', 'gaussian', 'param', 0.5);
        struct('nome', 'Wiener 3x3',          'tipo', 'wiener',   'param', [3 3]);
        ];

    best_models = grid_search_denoising(metodi_da_testare, ...
        cartella_immagini, cartella_labels, ...
        lista_immagini(idx_valid), lista_labels(idx_valid), ...
        num_valid, filtri_candidati);

    lista_valutazione_img = lista_immagini(idx_test);
    lista_valutazione_lbl = lista_labels(idx_test);
    num_valutazione       = length(idx_test);

else
    default_filtro = struct('nome', 'Median 3x3', 'tipo', 'median', 'param', [3 3]);
    num_metodi     = numel(metodi_da_testare);
    best_models    = struct();

    for m = 1:num_metodi
        best_models(m).nome    = metodi_da_testare{m};
        best_models(m).filtro  = default_filtro;
        best_models(m).F1_best = NaN;
        switch metodi_da_testare{m}
            case 'Otsu',      best_models(m).config = 3;
            case 'K-Means',   best_models(m).config = 3;
            case 'Watershed', best_models(m).config = 0.95;
        end
    end

    lista_valutazione_img = lista_immagini;
    lista_valutazione_lbl = lista_labels;
    num_valutazione       = num_pazienti;
end

% =========================================================================
% ORACLE SANITY CHECK — VALUTAZIONE SOLO SULLE SLICE CON TUMORE
%   solo_slice_positive = true: le slice prive di tumore nella GT vengono
%   saltate. Questo elimina il contributo dei FP "ciechi" sulle fette sane
%   e mostra il limite superiore del segmentatore (upper bound).
%   NOTA METODOLOGICA: usa la GT come oracolo → data leakage dichiarato
%   sul protocollo di valutazione (non sul modello di segmentazione).
% =========================================================================
fprintf('\n============================================================\n');
if ESEGUI_GRID_SEARCH
    fprintf('   ORACLE SANITY CHECK — TEST SET (%d pazienti)\n', num_valutazione);
else
    fprintf('   ORACLE SANITY CHECK — DATASET COMPLETO (%d pazienti)\n', num_valutazione);
end
fprintf('   [Solo slice con tumore nella GT | Grid Search con F1]\n');
fprintf('============================================================\n\n');

Risultati_Oracle = cell(1, numel(metodi_da_testare));

for m = 1:numel(metodi_da_testare)
    nome_metodo = best_models(m).nome;
    filtro_best = best_models(m).filtro;
    config_best = best_models(m).config;

    fprintf('\n======================================================\n');
    fprintf('   METODO:  %s\n', nome_metodo);
    fprintf('   Filtro:  %s\n', filtro_best.nome);
    fprintf('   Config:  %s\n', format_config_main(nome_metodo, config_best));
    if ESEGUI_GRID_SEARCH && ~isnan(best_models(m).F1_best)
        fprintf('   F1 val (validation):  %.4f\n', best_models(m).F1_best);
    end
    fprintf('======================================================\n');

    % Oracle: solo slice con tumore (solo_slice_positive = true)
    [Acc, Sens, Spec, Jacc, F1] = esegui_pipeline_metodo( ...
        nome_metodo, cartella_immagini, cartella_labels, ...
        lista_valutazione_img, lista_valutazione_lbl, num_valutazione, ...
        filtro_best, config_best, false, true);

    Risultati_Oracle{m}.nome        = nome_metodo;
    Risultati_Oracle{m}.filtro      = filtro_best.nome;
    Risultati_Oracle{m}.config      = config_best;
    Risultati_Oracle{m}.Accuracy    = Acc;
    Risultati_Oracle{m}.Sensitivity = Sens;
    Risultati_Oracle{m}.Specificity = Spec;
    Risultati_Oracle{m}.Jaccard     = Jacc;
    Risultati_Oracle{m}.F1          = F1;
end

% =========================================================================
% TABELLA ORACLE — RISULTATI FINALI
% =========================================================================
fprintf('\n\n******************************************************\n');
fprintf('   ORACLE SANITY CHECK — RISULTATI\n');
fprintf('   [Specificity esclude le slice senza tumore]\n');
fprintf('******************************************************\n\n');

fprintf('%-12s  %-22s  %-14s  %6s  %6s  %6s  %6s  %6s\n', ...
    'Metodo', 'Filtro', 'Config', 'Acc', 'Sens', 'Spec', 'Jacc', 'F1');
fprintf('%s\n', repmat('-', 1, 84));

for m = 1:numel(metodi_da_testare)
    Media_Acc  = mean(Risultati_Oracle{m}.Accuracy);
    Media_Sens = mean(Risultati_Oracle{m}.Sensitivity);
    Media_Spec = mean(Risultati_Oracle{m}.Specificity);
    Media_Jacc = mean(Risultati_Oracle{m}.Jaccard);
    Media_F1   = mean(Risultati_Oracle{m}.F1);
    cfg_label  = format_config_main(Risultati_Oracle{m}.nome, Risultati_Oracle{m}.config);
    fprintf('%-12s  %-22s  %-14s  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        Risultati_Oracle{m}.nome, Risultati_Oracle{m}.filtro, cfg_label, ...
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