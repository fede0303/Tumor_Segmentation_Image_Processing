function best_models = grid_search_denoising(metodi_test, cartella_immagini, ...
    cartella_labels, lista_immagini, lista_labels, num_pazienti, filtri)
% GRID_SEARCH_DENOISING — Trova la combinazione ottimale (filtro, config) per ogni metodo.
%
%   best_models = grid_search_denoising(metodi_test, cartella_immagini,
%       cartella_labels, lista_immagini, lista_labels, num_pazienti, filtri)
%
%   INPUT
%   metodi_test      : cell array di stringhe (es. {'Otsu','K-Means','Watershed'})
%   cartella_immagini, cartella_labels, lista_immagini, lista_labels, num_pazienti :
%                      dati del validation set (già verificati dal Main)
%   filtri           : array di strutture con campi 'nome', 'tipo', 'param'
%
%   OUTPUT
%   best_models      : struct array (1 × num_metodi) con campi:
%                        .nome    — nome del metodo
%                        .filtro  — struct del miglior filtro trovato
%                        .config  — valore della miglior configurazione trovata
%                        .Accuracy_best — Accuracy media sui pazienti di validazione per (filtro, config) best
%
%   LOGICA
%   I loop sono organizzati con il METODO come ciclo esterno e (FILTRO × CONFIG)
%   come cicli interni. Per ogni metodo si effettua una ricerca esaustiva su
%   tutte le combinazioni (filtro, config) e si seleziona quella con il massimo
%   Accuracy media sui pazienti di validazione. Il risultato è che ogni metodo viene
%   ottimizzato indipendentemente, con il proprio filtro e la propria configurazione
%   migliore. La valutazione finale nel Main confronterà quindi le pipeline
%   ottimizzate di ciascun metodo.

% =========================================================================
% DEFINIZIONE DELLE CONFIGURAZIONI PER OGNI METODO
% =========================================================================
% Le configurazioni da esplorare sono definite qui, in quanto la grid search
% deve iterarle esplicitamente. esegui_pipeline_metodo riceve ora la singola
% config come parametro esterno.
configs_map = struct();
configs_map.Otsu      = [2, 3, 4];
configs_map.KMeans    = [2, 3, 4];       % campo senza trattino per compatibilità MATLAB
configs_map.Watershed = [0.90, 0.95, 0.98];

num_metodi = length(metodi_test);
num_filtri = length(filtri);

% Preallocazione della struct array di output
best_models = struct();

fprintf('\n============================================================\n');
fprintf('   GRID SEARCH PER-METHOD (%d pazienti nel validation set)\n', num_pazienti);
fprintf('   Logica: max(Accuracy media sui pazienti) per ogni (metodo x filtro x config)\n');
fprintf('============================================================\n\n');

% =========================================================================
% LOOP ESTERNO: METODI
% =========================================================================
for m = 1:num_metodi
    metodo = metodi_test{m};

    % Recupera le configurazioni: 'K-Means' → field 'KMeans' (rimozione trattino)
    fieldname = strrep(metodo, '-', '');
    configs   = configs_map.(fieldname);
    num_conf  = length(configs);

    fprintf('╔══════════════════════════════════════════════════════════╗\n');
    fprintf('║  OTTIMIZZAZIONE METODO: %-34s║\n', metodo);
    fprintf('╚══════════════════════════════════════════════════════════╝\n');

    % Tracciamento del best per questo metodo
    Accuracy_best = -Inf;
    filtro_best_m = [];
    config_best_m = NaN;

    % =====================================================================
    % LOOP INTERNO: FILTRI × CONFIGURAZIONI
    % =====================================================================
    for f = 1:num_filtri
        fprintf('\n  Filtro %d/%d: %s\n', f, num_filtri, filtri(f).nome);

        for c = 1:num_conf
            cfg = configs(c);

            tic;
            % Esecuzione della pipeline con la singola (filtro, config)
            % mostra_sanity = false: nessuna figura durante la grid search
            [Dataset_Acc, ~, ~, ~, ~] = esegui_pipeline_metodo( ...
                metodo, cartella_immagini, cartella_labels, ...
                lista_immagini, lista_labels, num_pazienti, ...
                filtri(f), cfg, false);
            tempo = toc;

            Accuracy_media  = mean(Dataset_Acc);
            cfg_label = format_config(metodo, cfg);

            fprintf('    Config %-14s → Accuracy media = %.4f  (%.1f sec)\n', ...
                cfg_label, Accuracy_media, tempo);

            % Aggiorna il best se questo (filtro, config) supera il massimo corrente
            if Accuracy_media > Accuracy_best
                Accuracy_best       = Accuracy_media;
                filtro_best_m = filtri(f);
                config_best_m = cfg;
            end
        end
    end

    % Salva il best trovato per questo metodo
    best_models(m).nome    = metodo;
    best_models(m).filtro  = filtro_best_m;
    best_models(m).config  = config_best_m;
    best_models(m).Accuracy_best = Accuracy_best;

    fprintf('\n  >> BEST per %-10s: filtro="%s"  config=%-12s  Accuracy=%.4f\n\n', ...
        metodo, filtro_best_m.nome, format_config(metodo, config_best_m), Accuracy_best);
end

% =========================================================================
% RIEPILOGO FINALE
% =========================================================================
fprintf('============================================================\n');
fprintf('   RIEPILOGO GRID SEARCH — CONFIGURAZIONI OTTIMALI\n');
fprintf('============================================================\n');
fprintf('%-12s  %-22s  %-14s  %s\n', 'Metodo', 'Filtro Ottimale', 'Config Ottimale', 'Accuracy val');
fprintf('%s\n', repmat('-', 1, 66));
for m = 1:num_metodi
    fprintf('%-12s  %-22s  %-14s  %.4f\n', ...
        best_models(m).nome, ...
        best_models(m).filtro.nome, ...
        format_config(best_models(m).nome, best_models(m).config), ...
        best_models(m).Accuracy_best);
end
fprintf('\n');

end

% =========================================================================
% FUNZIONE DI SUPPORTO: Formatta il valore di config come stringa leggibile
% =========================================================================
function s = format_config(metodo, cfg)
    switch metodo
        case 'Otsu'
            s = sprintf('%d soglie', cfg);
        case 'K-Means'
            s = sprintf('k=%d', cfg);
        case 'Watershed'
            s = sprintf('p=%.2f', cfg);
        otherwise
            s = sprintf('%g', cfg);
    end
end