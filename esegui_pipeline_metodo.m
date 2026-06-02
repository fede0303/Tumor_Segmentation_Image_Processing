% ESEGUI_PIPELINE_METODO — Incapsula l'intera pipeline per un singolo metodo di segmentazione.
%
%   INPUT
%   - nome_metodo        : stringa, identifica l'algoritmo di segmentazione
%                           ammesso tra 'Otsu', 'K-Means', 'Watershed'.
%   - cartella_immagini  : percorso alla directory contenente i volumi MRI (.nii.gz).
%   - cartella_labels    : percorso alla directory delle corrispondenti Ground Truth.
%   - lista_immagini     : struct generata da dir() con l'elenco dei file delle immagini.
%   - lista_labels       : struct analoga per le label.
%   - num_pazienti       : numero di soggetti da elaborare.
%   - filtro             : struct con campi 'tipo' e 'param' per il filtro di denoising.
%                          Esempi: struct('tipo','median','param',[3 3])
%                                  struct('tipo','gaussian','param',0.5)
%                                  struct('tipo','wiener','param',[3 3])
%   - config_singola     : scalare. Configurazione ottimale del metodo fornita dall'esterno.
%                          Per Otsu e K-Means: intero in {2, 3, 4}.
%                          Per Watershed: reale in {0.90, 0.95, 0.98}.
%   - mostra_sanity      : (opzionale) logico, default true. Se false, disabilita
%                          la visualizzazione di controllo (utile durante grid search).
%
%   OUTPUT
%   - Dataset_Acc        : vettore (num_pazienti × 1) — Accuracy per paziente.
%   - Dataset_Sens       : idem per Sensitivity.
%   - Dataset_Spec       : idem per Specificity.
%   - Dataset_Jacc       : idem per Jaccard (IoU).
%   - Dataset_F1         : idem per Dice coefficient (F1-Score).
%
%   DESCRIZIONE
%   La funzione realizza il ciclo completo su tutti i pazienti:
%       1. Carica il volume FLAIR (canale 1 del 4D) e la Ground Truth;
%       2. Esegue la normalizzazione VOLUMETRICA dell'intensità: calcola
%          min e max globali sull'intero volume FLAIR (non slice per slice)
%          per garantire coerenza della scala tra le slice del paziente;
%       3. Definisce la ROI anatomica (40% centrale del volume);
%       4. Itera su tutte le slice della ROI effettuando pre-processing
%          (normalizzazione volumetrica, filtro di denoising, resize) e
%          chiama la funzione di segmentazione con la config_singola ricevuta;
%       5. Accumula i conteggi TP/TN/FP/FN e calcola le metriche a fine
%          volume (approccio volumetrico, anti-paradosso di Simpson);
%       6. Ripete per tutti i pazienti e restituisce i vettori finali.

function [Dataset_Acc, Dataset_Sens, Dataset_Spec, Dataset_Jacc, Dataset_F1] = ...
    esegui_pipeline_metodo(nome_metodo, cartella_immagini, cartella_labels, ...
    lista_immagini, lista_labels, num_pazienti, filtro, config_singola, mostra_sanity, ...
    solo_slice_positive)

    % Imposta il default per la visualizzazione di controllo
    if nargin < 9
        mostra_sanity = true;
    end

    % Se true, nel ciclo sulle slice vengono considerate SOLO le slice che
    % contengono almeno un voxel tumorale nella Ground Truth.
    % Questo permette di confrontare le metriche con e senza il contributo
    % dei falsi positivi generati sulle slice prive di lesione.
    % Default false: comportamento identico alla versione originale.
    if nargin < 10
        solo_slice_positive = false;
    end

% =========================================================================
% SELEZIONE DELLA FUNZIONE DI SEGMENTAZIONE
% =========================================================================
% Lo switch associa al nome del metodo il function handle seg_func.
% La config_singola viene ricevuta dall'esterno (dalla grid search o dal Main):
% non viene più definita internamente, separando la responsabilità di
% selezione della configurazione da quella di esecuzione della pipeline.

SE = strel('disk', 5);

switch nome_metodo
    case 'Otsu'
        seg_func = @segmenta_otsu;

    case 'K-Means'
        seg_func = @segmenta_kmeans;

    case 'Watershed'
        SE_grad   = strel('disk', 1);
        SE_smooth = strel('disk', 3);
        % Pattern Adapter tramite Closure: ignora SE standard e inietta i due SE specifici.
        seg_func = @(img, cfg, ~) segmenta_watershed(img, cfg, SE_grad, SE_smooth);

    otherwise
        error('Metodo sconosciuto: %s', nome_metodo);
end

% Pre-allocazione dei vettori di output (num_pazienti × 1).
% Con config_singola l'output è un vettore, non una matrice.
Dataset_Acc  = zeros(num_pazienti, 1);
Dataset_Sens = zeros(num_pazienti, 1);
Dataset_Spec = zeros(num_pazienti, 1);
Dataset_Jacc = zeros(num_pazienti, 1);
Dataset_F1   = zeros(num_pazienti, 1);

% Dimensione standard (BraTS Task01 è nativo 240×240; imresize è no-op ma garantisce uniformità).
dim_standard = [240 240];

% =========================================================================
% CICLO PRINCIPALE SUI PAZIENTI
% =========================================================================
for paziente_idx = 1:num_pazienti

    % -----------------------------------------------------------------
    % CARICAMENTO DEI VOLUMI
    % -----------------------------------------------------------------
    nome_file_img = fullfile(cartella_immagini, lista_immagini(paziente_idx).name);
    nome_file_lbl = fullfile(cartella_labels,   lista_labels(paziente_idx).name);

    fprintf('Paziente %d/%d: %s\n', paziente_idx, num_pazienti, lista_immagini(paziente_idx).name);

    % Caricamento volume 4D ed estrazione canale FLAIR (modalità 0 = indice 1 in MATLAB)
    info_img  = niftiinfo(nome_file_img);
    vol_4D    = niftiread(info_img);
    vol_FLAIR = vol_4D(:, :, :, 1);

    % Caricamento Ground Truth e binarizzazione Whole Tumor (etichette 1+2+4 → >0)
    info_lbl   = niftiinfo(nome_file_lbl);
    vol_GT_raw = niftiread(info_lbl);
    vol_GT     = vol_GT_raw > 0;

    % -----------------------------------------------------------------
    % NORMALIZZAZIONE VOLUMETRICA
    % -----------------------------------------------------------------
    % Calcolo di min e max UNA SOLA VOLTA sull'intero volume 3D.
    % Questo garantisce che ogni slice abbia la stessa scala assoluta:
    % un pixel con intensità X nella slice 40 e un pixel con la stessa
    % intensità X nella slice 80 producono lo stesso valore normalizzato.
    % Con mat2gray slice-by-slice questo non sarebbe garantito: ogni slice
    % avrebbe la propria scala indipendente, falsando le soglie di Otsu,
    % K-Means e il percentile del Watershed (i pixel sani nelle slice prive
    % di tumore verrebbero "stirati" verso 1.0 e classificati come tumore).
    min_vol = double(min(vol_FLAIR(:)));
    max_vol = double(max(vol_FLAIR(:)));
    denom   = max(max_vol - min_vol, 1);   % sicurezza su volumi costanti

    % -----------------------------------------------------------------
    % DEFINIZIONE DELLA ROI ANATOMICA (40% Centrale)
    % -----------------------------------------------------------------
    % Esclude il 30% inferiore (base cranica, bulbi oculari) e il 30%
    % superiore (vertice, artefatti) concentrandosi sulle slice dove il
    % tumore è più frequentemente localizzato. Scelta anatomica indipendente
    % dalla Ground Truth.
    profondita     = info_img.ImageSize(3);
    slice_iniziale = round(profondita * 0.3);
    slice_finale   = round(profondita * 0.7);

    % -----------------------------------------------------------------
    % VISUALIZZAZIONE DI CONTROLLO (SANITY CHECK)
    % -----------------------------------------------------------------
    % Per il primo paziente soltanto, se abilitato, viene generata una figura
    % con l'immagine FLAIR filtrata e la maschera della config_singola ottimale.
    if paziente_idx == 1 && mostra_sanity
        slice_equatore = round((slice_iniziale + slice_finale) / 2);
        visualizza_sanity_check(vol_FLAIR, slice_equatore, nome_metodo, ...
            seg_func, config_singola, SE, dim_standard, filtro, min_vol, denom);
    end

    % -----------------------------------------------------------------
    % INIZIALIZZAZIONE ACCUMULATORI VOLUMETRICI (scalari per config singola)
    % -----------------------------------------------------------------
    TP_tot = 0;
    TN_tot = 0;
    FP_tot = 0;
    FN_tot = 0;

    % Seed fisso per riproducibilità di imsegkmeans (K-Means)
    rng(42);

    % =================================================================
    % CICLO PRINCIPALE SULLE SLICE DELLA ROI
    % =================================================================
    for z = slice_iniziale:slice_finale

        % ---------------------------------------------------------
        % GROUND TRUTH DELLA SLICE
        % ---------------------------------------------------------
        gt_slice = vol_GT(:, :, z);
        gt_slice = imresize(gt_slice, dim_standard, 'nearest');

        % ---------------------------------------------------------
        % FILTRO SLICE: se abilitato, salta le slice senza tumore
        % ---------------------------------------------------------
        % Quando solo_slice_positive = true, viene valutata la segmentazione
        % SOLO sulle slice in cui la GT contiene almeno un voxel tumorale.
        % La presenza del tumore è derivata runtime dalla GT stessa;
        % non viene usata nessuna informazione sulla segmentazione prevista.
        if solo_slice_positive && ~any(gt_slice(:))
            continue;   % slice senza tumore: non contribuisce ad alcun contatore
        end

        % ---------------------------------------------------------
        % PRE-PROCESSING DELLA SLICE
        % ---------------------------------------------------------
        slice    = vol_FLAIR(:, :, z);
        % Normalizzazione volumetrica: usa min_vol e denom del volume intero
        img_gray = (double(slice) - min_vol) / denom;

        % Applicazione del filtro di denoising selezionato
        switch filtro.tipo
            case 'median'
                img_filtrata = medfilt2(img_gray, filtro.param);
            case 'gaussian'
                img_filtrata = imgaussfilt(img_gray, filtro.param);
            case 'wiener'
                img_filtrata = wiener2(img_gray, filtro.param);
            otherwise
                error('Filtro sconosciuto: %s', filtro.tipo);
        end

        img_processata = imresize(img_filtrata, dim_standard);

        % ---------------------------------------------------------
        % SEGMENTAZIONE CON LA CONFIG SINGOLA
        % ---------------------------------------------------------
        mask = seg_func(img_processata, config_singola, SE);

        % Accumulo TP/TN/FP/FN
        TP_tot = TP_tot + sum(mask(:) == 1 & gt_slice(:) == 1);
        TN_tot = TN_tot + sum(mask(:) == 0 & gt_slice(:) == 0);
        FP_tot = FP_tot + sum(mask(:) == 1 & gt_slice(:) == 0);
        FN_tot = FN_tot + sum(mask(:) == 0 & gt_slice(:) == 1);
    end

    % =================================================================
    % CALCOLO DELLE METRICHE PER IL PAZIENTE CORRENTE
    % =================================================================
    % Metriche calcolate sull'intero volume della ROI (approccio volumetrico).
    TP = TP_tot; TN = TN_tot; FP = FP_tot; FN = FN_tot;

    % Accuracy: sempre definita (denominatore = totale voxel della ROI).
    Dataset_Acc(paziente_idx) = (TP + TN) / (TP + TN + FP + FN);

    % Sensitivity: indefinita se la GT non contiene tumore nella ROI.
    if (TP + FN) > 0
        Dataset_Sens(paziente_idx) = TP / (TP + FN);
    else
        Dataset_Sens(paziente_idx) = 0;
    end

    % Specificity: TN+FP > 0 garantito dalla presenza di tessuto sano.
    Dataset_Spec(paziente_idx) = TN / (TN + FP);

    % Jaccard e F1: indefiniti se predizione e GT sono entrambe vuote.
    if (TP + FP + FN) > 0
        Dataset_Jacc(paziente_idx) = TP / (TP + FP + FN);
        Dataset_F1(paziente_idx)   = (2 * TP) / (2 * TP + FP + FN);
    else
        Dataset_Jacc(paziente_idx) = 0;
        Dataset_F1(paziente_idx)   = 0;
    end

end % Fine ciclo pazienti
end