% VISUALIZZA_SANITY_CHECK – Verifica qualitativa della maschera segmentata (primo paziente)
%
%   INPUT
%   vol_FLAIR           : volume 3D FLAIR del paziente (X × Y × Z).
%   slice_equatore      : indice della slice da visualizzare (slice centrale della ROI).
%   metodo_nome         : stringa con il nome del metodo di segmentazione (es. 'Otsu').
%   seg_func            : function handle alla funzione di segmentazione.
%   config_singola      : scalare. La singola configurazione ottimale del metodo.
%   SE                  : elemento strutturante morfologico (strel).
%   dim_standard        : dimensioni di output [240 240] a cui viene portata la slice.
%   filtro              : struct con campi 'tipo' e 'param' per il filtro di denoising.
%   min_vol             : minimo volumetrico pre-calcolato in esegui_pipeline_metodo.
%   denom               : denominatore volumetrico (max_vol - min_vol) pre-calcolato.
%
%   DESCRIZIONE
%   La funzione riproduce fedelmente la stessa pipeline di pre-processing usata
%   in esegui_pipeline_metodo (normalizzazione volumetrica → filtro → resize)
%   e chiama la funzione di segmentazione con la config_singola ottimale per
%   generare e visualizzare la maschera binaria risultante.
%   Mostra una figura con due subplot:
%       - l'immagine FLAIR filtrata (normalizzazione volumetrica + filtro selezionato);
%       - la maschera ottenuta con la configurazione ottimale del metodo.
%   Questo controllo qualitativo non influisce sulle metriche e viene eseguito
%   unicamente sul primo paziente del batch.

function visualizza_sanity_check(vol_FLAIR, slice_equatore, metodo_nome, ...
    seg_func, config_singola, SE, dim_standard, filtro, min_vol, denom)

% =========================================================================
% PRE-PROCESSING DELLA SLICE DI CONTROLLO
% =========================================================================
% Si riproduce fedelmente la stessa pipeline applicata in esegui_pipeline_metodo:
%   1. Normalizzazione VOLUMETRICA: usa min_vol e denom calcolati sull'intero
%      volume, coerentemente con la pipeline principale. Non si usa mat2gray
%      (che opererebbe slice per slice con scala indipendente).
%   2. Filtraggio con il filtro ottimale selezionato dalla grid search.
%   3. Ridimensionamento alle dimensioni standard (240×240).

slice_raw    = vol_FLAIR(:, :, slice_equatore);
img_gray     = (double(slice_raw) - min_vol) / denom;

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

img_finale = imresize(img_filtrata, dim_standard);

% =========================================================================
% GENERAZIONE DELL'ETICHETTA DI CONFIGURAZIONE
% =========================================================================
switch metodo_nome
    case 'Otsu'
        cfg_label = sprintf('Otsu %d soglie', config_singola);
    case 'K-Means'
        cfg_label = sprintf('K-Means k=%d', config_singola);
    case 'Watershed'
        cfg_label = sprintf('Watershed p=%.2f', config_singola);
    otherwise
        cfg_label = sprintf('Config %g', config_singola);
end

% =========================================================================
% CREAZIONE DELLA FIGURA
% =========================================================================
figure('Name', sprintf('Sanity Check: %s (Slice %d)', metodo_nome, slice_equatore), ...
    'Position', [100, 200, 900, 420]);

% Subplot 1: immagine FLAIR dopo pre-processing
subplot(1, 2, 1);
imshow(img_finale, []);
title(sprintf('FLAIR Filtrata\n%s — Filtro: %s', metodo_nome, filtro.nome), ...
    'FontSize', 10);

% Subplot 2: maschera della configurazione ottimale
maschera = seg_func(img_finale, config_singola, SE);
subplot(1, 2, 2);
imshow(maschera);
title(sprintf('Maschera Segmentata\n%s', cfg_label), 'FontSize', 10);

end