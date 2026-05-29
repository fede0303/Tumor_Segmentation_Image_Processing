% SEGMENTA_WATERSHED – Segmentazione tramite trasformata Watershed marker-controlled.
%
%   La funzione implementa il modulo di segmentazione basato sull'algoritmo Watershed.
%   Riceve una slice 2D già pre‑processata (normalizzata, filtrata e ridimensionata a
%   240×240) e restituisce una maschera binaria della regione tumorale (Whole Tumor),
%   pronta per l'accumulo nella matrice di confusione volumetrica.
%
%   INPUT
%   img_processata : matrice double [240×240] normalizzata in [0, 1].
%                    Slice FLAIR dopo pre‑processing (Fase 2 del Main).
%   soglia_marker  : scalare in (0, 1). Percentile superiore dell'intensità usato 
%                    per la definizione dinamica dei semi di foreground.
%   SE_grad        : elemento strutturante morfologico (disco, r=1) allocato 
%                    esternamente, impiegato per il calcolo del gradiente.
%   SE_smooth      : elemento strutturante morfologico (disco, r=3) allocato 
%                    esternamente, impiegato per la regolarizzazione dei massimi.
%
%   OUTPUT
%   maschera_finale : matrice logica [240×240] con valori {0, 1}.
%                     1 = pixel classificato come tumore, 0 = tessuto sano.
%
%   METODOLOGIA
%   1. Calcolo del gradiente morfologico (dilatazione − erosione) come superficie.
%   2. Marcatori di foreground (tumore): apertura/chiusura morfologica, estrazione
%      dei massimi regionali e filtraggio stringente tramite il percentile superiore,
%      calcolato ESCLUSIVAMENTE sul parenchima cerebrale.
%   3. Marcatori di background (sfondo): calcolo adattivo della mediana del parenchima 
%      sano, combinato topologicamente con il volume extra-cranico (aria) per
%      inibire il fenomeno di over-segmentation e arginare l'allagamento.
%   4. Modifica topografica imponendo i marcatori calcolati (imimposemin).
%   5. Trasformata Watershed sul gradiente modificato e isolamento del bacino bersaglio.
%   6. Riempimento dei vuoti interni tramite operazione morfologica di closing/fill.
%
%   CONFRONTO CON OTSU E K-MEANS
%   A differenza di Otsu e K-Means, che partizionano lo spazio puramente 
%   sulla base della densità globale delle intensità, il Watershed impiega 
%   vincoli topologici locali (minimi imposti sul gradiente). L'inclusione 
%   della regione extra-cranica nei marker fornisce un confine geometrico 
%   intrinseco, essenziale per un'equa comparazione delle performance.

function maschera_finale = segmenta_watershed(img_processata, soglia_marker, SE_grad, SE_smooth)

% =========================================================================
% 1. CALCOLO DEL GRADIENTE MORFOLOGICO
% =========================================================================
% Il gradiente morfologico definisce la mappa dei rilievi (bordi topografici). 
% L'operazione sfrutta l'elemento strutturante pre-allocato (SE_grad) e passato
% dall'orchestratore, annullando l'overhead di instanziazione durante 
% l'iterazione volumetrica sulle slice.
gradiente = imdilate(img_processata, SE_grad) - imerode(img_processata, SE_grad);

% =========================================================================
% MASCHERA CEREBRALE (esclude l'aria fuori dal cranio)
% =========================================================================
% Si definisce una maschera che isola il parenchima cerebrale, escludendo i voxel
% di background (aria) che altrimenti falserebbero i percentili. Questa maschera
% sarà utilizzata sia per il calcolo del percentile tumorale (foreground) sia per
% il calcolo del percentile di sfondo, garantendo coerenza statistica e robustezza.
soglia_aria = mean(img_processata(:)) * 0.5;
maschera_cervello = img_processata > soglia_aria; 

% =========================================================================
% 2. DEFINIZIONE DEI MARCATORI DI FOREGROUND (TUMORE)
% =========================================================================
% La regolarizzazione morfologica preliminare attenua le alte frequenze 
% non significative e il rumore residuo, avvalendosi di SE_smooth.
img_smooth = imclose(imopen(img_processata, SE_smooth), SE_smooth);

% Estrazione dei massimi regionali (potenziali centroidi tumorali).
marker_tumore = imregionalmax(img_smooth);

% I massimi locali vengono filtrati trattenendo esclusivamente le intensità
% superiori al percentile, operazione statisticamente idonea per isolare 
% la natura iper-intensa dei gliomi in sequenze FLAIR.
% Il percentile è calcolato ESCLUSIVAMENTE sui voxel cerebrali (maschera_cervello)
% per evitare che la quantità variabile di aria esterna distorca la soglia.
soglia_tumore = prctile(img_processata(maschera_cervello), soglia_marker * 100);
marker_tumore = marker_tumore & (img_processata >= soglia_tumore);

% Riduzione dei cluster connessi a seed-point discreti (singoli pixel).
marker_tumore = bwmorph(marker_tumore, 'shrink', Inf);

% =========================================================================
% 3. DEFINIZIONE DEI MARCATORI DI BACKGROUND (TESSUTO SANO E REGIONE ESTERNA)
% =========================================================================
% Per isolare spazialmente l'evoluzione della trasformata ed inibire
% espansioni extra-craniche del bacino tumorale, il background è un dominio composito.

% A. Volume Cerebrale Sano: Mediana adattiva sul solo parenchima.
intensita_cervello = img_processata(maschera_cervello);
soglia_sfondo = prctile(intensita_cervello, 50);
marker_sfondo_cervello = maschera_cervello & (img_processata <= soglia_sfondo);

% B. Volume Extra-Cranico (Aria): Vincolo topologico confinante.
maschera_aria = ~maschera_cervello;

% Combinazione logica dei marker di confinamento.
marker_sfondo_totale = marker_sfondo_cervello | maschera_aria;
marker_sfondo = bwmorph(marker_sfondo_totale, 'shrink', Inf);

% =========================================================================
% 4. IMPOSIZIONE DEI MARCATORI SUL GRADIENTE (imimposemin)
% =========================================================================
% Mappatura logica: 1 = Foreground (tumore), 2 = Background (sfondo+aria).
% In caso di sovrapposizione spaziale dei marker (evento eccezionale),
% l'assegnazione sequenziale privilegia il tumore (riga 108 sovrascrive riga 107).
% imimposemin riceve una maschera booleana unificata: la separazione
% foreground/background è garantita dalla loro disgiunzione spaziale, non
% dall'API di imimposemin.
marker_combined = zeros(size(img_processata), 'uint8');
marker_combined(marker_sfondo) = 2;
marker_combined(marker_tumore) = 1;

% Forzatura dei minimi sui marker stabiliti.
gradiente_marked = imimposemin(gradiente, marker_combined > 0);

% =========================================================================
% 5. TRASFORMATA WATERSHED
% =========================================================================
label_map = watershed(gradiente_marked);

% =========================================================================
% 6. ESTRAZIONE DELLA MASCHERA TUMORALE E RIEMPIMENTO BUCHI (Closing)
% =========================================================================
% Identificazione dell'etichetta correlata ai marker di foreground.
idx_tumore = unique(label_map(marker_tumore));
idx_tumore = idx_tumore(idx_tumore > 0); % Esclusione perimetro di displuvio

maschera_binaria = ismember(label_map, idx_tumore);

% Regolarizzazione delle lacune indotte da disomogeneità tissutali o necrosi.
maschera_finale = imfill(maschera_binaria, 'holes');

end