% SEGMENTA_KMEANS – Segmentazione tramite K-Means clustering e raffinamento morfologico.
%
% La funzione implementa il modulo di segmentazione basato sull'algoritmo K-Means.
% Riceve una slice 2D già pre‑processata (normalizzata, filtrata e ridimensionata a
% 240×240) e restituisce una maschera binaria della regione tumorale (Whole Tumor),
% pronta per l'accumulo nella matrice di confusione volumetrica.
%
% INPUT
% img_processata : matrice double [240×240] normalizzata in [0, 1].
%                  Slice FLAIR dopo pre‑processing (Fase 2 del Main).
% num_cluster    : intero scalare in {2, 3, 4}. Numero di cluster K in cui
%                  suddividere i pixel dell'immagine.
% SE             : elemento strutturante morfologico (output di strel).
%                  Deve essere lo stesso usato nel Main (strel('disk', 5)).
%
% OUTPUT
% maschera_finale : matrice logica [240×240] con valori {0, 1}.
%                   1 = pixel classificato come tumore, 0 = tessuto sano.
%
% METODOLOGIA
% 1. Clustering K-Means dell'immagine in num_cluster gruppi, basato sulla sola
%    intensità del pixel (distanza euclidea). L'algoritmo viene eseguito 3 volte
%    (NumAttempts = 3) con inizializzazioni diverse e si conserva la soluzione a
%    minima varianza intra‑cluster.
% 2. Calcolo dell'intensità media di ciascun cluster e ordinamento decrescente.
% 3. Selezione dei due cluster a media più alta (o del solo cluster più luminoso
%    se K = 2). Questa logica ricalca la fusione delle due classi superiori di Otsu.
% 4. Costruzione della maschera binaria assegnando 1 ai pixel dei cluster selezionati.
% 5. Apertura morfologica (erosione → dilatazione) per rimuovere falsi positivi isolati.
% 6. Chiusura morfologica (dilatazione → erosione) per sigillare lacune interne e
%      regolarizzare i contorni.
%
% CONFRONTO CON OTSU
% Otsu opera sull'istogramma globale e determina soglie che massimizzano la varianza
% inter‑classe. K-Means minimizza invece la distanza intra‑cluster di ogni pixel dal
% proprio centroide, risultando più flessibile in presenza di distribuzioni multimodali
% o asimmetriche. La selezione dei cluster tumorali segue la stessa filosofia di Otsu
% (prendere le classi più luminose), garantendo un confronto equo.

function maschera_finale = segmenta_kmeans(img_processata, num_cluster, SE)

% =========================================================================
% CLUSTERING K-MEANS DELL'IMMAGINE
% =========================================================================
% imsegkmeans accetta immagini single o uint8. Convertiamo per compatibilità.
img_single = single(img_processata);

% Esecuzione del clustering con K = num_cluster.
% 'NumAttempts', 3 esegue l'intero algoritmo 3 volte con centroidi iniziali
% diversi (inizializzazione K-Means++) e restituisce la soluzione con la
% minima varianza intra‑cluster, riducendo il rischio di minimi locali.
label_map = imsegkmeans(img_single, num_cluster, 'NumAttempts', 3);

% =========================================================================
% IDENTIFICAZIONE DEI CLUSTER TUMORALI (PIÙ INTENSI)
% =========================================================================
% Calcola l'intensità media di ciascun cluster basandosi sull'immagine originale.
% I cluster con media più alta corrispondono alle regioni iperintense in FLAIR
% (edema e core tumorale).
intensita_media_cluster = zeros(1, num_cluster);
for k = 1:num_cluster
    % Estrae i pixel del cluster k dall'immagine originale e ne calcola la media
    intensita_media_cluster(k) = mean(img_processata(label_map == k));
end

% Ordina i cluster in ordine decrescente di intensità media.
[~, idx_ordinati] = sort(intensita_media_cluster, 'descend');

% Seleziona i cluster corrispondenti al tumore, con logica coerente a Otsu:
%   K = 2 → si prende il solo cluster più luminoso (altrimenti maschera tutta 1).
%   K ≥ 3 → si prendono i due cluster più luminosi (edema + core).
num_cluster_tumore = min(2, num_cluster - 1);
cluster_tumore = idx_ordinati(1:num_cluster_tumore);

% Costruisce la maschera binaria: 1 per i pixel nei cluster tumorali, 0 altrove.
% ismember verifica pixel per pixel se il label appartiene ai cluster selezionati.
maschera_binaria = ismember(label_map, cluster_tumore);

% =========================================================================
% APERTURA MORFOLOGICA (Opening) - Rimozione dei falsi positivi isolati
% =========================================================================
% L'apertura (erosione → dilatazione) elimina piccoli cluster sparsi di pixel
% che non appartengono alla massa tumorale principale.
img_opening = imopen(maschera_binaria, SE);

% =========================================================================
% CHIUSURA MORFOLOGICA (Closing) - Riempimento delle lacune interne
% =========================================================================
% La chiusura (dilatazione → erosione) colma i vuoti interni alla lesione
% e regolarizza la maschera, rendendola compatta.
maschera_finale = imclose(img_opening, SE);

end