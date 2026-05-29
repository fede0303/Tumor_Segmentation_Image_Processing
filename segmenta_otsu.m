% SEGMENTA_OTSU – Segmentazione tramite Otsu Multi‑livello e raffinamento morfologico.
%
%   La funzione implementa il modulo di segmentazione basato sull'algoritmo di Otsu
%   multi‑livello. Riceve una slice 2D già pre‑processata (normalizzata, filtrata e
%   ridimensionata a 240×240) e restituisce una maschera binaria della regione tumorale
%   (Whole Tumor), pronta per l'accumulo nella matrice di confusione volumetrica.
%
%   INPUT
%   img_processata : matrice double [240×240] normalizzata in [0, 1].
%                    Slice FLAIR dopo pre‑processing (Fase 2 del Main).
%   num_soglie     : intero scalare in {2, 3, 4}. Numero di soglie per multithresh.
%                    Produce (num_soglie + 1) classi di intensità.
%   SE             : elemento strutturante morfologico (output di strel).
%                    Deve essere lo stesso usato nel Main (strel('disk', 5)).
%
%   OUTPUT
%   maschera_finale : matrice logica [240×240] con valori {0, 1}.
%                     1 = pixel classificato come tumore, 0 = tessuto sano.
%
%   METODOLOGIA
%   1. Calcolo delle soglie ottime tramite multithresh (Otsu multi‑livello).
%   2. Quantizzazione dell'immagine in (num_soglie + 1) livelli discreti.
%   3. Costruzione della maschera binaria selezionando le due classi a più alta
%      intensità (condizione img_quant >= num_soglie). Questa scelta consente di
%      catturare sia l'edema iperintenso sia il core tumorale, evitando lacune
%      interne che la sola classe apicale produrrebbe in presenza di necrosi.
%   4. Apertura morfologica (erosione → dilatazione) per rimuovere falsi positivi
%      isolati di dimensione inferiore al raggio dell'elemento strutturante.
%   5. Chiusura morfologica (dilatazione → erosione) per sigillare i vuoti interni
%      e regolarizzare i contorni della maschera.

function maschera_finale = segmenta_otsu(img_processata, num_soglie, SE)

% =========================================================================
% CALCOLO DELLE SOGLIE OTTIME E QUANTIZZAZIONE
% =========================================================================
% multithresh determina le N soglie che massimizzano la varianza inter‑classe
% sull'istogramma dell'immagine (criterio di Otsu).

soglie = multithresh(img_processata, num_soglie);

% imquantize assegna a ogni pixel un livello intero da 1 (più scuro) a
% (num_soglie + 1) (più luminoso), in base alle soglie calcolate.
img_quant = imquantize(img_processata, soglie);

% =========================================================================
% COSTRUZIONE DELLA MASCHERA BINARIA (Unione delle due classi con più alta intensità)
% =========================================================================
% La condizione (img_quant >= num_soglie) seleziona i livelli N e N+1, ovvero le
% due fasce più iperintense della sequenza FLAIR. Questa fusione include sia
% l'edema sia il core tumorale, garantendo la copertura del "Whole Tumor".
maschera_binaria = (img_quant >= num_soglie);

% =========================================================================
% APERTURA MORFOLOGICA (Opening) - Rimozione dei falsi positivi isolati
% =========================================================================
% L'apertura (erosione seguita da dilatazione) elimina i piccoli cluster di pixel
% sparsi (rumore), senza alterare la geometria della massa tumorale principale.
img_opening = imopen(maschera_binaria, SE);

% =========================================================================
% CHIUSURA MORFOLOGICA (Closing) - Riempimento delle lacune interne
% =========================================================================
% La chiusura (dilatazione seguita da erosione) colma i vuoti interni alla 
% lesione e rende la maschera compatta e coerente con la morfologia tumorale 
% attesa.
maschera_finale = imclose(img_opening, SE);
end