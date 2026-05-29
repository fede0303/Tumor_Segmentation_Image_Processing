# 🧠 Brain Tumor Segmentation (Image Processing)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Dataset](https://img.shields.io/badge/Dataset-Decathlon%20Task01-orange.svg)](http://medicaldecathlon.com/)

Questo progetto implementa e confronta diverse tecniche di **image processing** per la segmentazione automatica di tumori cerebrali (Whole Tumor) in risonanze magnetiche (MRI) volumetriche 3D (canale FLAIR), utilizzando il dataset **Task01_BrainTumour** (Medical Decathlon).

---

## 📌 Caratteristiche del Progetto

Il sistema è strutturato come un framework modulare per il confronto equo (*fair comparison*) di tre algoritmi classici di segmentazione:
1. **Multi-Otsu Thresholding**: Sogliatura globale ottimizzata per la separazione in classi d'intensità.
2. **K-Means Clustering**: Clustering non supervisionato dei pixel basato sulle caratteristiche di intensità.
3. **Marker-Controlled Watershed**: Segmentazione basata sulla trasformata spartiacque con controllo del gradiente e filtraggio dei marcatori.

### ⚙️ Pipeline di Elaborazione
Ogni paziente viene elaborato attraverso i seguenti step:
* **Normalizzazione Volumetrica**: Calcolo dei valori minimi e massimi sull'intero volume 3D (FLAIR) per preservare la coerenza di scala tra le slice, evitando le distorsioni tipiche della normalizzazione *slice-by-slice*.
* **Definizione della ROI (Region of Interest)**: Selezione automatica del 40% centrale delle slice del volume per escludere artefatti e zone non cerebrali.
* **Pre-processing (Denoising)**: Supporto a filtri Mediano, Gaussiano e Wiener con parametri ottimizzabili.
* **Grid Search (No Data Leakage)**: Ottimizzazione indipendente per ciascun metodo sul **Validation Set (80%)** per trovare la miglior combinazione di filtro e iperparametri.
* **Valutazione Finale**: Eseguita esclusivamente sul **Test Set (20%)** calcolando metriche volumetriche:
  $$\text{Accuracy, Sensitivity, Specificity, Jaccard Index (IoU), F1-Score (Dice Coefficient)}$$

---

## 📂 Struttura delle Cartelle

```bash
├── Main_Ottimizzato.m         # Script principale (orchestratore)
├── esegui_pipeline_metodo.m   # Incapsulamento della pipeline di segmentazione per paziente
├── grid_search_denoising.m    # Grid search per l'ottimizzazione degli iperparametri
├── segmenta_otsu.m            # Algoritmo di segmentazione Otsu
├── segmenta_kmeans.m          # Algoritmo di segmentazione K-Means
├── segmenta_watershed.m       # Algoritmo di segmentazione Watershed con marcatori
├── visualizza_sanity_check.m  # Funzione di visualizzazione e debug grafico
├── .gitignore                 # File per escludere file temporanei e dataset pesanti
└── README.md                  # Questo file
```

> [!NOTE]
> La cartella `data/` contenente il dataset originale `Task01_BrainTumour` è esclusa dal tracciamento Git tramite il file `.gitignore` per evitare il caricamento di file pesanti.

---

## 🚀 Come Iniziare

### 1. Prerequisiti
* **MATLAB** (consigliata versione R2022b o successiva)
* **Image Processing Toolbox**

### 2. Configurazione del Dataset
Scarica il dataset **Task01_BrainTumour** e posizionalo all'interno della cartella del progetto rispettando la seguente struttura:
```bash
data/
└── Task01_BrainTumour/
    ├── immaginiMRI/    # File .nii.gz delle immagini 4D (es. BRATS_001.nii.gz)
    └── GroundTruth/    # File .nii.gz delle label (es. BRATS_001.nii.gz)
```

### 3. Esecuzione
Apri MATLAB, naviga nella directory del progetto ed esegui lo script principale:
```matlab
run('Main_Ottimizzato.m')
```

Se `ESEGUI_GRID_SEARCH = true` lo script:
1. Eseguirà la ricerca a griglia sul validation set per trovare i migliori parametri di ciascun metodo.
2. Eseguirà la validazione finale sul test set stampando in console la tabella comparativa delle metriche.
3. Mostrerà un grafico di *sanity check* (controllo visivo) per il primo paziente del set.

---

## 📊 Metriche Valutate

Il modulo restituisce per ogni algoritmo le seguenti medie sui pazienti di test:
* **Accuracy (Acc)**: Frazione di voxel totali classificati correttamente.
* **Sensitivity (Sens)**: Capacità di identificare correttamente i voxel tumorali.
* **Specificity (Spec)**: Capacità di escludere i voxel di tessuto sano.
* **Jaccard (IoU)**: Indice di sovrapposizione tra la segmentazione e la Ground Truth.
* **F1-Score / Dice**: Coefficiente di similarità standard per la segmentazione medica.
