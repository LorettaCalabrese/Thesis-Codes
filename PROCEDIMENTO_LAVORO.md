Il lavoro si è suddiviso in due fasi: la prima fase in MATLAB di processamento dei segnali elettromiografici con conseguente creazione dei due dataset contenenti i dati utili
e una seconda fase in Python per l'implementazione di un algoritmo di deep learning per classificazione binaria e multiclasse di patologia (paraparesi spastica ereditaria). 

- I file nominati 'A' sono riferiti agli script MATLAB per andare a costruire il dataset per successiva classificazione binaria, costuito pertanto dalla suddivisione in sotto fasi del cammino:
    - A1:preprocessamento dei segnali elettromiografici con successiva creazione per ciascun paziente di singoli file excel per poter registrare gli inviluppi medi e il numero di campioni e singoli file csv contenenti i segnali EMG corrispondente a ciascuna sottofase del ciclo del passo;lo stesso codice, ma inserendo i dati dei soggetti paraparetici è stato ripetuto per quest'altra famiglia di dati di interesse. Si è deciso di creare due codici identici, ma separati, per comodità di gestione dei dati e dei file excel e csv creati;
    - A2:questo codice prevede la lettura di tutti i dati contenuti nei file excel relativi ai normativi (e un codice uguale per i paraparetici), con conseguente realizzaizone grafica degli inviluppi medi, con deviazione standard, per poter verificare il corretto valore assegnato alla soglia di attivazione muscolare e della fase di prepocessamento del segnale in generale, verificando inoltre il corretto andamento dei segnali muscolari;
    - A3:creazione del dataset finale contenente i segnali muscolari di ciascun paziente, nelle varie sottofasi del cammino, a partire dai singoli file csv creati al punto A1 e A2. 

- I file nominati 'B' sono riferiti invece agli script MATLAB per andare a costruire il dataset per successiva classificazione multiclasse, costutito pertanto dalle percentuali di attivaizoni muscolari per ogni   sottofase dle cammino:
    - B1:preprocessamento dei segnali elettromiografici con successiva creazione di singoli file csv per ciascun paziente, contente le percentuali di attvazione di ciascun muscolo, in ogni sottofase;
    - B2:apertura di ciascun file csv creato in precedenza con codice B1 e creazione del dataset finale.

- I file nominati 'C' sono riferiti agli script Pyhton utilizzati per :
  - C1: preprocessare entrambi i dataset;
  - C2: creare il modello LSTM per entrambi i dataset utilizzati; in particolare i file sono due ma prevedono la stessa struttura, perchè sono per i due diversi dataset; lo script relativo al dataset comprendente le sottofasi del passo è stato ripetuto per ogni sottofase specifica (i cui dati sono stati inclusi in dataset separati), per poter estrarre gli iperparametri necessari per la fase di ottimizzazioe del LSTM successiva. 
  - C3: modello LSTM con iperparametri ottimizzati (nel seguente codice è riportato il modello con gli iperparametri estratti per la sottofase di heel contact, solo come esempio di struttura);
  - C4: valutazione delle metriche di performance e algoritmo di permutation importance, ripetuto per entrambi i dataset. 
