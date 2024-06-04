# LinearSystemSolver

Progetto di Metodi del Calcolo Scientifico, A.A. 2023/2024.

## Descrizione generica

Si utilizzi un linguaggio di programmazione a vostra scelta: c++, fortran, java, pyton, etc. Lo scopo del progetto è quello di implementare una mini libreria che esegua i seguenti solutori iterativi, limitatamente al caso di matrici simmetriche e definite positive:

1. metodo di Jacobi;
2. metodo di Gauß-Seidel;
3. metodo del Gradiente;
4. metodo del Gradiente coniugato.

Per la gestione della struttura dati e le operazioni elementari fra matrici è richiesto di partire da una libreria _open-source_, come Eigen, Armadillo, blas/lapack. Oppure, qualora il linguaggio di programmazione lo permetta, utilizzare vettori e matrici già implementate al suo interno.

## Richieste sulla libreria

Sarà valutato positivamente che la libreria abbia un’architettura e sia ben strutturata, invece che sia una sequenza di funzioni indipendenti le une dalle altre. Inoltre, deve soddisfare i seguenti requisiti:

<ol type="a">
    <li>
        operare con ognuno dei metodo iterativi (1)-(4) sopra menzionati. A tal proposito la libreria scelta come base deve fornire <i>solamente</i> la struttura dati di matrici e vettori NON possono essere utilizzati i metodi relativi alla risoluzione dei sistemi lineari già implementati al suo interno. Per esempio la libreria Eigen contiene già al suo interno l’implementazione del metodo di Jacobi, tale metodo NON è valido ai fini della consegna (tuttavia questi metodi possono essere usati in fase di debug).
    </li>
    <li>
        i metodi iterativi devono partire dal vettore iniziale nullo (vettore con tutte le entrate pari a zero) e arrestarsi qualora la k-esima iterata $x^{(k)}$
        soddisfa
        $$\frac{||Ax^{(k)} − b||}{||b||} < tol$$
        con tol tolleranza assegnata dal utente. Oltre a questo controllo sulla soluzione $x^{(k)}$ i metodi iterativi dovranno avere un controllo sul numero massimo di iterazioni. Tutte le routine dovranno quindi arrestarsi (e segnalare di non essere giunte a convergenza) se
        $$k > maxIter$$
        dove <b>maxIter</b> è un <i>numero</i> molto elevato a vostra scelta (non inferiore
        a 20000).
    </li>
</ol>

La routine di cui sopra va poi applicata alle matrici sparse salvate nei file in allegato: `spa1.mtx`, `spa2.mtx`, `vem1.mtx` e `vem2.mtx` al variare del metodo iterativo di soluzione e considerando questi valori come tolleranza:
$$tol = [10−4, 10−6, 10−8, 10−10]$$

## Nota per il formato .mtx

Il formato “**.mtx**” è uno standard per lo scambio delle matrici in ambito numerico. È un file ASCII la cui descrizione è riportata nel seguente link: <https://math.nist.gov/MatrixMarket/formats.html>.  
Inoltre la numerazione delle entrate delle matrici parte da 1.
