# Recombine

Library Construction and Recombination Detection on coronaviruses, as presented in the paper "_Investigating Recombination Patterns in Coronaviruses Using Large-scale Genomic Data_".


**#Setting**

joblib==1.2.0

pandas==1.5.0

BioPython==1.78

matplotlib==2.2.2

blast-2.12.0


**#USAGE**

1. Run Recombine_preprocess.py for local library construction.
   
2. Run Recombine_fragment.py for recombination detection.
   
3. Run Recombine_visualize.py for visualized results.


**#NOTE**

0. BLAST should be installed, with blast-2.12.0 recommended.

1. For all .fasta files used in this project, the sequence must be in one line rather than in several lines, otherwise ERROR.
   
2. The recombine_lib/COV-reselect-Cluster7.fasta is the .fasta file for this project.

3. Recombine_spacetime_info_raw.csv includes the spatio-temporal information of filtered sequences for local library construction.

4. Recombine_FinalFrag.csv is the Example outputfile produced by Recombine_fragment.py.

5. Once the local library is constructed, it will no longer need re-constructing. However, changing the running platform requires library construction again.

6. Non-recombiant strains are not visualized.

7. Recombine_fragment.py has running modes in its .py code. 
Running_mode=1, fast detection; 2, regular detection (recommended); 3, precise detecion.
