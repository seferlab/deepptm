# DeepPTM: Protein Post-translational Modification Prediction from Protein Sequences by Combining BERT with Vision Transformers

DeepPTM: Protein Post-translational Modification Prediction from
Protein Sequences by Combining BERT with Vision Transformers

Prepared by Necla Nisa Soylu and Emre Sefer

We have obtained protein sequences and the corresponding post-translational modifications from CPLM database (http://cplm.biocuckoo.cn/) and focused on 4 important post-translational modifications over lysine in our experiments: Succinylation, Glutarylation, Crotonylation, and Glycation. Even though our main focus was to predict the modification sites on Homo sapiens, we also tested DeepPtm on Mus musculus and Saccharomyces cerevisiae. List of the datasets are given below:

- Succinylation (Homo Sapiens)
- Succinylation (Mus musculus)
- Succinylation (S. cerevisiae)
- Ubiquitination (Homo Sapiens)
- Crotonylation (Homo Sapiens)
- Glycation (Homo Sapiens)

Proteins were converted into peptides with a fixed length of 21, where lysine is at the central residue (modification either occurs or not), and 10 residues exist at each upstream and downstream. In order to remove the redundancy in our datasets, we have removed the homologous sequences exhibiting more than 40% similarity by using CD-HIT tool (https://sites.google.com/view/cd-hit). You can find our datasets in the following files:

- Suc_Human_cdhit40_Equal.xlsx
- Suc_Mouse_cdhit40_Equal.xlsx
- Suc_yeast_cdhit40_Equal.xlsx
- Ubi_Human_dhit40_Equal.xlsx
- Crot_Human_cdhit40_Equal.xlsx
- Gly_Human_dhit40_Equal.xlsx

We have trained the baseline machine learning approaches by converting each amino acid into numeric values via one-hot encoding strategy and each of the 20 amino acid is mapped from 20×1 vectors [1, 0, . . . , 0] to [0, 0, . . . , 1] respectively. As a result, we have obtained vectors with lengths 420, for each 21 nucleotide long protein sequence. We used Random Forest (RF) and Extreme Gradient Boosting (XGBoost) as our baseline models by obtaining optimal model parameters with using hyperparameter tuning. You can find our baseline models in the following files:

- RF ?
- XGb ?

In terms of deep learning approaches, each amino acids were considered as a single word in human language and corresponding vector embeddings were created by using one of the well-known natural language processing model BERT (Bidirectional Encoder Representations from Transformers). Following files contain the vector creation process for each peptides:

- Suc_Human_Equal_Vector_Creation.ipynb
- Suc_Mouse_Equal_Vector_Creation.ipynb
- Suc_Yeast_Equal_Vector_Creation.ipynb
- Ubi_Human_Vector_Creation.ipynb
- Crot_Human_Equal_Vector_Creation.ipynb
- Gly_Human_Equal_Vector_Creation.ipynb

Then we compared the performance of BERT + Random Forest, BERT + XGBoost, BERT + 1D CNN and BERT + 2D CNN with our proposed method DeepPTM: BERT + ViT. You can find our codes in the following files:

- BERT + RF ?
- BERT + XGB ?
- Suc_Human_1D_2D_CNN_and_ViT.ipynb
- Suc_Mouse_1D_2D_CNN_and_ViT.ipynb
- Suc_Yeast_1D_2D_CNN_and_ViT.ipynb
- Ubi_Human_1D_2D_CNN_and_ViT.ipynb 
- Crot_Human_1D_2D_CNN_and_ViT.ipynb
- Gly_Human_1D_2D_CNN_and_ViT.ipynb

Each file contains following models inside of it:

BERT + Random Forest
BERT + XGBoost
BERT + 1D CNN
BERT + 2D CNN
DeepPTM: BERT + ViT

Cross species and cross-modification were also implemented and can be found in the following file:

- Cross-Species for Succinylation and Cross-Modification.ipynb
