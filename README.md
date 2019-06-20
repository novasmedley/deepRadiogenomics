# Radiogenomic neural networks
`deepRadiogenomics` contains the source code to analyses in the paper entitled [TBA].

- [repo contents](###repo-contents)
- [data](###data)
- [system requirements](###system-requirements) 
- [usage](###usage)
- [citation](###citation)

### Repo contents
- [bash](./bash) cmd line scripts for model training, e.g., grid search
- [demo_data](./demo) toy data
- [R](./R): scripts for many post-modeling analyses, association testing, etc.

* general modeling functions:
    * `neuralnet.py`
    * `other_models_utils.py`
    * `bootstrap.py`
    * `custom_callbacks.py`
    * `sparse_x_generator.py`

* glioblastoma (gbm) specific functions:
    * training:
        * `setup.py`
        * `train_gene_ae.py` deep transcriptomic autoencoder
        * `train_nn.py` supervised radiogenomic neural network
        * `train_others.py` comparative models (logit, gbt, rf, svm)
    * extract radiogenomic associations
        * `gene_masking.py`
        * `get_masking.py`
        * `gene_saliency.py`
        * `get_saliency.py`
    * misc.
        * `parse_cv.py`
        * `demo_gene_masking.py`
        * `demo_gene_saliency.py`
   * all others, see [R](./R)

### Data

All data was originally taken from public repositories, where identifiable information was scrubbed.

* Transcriptomic data was downloaded from the legacy version of The Cancer Genome Archive (TCGA). 

* Imaging studies were download from from The Cancer Imaging Archive (TCIA). Vasarit traits were annotated by Dr. Suzie El-Saden and required pre-operative magnetic resonance imaging studies. 
    * Vasari guidelines
    * Annotation Google Form
    
Data files are available for download at XX in supplemental materials or at external links:

* training data
    * `gene_expression.txt` - gene expression profiles
    * `vasari_annotations.csv` - imaging traits
    * `nationwidechildrens.org_clinical_patient_gbm.txt` - TCGA-GBM clinical traits
* gene sets
    * `TCGA_unified_CORE_ClaNC840.txt` - Verhaak gene sets, [paper](https://doi.org/10.1016/j.ccr.2009.12.020)
    * `gene_sets_Puchalski` - Puchalski gene sets, [paper](https://doi.org/10.1126/science.aaf2666)
    * `msigdb_v6.2_GMTs` - MSigDB gene sets, available [here](http://software.broadinstitute.org/gsea/downloads.jsp)

For more details, see our paper.

### Install

* Neural networks were trained on Amazon Web Services using Deep Learning AMI with Ubuntu 16.04.4 LTS and the `tensorflow_p36` environment. All other classifiers were implemented on an Ubuntu 18.04.1 LST machine.

    * Check out AWS's environment documentation at XX.
    * Python 3.6 dependencies (training of comparative models, gene masking, gene saliency):
    
        ```
        keras 2.2.2
        keras-vis 0.4.1
        numpy 1.14.3
        pandas 0.23.0
        scikit-learn 0.20.0
        scipy 1.1.0
        seaborn 0.8.1
        tensorflow 1.10.0
        xgboost 0.80
        ```
    
    * R 3.4.4 dependencies (gene set enrichment analysis, but mostly figure generation):
        ```
        awtools 0.1.0
        broom 0.5.1
        cowplot 0.9.3
        data.table 1.11.8
        doParallel 1.0.14
        dplyr 0.7.6
        egg 0.4.0
        fgsea 1.4.1
        foreach 1.4.4
        ggrepel 0.8.0
        ggplot2 3.1.0
        grid 3.4.4
        gridExtra 2.3
        ggpubr 0.2
        pheatmap 1.0.12
        plyr 1.8.4
        qvalue 2.10.1
        rcartocolor 1.0.0
        RColorBrewer 1.1
        reshape2 1.4.3
        scales 1.0.0
        tidyr 0.8.1
        tidyverse 1.2.1
        viridis 0.5.1
        wesanderson 0.3.6
        ```
* install from Github using `git`:

    ```
    git clone https://github.com/novasmedley/deepRadiogenomics.git
    ```
### Usage

* wget data files
* 

#### demos
Demos were run using demo data, a small subset of the published dataset, on Ubuntu 18.04.1 LTS with 15.5 GB memory.

**Neural network pipeline**:

1. Train gene expression autoencoder (ae) - cross-validation(cv), 15 secs:

    ```
    $ python3 train_gene_ae.py --exp ae_test --dir ../demo --data ../demo \
    --label autoencoder --predType regression  --loss mae --opt Nadam --act tanh \
    --h1 200 --h2 100 --h3 50 --epoch 2 --folds 2 --patience 2
    ```

1. (optional) Parse cv results: `$ python3 parse_cv.py --dir ../demo/ae_cv --model nn`
        
1. Retrain ae - 15 secs: 
    
    run `train_gene_ae.py` using the same parameters as above but set `--retrain` to `1` 
 
1. Train radiogenomic model - cv, 21 secs:

    ```
    $ python3 train_nn.py ---exp nn_test --dir ../demo --data ../demo \
    --pretrain ../demo/ae_retrain_test/autoencoder/neuralnets/200_100_50_0_0_tanh_decay_0_drop_0_opt_Nadam_loss_mae_bat_10_eph_2 \
    --label f5 --opt Nadam --act tanh \
    --h1 200 --h2 100 --h3 50 --epoch 2 --folds 2 --patience 2 --freeze 0 --num_ae_layers 3
    ```

1. (optional) Parse cv results: `$ python3 parse_cv.py --dir ../demo/nn_cv --model nn`

1. Retrain radiogenomic model - 16 secs: 
    
    run `train_nn.py` using the same parameters as above but set `--retrain` to `1` 

1. Gene masking - 11 secs

    `$ python3 demo_gene_masking.py --label f5 --geneset verhaak --cpus 7`
    
1. Gene saliency - 19 secs
    
    `$ python3 demo_gene_saliency.py`


**Train other models**:

1. Fit logit with l1 regularization  - 1.5 mins, fitting 1000 hyperparameters
    ```
    $ train_others.py --exp other_cv --dir ../demo --data ../demo \
    --dataType vasari --predType binaryClass --label f5 --model logit1 --folds 2 --cpus 7
    ```

2. Parse cv results:

    `$ python3 parse_cv.py --dir ../demo/other_cv --model other`

### Citation

If you want to cite this work, please cite the paper:

```
```

and the repo:

```
@misc{smedleyRadiogenomics,
  title={deepRadiogenomics},
  author={Smedley, Nova F},
  year={2019},
  publisher={GitHub},
  howpublished={\url{https://github.com/novasmedley/deepRadiogenomics}},
}
```