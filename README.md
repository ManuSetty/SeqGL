# SeqGL

SeqGL is a  group lasso-based algorithm to extract multiple transcription factor (TF) binding signals from ChIP- and DNase-seq profiles. Benchmarked on over 100 ChIP-seq experiments, SeqGL outperformed traditional motif discovery tools in discriminative accuracy and cofactor detection. SeqGL successfully scales to DNase-seq data, identifying a large multiplicity of TF signals confirmed by ChIP, and can be used with multitask training to learn genomic-context and cell-type specific TF signals.

This repository is an R package to run SeqGL. Please refer to the vignette for installation and usage instructions. This package has been tested with R 3.2. The package and the vignette are available to download from the "Releases" tab. 
This package is also hosted on <a href="https://bitbucket.org/leslielab/seqgl/wiki/Home"> Bitbucket </a>.

#### Citation

SeqGL manuscript is available from [PLoS Computational Biology](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004271). If you use SeqGL for your work, please cite our paper.

        @article{SeqGL_2015,
                title = {SeqGL Identifies Context-Dependent Binding Signals in Genome-Wide Regulatory Element Maps},
                author = {Manu Setty and Christina S Leslie},
                journal = {PLoS Computational Biology},
                year = {2015},
                month = {may},
                url = {http://dx.doi.org/10.1371/journal.pcbi.1004271},
                doi = {dx.doi.org/10.1371/journal.pcbi.1004271}
        }
