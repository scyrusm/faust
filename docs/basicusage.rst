
===========
Basic Usage
===========

FAUST is mostly a python API. It is designed to compute effect sizes and p-values/q-values (see :doc:`nullhypothesis`) using gRNA-UMI counts.
These counts are most commonly enriched via polymerase chain reaction (PCR). We typically generate counts using `PoolQ3 <https://portals.broadinstitute.org/gpp/public/software/poolq>`_, from the `Broad Institute <https://www.broadinstitute.org/>`_.
It is also possible to generate these counts directly with FAUST.

The core functionality of FAUST may be found in the function :func:`faust.utilities.get_summary_df`.
This function expects a pandas dataframe `df` with the following example format:

======= === ======== ====== === === ====== ======
barcode UMI gene  input1 ln1 ln2 tumor1 tumor2
======= === ======== ====== === === ====== ======
TCGA... T.. gene1    532    300 251 10000  5238
ATCG... C.. gene2    102    10  11  12     23
GATC... G.. gene3    400    17  21  33     12
CGAT... A.. control1 310    177 212 120    242
======= === ======== ====== === === ====== ======

The next argument, `controls` should be a list of control targets.  
These generally correspond to gRNAs that target intergenic regions, or that target no site in the genome at all.

The next arguments, `inputs` and `outputs`, should be a list of columns of `df` that correspond to input and output sites, respectively.
FAUST will compute the ratio, for each gRNA-UMI, between the output and the input sites provided. 
The exact way this is done will depend on the argument ``input_type``.  
If ``input_type`` is 'single', FAUST will take the row-wise sum over all input columns in `inputs`; each entry in each column in `outputs` will then be divided (row-wise) by this sum to obtain the "factor of expansion" :math:`F_e` for each gRNA-UMI.
If ``input_type`` is 'matched', FAUST will compute this factor of expansion for matched elements in `inputs` and `outputs`. 
