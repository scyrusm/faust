
===========
Basic Usage
===========

FAUST is mostly a python API. It is designed to compute effect sizes and p-values/q-values (see :doc:`nullhypothesis`) using gRNA-UMI counts.
These counts are most commonly enriched via polymerase chain reaction (PCR). We typically generate counts using `PoolQ3 <https://portals.broadinstitute.org/gpp/public/software/poolq>`_, from the `Broad Institute <https://www.broadinstitute.org/>`_.
It is also possible to generate these counts directly with FAUST.

The core functionality of FAUST may be found in the function :func:`faust.utilities.get_summary_df`.
This function expects a pandas dataframe `df` with the following example format:

 .. csv-table:: example argument ``df`` for function :func:`faust.utilities.get_summary_df`
    :file: examplesummarydf.csv    
    :header-rows: 1

The next argument, `controls` should be a list of control targets.  
These generally correspond to gRNAs that target intergenic regions, or that target no site in the genome at all. In the table above, `controls` would take the value ``["control1","control2"]``

The next arguments, `inputs` and `outputs`, should be a list of columns of `df` that correspond to input and output sites, respectively.
FAUST will compute the ratio, for each gRNA-UMI, between the output and the input sites provided. 
The exact way this is done will depend on the argument ``input_type``.  
If ``input_type`` is 'single', FAUST will take the row-wise sum over all input columns in `inputs`; each entry in each column in `outputs` will then be divided (row-wise) by this sum to obtain the "factor of expansion" :math:`F_e` for each gRNA-UMI.
If ``input_type`` is 'matched', FAUST will compute this factor of expansion for matched elements in `inputs` and `outputs`. 

Let's suppose we want to test the null hypothesis :math:`H_0` that the factor of expansion :math:`F_e` between a common input aliquot and a particular lymph node is equally likely to be greater or lesser for gRNA-UMIs targeting gene1 vs. gRNA-UMIs targeting control loci.  
To do this, set ``input_type`` to be 'single', `controls` to be ``["control1","control2"]``, `inputs` to be ``["input1","input2"]``, and `outputs` to be ``["ln1","ln2"]``.  
FAUST will pool the counts for all the input aliquot measurements (we will test :math:`H_0` using a Mann-Whitney U test, so whether we sum or average these input counts won't affect our final result). 
FAUST will then evalute :math:`H_0` separately for ln1 and ln2.  

Let's now suppose we want to test the null hypothesis :math:`H_0` that the factor of expansion :math:`F_e` between a particular lymph node and a *matched* tumor is equally likely to be greater or lesser for gRNA-UMIs targeting gene1 vs. gRNA-UMIs targeting control loci.  
To do this, set ``input_type`` to be 'matched', `controls` to be ``["control1","control2"]``, `inputs` to be ``["ln1","ln2"]``, and `outputs` to be ``["tumor1","tumor2"]``.  
FAUST will then evalute :math:`H_0` between ln1 and tumor1, then ln2 and tumor2.
That is, the ordering of `inputs` and `outputs` matters, and they should be lists of the same length. 
Tried to run with ``input_type`` 'matched' with `inputs` and `outputs` of unequal lengths will raise an exception.
