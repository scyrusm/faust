
=========================
The FAUST null hypothesis
=========================

The FAUST algorithm wishes to test the following null hypothesis:

.. math::

    H_0: \mathrm{Prob(F_e > F_c) = Prob(F_c > F_e)}

where :math:`F_e` is a factor of expansion targeting an experimental guide (i.e. with a genuine genetic target), and :math:`F_c` is a factor of expansion targeting a control guide (i.e., targeting an intergenic region, or targeting nothing in the genome at all). 

Within FAUST, :math:`F` is calculated by computing a fold change of counts between an input sample and an output sample. 
