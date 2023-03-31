
========
Overview
========

.. figure:: https://upload.wikimedia.org/wikipedia/commons/4/40/Jean_Paul_Laurens_-_Dr._Fausto.jpg
   :alt: Hmm, surely one of these genes does SOMETHING... 

FAUST:  The FITS Algorithm for UMI-based Screens of T cells

FITS:  The Framework for In Vivo T cell Screens

FITS is a framework for designing in vivo T cell screens, and FAUST is a computational approach that goes with it. Briefly, if you are designing a CRISPR screen, you have to decide between:
  - How many genes you independently target (panel size)
  - How many biological replicates to use
  - What sort of type I and type II errors you are willing to tolerate
  - How many cells to transfer into recipient animals

These tradeoffs are considerably more pressing in the in vivo setting than in vitro. Conversely, the in vivo setting can help resolve context-specific gene perturbation function (for example, whether a gene perturbation plays one role in a lymph node, and another in a tumor). 
