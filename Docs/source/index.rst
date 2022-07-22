Sequence-based stability changes prediction
============================================================

.. image:: https://raw.githubusercontent.com/hurraygong/MU3DSP/master/pictures/Figure.2.png
   :width: 400px
   :align: left

MU3DSP, a residue level 3D structure-based prediction tool to assess single point mutation effects on protein thermodynamic stability and applying to dingle-domain monomeric proteins. Given protein sequence with single mutations as the input, the proposed model integrated both sequence level features of mutant residues and residue level mutation-based 3D structure features.


MU3DSP's highlights
^^^^^^^^^^^^^^^^^^^^^^^^^
- We propose a fast computational method, MU3DSP, based on LightGBM to predict stability changes upon single-residue mutations on proteins by fusing information from 3D structure profiles. The tertiary structure of the protein does not have to be available.
- MU3DSP uses a structure-based feature from either homology models of query variants if they are available or the annotated genomic variants database G2S (Genome to Structure).
- MU3DSP can achieve real-time prediction, and it only takes less than 1 minute on average to compute the same mutated position on one protein. A software tool is also provided to allow easy use of the tool.
- MU3DSP achieves state-of-the-art performance on two independent testing datasets. It is a reliable tool to assess both somatic and germline substitution mutations and assist in protein design.

Reference
^^^^^^^^^

ACKNOWLEDGEMENTS
^^^^^^^^^^^^^^^^

Support
^^^^^^^
Feel free to submit an `issue <https://github.com/hurraygong/MU3DSP/issues/new>`_
or send us an `email <gongjt057@nenu.edu.cn>`_.
Your help to improve MU3DSP is highly appreciated.


.. toctree::
   :caption: Main
   :maxdepth: 1
   :hidden:

   About
   Installation
   Datasets
   Quickstart
   Mutation structure preparation
   Results Analysis
   Casestudy
   Release
   References
   Contact


