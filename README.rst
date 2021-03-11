redGEM, redGEMX and lumpGEM
========

Reduction of human genome-scale models. Papers: 
* **redGEM**: Ataman, M., et al., "redGEM: Systematic reduction and analysis of genome-scale metabolic reconstructions for development of consistent core metabolic models". Plos Computational Biology, 2017. 13(7).
* **redGEMX**: Maria Masid, Meriç Ataman and Vassily Hatzimanikatis. "redHUMAN: analyzing human metabolism and growth media through systematic reductions of thermodynamically curated genome-scale models" 
* **lumpGEM**: Meric Ataman and Vassily Hatzimanikatis."lumpGEM: Systematic generation of subnetworks and elementally balanced lumped reactions for the biosynthesis of target metabolites". Plos Computational Biology, 2017. 13(7).



Requirements
------------

You will need to have `Git LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/redhgem.git /path/to/redgem
    cd /path/to/redgem
    git lfs install
    git lfs pull

The scripts have been developed with Matlab 2017b, and CPLEX 12.7 (freely downloadable with the `IBM Academic initiative <https://developer.ibm.com/academic/>`_), and successfully ran on several other versions of both softwares. However, it is important to respect the IBM compatibility specs sheets between Matlab, CPLEX, and the computer OS - available `on IBM's website <https://www.ibm.com/software/reports/compatibility/clarity/index.html>`_.

This module requires `matTFA <https://github.com/EPFL-LCSB/mattfa/>`_

Generating reduced models
-------------------------
1. Place the thermodynamic data for the corresponding orgnanism into the `matTFA thermoDatabases <https://github.com/EPFL-LCSB/matTFA/thermoDatabases>`_ folder.
2. Place the corresponding curated GEM into the `GEMs <https://github.com/EPFL-LCSB/redgem/GEMs>`_ folder.
3. Place the *get* file into the `runFileExample <https://github.com/EPFL-LCSB/redgem/runFileExample>`_  folder.
4. Run the *get* file


License
=======
The software in this repository is put under an APACHE licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/redgem/blob/master/LICENSE>`_ file for more details.

