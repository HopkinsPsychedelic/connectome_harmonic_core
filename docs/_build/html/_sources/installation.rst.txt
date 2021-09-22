.. include:: links.rst

------------
Installation
------------

To install CHAP, you can either:

* as recommended, install in a containerized envorinment 
* follow the guide to `Manual Installation Environment in Python (version 3.7 or above)`_ or

Containerized Environments (Docker and Singularity)
--------------------------------------------------
CHAP can be run through Docker or Singularity. For good information on running Docker and Singularity apps, `see this
<https://www.nipreps.org/apps/docker/#running-a-niprep-directly-interacting-with-the-docker-engine>`__. If using Singularity, pull the singularity image, using

singularity pull docker://winstonian3/connectome\_harmonic

The follow instructions in _Usage to run your singulari



Manual Installation Environment in Python (version 3.7 or above)
----------------------------------------------------------------

.. warning::
   This is not the recommended method for installing CHAP! Refer to section on Containerized Environments.

External Dependencies
---------------------

* ``CHAP`` is written using Python 3.7 or above

If you have set up your environment manually (bullet one), you are ready to execute ``CHAP`` via the command line using the `entrypoint script`_ within the directory associated with the CHAP `Github repository`_. The command-line options for ``CHAP`` are documented in the :ref:`usage` section. 

``CHAP`` conforms to the BIDS_ recommended user interface formatting, specifically

::$ python entrypoint_script.py <input_bids_path> <derivatives_path> <analysis_level> <named_options> (TODO: check if this is correct)

