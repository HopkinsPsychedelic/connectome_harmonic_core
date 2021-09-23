.. include:: links.rst


------------
Installation
------------

To install CHAP, you can either:

* as recommended, install in a containerized envorinment 
* follow the guide to `Manual Installation Environment in Python (version 3.7 or above)`_ or

Containerized Environments (Docker and Singularity)
---------------------------------------------------
CHAP can be run through Docker or Singularity. For good information on running Docker and Singularity apps, `see this
<https://www.nipreps.org/apps/docker/#running-a-niprep-directly-interacting-with-the-docker-engine>`__. If using Singularity, pull the singularity image, using

singularity pull docker://winstonian3/connectome\_harmonic

Then follow instructions in Usage Notes to run your Singularity or Docker image

Manual Installation Environment in Python (version 3.7 or above)
----------------------------------------------------------------

.. warning::
   This is not the recommended method for installing CHAP! Refer to section on Containerized Environments.

* ``CHAP`` is written using Python 3.7 or above

First ensure that the `External Dependencies`_ have been installed on your system. 

Specifications for configuring the proper python environment can be found in `Dockerfile <https://github.com/hptaylor/connectome_harmonic_core/blob/master/Dockerfile>`_. It is recommended to install the python dependencies via `Anaconda`_. 

If you have set up your environment manually, you are ready to execute ``CHAP`` via the command line using the `entrypoint script`_ within the directory associated with the CHAP `Github repository`_. The command-line options for ``CHAP`` are documented in the :ref:`usage` section. 

``CHAP`` conforms to the BIDS_ recommended user interface formatting, specifically ::

    $ python entrypoint_script.py <input_bids_path> <derivatives_path> <analysis_level> <named_options>


External Dependencies
---------------------
These external packages are present in the containerized build of CHAP, but are not managed by Python packaging (PyPi), and must be installed on your system sparately. 

- FSL_ (version 6.0.3)

- ANTs_ (version 2.2.0 - NeuroDocker build)

- FreeSurfer_ (version 7.1.0)

- fMRIPrep_ (version 20.2.3)

- `Ciftify`_ (version 2.2.3)

- `connectome-workbench <https://www.humanconnectome.org/software/connectome-workbench>`_ (version Debian-1.3.2)

- `bids-validator <https://github.com/bids-standard/bids-validator>`_ (version 1.4.0)



