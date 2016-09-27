.. _GPICommunicator:

GPICommunicator
===============

GPICommunicator is a derived Communicator class that implements
the communication routines via GPI.

Communicator Factory
^^^^^^^^^^^^^^^^^^^^

A GPI communicator is available via the factory as follows:

.. code-block:: c++

   CommunicatorPtr gpiComm = Communicator::getCommunicatorPtr( Communicator::GPI );

For portability, we recommend using in your applications the default communicator:

.. code-block:: c++

   CommunicatorPtr comm = Communicator::getCommunicatorPtr();

Running Parallel Programs with the GPI Communicator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``gaspi_run`` utility is used to start and run GPI-2 applications. 
A machine file with the hostnames of nodes where the
application will run, must be provided. A machine name can appear several
times in the file if more processes per node should be started.

The gaspi_run utility is invoked as follows:

.. code-block:: bash

    gaspi_run -m <machinefile> [OPTIONS] <path GASPI program>

The actual communicator can be specified at runtime.

.. code-block:: bash

    gaspi_run -m machines $PWD/lamaApplication.exe --SCAI_COMMUNICATOR=GPI

Some pitfalls regarding gaspi_run:

* The full path must be specified for the executable, i.e. the path to the program
  must exist on all nodes where the program should be started.
* If the executable is linked with dynamic libraries these libraries must be in the
  ``LD_LIBRARY_PATH``. This path must be set when an executable is started on other
  machines via ssh. Be careful: the ``gaspi_run`` utility gives no error message
  that it could not start the executable for this reason.
* In a script, you should add a ``sleep 1`` command between two immediate gaspi_run commands.

The most serious restriction for the current implementation of the utility is the fact
that it is not possible to pass environment variables with the start of the executable 
on other nodes, the values are only visible on the first processor.

* The use of command line arguments is absolutely recommned. In each application the
  command line arguments should be evaluated via the parseArg routine of the common::Settings
  class.

Environment variables used during the static initialization cannot be passed.

* ``SCAI_USE_MKL``, ``SCAI_USE_BLAS`` 

For the configuration of the logging please use the file ``${HOME}/.loggingrc``. For 
debugging of applicatins the following entries are very helpful:

.. code-block:: bash

  format = "#SCAI_COMM #date, #time #name @ #thread ( #func -> #file::#line ) #level #msg"
  flush  = 1

The format entry makes sure that the communicator appears in each logging entry. The flush
entry takes care that the logging entries from other processes are printed immeadiately and not
at the end of the run.

