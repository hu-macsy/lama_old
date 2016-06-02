.. _CommunicationPlan:

CommunicationPlan
=================

A communication plan is one central data structure used in the dmemo project.
It specifies for each processor the data that is required from other processors
and the data to be sent to other processors. The sent and received data is specified
by a set of indexes, therefore one communication plan can be applied to different 
arrays where the same communication pattern is required.

.. code-block:: c++

    0: needs { 1: { 3, 5 }, 2: { 8, 9 } }
    1: needs { 0: { 0, 2 }, 2: { 7, 9 } }
    2: needs { 0: { 1, 2 }, 1: { 4, 6 } }

By using an all-to-all communication scheme, all processors will know the
indexes for which data has to be sent to other processors.

.. code-block:: c++

    0: provides { 1: { 3, 5 }, 2: { 8, 9 } }
    1: provides { 0: { 0, 2 }, 2: { 7, 9 } }
    2: provides { 0: { 1, 2 }, 1: { 4, 6 } }

