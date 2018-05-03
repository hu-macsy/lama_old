.. _CyclicPartitioning:

Cyclic Partitioning
===================

Cyclic partitioning is another rather simple implementation of partitioning that
might achieve a better load balancing than the block distribution in cases where
the number of entries in the matrix differ much.

This class is mainly provided for comparison, but might also be used as a
fallback if no graph partitioning software is available.
