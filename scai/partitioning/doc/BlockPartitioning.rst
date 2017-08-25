.. _BlockPartitioning:

Block Partitioning
==================

Block partitioning is one rather trivial implementation of partitioning that
optimizes only for load distribution but not at all for communication.

This class is mainly provided for comparison, but might also be used as a
fallback if no graph partitioning software is available.
