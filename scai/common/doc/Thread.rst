.. _Thread:

Thread
======

Whenever threads are used in LAMA, we use the std::thread class of the C++11 standard.
The same is true for related classes like ``mutex``, ``recursive_mutex``, 
``condition_variable_any``, and ``unique_lock``. Previous releases provided a 
LAMA specific thread class, either based on pThreads or Boost thread library.

