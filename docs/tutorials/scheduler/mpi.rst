MPI simulation
==============

synopsis
--------

In this tutorial we make ouer simulation mpi capable:

.. literalinclude:: code/mpi/isingmpi.hpp
   :language: c++
   :lines: 14-34

the ``alps::check_schedule`` class is used to determin when to communicate next. ``alps::check_schedule`` has the following interface:

* ``schedule(t_min, t_max)`` creats a ``check_schedule`` object with intervals in [t_min, t_max]
* ``schedule.pending()`` return ``true`` if the timer has gone off else ``false``
* ``schedule.check_interval()`` return the number of seconds in the current interval
* ``schedule.update(fraction)`` compute the new interval length and reset the timer

constructor
-----------

.. literalinclude:: code/mpi/isingmpi.cpp
   :language: c++
   :lines: 5-11

.. literalinclude:: code/mpi/isingmpi.cpp
   :language: c++
   :lines: 13-15

.. literalinclude:: code/mpi/isingmpi.cpp
   :language: c++
   :lines: 17-29

.. literalinclude:: code/mpi/isingmpi.cpp
   :language: c++
   :lines: 31-41


the main function
-----------------

.. literalinclude:: code/mpi/main.cpp
   :language: c++
   :lines: 5-
