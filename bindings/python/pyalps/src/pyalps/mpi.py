"""MPI compatibility layer backed by :mod:`mpi4py`.

The historic module re-exported Boost.MPI's Boost.Python bindings. Rebuilding
that second Python binding stack would couple pyalps to Boost.Python again and
make ordinary wheels depend on one particular MPI implementation. Instead,
this module preserves the commonly used Boost.MPI Python spelling on top of
mpi4py. Install ``pyalps[mpi]`` to enable it.

The compatibility surface covers the world communicator, point-to-point
operations, collectives, status/request types, environment queries, and the
timer API. It intentionally does not reproduce Boost.MPI's C++/Python
serialization bridge or skeleton/content optimization; use mpi4py buffers for
that level of interoperability.
"""

from __future__ import annotations

import atexit as _atexit
from functools import reduce as _python_reduce
import sys
from typing import Any

try:
    import mpi4py as _mpi4py

    _mpi_module_was_loaded = "mpi4py.MPI" in sys.modules
    _previous_auto_initialize = _mpi4py.rc.initialize
    if not _mpi_module_was_loaded:
        # Delay mpi4py's automatic MPI_Init just long enough to distinguish an
        # externally initialized MPI process from an environment this module
        # must own. This reproduces Boost.MPI's finalize-only-what-we-created
        # behavior and leaves the global mpi4py setting as we found it.
        _mpi4py.rc.initialize = False
    try:
        from mpi4py import MPI as _MPI
    finally:
        if not _mpi_module_was_loaded:
            _mpi4py.rc.initialize = _previous_auto_initialize
except ImportError as error:  # pragma: no cover - depends on optional install
    raise ImportError(
        "pyalps.mpi requires mpi4py; install the optional dependency with "
        "'python -m pip install pyalps[mpi]'"
    ) from error


_initialized_here = False
if not _MPI.Is_initialized():
    _MPI.Init()
    _initialized_here = True


any_source = _MPI.ANY_SOURCE
any_tag = _MPI.ANY_TAG
Exception = _MPI.Exception
Status = _MPI.Status


class Request:
    """Non-value request with the Boost.MPI ``wait``/``test`` contract."""

    def __init__(self, request: Any):
        self._request = request

    def wait(self):
        status = Status()
        self._request.wait(status)
        return status

    def test(self):
        status = Status()
        flag, _value = self._request.test(status)
        return status if flag else None

    def cancel(self) -> None:
        self._request.cancel()


class RequestWithValue(Request):
    """Receive request whose completion returns ``(value, status)``."""

    def wait(self):
        status = Status()
        value = self._request.wait(status)
        return value, status

    def test(self):
        status = Status()
        flag, value = self._request.test(status)
        return (value, status) if flag else None


class RequestList(list):
    """Mutable request sequence used by the nonblocking helper functions."""


class Communicator:
    """Boost.MPI-compatible wrapper around an ``mpi4py.MPI.Comm``."""

    def __init__(self, comm: Any = None):
        if isinstance(comm, Communicator):
            comm = comm._comm
        self._comm = _MPI.COMM_WORLD if comm is None else comm

    @property
    def rank(self) -> int:
        return self._comm.rank

    @property
    def size(self) -> int:
        return self._comm.size

    def __bool__(self) -> bool:
        return self._comm != _MPI.COMM_NULL

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Communicator) and self._comm == other._comm

    def send(self, dest: int, tag: int = 0, value: Any = None) -> None:
        self._comm.send(value, dest=dest, tag=tag)

    def recv(
        self,
        source: int = any_source,
        tag: int = any_tag,
        return_status: bool = False,
    ) -> Any:
        status = _MPI.Status() if return_status else None
        value = self._comm.recv(source=source, tag=tag, status=status)
        return (value, status) if return_status else value

    def isend(self, dest: int, tag: int = 0, value: Any = None):
        return Request(self._comm.isend(value, dest=dest, tag=tag))

    def irecv(self, source: int = any_source, tag: int = any_tag):
        return RequestWithValue(self._comm.irecv(source=source, tag=tag))

    def probe(self, source: int = any_source, tag: int = any_tag):
        status = _MPI.Status()
        self._comm.probe(source=source, tag=tag, status=status)
        return status

    def iprobe(self, source: int = any_source, tag: int = any_tag):
        status = _MPI.Status()
        return status if self._comm.iprobe(source=source, tag=tag, status=status) else None

    def barrier(self) -> None:
        self._comm.barrier()

    def split(self, color: int, key: int = 0) -> "Communicator":
        return Communicator(self._comm.Split(color=color, key=key))

    def abort(self, errcode: int) -> None:
        self._comm.Abort(errcode)


world = Communicator(_MPI.COMM_WORLD)
rank = world.rank
size = world.size


def _unwrap(comm: Any):
    return comm._comm if isinstance(comm, Communicator) else comm


def _collective_values(comm: Any, value: Any) -> tuple[Any, ...]:
    return tuple(_unwrap(comm).allgather(value))


def all_gather(comm: Any = world, value: Any = None) -> tuple[Any, ...]:
    return _collective_values(comm, value)


def all_to_all(comm: Any = world, values: Any = None) -> tuple[Any, ...]:
    return tuple(_unwrap(comm).alltoall(values))


def broadcast(comm: Any = world, value: Any = None, root: int = 0) -> Any:
    return _unwrap(comm).bcast(value, root=root)


def gather(comm: Any = world, value: Any = None, root: int = 0):
    values = _unwrap(comm).gather(value, root=root)
    return tuple(values) if _unwrap(comm).rank == root else None


def scatter(comm: Any = world, values: Any = None, root: int = 0) -> Any:
    return _unwrap(comm).scatter(values, root=root)


def _apply_operation(values: tuple[Any, ...], op: Any) -> Any:
    if op is None:
        raise TypeError("an operation callable is required")
    return _python_reduce(op, values)


def reduce(comm: Any = world, value: Any = None, op: Any = None, root: int = 0):
    # Boost.MPI accepted arbitrary Python callables. Gathering before the
    # Python reduction preserves that behavior; users wanting native MPI
    # reductions can call the underlying ``world._comm`` directly.
    values = gather(comm, value, root)
    return _apply_operation(values, op) if _unwrap(comm).rank == root else None


def all_reduce(comm: Any = world, value: Any = None, op: Any = None) -> Any:
    return _apply_operation(_collective_values(comm, value), op)


def scan(comm: Any = world, value: Any = None, op: Any = None) -> Any:
    values = _collective_values(comm, value)
    return _apply_operation(values[: _unwrap(comm).rank + 1], op)


def _check_requests(requests) -> None:
    if not requests:
        raise ValueError("cannot wait on an empty request vector")
    if not all(isinstance(request, Request) for request in requests):
        raise TypeError("requests must contain pyalps.mpi Request objects")


def _raw_requests(requests):
    _check_requests(requests)
    return [request._request for request in requests]


def wait_any(requests):
    status = Status()
    index, value = _MPI.Request.waitany(_raw_requests(requests), status)
    return value, status, index


def test_any(requests):
    status = Status()
    index, flag, value = _MPI.Request.testany(_raw_requests(requests), status)
    return (value, status, index) if flag else None


def wait_all(requests, callable=None) -> None:
    statuses = [Status() for _ in requests]
    values = _MPI.Request.waitall(_raw_requests(requests), statuses)
    if callable is not None:
        for value, status in zip(values, statuses):
            callable(value, status)


def test_all(requests, callable=None) -> bool:
    statuses = [Status() for _ in requests]
    flag, values = _MPI.Request.testall(_raw_requests(requests), statuses)
    if flag and callable is not None and values is not None:
        for value, status in zip(values, statuses):
            callable(value, status)
    return bool(flag)


def wait_some(requests, callable=None) -> int:
    statuses = [Status() for _ in requests]
    indices, values = _MPI.Request.waitsome(_raw_requests(requests), statuses)
    return _finish_some(requests, indices, values, statuses, callable)


def test_some(requests, callable=None) -> int:
    statuses = [Status() for _ in requests]
    indices, values = _MPI.Request.testsome(_raw_requests(requests), statuses)
    return _finish_some(requests, indices, values, statuses, callable)


def _finish_some(requests, indices, values, statuses, callable) -> int:
    if not indices:
        return len(requests)
    if callable is not None:
        for value, status in zip(values, statuses):
            callable(value, status)

    # Boost.MPI partitions the mutable RequestList into pending requests
    # followed by completed requests and returns the first completed index.
    completed = set(indices)
    pending_requests = [r for i, r in enumerate(requests) if i not in completed]
    completed_requests = [requests[i] for i in indices]
    requests[:] = pending_requests + completed_requests
    return len(pending_requests)


class Timer:
    def __init__(self):
        self.restart()

    def restart(self) -> float:
        previous = getattr(self, "_start", _MPI.Wtime())
        self._start = _MPI.Wtime()
        return self._start - previous

    @property
    def elapsed(self) -> float:
        return _MPI.Wtime() - self._start

    @property
    def elapsed_min(self) -> float:
        return _MPI.Wtick()

    @property
    def elapsed_max(self) -> float:
        return sys.float_info.max

    @property
    def time_is_global(self) -> bool:
        return bool(_MPI.COMM_WORLD.Get_attr(_MPI.WTIME_IS_GLOBAL))


def init(argv=None, abort_on_exception: bool = True) -> bool:
    del argv, abort_on_exception
    global _initialized_here
    if _MPI.Is_initialized():
        return False
    _MPI.Init()
    _initialized_here = True
    return True


def finalize() -> None:
    global _initialized_here
    if _initialized_here and _MPI.Is_initialized() and not _MPI.Is_finalized():
        _MPI.Finalize()
    _initialized_here = False


if _initialized_here:
    _atexit.register(finalize)


def abort(errcode: int) -> None:
    _MPI.COMM_WORLD.Abort(errcode)


def initialized() -> bool:
    return _MPI.Is_initialized()


def finalized() -> bool:
    return _MPI.Is_finalized()


collectives_tag = _MPI.COMM_WORLD.Get_attr(_MPI.TAG_UB)
max_tag = collectives_tag - 1
processor_name = _MPI.Get_processor_name()
_host_key = getattr(_MPI, "HOST", None)
_io_key = getattr(_MPI, "IO", None)
host_rank = _MPI.COMM_WORLD.Get_attr(_host_key) if _host_key is not None else None
io_rank = _MPI.COMM_WORLD.Get_attr(_io_key) if _io_key is not None else None


__all__ = [
    "Communicator", "Exception", "Request", "RequestList", "RequestWithValue",
    "Status", "Timer", "abort", "all_gather", "all_reduce", "all_to_all",
    "any_source", "any_tag", "broadcast", "collectives_tag", "finalize",
    "finalized", "gather", "host_rank", "init", "initialized", "io_rank",
    "max_tag", "processor_name", "rank", "reduce", "scan", "scatter", "size",
    "test_all", "test_any", "test_some", "wait_all", "wait_any", "wait_some",
    "world",
]
