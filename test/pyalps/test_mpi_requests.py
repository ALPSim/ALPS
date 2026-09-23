"""Exercise real object transport; run both normally and under two-rank MPI."""

import time

import numpy as np
import pytest

pytest.importorskip("mpi4py")
from pyalps import mpi


def test_communicator_equality_and_hashability():
    from mpi4py import MPI

    wrappers = (mpi.world, mpi.Communicator(), mpi.Communicator(mpi.world))
    for comm in wrappers:
        assert comm == mpi.world
        assert not (comm != mpi.world)
        assert comm != mpi.Communicator(MPI.COMM_NULL)
        assert comm != object()
        with pytest.raises(TypeError):
            hash(comm)


def poll(function):
    deadline = time.monotonic() + 15
    while time.monotonic() < deadline:
        result = function()
        if result is not None and result is not False:
            return result
        time.sleep(0.0001)
    pytest.fail("MPI request did not complete")


@pytest.mark.parametrize("completion", ["wait", "test", "wait_all", "test_all", "wait_any", "test_any", "wait_some", "test_some"])
def test_large_object_requests(completion):
    destination = (mpi.rank + 1) % mpi.size
    source = (mpi.rank - 1) % mpi.size
    payload = {"source": mpi.rank, "data": b"x" * (2 * 1024 * 1024)}
    send = mpi.world.isend(destination, 810, payload)
    receive = mpi.world.irecv(source, 810)
    callbacks = []
    if completion == "wait":
        value, status = receive.wait()
    elif completion == "test":
        value, status = poll(receive.test)
    elif completion.endswith("any"):
        value, status, index = poll(lambda: getattr(mpi, completion)([receive]))
        assert index == 0
    else:
        # Put the rendezvous send first, so a sequential blocking wait would
        # deadlock unless the receive is progressed before the send completes.
        requests = mpi.RequestList([send, receive])
        callback = lambda value, status: callbacks.append((value, status))
        if completion == "wait_all":
            mpi.wait_all(requests, callback)
        elif completion == "test_all":
            poll(lambda: mpi.test_all(requests, callback))
        else:
            deadline = time.monotonic() + 15
            while requests and time.monotonic() < deadline:
                boundary = getattr(mpi, completion)(requests, callback)
                del requests[boundary:]
            assert not requests
        value, status = next(pair for pair in callbacks if pair[0] is not None)
        assert len(callbacks) == 2
    send.wait()
    assert value == {"source": source, "data": payload["data"]}
    assert status.source == source and status.tag == 810
    mpi.world.barrier()


def test_posted_receives_keep_order_when_waited_backwards():
    # Two wrappers of the same communicator must share posting order.
    first = mpi.world.irecv(mpi.any_source, mpi.any_tag)
    second = mpi.Communicator(mpi.world).irecv(mpi.any_source, mpi.any_tag)
    sends = [mpi.world.isend(mpi.rank, 811, i) for i in (1, 2)]
    assert second.wait()[0] == 2
    assert first.wait()[0] == 1
    mpi.wait_all(sends)
    mpi.world.barrier()


def test_partial_test_all_keeps_received_values():
    requests = [mpi.world.irecv(mpi.rank, tag) for tag in (812, 813)]
    first_send = mpi.world.isend(mpi.rank, 812, "first")
    assert poll(requests[0].test)[0] == "first"
    callbacks = []
    assert not mpi.test_all(requests, lambda *pair: callbacks.append(pair))
    assert callbacks == []
    second_send = mpi.world.isend(mpi.rank, 813, "second")
    mpi.wait_all(requests, lambda *pair: callbacks.append(pair))
    mpi.wait_all([first_send, second_send])
    assert [value for value, status in callbacks] == ["first", "second"]
    assert [status.tag for value, status in callbacks] == [812, 813]
    mpi.world.barrier()


def test_cancelling_unmatched_receive_does_not_consume_later_message():
    receive = mpi.world.irecv(mpi.rank, 814)
    assert receive.test() is None
    receive.cancel()
    value, status = receive.wait()
    assert value is None and status.Is_cancelled()
    send = mpi.world.isend(mpi.rank, 814, "after cancellation")
    assert mpi.world.recv(mpi.rank, 814) == "after cancellation"
    send.wait()
    mpi.world.barrier()


def test_blocking_large_send_progresses_posted_receive():
    receive = mpi.world.irecv(mpi.rank, 815)
    mpi.world.send(mpi.rank, 815, b"y" * (2 * 1024 * 1024))
    assert receive.wait()[0] == b"y" * (2 * 1024 * 1024)
    mpi.world.barrier()


def test_large_numpy_payload_with_wildcard_receive():
    values = (np.arange(256 * 1024).reshape(512, 512) + mpi.rank * 1j).astype(np.complex128)
    source = (mpi.rank - 1) % mpi.size
    send = mpi.world.isend((mpi.rank + 1) % mpi.size, 816, {"array": values, "source": mpi.rank})
    value, status = mpi.world.irecv(mpi.any_source, 816).wait()
    send.wait()
    assert status.source == source and value["source"] == source
    np.testing.assert_array_equal(value["array"], np.arange(256 * 1024).reshape(512, 512) + source * 1j)
    mpi.world.barrier()
