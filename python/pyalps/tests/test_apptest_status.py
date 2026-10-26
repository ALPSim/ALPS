"""pyalps.apptest must not report a crashed or incomplete run as passing."""

import os

import pytest

import pyalps
import pyalps.apptest as apptest

PROPS = {
    "SCRIPT": "script.py",
    "TESTNAME": "status",
    "TOLERANCE": "auto",
    "SAVE_OUT_IF_FAIL": "no",
    "WRITE_RESULTS_TO_FILE": "no",
}


@pytest.fixture
def testprop(monkeypatch, tmp_path):
    def configure(refs):
        monkeypatch.setattr(apptest, "read_testprop_xml",
                            lambda testinputfile: (dict(PROPS), None, list(refs), None))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "script.py").write_text("")
    monkeypatch.setattr(pyalps.tools, "copyStylesheet", lambda dirname: None)
    return configure


def write_outputs_and_exit(status):
    def execute(cmdline):
        open("result.h5", "w").close()
        return status
    return execute


def test_run_test_fails_when_script_exits_non_zero(monkeypatch, testprop):
    testprop(["result.h5"])
    monkeypatch.setattr(pyalps, "executeCommand", write_outputs_and_exit(1))
    monkeypatch.setattr(apptest, "compareTest", lambda *args, **kwargs: True)
    assert apptest.runTest("status.testin.xml") is False


def test_run_test_compares_when_script_succeeds(monkeypatch, testprop):
    testprop(["result.h5"])
    monkeypatch.setattr(pyalps, "executeCommand", write_outputs_and_exit(0))
    monkeypatch.setattr(apptest, "compareTest", lambda *args, **kwargs: True)
    assert apptest.runTest("status.testin.xml") is True


def test_compare_test_fails_on_missing_reference_pair(monkeypatch, testprop, tmp_path):
    testprop(["ref1.h5", "ref2.h5"])
    monkeypatch.setattr(apptest, "checkProperties", lambda tf, rf: True)
    monkeypatch.setattr(apptest, "detectDataType",
                        lambda tfile: lambda *args, **kwargs: (True, []))
    outputs = [str(tmp_path / "out1.h5")]
    assert apptest.compareTest("status.testin.xml", outputs, str(tmp_path), "start") is False


def test_create_test_raises_when_reference_script_fails(monkeypatch, testprop, tmp_path):
    monkeypatch.setattr(pyalps, "executeCommand", lambda cmdline: 1)
    with pytest.raises(Exception, match="exited with status 1"):
        apptest.createTest("script.py", outputs=["result.h5"], refdir=str(tmp_path / "ref"))
    assert os.getcwd() == str(tmp_path)
