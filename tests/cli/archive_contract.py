"""Check restored archive commands in independent processes and inspect their output."""
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def run(program, *args, success=True):
    result = subprocess.run([program, *map(str, args)], capture_output=True, text=True,
                            timeout=15)
    assert (result.returncode == 0) == success, result.stdout + result.stderr
    return result.stdout


def main():
    text_command, xml_command, *sqlite_command = sys.argv[1:]
    with tempfile.TemporaryDirectory(prefix="alps archive ") as directory:
        root = Path(directory)
        data = root / "data.txt"
        data.write_text("1 3 0.1\n\n2 4 0.2\n")
        result = ET.fromstring(run(text_command, "-x", "T", "-y", "Energy", "-e", data))
        assert result.tag == "ARCHIVE"
        simulations = result.findall("SIMULATION")
        assert len(simulations) == 2
        for item, x, y, error in zip(simulations, (1, 2), (3, 4), (0.1, 0.2)):
            parameter = item.find("PARAMETERS/PARAMETER")
            measurement = item.find("AVERAGES/SCALAR_AVERAGE")
            assert parameter.attrib["name"] == "T" and float(parameter.text) == x
            assert measurement.attrib["name"] == "Energy"
            assert float(measurement.findtext("MEAN")) == y
            assert float(measurement.findtext("ERROR")) == error
        run(text_command, "-x", "T", "-y", "Energy", root / "absent", success=False)
        data.write_text("bad row\n")
        run(text_command, "-x", "T", "-y", "Energy", data, success=False)

        simulation = root / "simulation.xml"
        simulation.write_text('''<SIMULATION><PARAMETERS>
<PARAMETER name="T">1.25</PARAMETER></PARAMETERS><AVERAGES>
<SCALAR_AVERAGE name="Energy"><COUNT>64</COUNT><MEAN>-2</MEAN>
<ERROR converged="yes">0.125</ERROR><VARIANCE>1</VARIANCE><AUTOCORR>0</AUTOCORR>
</SCALAR_AVERAGE></AVERAGES></SIMULATION>''')
        job = root / "job.in.xml"
        job.write_text('''<JOB name="archive regression"><OUTPUT file="job.out.xml"/>
<TASK status="finished"><INPUT file="simulation.xml"/>
<OUTPUT file="simulation.xml"/></TASK></JOB>''')
        result = ET.fromstring(run(xml_command, job))
        assert result.tag == "ARCHIVE" and result.attrib["name"] == "archive regression"
        assert float(result.findtext("SIMULATION/PARAMETERS/PARAMETER")) == 1.25
        assert float(result.findtext("SIMULATION/AVERAGES/SCALAR_AVERAGE/MEAN")) == -2
        run(xml_command, root / "absent.xml", success=False)

        if sqlite_command:
            command = sqlite_command[0]
            database = root / "results.db"
            run(command, "--command", "install", "--db-file", database)
            run(command, "--command", "append", "--db-file", database, "--xml-file", simulation)
            with sqlite3.connect(database) as db:
                assert db.execute("SELECT value FROM parameter WHERE name='T'").fetchone() == ('1.25',)
                assert db.execute("SELECT count, mean, error FROM measurement WHERE name='Energy'").fetchone() == ('64', '-2', '0.125')
            # A second process must reopen the same schema and retain its data.
            listing = run(command, "--command", "list", "--db-file", database)
            assert "Energy" in listing and "1.25" in listing
            run(command, "--command", "append", "--db-file", database, "--xml-file", simulation)
            with sqlite3.connect(database) as db:
                assert db.execute("SELECT count(*) FROM measurement").fetchone() == (1,)


if __name__ == "__main__":
    main()
