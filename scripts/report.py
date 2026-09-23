from datetime import datetime
from pathlib import Path

HEADER_PREFIX = "# Report started: "
TS_FORMAT = "%Y-%m-%d %H:%M:%S.%f"


def _format_elapsed(seconds: float) -> str:
    """Format elapsed seconds as e.g. 01h:05min:09.37s"""
    seconds = round(seconds, 2)  # round first so 59.999 -> 01min:00.00s, not 60.00s
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    return f"{int(hours):02d}h:{int(minutes):02d}min:{secs:05.2f}s"


def _read_start_time(report_file: Path):
    """Return the start time stored in the file header, or None if not found."""
    with open(report_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(HEADER_PREFIX):
                try:
                    return datetime.strptime(line[len(HEADER_PREFIX):].strip(), TS_FORMAT)
                except ValueError:
                    return None
    return None


def write_report(path, message, overwrite=False, filename="run_report.txt"):
    """
    Write a timestamped message to a report file.

    Args:
        path: Folder to put the report in, or a full file path (e.g. 'out/run1.txt').
        message: Text to log. Multi-line messages are indented under the timestamp.
        overwrite: True  -> start a fresh report (old content is deleted).
                   False -> append to the existing report (created if missing).
        filename: Report name used when `path` is a folder.

    Returns:
        Path to the report file.
    """
    path = Path(path)
    report_file = path / filename if (path.is_dir() or path.suffix == "") else path
    report_file.parent.mkdir(parents=True, exist_ok=True)

    now = datetime.now()
    start = None
    if report_file.exists() and not overwrite:
        start = _read_start_time(report_file)

    header = ""
    if overwrite or not report_file.exists():
        mode, start = "w", now
        header = f"{HEADER_PREFIX}{now.strftime(TS_FORMAT)}\n"
    else:
        mode = "a"
        if start is None:  # existing file without our header -> mark a start point
            start = now
            header = f"{HEADER_PREFIX}{now.strftime(TS_FORMAT)}\n"

    elapsed = _format_elapsed((now - start).total_seconds())
    stamp = now.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
    msg = str(message).replace("\n", "\n" + " " * 4)

    with open(report_file, mode, encoding="utf-8") as f:
        f.write(header)
        f.write(f"[{stamp}] [+{elapsed}] {msg}\n")

    return report_file


if __name__ == "__main__":
    import time

    write_report("results", "Pipeline started", overwrite=True)
    time.sleep(1.2)
    write_report("results", "Loaded 3 samples from data/raw")
    time.sleep(0.5)
    write_report("results", "Analysis done\nSaved plots to results/plots")
    print(open("results/report.txt").read())
