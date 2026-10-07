"""Build the public numerical paper with pdflatex and BibTeX."""
import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(command, cwd, log):
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True,
                            encoding="utf-8", errors="replace")
    log.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}); see {log}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=ROOT / "output/pdf/numerical-range-certificates.pdf")
    args = parser.parse_args()
    for tool in ("pdflatex", "bibtex"):
        if not shutil.which(tool):
            raise SystemExit(f"Required executable missing: {tool}")
    source = ROOT / "paper/main.tex"
    work = ROOT / "output/latex"
    work.mkdir(parents=True, exist_ok=True)
    command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
               "-file-line-error", f"-output-directory={work.as_posix()}", source.name]
    run(command, source.parent, work / "pass-1.txt")
    run(["bibtex", (work / source.stem).as_posix()], source.parent, work / "bibtex.txt")
    for number in (2, 3, 4):
        run(command, source.parent, work / f"pass-{number}.txt")
    log = (work / "main.log").read_text(encoding="utf-8", errors="replace")
    warnings = re.findall(r"(?:Overfull[^\n]*|[^\n]*Warning[^\n]*|[^\n]*undefined[^\n]*)", log)
    if any("Overfull" in row or "undefined" in row for row in warnings):
        raise RuntimeError("Unresolved paper layout/reference warnings: " + "\n".join(warnings))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(work / "main.pdf", args.output)
    print(json.dumps({"source": "paper/main.tex", "output": str(args.output),
                      "passes": 4, "warnings": warnings}))


if __name__ == "__main__":
    main()
