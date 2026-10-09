"""Execute notebook assertions and render a static, downloadable tutorial site."""

import argparse
import base64
import hashlib
import importlib.metadata
import json
import os
import shutil
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def build(output, write_executed=False):
    output.mkdir(parents=True, exist_ok=True)
    # Keep Jupyter/matplotlib runtime files inside the build; kernel uses this interpreter.
    runtime = output / ".runtime"
    kernel_dir = runtime / "jupyter" / "kernels" / "gesc"
    kernel_dir.mkdir(parents=True, exist_ok=True)
    (kernel_dir / "kernel.json").write_text(
        json.dumps(
            {
                "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
                "display_name": "GESC tutorial Python",
                "language": "python",
            }
        )
    )
    for name, path in {
        "JUPYTER_DATA_DIR": runtime / "jupyter",
        "JUPYTER_RUNTIME_DIR": runtime / "jupyter-runtime",
        "IPYTHONDIR": runtime / "ipython",
        "MPLCONFIGDIR": runtime / "matplotlib",
    }.items():
        path.mkdir(parents=True, exist_ok=True)
        os.environ[name] = str(path)
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"
    import nbformat
    from nbclient import NotebookClient
    from nbconvert import HTMLExporter

    assets = output / "assets"
    assets.mkdir(exist_ok=True)
    shutil.copy(ROOT / "site/assets/style.css", assets / "style.css")
    downloads = output / "notebooks"
    downloads.mkdir(exist_ok=True)
    shutil.copytree(
        ROOT / "tutorials",
        output / "tutorials",
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.ipynb"),
    )
    reports = []
    for path in sorted((ROOT / "tutorials").glob("*.ipynb")):
        start = time.monotonic()
        notebook = nbformat.read(path, as_version=4)
        # Never trust stored outputs: execute every cell again in a fresh kernel.
        for cell in notebook.cells:
            if cell.cell_type == "code":
                cell.outputs = []
                cell.execution_count = None
        print(f"Executing {path.name}", flush=True)
        client = NotebookClient(
            notebook,
            timeout=1800,
            kernel_name="gesc",
            allow_errors=False,
            resources={"metadata": {"path": str(ROOT)}},
        )
        try:
            client.execute()
        except Exception:
            failed = output.parent / "failed-notebooks"
            failed.mkdir(parents=True, exist_ok=True)
            nbformat.write(notebook, failed / path.name)
            raise
        metrics, figures = [], []
        for cell in notebook.cells:
            for item in cell.get("outputs", []):
                if item.output_type == "error":
                    raise RuntimeError(f"Notebook error in {path.name}")
                data = item.get("data", {})
                if "application/json" in data:
                    metrics.append(data["application/json"])
                if "image/png" in data:
                    figures.append(data["image/png"])
        if not metrics or not all(m.get("checks_passed") is True for m in metrics):
            raise RuntimeError(f"Missing passing scientific checks in {path.name}")
        if not figures:
            raise RuntimeError(f"No scientific figures in {path.name}")
        for index, encoded in enumerate(figures):
            (assets / f"{path.stem}-{index}.png").write_bytes(base64.b64decode(encoded))
        if write_executed:
            nbformat.write(notebook, path)
        nbformat.write(notebook, downloads / path.name)
        exporter = HTMLExporter(
            template_name="gesc",
            extra_template_basedirs=[str(ROOT / "site/templates")],
            exclude_input_prompt=True,
            exclude_output_prompt=True,
        )
        html, _ = exporter.from_notebook_node(
            notebook,
            resources={
                "metadata": {"name": path.stem},
                "gesc_title": notebook.metadata.gesc.title,
                "gesc_slug": path.stem,
            },
        )
        (output / f"{path.stem}.html").write_text(html)
        reports.append(
            {
                "notebook": path.name,
                "title": notebook.metadata.gesc.title,
                "elapsed_seconds": round(time.monotonic() - start, 2),
                "figure_count": len(figures),
                "metrics": metrics,
                "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
        print(
            f"Passed {path.name}: {len(figures)} figures, {reports[-1]['elapsed_seconds']} s",
            flush=True,
        )
    manifest = {
        "built_utc": datetime.now(UTC).isoformat(),
        "commit": os.environ.get("GITHUB_SHA", "local-build"),
        "python": sys.version.split()[0],
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "pydantic", "matplotlib", "nbclient", "nbconvert")
        },
        "tutorials": reports,
    }
    (output / "validation.json").write_text(json.dumps(manifest, indent=2) + "\n")
    index = (ROOT / "site/index.html").read_text()
    index = index.replace("{{BUILD_DATE}}", manifest["built_utc"][:10])
    index = index.replace("{{FIGURE_COUNT}}", str(sum(item["figure_count"] for item in reports)))
    (output / "index.html").write_text(index)
    (output / ".nojekyll").touch()
    shutil.copy(ROOT / "tutorials/README.md", output / "tutorials/README.md")
    shutil.rmtree(runtime)
    validate_site(output, reports)
    print(f"Validated {len(reports)} executed notebooks in {output}", flush=True)


def validate_site(output, reports):
    from bs4 import BeautifulSoup

    for file in output.glob("*.html"):
        soup = BeautifulSoup(file.read_text(), "html.parser")
        if not soup.title or not soup.select("h1"):
            raise ValueError(f"Missing page title/h1: {file.name}")
        for item in soup.select("a[href], img[src], link[href]"):
            reference = item.get("href", item.get("src", ""))
            if reference.startswith(("http:", "https:", "data:", "#", "mailto:", "//")):
                continue
            if not (file.parent / reference.split("#")[0]).exists():
                raise ValueError(f"Broken local reference {reference} in {file.name}")
    if len(reports) != 5:
        raise ValueError("All five requested tutorials must be executed and published.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "build/site")
    parser.add_argument(
        "--write-executed",
        action="store_true",
        help="Save fresh outputs back into tracked source notebooks.",
    )
    args = parser.parse_args()
    build(args.output.resolve(), args.write_executed)
