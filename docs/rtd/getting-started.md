# Getting Started

This page gets you from zero to a working Magrathea build and a first run.

## 1) Prerequisites

- A C++17 compiler (gcc or clang)
- CMake or Make (depending on project)
- GSL (GNU Scientific Library)
- Git
- Python 3 (only for building docs)

> **Tip:** If you have a working build already, skip to **First Run**.

## 2) Build from source (example)

```bash
git clone https://github.com/<YOUR-ORG>/Magrathea.git
cd Magrathea
# Option A: Makefile
make
# Option B: CMake (if your project uses CMake)
# cmake -S . -B build && cmake --build build --config Release
```

- The build should produce a `planet` (or similar) executable.
- If build fails, check that GSL is installed and discoverable.

## 3) First Run

Try one of the example configuration files (replace with your actual path/names):

```bash
./planet mode1.cfg
```

Expected outputs (example names; update to match the code):
- `out_profile.txt` — radial profiles (P, T, ρ, etc.)
- `out_summary.txt` — bulk properties (M, R, layer masses)
- Logs to stdout with convergence info

If this works, you're good to move on.

## 4) Build the docs locally (optional)

```bash
# from repo root
pip install -r docs/requirements.txt
doxygen Doxyfile
sphinx-build -b html docs docs/_build/html
# open docs/_build/html/index.html
```

If you see a page titled **Magrathea Documentation**, the pipeline is working.

## 5) Where to put your content

- Tutorials and guides → Markdown files in `docs/`
- C++ comments for API → in your headers, using Doxygen-style comments
- Configuration reference → add a new page in `docs/user-guide/` or `docs/reference/`
