# Welcome to SenID's Documentation

**SenID** — A Python package to identify cell-intrinsic and cell-extrinsic features
of cellular senescence and senescence-driven communication.

## Installation

The easiest way to install SenID is within a virtual environment.

```bash
# Create and activate a virtual environment
python3 -m venv senidenv
source senidenv/bin/activate
```

**Requires Python ≥ 3.10, < 3.11**

### Option 1: pip install

```bash
pip install senid
```

### Option 2: Install from GitHub

```bash
git clone https://github.com/axelalmet/SenID.git
cd SenID
pip install .
```

```{toctree}
:maxdepth: 2
:caption: Tutorials

tutorials/index
```

```{toctree}
:maxdepth: 3
:caption: API Reference

autoapi/senid/index
```

```{toctree}
:caption: Links

GitHub <https://github.com/axelalmet/SenID>
```