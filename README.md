## Prerequisites

- Ubuntu/Debian: `sudo apt install python3-pip build-essential curl`
- Fedora/RedHat: `sudo dnf install python3-pip gcc make`

## Running

Assuming `results.json` will be the file where experimental data is written to.

### Evaluation experiments

```sh
./main.out eval results.json
```

### Local Search experiments

```sh
./main.out ls results.json
```

### Path Relinking experiments

```sh
./main.out pr results.json
```

## Generating LaTeX figures and tables

Python 3 and the PyLatex package are required.

```sh
pip3 install --user pylatex
```

To generate figures:

```sh
./genfigures.py results.json
```

To generate tables:

```sh
./genfigures.py results.json
```
