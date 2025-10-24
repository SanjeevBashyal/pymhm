# PymHM - Python Mesoscale Hydrological Model

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![PyPI version](https://badge.fury.io/py/pymhm.svg)](https://badge.fury.io/py/pymhm)

A Python package for mesoscale hydrological modeling, designed to work both as a standalone library and as a QGIS plugin.

## Features

- **Hydrological Modeling**: Core functionality for mesoscale hydrological processes
- **QGIS Integration**: Seamless integration with QGIS 3.0+ as a plugin
- **Data Processing**: Tools for handling hydrological data and model inputs
- **Visualization**: Built-in plotting and visualization capabilities
- **Extensible**: Modular design for easy extension and customization

## Installation

### From PyPI (Recommended)

```bash
pip install pymhm
```

### From Source

```bash
git clone https://github.com/SanjeevBashyal/pymhm.git
cd pymhm
pip install -e .
```

### For QGIS Plugin

1. Download the plugin from the QGIS Plugin Manager
2. Or install manually by copying the `pymhm` folder to your QGIS plugins directory

## Quick Start

### As a Python Package

```python
import pymhm

# Initialize the model
model = pymhm.HydrologicalModel()

# Run simulation
results = model.run_simulation()
```

### As a QGIS Plugin

1. Open QGIS
2. Go to Plugins → Manage and Install Plugins
3. Search for "PymHM"
4. Install and activate the plugin
5. Access the plugin from the Plugins menu

## Documentation

- [User Guide](https://github.com/SanjeevBashyal/pymhm/wiki)
- [API Reference](https://github.com/SanjeevBashyal/pymhm/wiki/API-Reference)
- [Examples](https://github.com/SanjeevBashyal/pymhm/tree/main/examples)

## Requirements

- Python 3.8+
- NumPy >= 1.20.0
- SciPy >= 1.7.0
- Pandas >= 1.3.0
- For QGIS plugin: QGIS 3.0+

## Development

### Setting up Development Environment

```bash
git clone https://github.com/SanjeevBashyal/pymhm.git
cd pymhm
pip install -e .[dev]
```

### Running Tests

```bash
pytest
```

### Building Documentation

```bash
pip install -e .[docs]
sphinx-build -b html docs/ docs/_build/
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the GNU General Public License v2.0 - see the [LICENSE](LICENSE) file for details.

## Citation

If you use PymHM in your research, please cite:

```bibtex
@software{pymhm2025,
  title={PymHM: Python Mesoscale Hydrological Model},
  author={Sanjeev Bashyal},
  year={2025},
  url={https://github.com/SanjeevBashyal/pymhm},
  version={0.1.0}
}
```

## Support

- **Issues**: [GitHub Issues](https://github.com/SanjeevBashyal/pymhm/issues)
- **Discussions**: [GitHub Discussions](https://github.com/SanjeevBashyal/pymhm/discussions)
- **Email**: sanjeev.bashyal01@gmail.com

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a list of changes and version history.

## Acknowledgments

- QGIS Development Team for the excellent platform
- Python Scientific Community for the amazing tools
- Contributors and users of PymHM