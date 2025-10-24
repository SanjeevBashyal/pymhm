# Contributing to PymHM

Thank you for your interest in contributing to PymHM! This document provides guidelines and information for contributors.

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you agree to uphold this code.

## How to Contribute

### Reporting Issues

- Use the [GitHub issue tracker](https://github.com/SanjeevBashyal/pymhm/issues)
- Search existing issues before creating new ones
- Provide detailed information about the problem
- Include steps to reproduce the issue

### Suggesting Enhancements

- Use the [GitHub issue tracker](https://github.com/SanjeevBashyal/pymhm/issues) with the "enhancement" label
- Describe the feature and its benefits
- Consider the impact on existing functionality

### Contributing Code

1. **Fork the repository**
2. **Create a feature branch**
   ```bash
   git checkout -b feature/amazing-feature
   ```
3. **Make your changes**
   - Follow the coding style guidelines
   - Add tests for new functionality
   - Update documentation as needed
4. **Test your changes**
   ```bash
   pytest
   ```
5. **Commit your changes**
   ```bash
   git commit -m 'Add some amazing feature'
   ```
6. **Push to the branch**
   ```bash
   git push origin feature/amazing-feature
   ```
7. **Open a Pull Request**

## Development Setup

### Prerequisites

- Python 3.8+
- Git
- QGIS 3.0+ (for plugin development)

### Installation

```bash
# Clone the repository
git clone https://github.com/SanjeevBashyal/pymhm.git
cd pymhm

# Install in development mode
pip install -e .[dev]

# Install pre-commit hooks (optional)
pre-commit install
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=pymhm

# Run specific test file
pytest tests/test_specific.py
```

### Code Style

We use the following tools for code quality:

- **Black**: Code formatting
- **Flake8**: Linting
- **MyPy**: Type checking

```bash
# Format code
black pymhm/

# Check linting
flake8 pymhm/

# Type checking
mypy pymhm/
```

## Pull Request Guidelines

### Before Submitting

- [ ] Code follows the project's style guidelines
- [ ] Tests pass locally
- [ ] Documentation is updated
- [ ] Changes are focused and atomic
- [ ] Commit messages are clear and descriptive

### PR Description

- Describe what the PR does
- Reference any related issues
- Include screenshots for UI changes
- List any breaking changes

## Documentation

### Docstrings

Use Google-style docstrings:

```python
def example_function(param1: str, param2: int) -> bool:
    """Example function with proper docstring.
    
    Args:
        param1: Description of param1
        param2: Description of param2
        
    Returns:
        Description of return value
        
    Raises:
        ValueError: When param2 is negative
    """
    pass
```

### Documentation Updates

- Update README.md for user-facing changes
- Update docstrings for API changes
- Add examples for new features
- Update CHANGELOG.md for releases

## Release Process

1. Update version numbers
2. Update CHANGELOG.md
3. Create release branch
4. Tag the release
5. Build and upload to PyPI

## Questions?

- Open a [GitHub Discussion](https://github.com/SanjeevBashyal/pymhm/discussions)
- Contact: sanjeev.bashyal01@gmail.com

Thank you for contributing to PymHM!
