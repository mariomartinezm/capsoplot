# pycapso

This package contains some utilities used to analyse the behavior of the CAPSO model.

## Installation

```bash
$ pip install pycapso
```

## Usage

`pycapso` contains functions used to analyze the behavior of the CAPSO model,
e.g., mean field approximation. It also contains wrappers around matplotlib to
make the plotting of data easier.

```python
import matplotlib.pyplot as plt
from pycapso.models import get_mean_field

prey_pred = get_mean_field(num_iter=1000)

plt.plot(prey_pred)

plt.show()
```

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`pycapso` was created by Mario Martínez. It is licensed under the terms of the GNU General Public License v3.0 license.

## Credits

`pycapso` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
