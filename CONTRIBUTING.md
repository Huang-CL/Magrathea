# Contributing to MAGRATHEA

Contributions to MAGRATHEA are welcome. These may include bug fixes, documentation, equations of state, phase diagrams, and improvements to the numerical methods.

## Questions and bug reports

Use [GitHub Issues](https://github.com/Huang-CL/Magrathea/issues) to ask questions or report problems. Please search existing issues before opening a new one.

When reporting a bug, include:

* Your operating system, compiler version, and GSL version.
* The MAGRATHEA version or commit used.
* The command, configuration file, and input files needed to reproduce the problem.
* The full error message or relevant terminal output.
* A brief description of the expected and observed behavior.

For private matters or possible research collaborations, contact one of the maintainers directly.

## Suggesting a feature

Feature requests are also welcome. Open an issue describing the scientific or technical problem, the proposed change, and any relevant publications or similar implementations.

For substantial changes to the solver, input or output formats, or user interface, please open an issue before beginning development. This allows the approach to be discussed before extensive work is done. (We'd love someone to couple an atmosphere model!)

## Scientific contributions

New equations of state, phase diagrams, and tabulated data should include:

* The original scientific source.
* Units and parameter definitions.
* The pressure and temperature range over which the model or data are valid.
* Any assumptions or extrapolations.
* An example showing that the implementation works.

## Pull requests

To contribute code:

1. Fork the repository and create a branch for your changes.

2. Make the changes and update any affected documentation.

3. Compile MAGRATHEA from the repository root with:

   ```bash
   make -B
   ```

4. Run the automated tests with:

   ```bash
   make test
   ```

5. Run any additional modes or examples affected by the change.

6. Open a pull request explaining what changed, why it is needed, and how it was tested.

Please keep pull requests focused and avoid combining unrelated changes. Changes to numerical behavior, defaults, or output formats should be clearly described.

## Citation and license

Citation instructions are provided in [CITATION.md](CITATION.md), and published studies using MAGRATHEA are listed in [citations.md](citations.md).

By submitting a contribution, you agree that it may be distributed under the repository license.
