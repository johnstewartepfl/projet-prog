# Chemical Equation Equilibrator

*The Chemical Equation Equilibrator* is a Python program that balances and equilibrates chemical equations based on reactants and products provided as SMILES strings.

## Features

- *SMILES Input*: Accepts reactants and products as SMILES strings, providing a convenient way to input molecular formulas.
- *Chemical Equation Balancing*: Automatically balances chemical equations by solving an integer linear programming problem to find stoichiometric coefficients.
- *Display Reaction*: Formats and displays the balanced chemical reaction in a human-readable format.

## Requirements

- *NumPy*: Required for handling arrays and matrices.

```
pip install numpy
```

- *RDKit*: A cheminformatics toolkit used for parsing SMILES strings and generating 2D molecular structures.

```
pip install -c conda-forge rdkit
```

- *PuLP*: A linear programming library used for solving integer linear programming problems.

```
pip install pulp
```

- *Rxn-INSIGHT* relies on NumPy, Pandas, RDKit, RDChiral, and RXNMapper.

A virtual environment can be installed with Anaconda as follows:

```console
conda create -n rxn-insight python=3.10
conda activate rxn-insight
```

```
git clone https://github.com/schwallergroup/Rxn-INSIGHT.git
cd Rxn-INSIGHT
pip install .
```

Or, for developing with the optional dependencies, which are required to run the tests
and build the docs:
``` 
pip install -e ".[test,doc]"
```

All of the test environments can be run using the command `tox` from the top directory.
Alternatively, individual test environments can be run using the `-e` flag as 
in `tox -e env-name`. To run the tests, tests with coverage report, style checks, and
docs build, respectively:
```
tox -e py3
tox -e py3-coverage
tox -e style
tox -e docs
```

## Usage

1. Ensure you have the required dependencies installed: NumPy, RDKit, and PuLP.
2. Run the Python script chemical_equation_equilibrator.py.
3. Follow the prompts to input the number of reactants and products, as well as the SMILES strings for each.
4. The program will automatically balance the chemical equation and display the balanced reaction.


##  Reference

`M. R. Dobbelaere, I. Lengyel, C. V. Stevens, and K. M. Van Geem, 
‘Rxn-INSIGHT: fast chemical reaction analysis using bond-electron matrices’, J. Cheminform., vol. 16, no. 1, Mar. 2024.`

```
@ARTICLE{Dobbelaere2024-es,
  title     = "{Rxn-INSIGHT}: fast chemical reaction analysis using
               bond-electron matrices",
  author    = "Dobbelaere, Maarten R and Lengyel, Istv{\'a}n and Stevens,
               Christian V and Van Geem, Kevin M",
  journal   = "J. Cheminform.",
  publisher = "Springer Science and Business Media LLC",
  volume    =  16,
  number    =  1,
  month     =  mar,
  year      =  2024,
  copyright = "https://creativecommons.org/licenses/by/4.0",
  language  = "en"
}
```
