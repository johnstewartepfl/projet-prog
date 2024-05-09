# Chemical Equation Equilibrator

The Chemical Equation Equilibrator is a Python program that balances and equilibrates chemical equations based on reactants and products provided as SMILES (Simplified Molecular Input Line Entry System) strings.

## Features

- **SMILES Input**: Accepts reactants and products as SMILES strings, providing a convenient way to input molecular formulas.
- **Chemical Equation Balancing**: Automatically balances chemical equations by solving an integer linear programming problem to find stoichiometric coefficients.
- **Display Reaction**: Formats and displays the balanced chemical reaction in a human-readable format.

## Requirements

- **NumPy**: Required for handling arrays and matrices.
- **RDKit**: A cheminformatics toolkit used for parsing SMILES strings and generating 2D molecular structures.
- **PuLP**: A linear programming library used for solving integer linear programming problems.

## Usage

1. Ensure you have the required dependencies installed: NumPy, RDKit, and PuLP.
2. Run the Python script `chemical_equation_equilibrator.py`.
3. Follow the prompts to input the number of reactants and products, as well as the SMILES strings for each.
4. The program will automatically balance the chemical equation and display the balanced reaction.
