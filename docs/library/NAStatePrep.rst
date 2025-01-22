Neutral Atom Logical State Preparation
======================================

This module provides functionality to generate an optimal operation sequence for zoned neutral atom architectures.
The input must be a circuit of single-qubit gates, followed by a set of entangling gates (CZ) and finally single-qubit gates on selected qubits.
Those circuits arise in the realm of error correction when the initial logical state must be prepared.

The process is divided into three steps:
1. Extract a list of qubit pairs from the circuit that represents the entangling gates with :code:`get_ops_for_solver`

    .. currentmodule:: mqt.qmap.na.state_preparation
    .. autofunction:: get_ops_for_solver

2. Supply the list of entangling operations to the solver and generate the optimal operation sequence with the :code:`NAStatePreparationSolver`.
   For further details on the employed abstraction of the 2D plane in the solver, please refer to the corresponding article :cite:labelpar:`stadeOptimalStatePreparation2024`.

    .. currentmodule:: mqt.qmap.na.state_preparation
    .. autoclass:: NAStatePreparationSolver

3. Generate code from the solver's result with :code:`generate_code`

    .. currentmodule:: mqt.qmap.na.state_preparation
    .. autofunction:: generate_code
