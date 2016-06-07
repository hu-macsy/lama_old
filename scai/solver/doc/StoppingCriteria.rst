.. _stopping-criteria:

Stopping Criteria
-----------------

LAMA has three different stopping criteria that can be combined by conjunction (&&) and disjunction (||).
The solver stops if the criteria can be evaluated to true after an iteration.

Possible stopping criteria:

- IterationCount: the iteration count starts with zero with the last solve call and is incremented in every iteration

- ResidualStagnation: the residual stagnation is the absolute difference between the residual of the last two iterations

- ResidualThreshold: the residual threshold is the absolute residual

.. code-block:: c++

    Scalar eps = 0.00001;
    NormPtr norm = NormPtr( new L2Norm() );
    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );

    mySolver.setStoppingCriterion( rt );
