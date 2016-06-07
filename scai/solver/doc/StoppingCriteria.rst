.. _stopping-criteria:

Stopping Criteria
-----------------

For running the iterative solvers at least one stopping criteria should be set. Otherwise the default is running one iteration, so you can step through.

SCAI solver has three different stopping criteria types that can be combined by conjunction (&&) and disjunction (||).
The solver stops if the criteria can be evaluated to true after an iteration.

Available stopping criteria types are:

* ``IterationCount``: Evaluates to true, if the iteration count is reached. The iteration count starts with zero with the last solve call and is incremented in every iteration.
* ``ResidualThreshold``: Evaluates to true, if the residual value falls below the given epsilon.
* ``ResidualStagnation``: Evaluates to true, if the difference between the residual of the last two iterations falls below the given epsilon.

For the ``ResidualThreshold`` and ``ResidualStagnation`` criteria you have to define a norm for the evalutation of the residual and the border that should be passed. For the ``ResidualThreshold`` criteria you also have to define the check mode: it may be ``Absolute`` for absolute residue reduction, ``Relative`` for relative residue reduction or ``Divergence`` for relative residue increase.

.. code-block:: c++

    Scalar eps = 0.00001;
    NormPtr norm = NormPtr( new L2Norm() );
    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );

    mySolver.setStoppingCriterion( rt );

	CriterionPtr it( new IterationCount( 10000 ) );
	CriterionPtr oneOfBoth = rt || it;

	CriterionPtr rs( new ResidualStagnation( norm, eps ) );
	CriterionPtr both = rt && rs;
