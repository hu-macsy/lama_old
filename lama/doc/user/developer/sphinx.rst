.. no warning for not included in toc

:orphan:

Examples of Sphinx for User Guide
=================================

Lists
-----

- A bullet list item

- Second item

  - A sub item

- Third item


1) An enumerated list item

2) Second item

   a) Sub item

      i) Sub-sub item

3) Third item

#) Another enumerated list item

#) Second item  

Images
------

.. image:: ../_images/test.jpg

Links
-----

A sentence with links to Wikipedia_ and the `Linux kernel archive`_.

.. _Wikipedia: http://www.wikipedia.org/
.. _Linux kernel archive: http://www.kernel.org/

Another sentence with an `anonymous link to the Python website`__.

__ http://www.python.org/

A link to Lists_ with `Lists_`

Code
----

code block::

   int main()
   {
      std::cout << "Hallo" << std::endl;
   }

More *italic* **bold** Text or a small ``inline::code::snippet()``

Maths
-----

.. math::
   :nowrap:

   \begin{eqnarray}
      y    & = & ax^2 + bx + c \\
      f(x) & = & x^2 + 2xy + y^2
   \end{eqnarray}
   
And here :math:`Ax=b` we have some inline math.
