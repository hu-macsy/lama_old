Source Code Formatting
======================

The formatting rules of LAMA are described in the `styleGuide_`

Existing code might be adapted to our formatting rules in one
of the following ways:

- Source code formatting in Eclipse (Shift+Ctrl+F)
- Using `astyle`_: script available at tools/lama_format

.. _astyle: http://astyle.sourceforge.net

Source code formatting within Eclipse is more rudimentary while the astyle
tool helps much better adapting the code.

In the following the different options used of astyle are shortly explained:

- **--style=ansi** guarantees our proposed arrangement of the curly braces.
- **--indent=spaces=4** supports our indentation rules
- **--convert-tabs**
- **--indent-switches**
- **--indent-preprocessor**
- **--indent-switches**
- **--break-blocks**
- **--convert-tabs**

 * **--pad-header**    

 InsertS space padding after paren headers (e.g. 'if', 'for'...).

::

	if( a > 0)

 becomes

::

	if ( a > 0 )

The following opitions are especially useful as they are not supported by Eclipse:

- **--add-brackets**    

 AddS brackets to unbracketed one line conditional statements.

::

	if ( c > 0 ) x = 4;

 becomes:

::

	if ( c > 0 )
    {
        x = 4;
    }

 * {{{--pad-oper}}}

 InsertS space paddings around operators.

::

	x = 4*( a+b )/3;


 becomes:

 ::
 
 	x = 4 * ( a + b ) / 3;

 * {{{--pad-paren-in}}}

 Inserts space padding around parenthesis on the inside only.

::

	if (x > 0) a = (a + b) * 3

 becomes: 

 ::
 
	if ( x > 0 ) a = (a + b) * 3

Still to be discussed:

- **--delete-empty-lines**

 Deletes empty lines within a function or method.
 It will NOT delete lines added by the break-blocks options.


- **--align-pointer=type**


::

	int *x;
	int &a;

becomes

::

	int* x;
	int& a;

Unsupported Features
--------------------


