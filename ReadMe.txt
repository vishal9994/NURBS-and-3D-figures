Readme File NURBS
Project Under guidance of Dr. PK Jain Sir
Dated: 10/6/2020
-------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------

The folder consists of four MATLAB files -

-: NURBS.m
   Contains the main code and file to be run. It plots the NURBS curve based on the input parameters.

-: DeBoor.m
   input Degree,Knot Vector,Start point, End Point
   Contains the Cox De Boor algorithm that returns the Basis function and the 't' matrix

-: DeBoor1.m   
   input Degree,Knot Vector,Start point, End Point, Tvalue
   Contains the Cox De Boor algorithm that returns the coordinates at particular value of T.

-: plot2svg.m
   This file just save the plot in scalar vector graphic format.
   File imported from % http://www.mathworks.nl/matlabcentral/fileexchange/7401-scalable-vector-graphics-svg-export-of-figures

-------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------- 

Steps to follow-:

-Run the NURBS.m file
-Select the points(by default 5) on the figure screen using left click of mouse.
-Press 'X' on the Keyboard to display set of commands that can be used to maipulate the curve.

Using Keyboard and Mouse commands you can manipulate the curve. The commands are-:

Keyboard / Mouse commands-:

	'a'or'A'
	Show axis on plot 

*******************************************************************************************************************	

	'b'or'B'
	Plot Basis of the curve 
	'Shift + 'B''
	Plots full basis of the curve

*******************************************************************************************************************

	'c'or'C'
	Revolves the curve around X axis to form a 3D figure.	

*******************************************************************************************************************	

	'd'or'D'
	Changes degree
    	To change first change the knot vector and then add/ delete the points
    	and then perform the operation of changing degree,On command prompt 
    	enter the degree of the curve,which will add/ delete the points according
    	to the input conditions.

*******************************************************************************************************************

        'e'or'E'
	Select mode to Edges
	Using right click of mouse select the edge of the control polygon and using
	the scroll change the knots

*******************************************************************************************************************

        'h'or'H'
	Plot First Derivative of B Spline

*******************************************************************************************************************

	'm'or'M'
	Move the control point
	Change mode to vertices then using right click select the vertex of control 
	polygon and then press G button on keyboard to move the point and left click 
	to release the point.

*******************************************************************************************************************

	'i'or'I'
	Shows Information  about curve on the command window.

*******************************************************************************************************************

	'j'or'J'
	Show/Hide the joints on the curve

*******************************************************************************************************************

	'k'or'K'
	change the knot vector of the curve
	Enter the new knot vector on command window.    eg.  [0,1,1,1,2,3,4,5,5,5,6]

*******************************************************************************************************************

	'l'or'L'
	Show/Hide labels on the curve like weights,knots.

*******************************************************************************************************************

	'p'or'P'
	Assign the coordinates to the vertex
	Using the left click select the vertex and press P button on keyboard and on
	command window enter the coordinates.           eg.   [2,3]

*******************************************************************************************************************

        'q'or'Q'
	Quit and close the window

*******************************************************************************************************************

        's'or'S'
	Save plot in .svg format

*******************************************************************************************************************

	't'or'T'
	Get coordinates of point at particular value of t
	Press T on keyboard and on command window enter value of t.  eg.  1.3

*******************************************************************************************************************

	'v'or'V'
	Select mode to Vertices. Using the mouse scroll you can change the weights of
	selected control point.  		        eg.   1.5

*******************************************************************************************************************

	'w'or'W'
	Assign weight to the selected control point
	Change mode to vertices then select the control point using right mouse click
	and press W on keyboard and on command window enter the new weight.

*******************************************************************************************************************

	'x'or'X'
	Shows a list adjacent to a plot to show a commands that can be used to manipulate the curve.

*******************************************************************************************************************

        'z'or'Z'
	Extrude curve along any vector			eg. [0 0 10]


*******************************************************************************************************************

        
	
