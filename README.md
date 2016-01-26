##The GazeReader Toolbox##

GazeReader is a Matlab toolbox for the statistical analysis of eye-tracking data. It's built around a generalized linear point-process modeling framework. 
In addition to tools for statistical analysis, it provides a user interface from which to load data, specify trial events, and data and results
superimposed on back-ground representations of the scene and regions of interest. 

Please refer to

<dl>
<p>&emsp;Kovach, CK <i>A Generalized Linear Model for Eye Movements</i> PhD thesis, University of Iowa, (2008).
</p>
<p>or</p>
<p>&emsp;Kovach CK, Adolphs R (2015) Investigating attention in complex visual search. <i>Vision Research</i> 116, Part B:127-141.
</p></dl>

To get started with the graphical interface at the Matlab command line, use:

```>> h = GazeReader ```

####Note####
Much of GazeReader is organized around a graphical user interface developed in GUIDE (a GUI development toolbox in matlab).
Initial development also predated much of Matlab's current object-oriented capability, and the closest thing available at
the time to a handle class object was to use the "setappdata" and "getappdata" functions to store data in fields associated 
with a graphics object handle. As a result, the processing pipeline is somewhat inflexibly tied to the graphics objects. The 
graphical user interfaces are also defined in figure files, so the toolbox is not self-contained within the code. If I ever have 
a chance to work on an updated version, I'll reorganize things to use matlab classes and define the user interfaces in code rather
than through GUIDE figures. Also the documentation will need some fleshing out - I'm happy to answer questions about how to use 
the toolbox in the meantime. 
