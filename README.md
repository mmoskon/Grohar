# Grohar: automated visualisation of genome-scale metabolic models and their pathways #

Grohar represents a computational tool for the visualisation of metabolic networks and their analysis using Matplotlib, networkx and cobra. Grohar was written in python together with pyqt GUI toolkit to provide a user friendly experience.  Grohar is licensed under the Creative Commons Attribution license.


### How do I get set up? ###

* Use *Python 3.5+* with *Spyder* (get the whole package with [Anaconda](https://docs.continuum.io/anaconda/install/))
* main file = *top.py*
* graphical user interface is located within the file *top.ui* and can be modified easily with Qt Designer (included within the package *pyqt-distutils*)

You will need to install the following python packages - additional to the packages that are preinstalled with Anaconda:

- *cobra 0.5.4 (pip install cobra==0.5.4),*
- *bioservices*
- *PyQt5*
- *msgpack*
- *python-libsbml*
- *bioservices*
- *plotly*
- *cobrababel*
- *pydot*
- *pydot-ng*

You need to configure IPython in Spyder to open the Matplotlib graphs in a separate window. In Spyder go to *Tools* -- *Preferences* -- *IPython console* -- *Graphics*. Here you will find a setting entitled *Graphics Backend*, which should be set to *Automatic* (its default value is *Inline*).

Additional software that is also required:

- *Gurobi*
    - *conda config --add channels http://conda.anaconda.org/gurobi*
    - *conda install gurobi*
    - You should also install the gurobi license for the solvers to work. The license and instructions are available at http://www.gurobi.com/academia/for-universities.
- Graphviz: 
    - *http://www.graphviz.org/download/* 
    - In Windows, add C:\Program Files (x86)\Graphviz2.38\bin to PATH


You are ready to go!

### How do I use Grohar? ###

The user manual with some examples is avaible at the following [link](grohar_manual.pdf).

### Who do I talk to? ###

You can direct all your questions, comments and critiques to [miha.moskon@fri.uni-lj.si](mailto:miha.moskon@fri.uni-lj.si).

### Who was Grohar? ###

Ivan Grohar was a Slovenian Impressionist painter. You can read more about his work on [Wikipedia](https://en.wikipedia.org/wiki/Ivan_Grohar).

```
#!python


```