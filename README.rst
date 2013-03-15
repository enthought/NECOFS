This is a simple GUI wrapper for the NECOFS numerical model.  Includes a traits_ model and user interfaces in TraitsUI_ and enaml_.  The TraitsUI version should run on released versons of EPD_, while the enaml verion will require installing a more recent version of enaml available from GitHub_ (tested with 0.6.8).

.. _traits: http://docs.enthought.com/traits/index.html
.. _TraitsUI: http://docs.enthought.com/traitsui/index.html
.. _enaml: http://docs.enthought.com/enaml/index.html
.. _EPD: http://www.enthought.com/products/epd.php
.. _GitHub: https://github.com/enthought/enaml

This model code is based on an ipython notebook version at http://nbviewer.ipython.org/5092905

NECOFS (Northeastern Coastal Ocean Forecast System) is run by groups
at the University of Massachusetts Dartmouth and the Woods Hole
Oceanographic Institution, led by Drs. C. Chen, R. C. Beardsley,
G. Cowles and B. Rothschild. Funding is provided to run the model by
the NOAA-led Integrated Ocean Observing System and the State of
Massachusetts.

NECOFS is a coupled numerical model that uses nested weather models, a
coastal ocean circulation model, and a wave model. The ocean model is
a volume-mesh model with horizontal resolution that is finer in
complicated regions. It is layered (not depth-averaged) and includes
the effects of tides, winds, and varying water densities caused by
temperature and salinity changes.

* Model description: http://fvcom.smast.umassd.edu/research_projects/NECOFS/model_system.html
* THREDDS server with other forecast and archive products: http://www.smast.umassd.edu:8080/thredds/catalog.html

FILES
=====

`necofs_model.py` : Traits model to express the code.  Running this file as a script will run the TraitsUI version of the GUI.
`enaml_gui.py` : Runs the enaml version of the GUI
`ocean_view.enaml` : defines the enaml GUI (used by `enaml_gui.py`)
`necofs_velocity.py` : export of original iPython notebook to a python script
