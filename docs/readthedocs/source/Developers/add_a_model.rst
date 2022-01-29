##################
Add a New Model
##################

New models can be added to SELF by creating classes that are type extensions of the :code:`Model1D`, :code:`Model2D`, or :code:`Model3D` classes that are defined in :code:`src/SELF_Model.f90`. Each of these base classes are abstract, but they provide a lot of the "boiler-plate" type-bound procedures that help with time integration, tendency calculation, and file I/O. However, these abstract classes have a set of deferred type-bound procedures that you are required to define in order to create a model for a specific conservation law. 


To help you get started, template source code is provided for 1-D, 2-D, and 3-D models in the :code:`src/templates` directory.


Create a New Module
####################

Copy a Template
----------------

Add any Necessary Variables
-----------------------------

Concretize the Flux Method
-----------------------------

Concretize the Riemann Flux Method
------------------------------------

Concretize the Source Method
-----------------------------


Override Model Procedures
######################################
In some cases, the provided type-bound procedures from the parent :code:`Model` class are insufficient for your problem. For example, if you want to implement a split-form flux divergence calculation, rather than a conservative-form flux divergencce, you will need to override the :code:`FluxDivergence` method. This section walks through overriding existing type-bound procedures so that you can further customize your model.

Add GPU Acceleration to your Model
=====================================
If you are starting from one of the provided templates, you may have noticed that the :code:`EnableGPUAccel` procedure always forces GPU acceleration to be disabled. This is done to help enforce correctness in your model, regardless of hardware, until you are ready to fully implement GPU acceleration for your model. To accomplish this, you will need to define HIP kernels for each of the type-bound procedures that you define in your model and link these in using :code:`ISO_C_BINDING` and external interface declarations.


