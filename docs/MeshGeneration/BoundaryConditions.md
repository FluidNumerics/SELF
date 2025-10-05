# Boundary Conditions

In most mesh generation software, you are given the ability to align a name (string/character) and an integer id with each face, node, or edge that lies on a physical boundary. This information is often used in physical modeling software, like SELF, to implement the appropriate boundary conditions for your model. SELF provides a flexible framework for mapping boundary condition names and integer ids to procedures that you can implement for specific boundary conditions. This section of the documentation provides an overview of how you can register custom boundary conditions for a model built with SELF. Additionally we'll cover how you will want to register boundary conditions in a way that is consistent with the meshes that you create in other external software (e.g. HOHQMesh).

## How boundary conditions are registered in SELF
