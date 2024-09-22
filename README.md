# Texture Parameterization
## Introduction
This application uses Floater's algorithm to create a 2D representation of a 3D mesh. This process is used to create texture coordinates.

## Compilation

This application was made using Linux and QT 5.15.3, which can be downloaded at  
https://www.qt.io/download-dev  

If using a Windows machine, the Linux terminal can be accessed by using WSL
Learn more here: https://learn.microsoft.com/en-us/windows/wsl/

To compile the program, enter the following commands in a terminal:  
    
    qmake -project QT+=opengl LIBS+=-lGLU
    qmake
    make

## Usage
Run the program using
`./Texture-Parameterization filename orientation`

### Arguments
`filename` 
- The mesh to import (must be **.obj** file)
- Whilst not essential, vertex colour values (vc) are recommended. These can be used for creating texture, but will also make visualising the mesh easier

 `orientation`
 - The 2 axes from the 3D mesh that will carry over to the 2D space
 - Can be either `xy`, `xz` or `zy` (any other input is rejected)

### Interface
The interface contains settings to change how the mesh is viewed including:
- An arcball to rotate the model
- Sliders that move the mesh in the X and Y directions
- Zoom slider

There are also many checkboxes which change the rendering mode. These are:
- `Wireframe` - Render the wireframe
- `Nrm -> RGB` - Render the vertex normals as colours
- `UVW -> RGB` - Render the UVW values (texture coordinates) as colours
- `Texture` - Render the texture map 
- `Normal Map` - Renders the normal map

![Image](/assets/texture%20parameterization.PNG)

### Exporting
`Export OBJ` exports the mesh to a new **.obj** file, called `exported.obj`. This file includes the newly generated texture coordinates (**vt** values).

`Export Texture` exports the texture map (i.e. the vertex colours placed on the 2D map) as a PNG called `texture.png`