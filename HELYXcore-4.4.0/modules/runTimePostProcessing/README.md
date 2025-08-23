# README #

### What is this repository for? ###

This repository contains the runtime postprocessing core module, 
which is responsible for writing image files during the execution of the case.

### How do I get set up? ###

1. Clone this repository into the modules folder of the core
2. Compile and pack the EVTK (use the branch with the same name as the branch you are using for the RTPP)
3. Open the userSettings.cmake of the core (etc folder) and add the VTK settings if they are not already there:
   ```
   set(VTK_REQUIRED
       ON
       CACHE STRING
       "Compile VTK"
       FORCE
   )
   set(VTK_ARCH_PATH
       "path to the folder where you packed the EVTK, for example the customVTK folder"
       CACHE PATH
       "Default path on which to look for the EVTK library"
       FORCE
   )
   set(VTK_RENDERING_BACKEND
       OSMesa
       CACHE STRING
       "The VTK rendering backend (one of OSMesa, OpenGL1, OpenGL2)"
       FORCE
   )
   ```
4. Replace the path to the EVTK folder in the settings above. For example: 
```
"~/git/EVTK-d.e.v/SRC/COMPILE/LINUX/customEVTK"
```
5. Compile the core.

To use it, create a Write Image function object in the GUI, and then add the desired scenes to it. 

### Contribution guidelines ###

* CLion is recommended
* There are tests implemented with HelyxCI, with the rtpp repo.

### Who do I talk to? ###

Emanuel Cesconeto - e.cesconeto@engys.com