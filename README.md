# project-fa19
CS 170 Fall 2019 Project
To run the code, pip install the following:
  
  NetworkX
  
  Google's OR-Tools

We don't think anything else needs to be installed that would not already be.

The code can be run normally by using "python3 solver.py inputs/ outputs/".
We had trouble with figuring out what flag to use to run solve_all, so solver.py will default to attempting to run on folders (we changed a line of code to say "if not args.all:". To fix this to run on single files, we just removed the "not". So, to run the code on multiple files don't add a flag.

We realized at the very end that we didn't verify that the path began at the starting location, so we have a lot of validation errors. If we had more time we could fix this, but we realized too late. So, that issue will present itself when the code runs.
