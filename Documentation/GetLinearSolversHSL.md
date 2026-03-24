
Get a library with HSL linear solvers for use with ipopt
========================================================

Before getting the source code make sure to read the terms of the licence, especially regarding sharing of the software.

1.	Download the solver
Go to : https://licences.stfc.ac.uk/product/coin-hsl
Click on HSL Academic Licence -> ORDER NOW
Create an account or Sign in if you already have an account.
You will receive a download link via email 1 or 2 days later.
Download CoinHSL 2024.05.15 (windows binaries) via de link received via email:
    
  	<img width="945" height="86" alt="image" src="https://github.com/user-attachments/assets/25c327ab-95e6-4a23-873f-8ea33b16d4f2" />
 
3.	Add the path of the solver to the system path of your computer
Once you unzipped the folder on your computer, you will see a subfolder named “bin”. It is important that you add this to the system path of your pc. If you don’t than windows will not know where the solver is and then you will run into problems.
Search in windows for “Edit the system environment variables” 
Click on “Path” in “System variables” and then “Edit” or Double click on “Path”
Add the path to you’re the subfolder “bin” in the list like this: C:\...\CoinHSL.v2024.5.15.x86_64-w64-mingw32-libgfortran5\bin (de path is depending on where you unzipped the              downloaded folder)
  	<img width="945" height="101" alt="image" src="https://github.com/user-attachments/assets/ffae9b3e-58a3-49ba-80af-f538d5ba3fbb" />

5.	Use the new solver.
In PredSim/main.m, add:

```
    S.solver.ipopt_options.hsllib = 'C:/…./bin/libhsl.dll'; (depending on where you unzipped the solver)
    S.solver.linear_solver = 'ma97';
```

Imporant! Use the most recent version of CASADI.  



