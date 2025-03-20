# Method
This solver imports geometry from an ABAQUS inp file that has the boundary condition and geometry information. This is a quadrilateral solver using a 2x2 gauss quadrature similar to the CPS4 solver in ABAQUS. The mesh spacing can be adjusted to the accuracy required in the inp file. To run this solver, open all files in MATLAB and run the triangle.m file.
# Results
![image](https://github.com/user-attachments/assets/6f211ed9-b037-404a-9d67-5122cf542ba9)

The MATLAB results were as predicted, with the relative error of the numerical number decreasing as the mesh spacing decreased. There was a slight difference in the MATLAB results 
with the corresponding CPS4 results. The stress value, like the other methods, demonstrated a higher accuracy than the displacement. The computation time of the method was to note as the medium mesh took 7 hours to compile while the fine mesh needed to compile and execute overnight. The computational expense was vast compared to the ABAQUS solver while 
demonstrating a lower accuracy. The trend of convergence was consistent.

![image](https://github.com/user-attachments/assets/14609b64-4b3a-48db-b9f0-46b59b36d158)


