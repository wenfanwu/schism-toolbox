This HYCOM database was generated using the function `get_hycom_online.m`.

## Database Format

To ensure compatibility with post-processing and interpolation routines, each data file must follow the structure below:

1. **Dimensional Variables**  
   
   - `lon`, `lat`, and `depth` must all be **monotonically increasing**.  
   - `depth` values must be **positive** and consistent across all files.

2. **Variable Dimensions**  
   
   - Variables must be stored as either:  
     - **2-D** arrays of size `lon × lat`, or  
     - **3-D** arrays of size `lon × lat × depth`.

3. **File Naming**  
   
   - Each file must include a **timestamp** in its filename to indicate its corresponding time step (e.g., `HYCOM_20200101.mat`).

> ✅ These formatting rules ensure efficient data loading, seamless integration with interpolation routines, and support for parallel extraction.