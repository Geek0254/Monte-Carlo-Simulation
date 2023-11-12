# Monte Carlo Land Area Estimation

This Python script utilizes Monte Carlo simulation to estimate the land area of various countries using their latitude and longitude coordinates. The script uses geospatial data in the form of shapefiles, specifically the GADM dataset. The estimation is performed by generating random points within the bounds of each country and determining whether these points fall within the country's polygons.

## Prerequisites
- Python 3.x
- Required Python packages: `geopandas`, `shapely`, `numpy`, `scipy`

## Usage
1. Clone the repository or download the Python script.
2. Ensure that the required Python packages are installed using the following:
    ```bash
    pip install geopandas shapely numpy scipy
    ```
3. Download the GADM shapefiles for the countries you want to analyze from [GADM website](https://gadm.org/download_world.html).
4. Update the script with the correct file paths for the shapefiles.

## How to Run
Run the script using a Python interpreter. It will prompt you to enter the name of a country, and it will display the estimated land area along with a visual representation of the country's polygons.

```bash
python monte_carlo_land_area.py
```

## Description
The script is organized into sections for different countries. Each section follows a similar structure:

1. Load Shapefile: Load the GADM shapefile for the specific country.

2. Generate Random Points: Using Monte Carlo simulation, generate random points within the bounding box of the country.

3. Check Intersection: Check if each random point falls within any of the country's polygons.

4. Estimate Land Area: Calculate the estimated land area based on the ratio of points inside the polygons.

5. Print Results: Display the estimated land area and a 95% confidence interval.

6. Hypothesis Test: Conduct a hypothesis test (t-test) between two countries to compare their estimated land areas.

7. Function for Generalization: The script includes a function to generalize the process for any country using its GADM ISO naming code.

## Note
- Ensure that the GADM shapefiles are in the correct format and have the necessary permissions for reading.
- The accuracy of the estimation depends on the number of simulations; adjust `num_simulations` for more accurate results.

Feel free to explore and modify the script for your specific use case or additional functionalities.
