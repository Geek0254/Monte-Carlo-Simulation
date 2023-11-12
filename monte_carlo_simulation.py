
import geopandas as gpd
import random
from shapely.geometry import Point, Polygon
import numpy as np
import scipy.stats as st
from scipy.stats import t
from scipy.stats import ttest_ind

# Load shapefile for China
china_shapefile = gpd.read_file("gadm41_CHN_3.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
china_polygons = china_shapefile['geometry']

# Define the total area of the country
total_area = china_polygons.area.sum() / 10**6 # Convert to square kilometers

# Set the number of simulations to run
num_simulations = 100

# Loop through each simulation
count = 0
for i in range(num_simulations):
    # Generate a random longitude and latitude within the bounds of the country
    random_point = Point(random.uniform(china_shapefile.total_bounds[0], china_shapefile.total_bounds[2]),
                         random.uniform(china_shapefile.total_bounds[1], china_shapefile.total_bounds[3]))

    # Check if the point is within any of the polygons of the country
    for polygon in china_polygons:
        if polygon.intersects(random_point):
            count += 1
            break

# Calculate the estimated total land area of China
estimated_total_area = (count / num_simulations) * total_area
print("The estimated land area of China is: {:.4f} square kilometers.".format(estimated_total_area))

#calculating the CONFIDENCE INTERVAL for the estimated area

# Estimated land area for China
estimated_total_area = 8217456.5907

# Generate 1000 random samples from the estimated land area
num_samples = 1000
samples = np.random.normal(estimated_total_area, size=num_samples)

# Compute the mean and standard error of the samples
mean = np.mean(samples)
stderr = st.sem(samples)

# Compute the 95% confidence interval using the t-distribution with the upper and lower bounds
lower_bound, upper_bound = t.interval(0.95, len(samples)-1, loc=mean, scale=stderr)

#print confidence interval
print("The 95% confidence interval for the estimated land area of China is: [{:.4f}, {:.4f}] square kilometers.".format(lower_bound, upper_bound))



# Load shapefile for CANADA
canada_shapefile = gpd.read_file("gadm41_CAN_3.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
canada_polygons = canada_shapefile['geometry']

# Define the total area of the country
total_area = canada_polygons.area.sum() / 10**6 # Convert to square kilometers

# Set the number of simulations to run
num_simulations = 100

# Loop through each simulation
count = 0
for i in range(num_simulations):
    # Generate a random longitude and latitude within the bounds of the country
    random_point = Point(random.uniform(canada_shapefile.total_bounds[0], canada_shapefile.total_bounds[2]),
                         random.uniform(canada_shapefile.total_bounds[1], canada_shapefile.total_bounds[3]))

    # Check if the point is within any of the polygons of the country
    for polygon in canada_polygons:
        if polygon.intersects(random_point):
            count += 1
            break

# Calculate the estimated total land area of Canada
estimated_total_area = (count / num_simulations) * total_area
print("The estimated land area of Canada is: {:.4f} square kilometers.".format(estimated_total_area))

#calculating the CONFIDENCE INTERVAL for the estimated area
estimated_total_area = 18901666.6912

# Generate 1000 random samples from the estimated land area
num_samples = 1000
samples = np.random.normal(estimated_total_area, size=num_samples)

# Compute the mean and standard error of the samples
mean = np.mean(samples)
stderr = st.sem(samples)

# Compute the 95% confidence interval using the t-distribution with the upper and lower bounds
lower_bound, upper_bound = t.interval(0.95, len(samples)-1, loc=mean, scale=stderr)

#print confidence interval
print("The 95% confidence interval for the estimated land area of Canada is: [{:.4f}, {:.4f}] square kilometers.".format(lower_bound, upper_bound))



#HYPOTHESIS TEST OF CHINA AND CANADA USING BOOTSTRAP

# Estimated land area for Canada and China
estimated_total_area_canada = 18901666.6912
estimated_total_area_china = 8217456.5907

# Generate 1000 random samples for each country
num_samples = 1000
samples_canada = np.random.normal(estimated_total_area_canada, size=num_samples)
samples_china = np.random.normal(estimated_total_area_china, size=num_samples)

# Compute the t-statistic and p-value
t_stat, p_val = ttest_ind(samples_canada, samples_china)

#print results
print("The t-statistic is: {:.4f}.".format(t_stat))
print("The p-value is: {:.4f}.".format(p_val))

#print results USING IF STATEMENT
if p_val < 0.05:
    print("The null hypothesis can be rejected.")
else:
    print("The null hypothesis cannot be rejected.")



# Load shapefile for Columbia
columbia_shapefile = gpd.read_file("gadm41_COL_2.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
columbia_polygons = columbia_shapefile['geometry']

# Define the total area of the country
total_area = columbia_polygons.area.sum() / 10**6 # Convert to square kilometers

# Set the number of simulations to run
num_simulations = 100

# Loop through each simulation
count = 0
for i in range(num_simulations):
    # Generate a random longitude and latitude within the bounds of the country
    random_point = Point(random.uniform(columbia_shapefile.total_bounds[0], columbia_shapefile.total_bounds[2]),
                         random.uniform(columbia_shapefile.total_bounds[1], columbia_shapefile.total_bounds[3]))

    # Check if the point is within any of the polygons of the country
    for polygon in columbia_polygons:
        if polygon.intersects(random_point):
            count += 1
            break

# Calculate the estimated total land area of Colombia
estimated_total_area = (count / num_simulations) * total_area
print("The estimated land area of Columbia is: {:.4f} square kilometers.".format(estimated_total_area))

#calculating the CONFIDENCE INTERVAL for the estimated area
estimated_total_area = 415419.7510

# Generate 1000 random samples from the estimated land area
num_samples = 1000
samples = np.random.normal(estimated_total_area, size=num_samples)

# Compute the mean and standard error of the samples
mean = np.mean(samples)
stderr = st.sem(samples)

# Compute the 95% confidence interval using the t-distribution with the upper and lower bounds
lower_bound, upper_bound = t.interval(0.95, len(samples)-1, loc=mean, scale=stderr)

#print confidence interval
print("The 95% confidence interval for the estimated land area of Colombia is: [{:.4f}, {:.4f}] square kilometers.".format(lower_bound, upper_bound))



# Load shapefile for Egypt
egypt_shapefile = gpd.read_file("gadm41_EGY_2.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
egypt_polygons = egypt_shapefile['geometry']

# Define the total area of the country
total_area = egypt_polygons.area.sum() / 10**6 # Convert to square kilometers

# Set the number of simulations to run
num_simulations = 100

# Loop through each simulation
count = 0
for i in range(num_simulations):
    # Generate a random longitude and latitude within the bounds of the country
    random_point = Point(random.uniform(egypt_shapefile.total_bounds[0], egypt_shapefile.total_bounds[2]),
                         random.uniform(egypt_shapefile.total_bounds[1], egypt_shapefile.total_bounds[3]))

    # Check if the point is within any of the polygons of the country
    for polygon in egypt_polygons:
        if polygon.intersects(random_point):
            count += 1
            break

# Calculate the estimated total land area of Egypt
estimated_total_area = (count / num_simulations) * total_area
print("The estimated land area of Egypt is: {:.4f} square kilometers.".format(estimated_total_area))

#calculating the CONFIDENCE INTERVAL for the estimated area
estimated_total_area = 990600.1459

# Generate 1000 random samples from the estimated land area
num_samples = 1000
samples = np.random.normal(estimated_total_area, size=num_samples)

# Compute the mean and standard error of the samples
mean = np.mean(samples)
stderr = st.sem(samples)

# Compute the 95% confidence interval using the t-distribution with the upper and lower bounds
lower_bound, upper_bound = t.interval(0.95, len(samples)-1, loc=mean, scale=stderr)

#print confidence interval
print("The 95% confidence interval for the estimated land area of Egypt is: [{:.4f}, {:.4f}] square kilometers.".format(lower_bound, upper_bound))



#HYPOTHESIS TEST OF CHINA AND CANADA USING BOOTSTRAP

# Estimated land area for Canada and China
estimated_total_area_colombia = 415419.7510
estimated_total_area_egypt = 990600.1459

# Generate 1000 random samples for each country
num_samples = 1000
samples_colombia = np.random.normal(estimated_total_area_colombia, size=num_samples)
samples_egypt = np.random.normal(estimated_total_area_egypt, size=num_samples)

# Compute the t-statistic and p-value
t_stat, p_val = ttest_ind(samples_colombia, samples_egypt)

#print results
print("The t-statistic is: {:.4f}.".format(t_stat))
print("The p-value is: {:.4f}.".format(p_val))

#print results USING IF STATEMENT
if p_val < 0.05:
    print("The null hypothesis can be rejected.")
else:
    print("The null hypothesis cannot be rejected.")



# Load shapefile for Cuba
cuba_shapefile = gpd.read_file("gadm41_CUB_2.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
cuba_polygons = cuba_shapefile['geometry']

# Define the total area of the country
total_area = cuba_polygons.area.sum() / 10**6 # Convert to square kilometers

# Set the number of simulations to run
num_simulations = 10

# Loop through each simulation
count = 0
for i in range(num_simulations):
    # Generate a random longitude and latitude within the bounds of the country
    random_point = Point(random.uniform(cuba_shapefile.total_bounds[0], cuba_shapefile.total_bounds[2]),
                         random.uniform(cuba_shapefile.total_bounds[1], cuba_shapefile.total_bounds[3]))

    # Check if the point is within any of the polygons of the country
    for polygon in cuba_polygons:
        if polygon.intersects(random_point):
            count += 1
            break

# Calculate the estimated total land area of Cuba
estimated_total_area = (count / num_simulations) * total_area
print("The estimated land area of Cuba is: {:.4f} square kilometers.".format(estimated_total_area))

#calculating the CONFIDENCE INTERVAL for the estimated area
estimated_total_area = 51520.9337

# Generate 1000 random samples from the estimated land area
num_samples = 1000
samples = np.random.normal(estimated_total_area, size=num_samples)

# Compute the mean and standard error of the samples
mean = np.mean(samples)
stderr = st.sem(samples)

# Compute the 95% confidence interval using the t-distribution with the upper and lower bounds
lower_bound, upper_bound = t.interval(0.95, len(samples)-1, loc=mean, scale=stderr)

#print confidence interval
print("The 95% confidence interval for the estimated land area of Cuba is: [{:.4f}, {:.4f}] square kilometers.".format(lower_bound, upper_bound))



# Load shapefile for Iceland with higher resolution
iceland_shapefile = gpd.read_file("gadm41_ISL_0.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
iceland_polygons = iceland_shapefile['geometry']

# Define the total area of the country
total_area = iceland_polygons.area.sum() / 10**6 # Convert to square kilometers

# Set the number of simulations to run
num_simulations = 10

# Loop through each simulation
count = 0
for i in range(num_simulations):
    # Generate a random longitude and latitude within the bounds of the country
    random_point = Point(random.uniform(iceland_shapefile.total_bounds[0], iceland_shapefile.total_bounds[2]),
                         random.uniform(iceland_shapefile.total_bounds[1], iceland_shapefile.total_bounds[3]))

    # Check if the point is within any of the polygons of the country
    for polygon in iceland_polygons:
        if polygon.intersects(random_point):
            count += 1
            break

# Calculate the estimated total land area of Iceland
estimated_total_area = (count / num_simulations) * total_area
print("The estimated land area of Iceland is: {:.4f} square kilometers.".format(estimated_total_area))

#calculating the CONFIDENCE INTERVAL for the estimated area
estimated_total_area = 113931.3687

# Generate 1000 random samples from the estimated land area
num_samples = 1000
samples = np.random.normal(estimated_total_area, size=num_samples)

# Compute the mean and standard error of the samples
mean = np.mean(samples)
stderr = st.sem(samples)

# Compute the 95% confidence interval using the t-distribution with the upper and lower bounds
lower_bound, upper_bound = t.interval(0.95, len(samples)-1, loc=mean, scale=stderr)

#print confidence interval
print("The 95% confidence interval for the estimated land area of Iceland is: [{:.4f}, {:.4f}] square kilometers.".format(lower_bound, upper_bound))



#HYPOTHESIS TEST OF CUBA AND ICELAND USING BOOTSTRAP

# Estimated land area for Cuba and Iceland
estimated_total_area_cuba = 51520.9337
estimated_total_area_iceland = 113931.3687

# Generate 1000 random samples for each country
num_samples = 1000
samples_cuba = np.random.normal(estimated_total_area_cuba, size=num_samples)
samples_iceland = np.random.normal(estimated_total_area_iceland, size=num_samples)

# Compute the t-statistic and p-value
t_stat, p_val = ttest_ind(samples_cuba, samples_iceland)

#print results
print("The t-statistic is: {:.4f}.".format(t_stat))
print("The p-value is: {:.4f}.".format(p_val))

#print results USING IF STATEMENT
if p_val < 0.05:
    print("The null hypothesis can be rejected.")
else:
    print("The null hypothesis cannot be rejected.")





#OTHER SMALL COUNTRIES EVALUATION

# Load shapefile for Cuba with higher resolution
lebanon_shapefile = gpd.read_file("gadm41_LBN_1.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
lebanon_polygons = lebanon_shapefile['geometry']

# Define the total area of the country
total_area = lebanon_polygons.area.sum() / 10**6 # Convert to square kilometers

# Set the number of simulations to run
num_simulations = 5

# Loop through each simulation
count = 0
for i in range(num_simulations):
    # Generate a random longitude and latitude within the bounds of the country
    random_point = Point(random.uniform(lebanon_shapefile.total_bounds[0], lebanon_shapefile.total_bounds[2]),
                         random.uniform(lebanon_shapefile.total_bounds[1], lebanon_shapefile.total_bounds[3]))

    # Check if the point is within any of the polygons of the country
    for polygon in lebanon_polygons:
        if polygon.intersects(random_point):
            count += 1
            break

# Calculate the estimated total land area of lebanon
estimated_total_area = (count / num_simulations) * total_area
print("The estimated land area of lebanon is: {:.4f} square kilometers.".format(estimated_total_area))

# Load shapefile for Jamaica
jamaica_shapefile = gpd.read_file("gadm41_JAM_1.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
jamaica_polygons = jamaica_shapefile['geometry']

# Define the total area of the country
total_area = jamaica_polygons.area.sum() / 10**6 # Convert to square kilometers

# Set the number of simulations to run
num_simulations = 4

# Loop through each simulation
count = 0
for i in range(num_simulations):
    # Generate a random longitude and latitude within the bounds of the country
    random_point = Point(random.uniform(jamaica_shapefile.total_bounds[0], jamaica_shapefile.total_bounds[2]),
                         random.uniform(jamaica_shapefile.total_bounds[1], jamaica_shapefile.total_bounds[3]))

    # Check if the point is within any of the polygons of the country
    for polygon in jamaica_polygons:
        if polygon.intersects(random_point):
            count += 1
            break

# Calculate the estimated total land area of jamaica
estimated_total_area = (count / num_simulations) * total_area
print("The estimated land area of Jamaica is: {:.4f} square kilometers.".format(estimated_total_area))

# Load shapefile for Seychelles
Seychelles_shapefile = gpd.read_file("gadm41_syc_0.shp").to_crs("EPSG:3857")

# Convert the shapefile into polygons
Seychelles_polygons = Seychelles_shapefile['geometry']

# Define the total area of the country
total_area = Seychelles_polygons.area.sum() / 10**6 # Convert to square kilometers

# Calculate the estimated total land area of Seychelles
print("The estimated land area of Seychelles is: {:.4f} square kilometers.".format(total_area))





"""FUNCTION TO CALCULATE THE AREA OF ANY COUNTRY USING THE GADM ISO NAMING CODES"""

import os
import random
import zipfile
from math import pi, sin, cos, sqrt
from urllib.request import urlretrieve

import geopandas as gpd
import matplotlib.pyplot as plt


# Define constants
MIN_SAMPLE_POINTS = 10
MAX_SAMPLE_POINTS = 100
THRESHOLD_AREA = 1e5  # km^2

GADM_URL = 'https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_{}_shp.zip'


def get_country_shapefile(country_name):
    """Download and extract the shapefile for a given country."""
    country_code = country_name.upper()[:3]  # ISO 3166-1 alpha-3 code
    url = GADM_URL.format(country_code)
    filename = os.path.basename(url)
    foldername = os.path.splitext(filename)[0]
    if not os.path.exists(foldername):
        print(f'Downloading {filename}...')
        urlretrieve(url, filename)
        print(f'Extracting {filename}...')
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(foldername)
    return foldername


def get_land_area(polygons):
    """Calculate the land area of a GeoDataFrame of polygons."""
    projected = polygons.to_crs(epsg=3857)
    return projected.geometry.area.sum() / 1e6  # km^2



def generate_sample_points(polygons, num_points):
    """Generate random sample points within the bounds of a GeoDataFrame of polygons."""
    min_lon, min_lat, max_lon, max_lat = polygons.total_bounds
    points = []
    while len(points) < num_points:
        lon, lat = random.uniform(min_lon, max_lon), random.uniform(min_lat, max_lat)
        point = (lon, lat)
        if polygons.contains(gpd.points_from_xy([lon], [lat])).any():
            points.append(point)
    return points


def estimate_land_area(polygons):
    """Estimate the land area of a GeoDataFrame of polygons using the Monte Carlo method."""
    land_area = get_land_area(polygons)
    if land_area < THRESHOLD_AREA:
        return land_area  # Direct calculation for small areas
    elif land_area < 10 * THRESHOLD_AREA:
        num_points = MIN_SAMPLE_POINTS
    else:
        num_points = MAX_SAMPLE_POINTS
    points = generate_sample_points(polygons, num_points)
    land_points = sum(1 for lon, lat in points if polygons.contains(gpd.points_from_xy([lon], [lat])).any())
    area_factor = land_points / num_points
    return land_area / area_factor


def main():
    # Get country name input from user
    country_name = input('Enter the name of a country: ')

    # Download and extract shapefile for the country
    foldername = get_country_shapefile(country_name)

    # Read polygons from shapefile and calculate land area
    polygons = gpd.read_file(foldername)
    land_area = estimate_land_area(polygons)

    # Display results
    print(f'The estimated land area of {country_name} is {land_area:.2f} km^2')
    ax = polygons.plot(alpha=0.5, edgecolor='black')
    ax.set_title(country_name)
    plt.show()


if __name__ == '__main__':
    main()

import os
import random
import zipfile
from math import pi, sin, cos, sqrt
from urllib.request import urlretrieve

import geopandas as gpd
import matplotlib.pyplot as plt


# Define constants
MIN_SAMPLE_POINTS = 10
MAX_SAMPLE_POINTS = 100
THRESHOLD_AREA = 1e5  # km^2

GADM_URL = 'https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_{}_shp.zip'


def get_country_shapefile(country_name):
    """Download and extract the shapefile for a given country."""
    country_code = country_name.upper()[:3]  # ISO 3166-1 alpha-3 code
    url = GADM_URL.format(country_code)
    filename = os.path.basename(url)
    foldername = os.path.splitext(filename)[0]
    if not os.path.exists(foldername):
        print(f'Downloading {filename}...')
        urlretrieve(url, filename)
        print(f'Extracting {filename}...')
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(foldername)
    return foldername


def get_land_area(polygons):
    """Calculate the land area of a GeoDataFrame of polygons."""
    projected = polygons.to_crs(epsg=3857)
    return projected.geometry.area.sum() / 1e6  # km^2



def generate_sample_points(polygons, num_points):
    """Generate random sample points within the bounds of a GeoDataFrame of polygons."""
    min_lon, min_lat, max_lon, max_lat = polygons.total_bounds
    points = []
    while len(points) < num_points:
        lon, lat = random.uniform(min_lon, max_lon), random.uniform(min_lat, max_lat)
        point = (lon, lat)
        if polygons.contains(gpd.points_from_xy([lon], [lat])).any():
            points.append(point)
    return points


def estimate_land_area(polygons):
    """Estimate the land area of a GeoDataFrame of polygons using the Monte Carlo method."""
    land_area = get_land_area(polygons)
    if land_area < THRESHOLD_AREA:
        return land_area  # Direct calculation for small areas
    elif land_area < 10 * THRESHOLD_AREA:
        num_points = MIN_SAMPLE_POINTS
    else:
        num_points = MAX_SAMPLE_POINTS
    points = generate_sample_points(polygons, num_points)
    land_points = sum(1 for lon, lat in points if polygons.contains(gpd.points_from_xy([lon], [lat])).any())
    area_factor = land_points / num_points
    return land_area / area_factor


def main():
    # Get country name input from user
    country_name = input('Enter the name of a country: ')

    # Download and extract shapefile for the country
    foldername = get_country_shapefile(country_name)

    # Read polygons from shapefile and calculate land area
    polygons = gpd.read_file(foldername)
    land_area = estimate_land_area(polygons)

    # Display results
    print(f'The estimated land area of {country_name} is {land_area:.2f} km^2')
    ax = polygons.plot(alpha=0.5, edgecolor='black')
    ax.set_title(country_name)
    plt.show()


if __name__ == '__main__':
    main()

import os
import random
import zipfile
from math import pi, sin, cos, sqrt
from urllib.request import urlretrieve

import geopandas as gpd
import matplotlib.pyplot as plt


# Define constants
MIN_SAMPLE_POINTS = 10
MAX_SAMPLE_POINTS = 100
THRESHOLD_AREA = 1e5  # km^2

GADM_URL = 'https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_{}_shp.zip'


def get_country_shapefile(country_name):
    """Download and extract the shapefile for a given country."""
    country_code = country_name.upper()[:3]  # ISO 3166-1 alpha-3 code
    url = GADM_URL.format(country_code)
    filename = os.path.basename(url)
    foldername = os.path.splitext(filename)[0]
    if not os.path.exists(foldername):
        print(f'Downloading {filename}...')
        urlretrieve(url, filename)
        print(f'Extracting {filename}...')
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(foldername)
    return foldername


def get_land_area(polygons):
    """Calculate the land area of a GeoDataFrame of polygons."""
    projected = polygons.to_crs(epsg=3857)
    return projected.geometry.area.sum() / 1e6  # km^2



def generate_sample_points(polygons, num_points):
    """Generate random sample points within the bounds of a GeoDataFrame of polygons."""
    min_lon, min_lat, max_lon, max_lat = polygons.total_bounds
    points = []
    while len(points) < num_points:
        lon, lat = random.uniform(min_lon, max_lon), random.uniform(min_lat, max_lat)
        point = (lon, lat)
        if polygons.contains(gpd.points_from_xy([lon], [lat])).any():
            points.append(point)
    return points


def estimate_land_area(polygons):
    """Estimate the land area of a GeoDataFrame of polygons using the Monte Carlo method."""
    land_area = get_land_area(polygons)
    if land_area < THRESHOLD_AREA:
        return land_area  # Direct calculation for small areas
    elif land_area < 10 * THRESHOLD_AREA:
        num_points = MIN_SAMPLE_POINTS
    else:
        num_points = MAX_SAMPLE_POINTS
    points = generate_sample_points(polygons, num_points)
    land_points = sum(1 for lon, lat in points if polygons.contains(gpd.points_from_xy([lon], [lat])).any())
    area_factor = land_points / num_points
    return land_area / area_factor


def main():
    # Get country name input from user
    country_name = input('Enter the name of a country: ')

    # Download and extract shapefile for the country
    foldername = get_country_shapefile(country_name)

    # Read polygons from shapefile and calculate land area
    polygons = gpd.read_file(foldername)
    land_area = estimate_land_area(polygons)

    # Display results
    print(f'The estimated land area of {country_name} is {land_area:.2f} km^2')
    ax = polygons.plot(alpha=0.5, edgecolor='black')
    ax.set_title(country_name)
    plt.show()


if __name__ == '__main__':
    main()

