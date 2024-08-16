import pandas as pd
import numpy as np

def interpolate_points(p1, p2, max_dist):
    """Interpolates points between p1 and p2 if the distance is greater than max_dist."""
    dist = np.linalg.norm(np.array(p2) - np.array(p1))
    if dist <= max_dist:
        return [p1, p2]
    
    num_interpolations = int(np.ceil(dist / max_dist))
    interpolated_points = [p1]
    
    for i in range(1, num_interpolations):
        ratio = i / num_interpolations
        new_point = np.array(p1) + ratio * (np.array(p2) - np.array(p1))
        interpolated_points.append(new_point.tolist())
    
    interpolated_points.append(p2)
    return interpolated_points

def remove_duplicates(points, tol=1e-8):
    """Removes duplicate points from a list, within a tolerance."""
    unique_points = []
    for point in points:
        if not unique_points or np.linalg.norm(np.array(point) - np.array(unique_points[-1])) > tol:
            unique_points.append(point)
    return unique_points

def process_csv(input_file, output_file, max_dist, scaling):
    # Read CSV file
    df = pd.read_csv(input_file)
    points = df[['x', 'y']].values.tolist()

    # Apply scaling to each point individually
    points = [[scaling * coord for coord in point] for point in points]

    # Interpolate points
    new_points = []
    for i in range(len(points) - 1):
        p1 = points[i]
        p2 = points[i + 1]
        interpolated = interpolate_points(p1, p2, max_dist)
        new_points.extend(interpolated)
    
    # Remove duplicate points
    new_points = remove_duplicates(new_points)

    # Create a new DataFrame
    new_df = pd.DataFrame(new_points, columns=['x', 'y'])
    new_df.to_csv(output_file, index=False, header=False)

# Example usage
process_csv('interlagos.csv', 'track.csv', max_dist=0.5, scaling=1.25)
